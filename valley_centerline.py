import argparse
import shapely
import pandas as pd
import geopandas as gpd
import networkx as nx
import momepy
import os

DEBUG = True
BUFFER_DISTANCE = 0.0001
SUBDIVISION_AMOUNT = 200

def main():
    # Get arguments for input and output files
    parser = argparse.ArgumentParser(description='Find the centerline between two walls.')
    parser.add_argument("input", help="the filepath of the feature containing the two walls between which the centerline will be found")
    parser.add_argument("output", help="the filepath where the centerline will be saved")


    args = parser.parse_args()
    input = gpd.read_file(args.input)
    output = args.output

    # Import wall vector file

    #Explode function ensures each wall is a LineString rather than a MultiLineString
    walls = input.explode()
    extract_centerline(walls, output)

def extract_centerline(walls, output):
    if not os.path.exists(os.path.join(os.getcwd(), "Output")):
        os.makedirs(os.path.join(os.getcwd(), "Output"))
    crs = walls.crs

    # Convert valley walls into a collection of points for voronoi algorithm

    coords = walls.iloc[0].geometry.coords[:] + walls.iloc[1].geometry.coords[::-1]
    valleypoly = shapely.geometry.Polygon(coords)
    start_boundary = shapely.geometry.LineString([walls.iloc[0].geometry.coords[0], walls.iloc[1].geometry.coords[0]])
    end_boundary = shapely.geometry.LineString([walls.iloc[0].geometry.coords[-1], walls.iloc[1].geometry.coords[-1]])

    # Handle the case of both walls being drawn in opposite directions
    if not valleypoly.is_valid:
        coords = walls.iloc[0].geometry.coords[:] + walls.iloc[1].geometry.coords[:]
        valleypoly = shapely.geometry.Polygon(coords)
        start_boundary = shapely.geometry.LineString([walls.iloc[0].geometry.coords[0], walls.iloc[1].geometry.coords[-1]])
        end_boundary = shapely.geometry.LineString([walls.iloc[0].geometry.coords[-1], walls.iloc[1].geometry.coords[0]])


    start_pole = start_boundary.interpolate(0.5, normalized=True)
    end_pole = end_boundary.interpolate(0.5, normalized=True)

    #coords = coords + extra_coords
    points = shapely.MultiPoint([shapely.Point(i[0], i[1]) for i in coords])

    # Generate polygon and buffer from the input walls
    buffer = valleypoly.buffer(BUFFER_DISTANCE)

    # Create Voronoi Polygons

    voronoi_edges = {'features': 0, 'geometry': [shapely.voronoi_polygons(points, extend_to=buffer, only_edges = True)]}
    voronoi_edges = gpd.GeoDataFrame(voronoi_edges, crs=crs).explode()
    voronoi_edges = voronoi_edges.sjoin(walls, how='left', predicate='intersects')
    voronoi_edges = voronoi_edges.loc[voronoi_edges['index_right'].isna()].reset_index().drop(columns=['features', 'index_right', 'begin', 'end', 'begin_2', 'end_2', 'index'])

    print("Voronoi Edges created.")


    ########################################
    # Find potential start/end edges       #
    ########################################

    features = [0, 1]
    geometry = [start_pole, end_pole]
    df = {'features': features, 'geometry': geometry}
    poles_gdf = gpd.GeoDataFrame(df, crs=crs)

    ########################################
    # Find nearest edges to poles          #
    ########################################
    nearest_edges = voronoi_edges.iloc[gpd.sjoin_nearest(poles_gdf, voronoi_edges)['index_right']].reset_index()
    cols = list(nearest_edges)
    print(cols)
    print("Found nearest edges!")

    ########################################
    # Find shortest path from start to end #
    ########################################

    graph = momepy.gdf_to_nx(voronoi_edges)

    start = nearest_edges['geometry'][0].coords[0]
    end = nearest_edges['geometry'][1].coords[0]
    path = shapely.LineString(nx.shortest_path(graph, source=start, target=end))
    print('Path found')

    ########################################
    # Export centerline                    #
    ########################################

    df = {'features': [0], 'geometry': path}
    centerline_gdf = gpd.GeoDataFrame(df, crs=crs)
    centerline_gdf.to_file(output)
    print("Exported")

    return centerline_gdf

if __name__ == "__main__":
    main()