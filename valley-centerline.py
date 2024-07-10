import argparse
import shapely
import pandas as pd
import geopandas as gpd
import networkx as nx
import momepy

DEBUG = True
BUFFER_DISTANCE = 0.0001
SUBDIVISION_AMOUNT = 200

def main():
    '''
    parser = argparse.ArgumentParser(description='Find the centerline between two lines.')
    parser.add_argument("input_walls.iloc[0].geometry", help="the filepath of the valley wall ")
    parser.add_argument("output", help="the filepath for the output shapefile")


    args = parser.parse_args()
    input_1 = gpd.read_file(args.input_walls.iloc[0].geometry)
    input_2 = gpd.read_file(args.input_walls.iloc[1].geometry)
    output = args.output
    '''

    # Import line shapefiles

    #Explode function ensures each wall is a LineString rather than a MultiLineString
    walls = gpd.read_file("Input\SampleData\Whitewater\WhitewaterWalls.gpkg").explode()
    find_centerline(walls)

def find_centerline(walls):
    crs = walls.crs

    '''
    def subdivide_wall(wall, subdivision_distance):
        subdivided_coords = []
        dist = subdivision_distance
        while dist < wall.length:
            subdivided_coords.append(wall.interpolate(dist))
            dist += subdivision_distance
        return subdivided_coords
    '''

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
    if DEBUG:
        df = {'features': [0], 'geometry': [valleypoly]}
        gdf = gpd.GeoDataFrame(df, crs=crs)
        gdf.to_file('Output\Valley_Polygon.gpkg')


        df = {'features': [0], 'geometry': [buffer]}
        buff = gpd.GeoDataFrame(df, crs=crs)
        buff.to_file('Output\ValleyBuffer.gpkg')

    # Create Voronoi Polygons

    voronoi_edges = shapely.voronoi_polygons(points, extend_to=buffer, only_edges = True)

    ## Clean this up later
    geom = [voronoi_edges]
    print(len(geom))

    voronoi_edges_df = {'features': 0, 'geometry': [voronoi_edges]}
    voronoi_edges = gpd.GeoDataFrame(voronoi_edges_df, crs=crs)
    voronoi_edges = voronoi_edges.explode()
    trimmed = voronoi_edges.sjoin(walls, how='left', predicate='intersects')
    trimmed = trimmed.loc[trimmed['index_right0'].isna()]
    trimmed = trimmed[['features', 'geometry']]
    trimmed = trimmed.reset_index()
    trimmed.to_file('Output\\trimmed.gpkg')
    voronoi_edges.to_file('Output\Voronoi.gpkg')



    print("Voronoi Edges created.")


    ########################################
    # Find potential start/end edges       #
    ########################################

    features = [0, 1]
    geometry = [start_pole, end_pole]
    df = {'features': features, 'geometry': geometry}
    poles_gdf = gpd.GeoDataFrame(df, crs=crs)
    if DEBUG:
        poles_gdf.to_file('Output\poles.gpkg')
        print("Exported")

    ########################################
    # Find nearest edges to poles          #
    ########################################
    nearest_edges = trimmed.iloc[gpd.sjoin_nearest(poles_gdf, trimmed)['index_right']].reset_index()
    cols = list(nearest_edges)
    print(cols)
    print("Found nearest edges!")

    ########################################
    # Find shortest path from start to end #
    ########################################

    graph = momepy.gdf_to_nx(trimmed)

    start = nearest_edges['geometry'][0].coords[0]
    end = nearest_edges['geometry'][1].coords[0]
    path = shapely.LineString(nx.shortest_path(graph, source=start, target=end))
    print('Path found')

    ########################################
    # Export centerline                    #
    ########################################

    df = {'features': [0], 'geometry': path}
    centerline_gdf = gpd.GeoDataFrame(df, crs=crs)
    centerline_gdf.to_file('Output\centerline.gpkg')
    print("Exported")
    return centerline_gdf

if __name__ == "__main__":
    main()