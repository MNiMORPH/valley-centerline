import argparse
import shapely
import pandas as pd
import geopandas as gpd
from geovoronoi import voronoi_regions_from_coords
import networkx as nx
import momepy
import numpy as np

DEBUG = False
BUFFER_DISTANCE = 0.0001
SUBDIVISION_AMOUNT = 200
'''
parser = argparse.ArgumentParser(description='Find the centerline between two lines.')
parser.add_argument("input_walls.iloc[0].geometry", help="the first valley wall shapefile")
parser.add_argument("input_walls.iloc[1].geometry", help="the second valley wall shapefile")
parser.add_argument("output", help="the filepath for the output shapefile")


args = parser.parse_args()
input_1 = gpd.read_file(args.input_walls.iloc[0].geometry)
input_2 = gpd.read_file(args.input_walls.iloc[1].geometry)
output = args.output
'''

# Import line shapefiles

#Explode function ensures each wall is a LineString rather than a MultiLineString
walls = gpd.read_file("Input\SampleData\Whitewater\WhitewaterWalls.shp").explode()
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
points = [shapely.Point(i[0], i[1]) for i in coords]

# Generate polygon and buffer from the input walls
buffer = valleypoly.buffer(BUFFER_DISTANCE)
if DEBUG:
    df = {'features': [0], 'geometry': [valleypoly]}
    gdf = gpd.GeoDataFrame(df)
    gdf = gdf.set_crs(crs)
    gdf.to_file('Output\Valley_Polygon.shp')


    df = {'features': [0], 'geometry': [buffer]}
    gdf = gpd.GeoDataFrame(df)
    gdf = gdf.set_crs(crs)
    gdf.to_file('Output\ValleyBuffer.shp')

# Create Voronoi Polygons
region_polys, region_pts = voronoi_regions_from_coords(points, buffer)
problem_polys = []
for key, value in region_polys.items():
    if value.geom_type == 'MultiPolygon':
        problem_polys.append(key)

#Export the voronoi polygons
for key in problem_polys:
    region_polys.pop(key)

if DEBUG:
    features = [i for i in range(len(region_polys))]
    geometry = [geom for geom in region_polys.values()]
    df = {'features': features, 'geometry': geometry}
    gdf = gpd.GeoDataFrame(df)
    gdf = gdf.set_crs(crs)

    gdf.to_file('Output\Voronoi.shp')


print("Voronoi Polygons created.")


########################################
# Find potential start/end edges       #
########################################
edges = np.array([])
for poly in region_polys.values():
    for i in range(len(poly.exterior.coords) - 1):
        edge = shapely.geometry.LineString((poly.exterior.coords[i], poly.exterior.coords[i+1]))
        if not(edge.intersects(walls.iloc[0].geometry) or edge.intersects(walls.iloc[1].geometry)):
            edges = np.append(edges, edge)


features = [i for i in range(len(edges))]
geometry = edges
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
edges_gdf = gdf.set_crs(crs)
if DEBUG:
    edges_gdf.to_file('Output\edges.shp')
    print("Exported edges")

features = [0, 1]
geometry = [start_pole, end_pole]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
poles_gdf = gdf.set_crs(crs)
if DEBUG:
    poles_gdf.to_file('Output\poles.shp')
    print("Exported")



########################################
# Find nearest edges to poles          #
########################################
nearest_edges = gpd.sjoin_nearest(poles_gdf, edges_gdf).merge(edges_gdf, left_on="index_right", right_index=True)
cols = list(nearest_edges)
print(cols)
print("Found nearest edges!")


features = [i for i in range(len(nearest_edges['geometry_y']))]
geometry = nearest_edges['geometry_y']
df = {'features': features, 'geometry': geometry}
nearest_edges = gpd.GeoDataFrame(df)
if DEBUG:
    nearest_edges.to_file('Output\start_end.shp')

########################################
# Find shortest path from start to end #
########################################

graph = momepy.gdf_to_nx(edges_gdf)

start = nearest_edges['geometry'][0].coords[0]
end = nearest_edges['geometry'][1].coords[0]
path = nx.shortest_path(graph, source=start, target=end)
print('Path found')
print(path[0:5])
path = shapely.LineString(path)

########################################
# Export centerline                    #
########################################

features = [0]
df = {'features': features, 'geometry': path}
centerline_gdf = gpd.GeoDataFrame(df)
centerline_gdf = centerline_gdf.set_crs(crs)
centerline_gdf.to_file('Output\centerline.shp')
print("Exported")