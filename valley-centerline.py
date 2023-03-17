import argparse
import shapely
import geopandas as gpd
from geovoronoi import voronoi_regions_from_coords
import numpy as np

BUFFER_DISTANCE = 0.0001
SUBDIVISION_AMOUNT = 200
'''
parser = argparse.ArgumentParser(description='Find the centerline between two lines.')
parser.add_argument("input_wall_1", help="the first valley wall shapefile")
parser.add_argument("input_wall_2", help="the second valley wall shapefile")
parser.add_argument("output", help="the filepath for the output shapefile")


args = parser.parse_args()
input_1 = gpd.read_file(args.input_wall_1)
input_2 = gpd.read_file(args.input_wall_2)
output = args.output
'''

# Import line shapefiles

input_1 = gpd.read_file("C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\WW_ValleyWall_N_10mpoints.shp")
input_2 = gpd.read_file("C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\WW_ValleyWall_S_10mpoints.shp")

wall_1 = shapely.geometry.shape(input_1['geometry'][0])
wall_2 = shapely.geometry.shape(input_2['geometry'][0])

def subdivide_wall(wall, subdivision_distance):
    subdivided_coords = []
    dist = subdivision_distance
    while dist < wall.length:
        subdivided_coords.append(wall.interpolate(dist))
        dist += subdivision_distance
    return subdivided_coords

#extra_coords = subdivide_wall(wall_1, SUBDIVISION_AMOUNT) + subdivide_wall(wall_2, SUBDIVISION_AMOUNT)

# Convert valley walls into a collection of points for voronoi algorithm

coords = wall_1.coords[:] + wall_2.coords[::-1]
valleypoly = shapely.geometry.Polygon(coords)
start_boundary = shapely.geometry.LineString([wall_1.coords[0], wall_2.coords[0]])
end_boundary = shapely.geometry.LineString([wall_1.coords[-1], wall_2.coords[-1]])

# Handle the case of both walls being drawn in opposite directions
if not valleypoly.is_valid:
    coords = wall_1.coords[:] + wall_2.coords[:]
    valleypoly = shapely.geometry.Polygon(coords)
    start_boundary = shapely.geometry.LineString([wall_1.coords[0], wall_2.coords[-1]])
    end_boundary = shapely.geometry.LineString([wall_1.coords[-1], wall_2.coords[0]])

start_pole = start_boundary.interpolate(0.5, normalized=True)
end_pole = end_boundary.interpolate(0.5, normalized=True)

#coords = coords + extra_coords
points = [shapely.Point(i[0], i[1]) for i in coords]

# Generate polygon and buffer from the input walls
buffer = valleypoly.buffer(BUFFER_DISTANCE)

features = [0]
geometry = [valleypoly]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\ValleyPoly.shp')

features = [0]
geometry = [buffer]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\ValleyBuffer.shp')
print(coords)
# Create Voronoi Polygons
region_polys, region_pts = voronoi_regions_from_coords(points, buffer)
print(type(region_polys))
problem_polys = []
for key, value in region_polys.items():
    if value.geom_type == 'MultiPolygon':
        problem_polys.append(key)
#Export the voronoi polygons
for key in problem_polys:
    region_polys.pop(key)

features = [i for i in range(len(region_polys))]
geometry = [geom for geom in region_polys.values()]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\Voronoi.shp')


print("Voronoi Polygons created.")


########################################
# Find potential start/end edges       #
########################################
edges = np.array()
for poly in region_polys.values():
    for i in range(len(poly.exterior.coords) - 1):
        edge = shapely.geometry.LineString((poly.exterior.coords[i], poly.exterior.coords[i+1]))
        edges.append(edge)

edges_geo = gpd.GeoSeries(edges)
print(edges_geo[:5])

'''
for i in range(edges_geo.size):
    duplicates = edges_geo.geom_equals(edges_geo[i])
    print(duplicates.value_counts())
'''
features = [0, 1]
geometry = [start_pole, end_pole]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf = gdf.set_crs("EPSG:32615")
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\poles.shp')
print("Exported")


########################################
# Find nearest edges to poles          #
########################################

########################################
# Find shortest path from start to end #
########################################

########################################
# Export centerline                    #
########################################
'''
voronoi_edges = []
edges = []
for poly in region_polys.values():
    if poly.geom_type == 'MultiPolygon':
        # If the buffer distance is too small (<0.0006 in this case), some of the voronoi outputs are MultiPolygons.
        # I'm just discarding them here, which seems to work alright, but I'm wondering if there are cases when this would give a bad result.
        pass
    else:
        for i in range(len(poly.exterior.coords) - 1):
            edge = shapely.geometry.LineString((poly.exterior.coords[i], poly.exterior.coords[i+1]))
            repeat = False

            # Edges will have duplicates if they're the border between polygons. If this is the case, we don't want to add the duplicates.
            for test in voronoi_edges:
                if edge.coords[:] == test.coords[:] or edge.coords[:] == test.coords[::-1]:
                    print("duplicate i dont think this ever does anything")
                    # In all the cases so far, it seems that the coordinates have been flipped, so I'm not sure if this will ever be the case.
                    repeat = True
            # Get rid of edges that are duplicates or outside of the valley
            if edge.within(valleypoly) and not repeat:
                voronoi_edges.append(edge)
            # Get the first and last segments, which should intersect the lines connecting the two valley walls
            if not(edge.intersects(wall_1) or edge.intersects(wall_2)):
                if edge.intersects(start_boundary):
                    start_edge = edge
                    print('Found start edge!')
                    edges.append(edge)
                    if len(edges) > 5:
                        features = [i for i in range(len(edges))]
                        geometry = [geom for geom in edges]
                        df = {'features': features, 'geometry': geometry}
                        gdf = gpd.GeoDataFrame(df)
                        gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\V_Edges.shp')
                        break
                if edge.intersects(end_boundary):
                    end_edge = edge
                    print('Found end edge!')
                    edges.append(edge)
                    if len(edges) > 5:
                        features = [i for i in range(len(edges))]
                        geometry = [geom for geom in edges]
                        df = {'features': features, 'geometry': geometry}
                        gdf = gpd.GeoDataFrame(df)
                        gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\V_Edges.shp')
                        break
print('edges found!')


# Export voronoi edges (for testing purposes)
features = [i for i in range(len(voronoi_edges))]
geometry = [geom for geom in voronoi_edges]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\Voronoi_Edges.shp')


# Recursive solution for getting only the centerline from beginning to end
def get_centerline_path(edge, from_edge):
    viable_edges = []
    global centerlineedges
    global tried_edges
    # Add self to the list of tried edges
    tried_edges.append(edge)
    if edge.touches(end_edge):
        # If we've found the end, add ourselves to the list and let the edge that found us know
        centerlineedges.append(edge)
        centerlineedges.append(end_edge)
        print("We've reached the end!")
        return True
    # See which neighbors haven't been checked yet
    for other_edge in voronoi_edges:
        if edge.touches(other_edge) and other_edge not in tried_edges:
            # Don't add the other edge if it's one connected to the parent edge
            if from_edge is None or not(from_edge.touches(other_edge)):
                viable_edges.append(other_edge)
    for other_edge in viable_edges:
        if get_centerline_path(other_edge, edge):
            # If one of our viable edges is on the path toward the end, add ourselves to the list and let the edge that found us know
            centerlineedges.append(edge)
            return True
    return False

centerlineedges = []
get_centerline_path(start_edge, None)

# Export centerline
centerline = shapely.geometry.MultiLineString(centerlineedges)
centerline_merged = shapely.ops.linemerge(centerline)
print(centerline)
features = [0]
geometry = [centerline_merged]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf = gdf.set_crs("EPSG:32615")
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\Centerline.shp')
print("Exported")
'''
