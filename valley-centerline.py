import shapely
import geopandas as gpd
from geovoronoi import voronoi_regions_from_coords

BUFFER_DISTANCE = 0.0001

# Import line shapefiles
valleywall_01 = gpd.read_file("C:\LSDTopoTools\Github\Valley-Centerline\Wall2.shp")
valleywall_02 = gpd.read_file("C:\LSDTopoTools\Github\Valley-Centerline\Wall1.shp")

valleywallshape_01 = shapely.geometry.shape(valleywall_01['geometry'][0])
valleywallshape_02 = shapely.geometry.shape(valleywall_02['geometry'][0])

# Convert valley walls into a collection of points for voronoi algorithm

# This bit assumes both lines were drawn in the same direction
# @TODO: make it so it works when lines are drawn in opposite directions
coords = valleywallshape_01.coords[:] + valleywallshape_02.coords[::-1]
start_boundary = shapely.geometry.LineString([valleywallshape_01.coords[0], valleywallshape_02.coords[0]])
end_boundary = shapely.geometry.LineString([valleywallshape_01.coords[-1], valleywallshape_02.coords[-1]])

points = shapely.geometry.MultiPoint(coords)

# Generate polygon and buffer from the input walls
valleypoly = shapely.geometry.Polygon(coords)
buffer = valleypoly.buffer(BUFFER_DISTANCE)

features = [0]
geometry = [valleypoly]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\\valleypoly.shp')

features = [0]
geometry = [buffer]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\\valleybuffer.shp')

# Create Voronoi Polygons
region_polys, region_pts = voronoi_regions_from_coords(points, buffer)

#Export the voronoi polygons

features = [i for i in range(len(region_polys))]
geometry = [geom for geom in region_polys.values()]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\\voronoi.shp')


voronoi_edges = []

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
                if edge.coords[:] == test.coords[:]:
                    # In all the cases so far, it seems that the coordinates have been flipped, so I'm not sure if this will ever be the case.
                    repeat = True
                if edge.coords[:] == test.coords[::-1]:
                    repeat = True
            # Get rid of edges that are duplicates or outside of the valley
            if edge.within(valleypoly) and not repeat:
                voronoi_edges.append(edge)
            # Get the first and last segments, which should intersect the lines connecting the two valley walls
            if edge.intersects(start_boundary) and not(edge.intersects(valleywallshape_01) or edge.intersects(valleywallshape_02)):
                start_edge = edge
            if edge.intersects(end_boundary) and not(edge.intersects(valleywallshape_01) or edge.intersects(valleywallshape_02)):
                end_edge = edge

# Export voronoi edges (for testing purposes)
features = [i for i in range(len(voronoi_edges))]
geometry = [geom for geom in voronoi_edges]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\\voronoi_edges.shp')

tried_edges = []

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
        return True
    # See which neighbors haven't been checked yet
    for other_edge in voronoi_edges:
        if edge.touches(other_edge) and other_edge not in tried_edges:
            # Check if the edge is one of from_edge's viable edges:
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
features = [0]
geometry = [centerline_merged]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\\centerlineedges.shp')
