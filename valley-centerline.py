import shapely
import geopandas as gpd
from geovoronoi import voronoi_regions_from_coords

# Import line shapefiles
valleywall_01 = gpd.read_file("C:\LSDTopoTools\Github\Valley-Centerline\Wall2.shp")
valleywall_02 = gpd.read_file("C:\LSDTopoTools\Github\Valley-Centerline\Wall1.shp")

valleywallshape_01 = shapely.geometry.shape(valleywall_01['geometry'][0])
valleywallshape_02 = shapely.geometry.shape(valleywall_02['geometry'][0])

# Convert valley walls into a collection of points for voronoi algorithm
coords = valleywallshape_01.coords[:] + valleywallshape_02.coords[:]
points = shapely.geometry.MultiPoint(coords)
print(points)

# Generate bounding box around valley walls for voronoi algorithm
x_min = None
x_max = None
y_min = None
y_max = None
for coord in coords:
    if x_min == None or coord[0] < x_min:
        x_min = coord[0]
    if x_max == None or coord[0] > x_max:
        x_max = coord[0]
    if y_min == None or coord[1] < y_min:
        y_min = coord[1]
    if y_max == None or coord[1] < y_max:
        y_max = coord[1]
boundingbox = shapely.geometry.box(x_min - 10, y_min - 10, x_max + 10, y_max + 10)

# Create Voronoi Polygons
region_polys, region_pts = voronoi_regions_from_coords(points, boundingbox)

#Export the voronoi polygons
features = [i for i in range(len(region_polys))]
geometry = [geom for geom in region_polys.values()]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\\voronoi.shp')

# Get necessary voronoi edges
voronoi_edges = []
for poly in region_polys.values():
    for i in range(len(poly.exterior.coords) - 1):
        edge = shapely.geometry.LineString((poly.exterior.coords[i], poly.exterior.coords[i+1]))
        # Get rid of edges that touch the bounding box or cross the valley walls
        # FIXME - Some edges inside and outside of the walls still slip through
        if not(edge.crosses(valleywallshape_01) or edge.crosses(valleywallshape_02) or edge.intersects(boundingbox.exterior)):
            voronoi_edges.append(edge)
# Export voronoi edges
features = [i for i in range(len(voronoi_edges))]
geometry = [geom for geom in voronoi_edges]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\\voronoi_edges.shp')
# TODO - Connect edges into complete centerline
# TODO - Export centerline

'''
# Project points along lines
midpoints = []
for point in valleywallshape_01.coords[:]:
    # Get distance along line
    project_point = shapely.geometry.Point(point)
    projection_distance = valleywallshape_02.project(project_point)
    
    # Get coordinates of that distance
    interpolated_point = valleywallshape_02.interpolate(projection_distance)
    print(interpolated_point)

    # Make a line between that point and the projected point
    line = shapely.geometry.LineString([interpolated_point, project_point])
    # Find the halfway point
    if(not(line.crosses(valleywallshape_01) or line.crosses(valleywallshape_02))):
        midpoint = line.interpolate(0.5, normalized=True)
        midpoints.append(midpoint.coords[:][0])

midpoint_line = shapely.geometry.LineString(midpoints)
features = [0]
df = {'features': features, 'geometry': [midpoint_line]}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\Zumbro_centerline_3.shp')
print(gdf)
'''