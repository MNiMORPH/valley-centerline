import shapely
import geopandas as gpd


# Import line shapefiles
valleywall_01 = gpd.read_file("Wall1.shp")
valleywall_02 = gpd.read_file("Wall2.shp")

valleywallshape_01 = shapely.geometry.shape(valleywall_01['geometry'][0])
valleywallshape_02 = shapely.geometry.shape(valleywall_02['geometry'][0])


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
    midpoint = line.interpolate(0.5, normalized=True)
    midpoints.append(midpoint.coords[:][0])

midpoint_line = shapely.geometry.LineString(midpoints)
features = [0]
df = {'features': features, 'geometry': [midpoint_line]}
gdf = gpd.GeoDataFrame(df)
gdf.to_file('zumbro_centerline_01.shp')
print(gdf)