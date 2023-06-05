import argparse
import shapely
import geopandas as gpd

SUBDIVISION_AMOUNT = 200

parser = argparse.ArgumentParser(description='Find the centerline between two lines.')
parser.add_argument("input_wall_1", help="the first valley wall shapefile")
parser.add_argument("input_wall_2", help="the second valley wall shapefile")
parser.add_argument("output", help="the filepath for the output shapefile")


args = parser.parse_args()
input_1 = gpd.read_file(args.input_wall_1)
input_2 = gpd.read_file(args.input_wall_2)
output = args.output


# Import line shapefiles
'''

input_1 = gpd.read_file("C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\WW_ValleyWall_N_10mpoints.shp")
input_2 = gpd.read_file("C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\WW_ValleyWall_S_10mpoints.shp")
'''


wall_1 = shapely.geometry.shape(input_1['geometry'][0])
wall_2 = shapely.geometry.shape(input_2['geometry'][0])

wall_1_points = []
wall_2_points = []
centerline_points = []

for i in range(SUBDIVISION_AMOUNT + 1):
    wall_1_points.append(wall_1.interpolate(i/SUBDIVISION_AMOUNT, normalized=True))
    wall_2_points.append(wall_2.interpolate(i/SUBDIVISION_AMOUNT, normalized=True))
    connecting_line = shapely.geometry.LineString([wall_1_points[i], wall_2_points[i]])
    centerline_points.append(connecting_line.interpolate(0.5, normalized=True))
centerline = shapely.geometry.LineString(centerline_points)
print(centerline)

features = [0]
geometry = [centerline]
df = {'features': features, 'geometry': geometry}
gdf = gpd.GeoDataFrame(df)
gdf = gdf.set_crs("EPSG:32615")
gdf.to_file('C:\LSDTopoTools\Github\Valley-Centerline\WW\Input\Centerline_interpolation.shp')
print("Exported")
