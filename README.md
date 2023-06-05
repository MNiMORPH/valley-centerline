# valley-centerline
Tool to pick the centerline of a valley (or other feature) based on the valley walls (or other edges) using Dijkstra's algorithm. Required libraries: geopandas, geovoronoi, momepy, networkx, numpy, shapely

## Ways to pick valley margins

In your GIS software of choice, create polylines on either side of the valley. The valley width that is chosen depends on the problem to be solved. For mapping terraces within the valley, it would be best to map at the break between the bedrock valley walls and the flat valley bottom and terraces, such that the defined valley encompasses all of the features of concern. Dealing with tributaries might be more difficult, but the idea that I would suggest is just to start with the simplest: snap a straight line across each tributary-valley mouth between the closest points before the valley starts to curve in towards the tributary.

## Using valley-centerline
Change the file paths for Input 1 and Input 2 in valley-centerline.py to point to your input valley walls. Run valley-centerline.py. The output centerline, along with some additional files, is found in the Output folder with the file name 'Centerline.shp'. This file will have the same CRS as the CRS for Input 1.