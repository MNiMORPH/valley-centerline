[![DOI](https://zenodo.org/badge/331517823.svg)](https://zenodo.org/badge/latestdoi/331517823)

# valley-centerline
Tool to pick the centerline of a valley (or other feature) based on the valley walls (or other edges) using Dijkstra's algorithm. Required libraries: os, argparse, shapely, pandas, geopandas, networkx, momepy

## Ways to pick valley margins

In your GIS software of choice, create a single vector line layer consisting of 2 lines - one on either side of the valley. The valley width that is chosen depends on the problem to be solved. For mapping terraces within the valley, it would be best to map at the break between the bedrock valley walls and the flat valley bottom and terraces, such that the defined valley encompasses all of the features of concern. Dealing with tributaries might be more difficult, but the idea that I would suggest is just to start with the simplest: snap a straight line across each tributary-valley mouth between the closest points before the valley starts to curve in towards the tributary.

## Using valley-centerline
Run valley-centerline.py in the terminal with the following command:

    python valley-centerline.py input output

Where input is the path to the file containing the walls between which the centerline will be found, and output is the path where the centerline will be saved. The output centerline will have the same coordinate reference system as the input file.