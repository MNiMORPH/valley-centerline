# valley-centerline
Tool to pick the centerline of a valley (or other feature) based on the valley walls (or other edges)

## Ways to pick valley margins

In a GIS, create polylines on either side of the valley. The valley width that is chosen depends on the problem to be solved. For mapping terraces within the valley, it would be best to map at the break between the bedrock valley walls and the flat valley bottom and terraces, such that the defined valley encompasses all of the features of concern. Dealing with tributaries might be more difficult, but the idea that I would suggest is just to start with the simplest: snap a straight line across each tributary-valley mouth between the closest points before the valley starts to curve in towards the tributary. You might also start with a shorter portion of the river to build out a smaller valley-wall data set and see how it works with the algorithm that you come up with (which will be stored here), and then modify and/or extend it depending on how the results look.
