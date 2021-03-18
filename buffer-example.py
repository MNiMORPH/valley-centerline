from shapely.geometry import LineString
from matplotlib import pyplot as plt

line = LineString([(0, 0), (1, 1), (1, 2), (2, 2), (3, 3), (4, 2), (5, 2)])
buff = line.buffer(0.5)

xl, yl = line.coords.xy
xb, yb = buff.exterior.coords.xy

plt.ion()
plt.figure()
plt.plot(xb, yb, 'ko-')
plt.plot(xl, yl, 'ro-')

