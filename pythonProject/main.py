import matplotlib
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import traja
import pandas as pd
import numpy as np
from traja import TrajaCollection
from geovoronoi import voronoi_regions_from_coords
from shapely.geometry import MultiPolygon, Polygon, Point
from matplotlib.patches import Polygon as mPolygon
from shapely.geometry.polygon import LinearRing
from matplotlib.patches import PathPatch
from matplotlib.path import Path

num_traj = 2

# GET DATA
# create dictionary to initialize traja collection
# can also initiate with df if given id column
trjs = {ind: traja.generate(n=1000, random=True, seed=ind) for ind in range(num_traj)}
coll = TrajaCollection(trjs)
# feed dictionary into traja collection to plot all trajectories at once
coll.plot()
plt.show()
print(trjs)



#frames
t = list(range(0, 1000, 1))

#find x and y bounds by searching for the min and max x and y vals
x_max = max([trjs[i]['x'].max() for i in range(num_traj)]) + 200
x_min = min([trjs[i]['x'].min() for i in range(num_traj)]) - 200
y_max = max([trjs[i]['y'].max() for i in range(num_traj)]) + 200
y_min = min([trjs[i]['y'].min() for i in range(num_traj)]) - 200

#initialize plot, axis, and collection of lines
fig, ax = plt.subplots()
lns = [plt.plot([], [], '.')[0] for i in range(num_traj)]

def init():
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    return lns

#update frame by reading the dfs
def update(frame):
    for i in range(num_traj):
        lns[i].set_data(trjs[i].loc[:frame]['x'], trjs[i].loc[:frame]['y'])
    return lns

ani = FuncAnimation(fig, update, frames=t,
                    init_func=init, blit=True, interval = 25, repeat=False)
plt.show()

#ani.save("movie.mp4")
FFwriter = animation.FFMpegWriter()
ani.save('animation.mp4', writer = FFwriter)

