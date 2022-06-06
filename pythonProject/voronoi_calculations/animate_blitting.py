import geovoronoi
from geovoronoi import voronoi_regions_from_coords
import numpy as np
from shapely.geometry import MultiPolygon, Polygon, Point
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as mPolygon
from shapely.geometry.polygon import LinearRing
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import traja
from traja import TrajaCollection
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
from random import random
from IPython.display import display
import pandas as pd
import matplotlib.collections as clt



def generate_traj(num_traj):
    # generates num_trial amount of random trajectories
    trjs = {ind: traja.generate(n=1000, random=True, seed=ind) for ind in range(num_traj)}
    for trj in trjs:
        derivs = traja.get_derivatives(trjs[trj])
        ang = traja.calc_angle(trjs[trj])
        trjs[trj] = trjs[trj].join(derivs).join(ang.rename('angles'))
        trjs[trj].drop(['displacement_time', 'speed_times', 'acceleration_times'], inplace=True, axis=1)
        trjs[trj].to_csv('test_csv.csv')
    coll = TrajaCollection(trjs)
    coll.plot()
    plt.show()
    return trjs

def get_bounds(trjs, num_traj):
    # find x and y bounds by searching for the min and max x and y vals
    x_max = max([trjs[i]['x'].max() for i in range(num_traj)]) + 200
    x_min = min([trjs[i]['x'].min() for i in range(num_traj)]) - 200
    y_max = max([trjs[i]['y'].max() for i in range(num_traj)]) + 200
    y_min = min([trjs[i]['y'].min() for i in range(num_traj)]) - 200
    return (x_min, x_max, y_min, y_max)


def RingCoding(ob):
    n = len(ob.coords)
    codes = np.ones(n, dtype=Path.code_type) * Path.LINETO
    codes[0] = Path.MOVETO
    return codes

def Pathify(polygon):
    vertices = np.concatenate(
                    [np.asarray(polygon.exterior)]
                    + [np.asarray(r) for r in polygon.interiors])
    codes = np.concatenate(
                [RingCoding(polygon.exterior)]
                + [RingCoding(r) for r in polygon.interiors])
    return Path(vertices, codes)

def CreatePatch(poly, area_override=None):
    MAX_DENSITY = 0.75
    area = poly.area
    if area_override is not None:
        area = area_override
    density = 1 / area
    color = (min(1, density / MAX_DENSITY), max(0, (MAX_DENSITY - density) / MAX_DENSITY), 0, 0.5)
    region_external_coords = list(poly.exterior.coords)

    if len(poly.interiors) > 0:
        path = Pathify(poly)
        patch = PathPatch(path, facecolor=color, edgecolor=color)
    else:
        patch = mPolygon(region_external_coords, True)
    patch.set_color(color)
    return patch


def plot_voronoi(frame, shapely_poly, trjs, num_traj):
    patches = []

    # points
    coords = np.array([[trjs[i].loc[frame]['x'], trjs[i].loc[frame]['y']] for i in range(num_traj)])


    #min_x, min_y = np.inf, np.inf
    #max_x, max_y = -np.inf, -np.inf
    #for poly in shapely_poly:
    #    b=poly.bounds
    #    min_x=min(b[0], min_x)
    #    max_x=max(b[2], max_x)
    #    min_y=min(b[1], min_y)
    #    max_y=max(b[3], max_y)
    region_polys, region_pts = voronoi_regions_from_coords(coords, shapely_poly)
    for i in region_polys:
        if type(region_polys[i]) is MultiPolygon:
            point=region_pts[i][0]
            temp_point=Point(coords[point])
            for poly in region_polys[i]:
                if poly.contains(temp_point):
                    patch=CreatePatch(poly)
                    patches.append(patch)
                    temp_area=poly.area
            for poly in region_polys[i]:
                if not poly.contains(temp_point):
                    patch=CreatePatch(poly, temp_area)
                    patches.append(patch)
        else:
            patch=CreatePatch(region_polys[i])
            patches.append(patch)
    points=list(zip(*coords))
    return patches, coords











if __name__=="__main__":
    num_traj = 10
    trjs = generate_traj(num_traj)
    bounds = get_bounds(trjs, num_traj)
    #colors = [(random(), random(), random()) for _ in range(num_traj)]
    colors = (47/255, 85/255, 138/255)
    # plot_voronoi(bounds, trjs, 700, num_traj)
    # print(colors)

    # DEFINE EXTERIOR POLYGONS HERE
    a = [(bounds[0], bounds[2]), (bounds[0], bounds[3]), (bounds[1], bounds[3]), (bounds[1], bounds[2])]

    # DEFINE INTERNAL HOLES HERE
    b = LinearRing([(0, -300), (0, -200), (100, -200), (100, -300), (0, -300)])
    c = LinearRing([(-400, 200), (-200, 300), (-300, 400)])

    shapely_poly = MultiPolygon([[a, [b, c]]])

    fig, ax = plt.subplots()

    scat = ax.scatter([0 for _ in range(num_traj)], [0 for _ in range(num_traj)], c = colors, marker='.', zorder=3, alpha=1.0)
    scat.set(animated=True)
    x = [-600, 400]
    y = [-400, 400]
    col = clt.PatchCollection([], animated=True)
    #ax.set_xlim(bounds[0] - 200, bounds[1] + 200)
    #ax.set_ylim(bounds[2] - 200, bounds[3] + 200)
    ax.set_xlim(-200, 200)
    ax.set_ylim(-200, 200)
    ax.add_collection(col)
    col.set(edgecolor='b', facecolor=colors, alpha=0.5)

    def init():
        return col, scat,


    def update(frame, shapely_poly, trjs, num_traj):
        patches, coords = plot_voronoi(frame+1, shapely_poly, trjs, num_traj)
        scat.set(offsets=coords)
        col.set_paths(patches)
        return scat, col

    anim = animation.FuncAnimation(fig, update, blit=True, interval=100, repeat=False, frames=100, fargs=(shapely_poly, trjs,
                                    num_traj))
    plt.show()
    FFwriter = animation.FFMpegWriter()
    anim.save('voronoi_animation_blitting.mp4', writer=FFwriter)

    #print(trjs)
    display(trjs[1])


