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
import random
import matplotlib.collections as clt



def generate_traj(num_traj):
    # generates num_trial amount of random trajectories
    trjs = {ind: traja.generate(n=1000, random=True, seed=ind) for ind in range(num_traj)}
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
    # Draws the polygon passed in
    n = len(ob.coords)
    codes = np.ones(n, dtype=Path.code_type) * Path.LINETO
    codes[0] = Path.MOVETO
    return codes

def Pathify(polygon):
    # Convert coordinates to path vertices. Objects produced by Shapely's
    # analytic methods have the proper coordinate order, no need to sort.
    # Takes a Shapely polygon and converts to vertices and path codes
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
    #color = c
    color = (min(1, density / MAX_DENSITY), max(0, (MAX_DENSITY - density) / MAX_DENSITY), 0, 0.5)
    region_external_coords = list(poly.exterior.coords)

    if len(poly.interiors) > 0:
        # if there are holes in the shape
        path = Pathify(poly)
        patch = PathPatch(path, facecolor=color, edgecolor=color)
    else:
        patch = mPolygon(region_external_coords, True)
    patch.set_color(color)
    return patch

#def init(scatter, collection):
   # return scatter, collection,





def update(frame, bounds, trjs, num_traj):
    #print(trjs)
    ax.clear()
    patches, points = plot_voronoi(frame+1, bounds, trjs, num_traj)
    for patch in patches:
        ax.add_patch(patch)
    plt.scatter(points[0], points[1])









def plot_voronoi(frame, bounds, trjs, num_traj):



    patches = []
    print(f"trjs {trjs}")

    # points
    coords = np.array([[trjs[i].loc[frame]['x'], trjs[i].loc[frame]['y']] for i in range(num_traj)])

    #DEFINE EXTERIOR POLYGONS HERE
    a = [(bounds[0], bounds[2]), (bounds[0], bounds[3]), (bounds[1],bounds[3]), (bounds[1], bounds[2])]


    #DEFINE INTERNAL HOLES HERE
    b = LinearRing([(0, -300), (0, -200), (100, -200), (100, -300), (0, -300)])
    c = LinearRing([(-400, 200), (-200, 300), (-300, 400)])

    shapely_poly = MultiPolygon([[a, [b, c]]])

    min_x, min_y = np.inf, np.inf
    max_x, max_y = -np.inf, -np.inf
    for poly in shapely_poly:
        b=poly.bounds
        min_x=min(b[0], min_x)
        max_x=max(b[2], max_x)
        min_y=min(b[1], min_y)
        max_y=max(b[3], max_y)

    #fig, ax = plt.subplots()
    #ax.set_xlim(min_x-200, max_x+200)
    #ax.set_ylim(min_y-200, max_y+200)

    # this creates a dictionary of polygons/multipolygons
    # and a dictionary of lists, indicating which point is in those polygons
    # (if there are identical points, those lists might have 2+ numbers in them)
    # region_polys is a dict that maps Voronoi region IDs to shapely Polygon objects
    # that represent the shape of the respective Voronoi region.
    # region_pts is a dict that maps the same Voronoi region IDs as in region_polys to
    # a list of indices into coords, i.e. these indices represent the points that belong to this
    # Voronoi region. Usually, this is only a single point.
    # shapely_poly is the constraints of the voronoid diagram
    region_polys, region_pts = voronoi_regions_from_coords(coords, shapely_poly)
    for i in region_polys:
        if type(region_polys[i]) is MultiPolygon:
            # this means that the voronoi cell is technically a multipolygon.
            # while you could argue whether this should ever occur, the current implementation
            # does this.
            # so, we should probably check which polygon actually contains the point.
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
    #plt.scatter(points[0], points[1])
    #plt.title(f"timestep_{frame}")
    #plt.savefig(f'timestep_{frame}.png')
    #plt.show()

    return patches, points






if __name__=="__main__":
    num_traj = 10
    #colors = [(random.randint(0,255)/255, random.randint(0,255)/255, random.randint(0,255)/255, 0.5) for _ in range(num_traj)]
    trjs = generate_traj(num_traj)
    bounds = get_bounds(trjs, num_traj)
    # plot_voronoi(bounds, trjs, 700, num_traj)
    # print(colors)

    #for i in range(1, 1000, 100):
        #plot_voronoi(bounds, trjs, i, num_traj)

    #plot_voronoi(bounds, trjs, 700, num_traj, colors)

    fig, ax = plt.subplots()
    ax.set_xlim(bounds[0]-200, bounds[1]+200)
    ax.set_ylim(bounds[2]-200, bounds[3]+200)
    anim = animation.FuncAnimation(fig, update, blit=False, interval=25, repeat=False, frames=1000, fargs=(bounds, trjs,
                                    num_traj))
    plt.show()








