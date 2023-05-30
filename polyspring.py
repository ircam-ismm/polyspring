import numpy as np
import shapely as sh
from shapely.ops import nearest_points
from scipy.spatial import Delaunay
from scipy.interpolate import griddata


def polygon_distance_function(region, vertices):
    multi = sh.MultiPoint(vertices)
    print(type(multi.geoms))
    return [p.distance(region) for p in multi.geoms]

def gauss2D(x, y, mx, my, sigx, sigy, theta):
    a = np.cos(theta)**2/(2*sigx**2) + np.sin(theta)**2/(2*sigy**2)
    b = -np.sin(2*theta)/(4*sigx**2) + np.sin(2*theta)/(4*sigy**2)
    c = np.sin(theta)**2/(2*sigx**2) + np.cos(theta)**2/(2*sigy**2)
    gauss = np.exp(- a*(x-mx)**2 - 2*b*(x-mx)*(y-my) - c*(y-my)**2)
    return gauss


class Corpus():

    def __init__(self, track, cols=(0,1)):
        self.dist_func = lambda points : polygon_distance_function(self.region, points)
        self.track = track
        self.h_dist = lambda x, y : 1
        self.interp = 0
        self.setCols(self, cols)

    def setInterp(self, value):
        self.interp = value

    def boundingRegion(self):
        # Region to bounding box
        xmin = min(self.points, key=Point.getX)
        xmax = max(self.points, key=Point.getX)
        ymin = min(self.points, key=Point.getY)
        ymax = max(self.points, key=Point.getY)
        vertices = ((xmin,ymin), (xmin,ymax), (xmax,ymax), (xmax,ymin))
        self.setRegion(sh.MultiPoint(vertices))

    def setRegion(self, region):
        self.region = region
        self.dist_func = lambda points : polygon_distance_function(self.region, points)
        if len(self.points) !=  0:
            self.l0_uni = np.sqrt(2 / (np.sqrt(3) * len(self.points) / self.region.area))

    def setCols(self, cols):
        self.buffers_md = {}
        all_buffer = []
        # concatenate all buffers and store the length of each buffer separately
        for key,buffer in self.track.items():
            all_buffer += buffer
            self.buffers_md[key] = len(buffer)
        self.points = tuple(Point(pt[cols[0]], pt[cols[1]]) for pt in all_buffer)
        self.l0_uni = np.sqrt(2 / (np.sqrt(3) * len(self.points) / self.region.area))

    def getScalingFactor(self):
        targetArea = 0
        nPair = 0
        average_dist = 0
        for point in self.points:
            for near in point.near:
                nPair += 1
                midX ,midY = point.midTo(near)
                average_dist += point.distTo(near)**2
                targetArea += 1 / self.h_dist(midX, midY)**2
        evol = average_dist / (nPair * self.l0_uni**2)
        #return np.sqrt(average_dist / targetArea)
        return self.l0_uni * np.sqrt(nPair / targetArea)

    def preUniformization(self, init=True):
        c = self.region.centroid.coords[0]
        s = np.sqrt(self.region.area) / 3
        x1, y1, x2, y2 = c[0]-s, c[1]-s, c[0]+s, c[1]+s
        allPoints = list(self.points[:]) # copy to preserve initial sorting of self.points
        nbPoints = len(allPoints)
        allPoints.sort(key=Point.getX)
        for i in range(nbPoints):
            allPoints[i].x = (i / (nbPoints - 1)) * (x2 - x1) + x1
        allPoints.sort(key=Point.getY)
        for i in range(nbPoints):
            allPoints[i].y = (i / (nbPoints - 1)) * (y2 - y1) + y1

    def delaunayTriangulation(self):
        allCoord = [[pt.x, pt.y] for pt in self.points]
        allCoord = np.asarray(allCoord)
        triangulation = Delaunay(allCoord)
        self.updateNearPoints(triangulation)
        return triangulation
    
    def updateNearPoints(self, triangulation):
        for point in self.points:
            point.resetNear()
            point.updateOrigin()
        for tri in triangulation.simplices:
            p1 = self.points[tri[0]]
            p2 = self.points[tri[1]]
            p3 = self.points[tri[2]]
            if p1 not in p2.near:
                p2.near.append(p1)
                p1.near.append(p2)
            if p1 not in p3.near:
                p3.near.append(p1)
                p1.near.append(p3)
            if p2 not in p3.near:
                p3.near.append(p2)
                p2.near.append(p3)

    def init_point_in_region(self):
        for point in self.points:
            if not point.shap.within(self.region):
                point.moveTo(nearest_points(self.region, point.shap)[0].coords[0])

    def distribute(self, exportPeriod=0, uni=False, init=True):
        for point in self.points:
            point.recallOg()
        # pre-uniformization
        self.preUniformization(init=init)
        # simulation parameters
        dt = 0.2
        stop_tol = 0.005
        tri_tol = 0.1
        int_pres = 1.2
        k = 1
        # variable initialization
        hScale = self.l0_uni
        tot_count = 0
        tri_count = 0
        updateTri = True
        exit = False
        # main loop
        while not exit:
            exit = True
            # update triangulation if necessary
            if updateTri:
                tri_count += 1
                self.delaunayTriangulation()
                updateTri = False
            # compute rest length scaling factor
            hScale = self.getScalingFactor()
            # sum repulsive actions for each point
            for point in self.points: 
                for near in point.near:
                    midX ,midY = point.midTo(near)
                    f = k * (int_pres * hScale / self.h_dist(midX, midY) - point.distTo(near))
                    if f > 0:
                        near.repulsiveForce(dt * f, point)
            # second loop after all forces computation
            for point in self.points: 
                # check stop condition if inside region, else move it back inside
                if point.shap.within(self.region):
                    if exit and point.moveDist() / self.l0_uni > stop_tol: 
                        exit = False
                else:
                    point.moveTo(nearest_points(self.region, point.shap)[0].coords[0])
                # update point positions
                point.update()
                # check if triangulation needs to be updated
                if not updateTri and point.distFromOrigin() / self.l0_uni > tri_tol:
                    updateTri = True
            # Increment step counter
            tot_count += 1
            # intermediary export to max
            if exportPeriod != 0  and tot_count%exportPeriod == 0:
                self.export()
        # reset triangulation and return various counts
        for point in self.points:
            point.resetNear()
            point.storeUni()
        return tot_count, tri_count

    def simple_attractors(self, gaussians_param, reset=False):
        for point in self.points:
            point.recallUni()
        if reset:
            self.export()
            return
        N = 2 * int(np.ceil(np.sqrt(len(self.points))))
        coord = np.linspace(0, 1, N)
        x, y = np.meshgrid(coord, coord)
        density = np.zeros_like(x)
        for param in gaussians_param:
            gaussian = gauss2D(x, y, *param)
            gaussian /= gaussian.max()
            density += gaussian
        density = self.l0_uni * (density - density.min()) / (density.max() - density.min())
        grad = np.gradient(density)
        grid = np.array( (x.flatten(), y.flatten()) ).T
        coords_x = np.array([pt.x for pt in self.points])
        coords_y = np.array([pt.y for pt in self.points])
        grad_x = griddata(grid, grad[0].flatten(), (coords_x, coords_y))
        grad_y = griddata(grid, grad[1].flatten(), (coords_x, coords_y))
        interp_density = griddata(grid, density.flatten(), (coords_x, coords_y))
        grad_norm = np.sqrt(grad_x**2 + grad_y**2)
        disp_x = interp_density * grad_y / (grad_norm + grad_norm.max()/1000)
        disp_y = interp_density * grad_x / (grad_norm + grad_norm.max()/1000)
        for i, point in enumerate(self.points):
            point.moveTo((point.x + disp_x[i], point.y + disp_y[i]))
            point.update()
        self.export()

    def export(self):
        pass

class Point():

    def __init__(self, x, y):
        self.og_x = x # original position
        self.og_y = y
        self.x = x # current position
        self.y = y
        self.shap = sh.Point(x, y)
        self.uni_x = x # position after uniformisation
        self.uni_y = y
        self.prev_x = x # position at previous triangulation
        self.prev_y = y
        self.push_x = 0.0 # amount of pushing for next step
        self.push_y = 0.0
        self.near = []
    
    def midTo(self, point):
        midX = (self.x + point.x)/2
        midY = (self.y + point.y)/2
        return midX, midY

    def getX(self):
        return self.x
    
    def getY(self):
        return self.y
    
    def distTo(self, point):
        return np.sqrt((self.x-point.x)**2 + (self.y-point.y)**2)

    def repulsiveForce(self, f, point):
        angle = np.arctan2(self.y - point.y, self.x - point.x)
        self.push_x += f * np.cos(angle)
        self.push_y += f * np.sin(angle)
        self.shap = sh.Point(
            self.x + self.push_x,
            self.y + self.push_y) # update the shapely point now for outside observation

    def update(self):
        self.x += self.push_x
        self.y += self.push_y
        self.push_x = 0.0
        self.push_y = 0.0
        self.shap = sh.Point(self.x, self.y)
        
    def updateOrigin(self):
        self.prev_x = self.x
        self.prev_y = self.y
        
    def distFromOrigin(self):
        return np.sqrt((self.x-self.prev_x)**2 + (self.y-self.prev_y)**2)
    
    def resetNear(self):
        self.near = []
        
    def moveTo(self, coords):
        nextX = coords[0]
        nextY = coords[1]
        self.push_x = nextX - self.x
        self.push_y = nextY - self.y

    def moveDist(self):
        return np.sqrt(self.push_x ** 2 + self.push_y ** 2)

    def recallOg(self):
        self.x = self.og_x
        self.y = self.og_y
        self.shap = sh.Point(self.x, self.y)

    def storeUni(self):
        self.uni_x = self.x
        self.uni_y = self.y

    def recallUni(self):
        self.x = self.uni_x
        self.y = self.uni_y
        self.shap = sh.Point(self.x, self.y)

    def __str__(self):
        return str(self.x) + ' ' + str(self.y)

    def __repr__(self):
        return str(self.x) + ' ' + str(self.y)


if __name__ == '__main__':
    print('no main process')