#!/usr/bin/env python3

import numpy as np
import sys
import logging
import yaplotlib as yp
import itertools as it

__version__ = "0.1.2"

class Grid():
    def __init__(self, grid=None, file=None, ngrid=None, center=False, pbc=False):
        logger = logging.getLogger()
        if grid is not None:
            self.values = grid.copy()
        elif file is not None:
            self.load(file)
        elif ngrid is not None:
            self.values = np.zeros(ngrid)
        if center:
            x,y,z = self.values.shape
            self.values = np.roll(self.values, x//2, axis=0)
            self.values = np.roll(self.values, y//2, axis=1)
            self.values = np.roll(self.values, z//2, axis=2)
        if pbc:
            self.pbc()
        
    def load(self, file):
        line = file.readline()
        nx,ny,nz = [int(x) for x in line.split()]
        self.values = np.zeros((nx,ny,nz))
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    value = float(file.readline())
                    self.values[x,y,z] = value

    def serialize(self):
        s = "@GRID\n"
        nx,ny,nz = self.values.shape
        s += "{0} {1} {2}\n".format(nx,ny,nz)
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    s+= "{0}\n".format(self.values[x,y,z])
        return s

    def pbcx(self):
        """
        extend the array
        """
        values = self.values
        x,y,z = values.shape
        newvalues = np.zeros((x+1,y,z))
        newvalues[:x,:,:] = values
        newvalues[x,:,:] = values[0,:,:]
        self.values = newvalues

    def pbcy(self):
        """
        extend the array
        """
        values = self.values
        x,y,z = values.shape
        newvalues = np.zeros((x,y+1,z))
        newvalues[:,:y,:] = values
        newvalues[:,y,:] = values[:,0,:]
        self.values = newvalues

    def pbcz(self):
        """
        extend the array
        """
        values = self.values
        x,y,z = values.shape
        newvalues = np.zeros((x,y,z+1))
        newvalues[:,:,:z] = values
        newvalues[:,:,z] = values[:,:,0]
        self.values = newvalues

    def pbc(self):
        """
        extend the array
        """
        self.pbcx()
        self.pbcy()
        self.pbcz()

    def doublex(self):
        """
        double the lattice by linear interpolation
        """
        x,y,z = self.values.shape
        newvalues = np.zeros((x*2-1,y,z))
        newvalues[::2,:,:] = self.values[:,:,:]
        newvalues[1::2,:,:] = (self.values[0:x-1,:,:] + self.values[1:x,:,:])/2
        self.values = newvalues
    
    def doubley(self):
        """
        double the lattice by linear interpolation
        """
        x,y,z = self.values.shape
        newvalues = np.zeros((x,y*2-1,z))
        newvalues[:,::2,:] = self.values[:,:,:]
        newvalues[:,1::2,:] = (self.values[:,0:y-1,:] + self.values[:,1:y,:])/2
        self.values = newvalues
    
    def doublez(self):
        """
        double the lattice by linear interpolation
        """
        x,y,z = self.values.shape
        newvalues = np.zeros((x,y,z*2-1))
        newvalues[:,:,::2] = self.values[:,:,:]
        newvalues[:,:,1::2] = (self.values[:,:,0:z-1] + self.values[:,:,1:z])/2
        self.values = newvalues
    

    def double(self):
        """
        double the lattice by linear interpolation
        """
        self.doublex()
        self.doubley()
        self.doublez()



class Contour(Grid):
    """
    Draw contour surface at value zero.
    """
    #vertex ID: x*4+y*2+z for x,y,z in (0,1)
    vertices = np.array([(0.0,0.0,0.0),(0.0,0.0,1.0),(0.0,1.0,0.0),(0.0,1.0,1.0),
                         (1.0,0.0,0.0),(1.0,0.0,1.0),(1.0,1.0,0.0),(1.0,1.0,1.0)])
    
    edges = [(0,4), (0,2), (0,1),
             (1,5), (1,3), (2,3),
             (3,7), (5,7), (6,7),
             (2,6), (4,6), (4,5)]

    neibor = [[None,      [11, 3, 2], [1, 9,10]],
              [[2, 4, 5], None,       [9,10, 0]],
              [[4, 5, 1], [0,11, 3],  None],
              [None,      [2, 0,11],  [7, 6, 4]],
              [[5, 1, 2], None,       [3, 7, 6]],
              [[1, 2, 4], [6, 8, 9],  None],
              [None,      [8, 9, 5],  [4, 3, 7]],
              [[11,10, 8], None,      [6, 4, 3]],
              [[7,11,10], [9, 5, 6],  None],
              [None,      [5, 6, 8],  [10, 0, 1]],
              [[8, 7,11], None,       [0, 1, 9]],
              [[10, 8, 7],[3, 2, 0],  None]]

    
    def cube_next(self, fi, edge, eorder):
        for fj in range(3):
            if fj != fi and self.neibor[edge][fj] is not None:
                for nextedge in self.neibor[edge][fj]:
                    if nextedge not in self.emark:
                        self.emark.add(nextedge)
                        if 0 <= self.ep[nextedge] <= 1:
                            neworder = self.cube_next(fj,nextedge, eorder+[edge])
                            return neworder
        return eorder+[edge]

    def cube_contour_edge_order(self):
        eorders = []
        for e in range(12):
            if e not in self.emark:
                self.emark.add(e)
                if 0 <= self.ep[e] <= 1:
                    for fi in range(3):
                        edges = self.neibor[e][fi]
                        if edges is not None:
                            for nextedge in edges:
                                if nextedge not in self.emark:
                                    self.emark.add(nextedge)
                                    if 0 <= self.ep[nextedge] <= 1:
                                        eorders.append(self.cube_next(fi,nextedge, [e,]))
        return eorders
                            
                
    def contour_surface_in_a_cube(self, cube):
        """
        generates the contour magically
        """
        self.ep = []
        self.emark = set()
        for k in range(len(self.edges)):
            i,j = self.edges[k]
            if cube[i] == cube[j]:
                if cube[i] == 0.0:
                    p = 0.5
                else:
                    p = -1
            else:
                p = (- cube[i])/(cube[j] - cube[i])
            self.ep.append(p)
            if not (0 <= p <= 1):
                self.emark.add(k)
                
        edge_orders = self.cube_contour_edge_order()
        #from edge orders to vertex positions
        polys = []
        for edge_order in edge_orders:
            poly = []
            for edge in edge_order:
                i,j = self.edges[edge]
                p = self.ep[edge]
                v = self.vertices[i]*(1-p) + self.vertices[j]*p
                poly.append(v)
            polys.append(poly)
        return polys
                                
        
    def contour_flakes(self, value):
        """
        Iterator
        """
        xw,yw,zw = self.values.shape
        for x in range(xw-1):
            for y in range(yw-1):
                for z in range(zw-1):
                    subgrid = self.values[x:x+2, y:y+2, z:z+2].reshape([8]) - value
                    flakes = self.contour_surface_in_a_cube(subgrid)
                    for flake in flakes:
                        yield flake + np.array([x, y, z])

    
    def tetrahedron(self, points, values):
        parity = values > 0.0
        n = np.sum(parity)
        if n in (0,4):
            return
        elif n == 1:
            apex = np.where(values > 0.0)[0][0] #single point
            base = np.where(values <= 0.0)[0]   #three points
            return [(values[apex]*points[v]-values[v]*points[apex])/(values[apex]-values[v]) for v in base]
        elif n == 3:
            apex = np.where(values <= 0.0)[0][0] #single point
            base = np.where(values > 0.0)[0]   #three points
            return [(values[apex]*points[v]-values[v]*points[apex])/(values[apex]-values[v]) for v in base]
        else: # n==2
            top = np.where(values > 0.0)[0]  #two points
            bot = np.where(values <= 0.0)[0] #two points
            combi = ((top[0],bot[0]),(top[1],bot[0]),(top[1],bot[1]),(top[0],bot[1]))
            return [(values[i]*points[j]-values[j]*points[i])/(values[i]-values[j]) for i,j in combi]
    def facets(self, value):
        """
        another algorithm; iterator
        """
        values = self.values
        sx,sy,sz = values.shape
        # values: 3-dim values on the grid
        for ix,iy,iz in it.product(range(sx-1),range(sy-1),range(sz-1)):
            # confirmed the order
            # The whole grid ranges in -0.5 .. +0.5
            # Eight corner positions of a cube
            qp = np.array([(dx,dy,dz) for dx in range(2) for dy in range(2) for dz in range(2)])
            # Eight corner values of a cube
            qv = values[ix:ix+2, iy:iy+2, iz:iz+2].reshape([8]) - value
            # divide a cube into face-sharing five tetrahedra and
            # render the contour surfaces individually.
            if (ix+iy+iz)%2 == 0:
                for i,j,k,l in ((0,1,2,4), (1,4,5,7), (1,2,3,7), (2,4,6,7), (1,2,4,7)):
                    face = self.tetrahedron(qp[[i,j,k,l]], qv[[i,j,k,l]])
                    if face is not None:
                        yield face + np.array([ix,iy,iz])
            else:
                for i,j,k,l in ((0,1,3,5), (3,5,6,7), (0,4,5,6), (0,2,3,6), (0,3,5,6)):
                    face = self.tetrahedron(qp[[i,j,k,l]], qv[[i,j,k,l]])
                    if face is not None:
                        yield face + np.array([ix,iy,iz])


def test3():
    file = open(sys.argv[1])
    while True:
        line = file.readline()
        if "@GRID" in line:
            g = Contour(file=file)
            g.double()
            print("@ 3")
            print(g.contour_yaplot(g.contour_flakes(0.3)), end="")
            g.grid = np.roll(g.grid[::-1, :, :], 1, axis=0)
            print("@ 4")
            print(g.contour_yaplot(g.contour_flakes(0.1)))
            sys.exit(0)

def main():
    contours = [float(x) for x in sys.argv[1:]]
    while True:
        line = sys.stdin.readline()
        if len(line) == 0:
            break
        if "@GRID" in line:
            g = Contour(file=sys.stdin)
            g.double()
            # two algorithms
            s = ""
            for i,c in enumerate(contours):
                s += yp.Color(3+i)
                s += yp.Layer(1+i)
                for flake in g.facets(c):
                    s += yp.Polygon(flake)
            for i,c in enumerate(contours):
                s += yp.Color(3+i)
                s += yp.Layer(1+i)
                for flake in g.contour_flakes(c):
                    s += yp.Polygon(flake)
            print(s)

if __name__ == "__main__":
    main()
