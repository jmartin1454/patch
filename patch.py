#!/usr/bin/python

from scipy.constants import mu_0, pi
import numpy as np

def b_segment(i,p0,p1,r):
    # p0 is one end (vector in m)
    # p1 is the other (m)
    # r is the position of interest (m)
    # i is the current (A)
    d0 = r - p0
    d1 = r - p1
    ell = p1 - p0
    lend0 = np.sqrt(d0.dot(d0))
    lend1 = np.sqrt(d1.dot(d1))
    lenell = np.sqrt(ell.dot(ell))
    costheta0 = np.inner(ell,d0)/lenell/lend0
    costheta1 = -np.inner(ell,d1)/lenell/lend1
    ellcrossd0 = np.cross(ell,d0)
    lenellcrossd0 = np.sqrt(ellcrossd0.dot(ellcrossd0))
    modsintheta0 = lenellcrossd0/lenell/lend0
    a = lend0 * modsintheta0
    if(lenellcrossd0>0):
        nhat = ellcrossd0/lenellcrossd0
    else:
        nhat = np.array([0,0,0])

    if(a>0):
        b_total=mu_0*i/4.0/pi/a*(costheta0+costheta1)*nhat
    else:
        b_total = np.array([0,0,0])

    return b_total

def b_loop(i,points,r):
    # i is the current (A)
    # points is a list of numpy 3-arrays defining the loop (m)
    # r is the position of interest (m)
    b_total = np.array([0,0,0])
    for j in range(len(points)):
        b_total = b_total + b_segment(i,points[j-1],points[j],r)
    return b_total

# Halliday & Resnick, 10th ed., question 29.13
p0 = np.array([0,0,0])
p1 = np.array([0.18,0,0])
r = np.array([0.09,0.131,0])
i = 0.0582
print(b_segment(i,p0,p1,r))

# Halliday & Resnick, 10th ed., question 29.17
p0 = np.array([0,0,0])
p1 = np.array([0.136,0,0])
r = np.array([0.136,0.251,0])
i = 0.693
print(b_segment(i,p0,p1,r))

# Halliday & Resnick, 10th ed., question 29.31
a = 0.047
i = 13.0
p0 = np.array([0,0,0])
p1 = np.array([2*a,0,0])
p2 = np.array([2*a,a,0])
p3 = np.array([a,a,0])
p4 = np.array([a,2*a,0])
p5 = np.array([0,2*a,0])
points = (p0,p1,p2,p3,p4,p5)
r = np.array([2*a,2*a,0])
print(b_loop(i,points,r))

# Halliday & Resnick, 10th ed., question 29.83
a = 0.08
i = 10.0
p0 = np.array([0,0,0])
p1 = np.array([a,0,0])
p2 = np.array([a,-a,0])
p3 = np.array([0,-a,0])
points = (p0,p1,p2,p3)
r = np.array([a/4,-a/4,0])
print(b_loop(i,points,r))

class coilcube:
    def __init__(self,xdim,ydim,zdim,corners):
        self.xdim = xdim
        self.ydim = ydim
        self.zdim = zdim
        self.corners = corners
        x = corners[1]-corners[0]
        y = corners[2]-corners[0]
        z = corners[3]-corners[0]
        self.face = []
        thesecorners=(corners[0],corners[1],corners[2])
        self.face.append(face(xdim,ydim,thesecorners))
        thesecorners=thesecorners+z
        self.face.append(face(xdim,ydim,thesecorners))
        thesecorners=(corners[0],corners[1],corners[3])
        self.face.append(face(xdim,zdim,thesecorners))
        thesecorners=thesecorners+y
        self.face.append(face(xdim,zdim,thesecorners))
        thesecorners=(corners[0],corners[2],corners[3])
        self.face.append(face(ydim,zdim,thesecorners))
        thesecorners=thesecorners+x
        self.face.append(face(ydim,zdim,thesecorners))
        self.numcoils = (xdim*ydim + xdim*zdim + ydim*zdim)*2
    def coil(self,number):
        xdim=self.xdim
        ydim=self.ydim
        zdim=self.zdim
        if(number<xdim*ydim):
            return self.face[0].coil[number]
        elif(number<xdim*ydim*2):
            return self.face[1].coil[number-xdim*ydim]
        elif(number<xdim*ydim*2+xdim*zdim):
            return self.face[2].coil[number-xdim*ydim*2]
        elif(number<xdim*ydim*2+xdim*zdim*2):
            return self.face[3].coil[number-xdim*ydim*2-xdim*zdim]
        elif(number<xdim*ydim*2+xdim*zdim*2+ydim*zdim):
            return self.face[4].coil[number-xdim*ydim*2-xdim*zdim*2]
        else:
            return self.face[5].coil[number-xdim*ydim*2-xdim*zdim*2-ydim*zdim]
    def draw_coil(self,number,ax):
        coil = self.coil(number)
        points = coil.corners + (coil.corners[0],)
        x = ([p[0] for p in points])
        y = ([p[1] for p in points])
        z = ([p[2] for p in points])
        ax.plot(x,y,z,label='coil')
    def draw_coils(self,ax):
        for number in range(self.numcoils):
            self.draw_coil(number,ax)
    def b(self,r):
        b_total = 0.0
        for number in range(self.numcoils):
            b_total = b_total + self.coil(number).b(r)
        return b_total

class face:
    def __init__(self,xdim,ydim,corners):
        self.xdim = xdim
        self.ydim = ydim
        self.corners = corners
        x = corners[1]-corners[0]
        xstep = x/xdim
        y = corners[2]-corners[0]
        ystep = y/ydim
        coilnum = 0
        self.coil = []
        for i in range(xdim):
            for j in range(ydim):
                p0 = corners[0]+xstep*i+ystep*j
                p1 = p0+xstep
                p2 = p1+ystep
                p3 = p2-xstep
                points = (p0,p1,p2,p3)
                self.coil.append(coil(points,0.0))
                coilnum = coilnum + 1
        self.coilnum = coilnum

class coil:
    def __init__(self,corners,current):
        self.corners = corners
        self.current = current
    def set_current(self,current):
        self.current = current
    def b(self,r):
        return b_loop(self.current,self.corners,r)

# repeat of 29.83:
thiscoil = coil(points,i)
print(thiscoil.b(r))
thiscoil.set_current(i/2.0)
print(thiscoil.b(r))

# test of face class -- made a slight change after this test
points = (p0,p1,p3)
thisface = face(2,2,points)
print(thisface.coilnum)
print(thisface.coil[0].corners)
print(thisface.coil[1].corners)
print(thisface.coil[2].corners)
print(thisface.coil[3].corners)

# test of coilcube class
print("Coilcube test")
a = 1.0
p0 = np.array([-a/2,-a/2,-a/2])
p1 = p0 + np.array([a,0,0])
p2 = p0 + np.array([0,a,0])
p3 = p0 + np.array([0,0,a])
points = (p0,p1,p2,p3)
print('hello')
print(points)
mycube = coilcube(2,2,2,points)
#print(mycube.face[1].coil[2].corners)
#print(mycube.coil(6).corners)
#print(mycube.coil(6).current)
print(mycube.numcoils)

class sensor:
    def __init__(self,pos):
        self.pos = pos

class sensorarray:
    def __init__(self,xdim,ydim,zdim,corners):
        x = corners[1]-corners[0]
        y = corners[2]-corners[0]
        z = corners[3]-corners[0]
        self.sensors = []
        for i in range(xdim):
            for j in range(ydim):
                for k in range(zdim):
                    pos = corners[0]+x*i/(xdim-1)+y*j/(ydim-1)+z*k/(zdim-1)
                    self.sensors.append(sensor(pos))
        self.numsensors = len(self.sensors)
    def draw_sensor(self,number,ax):
        x = self.sensors[number].pos[0]
        y = self.sensors[number].pos[1]
        z = self.sensors[number].pos[2]
        c = 'r'
        m = 'o'
        ax.scatter(x,y,z,c=c,marker=m)
    def draw_sensors(self,ax):
        for number in range(self.numsensors):
            self.draw_sensor(number,ax)

# test of sensorarray class
a = 0.8
p0 = np.array([-a/2,-a/2,-a/2])
p1 = p0 + np.array([a,0,0])
p2 = p0 + np.array([0,a,0])
p3 = p0 + np.array([0,0,a])
points = (p0,p1,p2,p3)
myarray = sensorarray(11,11,11,points)
print(myarray.sensors[0].pos)
print(myarray.numsensors)
print(myarray.sensors[myarray.numsensors-1].pos)
print(myarray.sensors[myarray.numsensors-2].pos)

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
mycube.draw_coils(ax)
myarray.draw_sensors(ax)
ax.legend()
plt.show()



print(mycube.b(myarray.sensors[0].pos))
mycube.coil(0).set_current(1.0)
print(mycube.b(myarray.sensors[0].pos))
mycube.coil(0).set_current(0.0)


m = np.zeros((mycube.numcoils,myarray.numsensors*3))
print(m)

# test each coil by graphing field at each sensor
#for i in range(mycube.numcoils):
#    fig = plt.figure()
#    ax = fig.gca(projection='3d')
#    mycube.draw_coil(i,ax)
#    mycube.coil(i).set_current(1.0)
#    for j in range(myarray.numsensors):
#        r = myarray.sensors[j].pos
#        b=mycube.b(r)
#        bhat=b*5.e4
#        points = []
#        points.append(r)
#        points.append(r+bhat)
#        xs = ([p[0] for p in points])
#        ys = ([p[1] for p in points])
#        zs = ([p[2] for p in points])
#        ax.plot(xs,ys,zs)
#    mycube.coil(i).set_current(0.0)
#    ax.legend()
#    plt.show()

# fill m
for i in range(mycube.numcoils):
    if(i%2 == 0):
        print('Even')
        mycube.coil(i).set_current(1.0)
    else:
        print('Odd')
        mycube.coil(i).set_current(-1.0)
    for j in range(myarray.numsensors):
        r = myarray.sensors[j].pos
        b = mycube.b(r)
        for k in range(3):
            m[i,j*3+k]=b[k]
    mycube.coil(i).set_current(0.0)
plt.imshow(m,interpolation='none')
plt.colorbar()
plt.show()

minv = np.linalg.pinv(m)
plt.imshow(minv,interpolation='none')
plt.colorbar()
plt.show()

print(np.linalg.cond(m))
