#!/usr/bin/python

from scipy.constants import mu_0, pi
import numpy as np
import time
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import math
from decimal import Decimal

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
#        b_total=mu_0*i/4.0/pi/a*(costheta0+costheta1)*nhat
        b_total=(mu_0*i/4.0/pi/a*(costheta0+costheta1)*nhat)*1000000000 #convert to nT
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
#print(b_segment(i,p0,p1,r))

# Halliday & Resnick, 10th ed., question 29.17
p0 = np.array([0,0,0])
p1 = np.array([0.136,0,0])
r = np.array([0.136,0.251,0])
i = 0.693
#print(b_segment(i,p0,p1,r))

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
#print(b_loop(i,points,r))

# Halliday & Resnick, 10th ed., question 29.83
a = 0.08
i = 10.0
p0 = np.array([0,0,0])
p1 = np.array([a,0,0])
p2 = np.array([a,-a,0])
p3 = np.array([0,-a,0])
points = (p0,p1,p2,p3)
r = np.array([a/4,-a/4,0])
#print(b_loop(i,points,r))

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
        ax.plot(x,y,z,label='coil'+str(number+1))
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
#print(thiscoil.b(r))
#quit()

thiscoil.set_current(i/2.0)
#print(thiscoil.b(r))
#quit()

# test of face class -- made a slight change after this test
points = (p0,p1,p3)
thisface = face(2,2,points)
#print(thisface.coilnum)
#print(thisface.coil[0].corners)
#print(thisface.coil[1].corners)
#print(thisface.coil[2].corners)
#print(thisface.coil[3].corners)
#quit()

# test of coilcube class
#print("Coilcube test")
#a = 1.17
a = 1.24

p0 = np.array([-a/2,-a/2,-a/2])
p1 = p0 + np.array([a,0,0])
p2 = p0 + np.array([0,a,0])
p3 = p0 + np.array([0,0,a])
points = (p0,p1,p2,p3)
#print('hello')
#print(points)
#quit()
mycube = coilcube(1,1,1,points)
#print(mycube.face[1].coil[2].corners)
#print(mycube.coil(6).corners)
#print(mycube.coil(6).current)
print "No of coils : ", mycube.numcoils
#quit()
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
a = 1.0
p0 = np.array([-a/2,-a/2,-a/2])
p1 = p0 + np.array([a,0,0])
p2 = p0 + np.array([0,a,0])
p3 = p0 + np.array([0,0,a])
points = (p0,p1,p2,p3)
print points
myarray = sensorarray(2,2,2,points)
print 'No. of sensors : ', myarray.numsensors
########################################
#Define new positions of the sensors

#print(myarray.sensors[0].pos)
#myarray.sensors[0].pos=[0.46,0.53,-0.56]
#print(myarray.sensors[0].pos)
#myarray.sensors[myarray.numsensors-1].pos=[ -0.52,-0.56,-0.56]
#myarray.sensors[myarray.numsensors-1].pos=[ -0.62,0,0]
#myarray.sensors[myarray.numsensors-2].pos=[ 0.38,-0.56,0.61]
#myarray.sensors[myarray.numsensors-2].pos=[ 0,-0.62,0]
#myarray.sensors[myarray.numsensors-3].pos=[ -0.51,0.53,0.61,]
#myarray.sensors[myarray.numsensors-4].pos=[ 0.4,-0.42,-0.39]
##myarray.sensors[myarray.numsensors-5].pos=[-0.4,0.42,0.39]
#myarray.sensors[myarray.numsensors-6].pos=[-0.47,0.42,-0.39]
##myarray.sensors[myarray.numsensors-7].pos=[-0.4,-0.42,0.39]
##############################################################
#Draw coil
mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
mycube.draw_coils(ax)
myarray.draw_sensors(ax)
ax.legend(loc='upper left')
plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0,wspace=0.0, hspace=0.0)
plt.show()

############################################################
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
####################################################################
m = np.zeros((mycube.numcoils,myarray.numsensors*3))

# fill m
for i in range(mycube.numcoils):
#    mycube.coil(i).set_current(1.0)
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
    print m
    time.sleep(1)
    mycube.coil(i).set_current(0.0)
print(np.linalg.cond(m))
##########################################
#Plot M and other components of M
mp=1000000
mp4=1000
mp6=100
mp4s='%.1E' % Decimal(1000)
mp6s='%.1E' % Decimal(100)
mps='%.1E' % Decimal(1000000)

column_labels = ['1x', '1y', '1z','2x', '2y', '2z','3x', '3y', '3z','4x', '4y', '4z','6x', '6y', '6z', '8x', '8y', '8z']
row_labels_all = ['X-(1)', 'X+(2)', 'Y-(3)', 'Y+(4)', 'Z-(5)', 'Z+(6)']
row_labels=[]
for d in range (0,len(np.arange(m.shape[0])),1):
	row_labels.append(row_labels_all[d])

M=m.T #(Here, M=s*c=sensors*coils Matrix  )
M1=np.linalg.pinv(M)

U, V, Wt = np.linalg.svd(M, full_matrices=True)
W, Vinv, Ut = np.linalg.svd(M1, full_matrices=True)
Vmat=np.zeros([len(np.arange(m.shape[1])),len(np.arange(m.shape[0]))]) 
for d in range (0,len(np.arange(m.shape[0])),1):
	Vmat[d][d]=V[d]
Vmat_T=Vmat.T
L_V=[]
for i in range (0,len(np.arange(m.shape[0])),1):
	L_V.append(round(math.log10(V[i]),1))
print "The Diagonal Matrix in log form is : ", L_V

fig1, ax1 = plt.subplots(figsize=(9, 4.5))
plt.subplots_adjust(left=0.10, bottom=0.10, right=0.99, top=0.94,wspace=0.0, hspace=0.0)
fig2, ax2 = plt.subplots(figsize=(9, 4.5))
plt.subplots_adjust(left=0.10, bottom=0.10, right=0.99, top=0.94,wspace=0.0, hspace=0.0)
fig3, ax3 = plt.subplots(figsize=(9, 4.5))
plt.subplots_adjust(left=0.03, bottom=0.04, right=0.99, top=0.94,wspace=0.0, hspace=0.0)
fig4, ax4 = plt.subplots(figsize=(9, 4.5))
plt.subplots_adjust(left=0.03, bottom=0.04, right=0.99, top=0.94,wspace=0.0, hspace=0.0)
fig6, ax6 = plt.subplots(figsize=(9, 4.5))
plt.subplots_adjust(left=0.03, bottom=0.04, right=0.99, top=0.94,wspace=0.0, hspace=0.0)



ax1.imshow(m, interpolation='none', cmap=cm.bwr, aspect='auto' )
ax2.imshow(M1, interpolation='nearest', cmap=cm.bwr, aspect='auto' )
ax3.imshow(Vmat, interpolation='nearest', cmap=cm.bwr, aspect='auto' )
ax4.imshow(U, interpolation='nearest', cmap=cm.bwr, aspect='auto' )
ax6.imshow(Wt, interpolation='nearest', cmap=cm.bwr, aspect='auto' )

#ax1.set_xticklabels(column_labels)
#ax2.set_xticklabels(column_labels)

ax1.set_xticks(np.arange(m.shape[1]))
ax1.set_yticks(np.arange(m.shape[0]))
ax1.set_xlabel('Fluxgate positions')
ax1.set_ylabel('Coils')
ax1.set_yticklabels(row_labels)
ax1.set_title('Matrix M* (nT/A) ('+str(len(np.arange(m.shape[0])))+'coils * '+str(len(np.arange(m.shape[1])))+'sensors)')

ax2.set_xticks(np.arange(M1.shape[1]))
ax2.set_yticks(np.arange(M1.shape[0]))
ax2.set_xlabel('Fluxgate positions')
ax2.set_ylabel('Coils')
ax2.set_yticklabels(row_labels)
ax2.set_title('Pseudoinverse of M (*'+str(mps)+') (A/nT) ('+str(len(np.arange(m.shape[0])))+'coils * '+str(len(np.arange(m.shape[1])))+'sensors)')

ax3.set_xticks(np.arange(Vmat.shape[1]))
ax3.set_yticks(np.arange(Vmat.shape[0]))
ax3.set_title('V-Sqrt of eigenvalues of M*M & MM* ('+str(len(np.arange(m.shape[1])))+'* '+str(len(np.arange(m.shape[0])))+')')

ax4.set_xticks(np.arange(U.shape[1]))
ax4.set_yticks(np.arange(U.shape[0]))
ax4.set_title('U-Orthonormal eigenvectors(*'+str(mp4s)+') of MM* ('+str(len(np.arange(m.shape[1])))+'*'+str(len(np.arange(m.shape[1])))+')')

ax6.set_xticks(np.arange(Wt.shape[1]))
ax6.set_yticks(np.arange(Wt.shape[0]))
ax6.set_title('W*-Orthonormal eigenvectors(*'+str(mp6s)+') of M*M ('+str(len(np.arange(m.shape[0])))+'*'+str(len(np.arange(m.shape[0])))+')')

for c in range (0,len(np.arange(m.shape[0]))):
	for s in range (0,len(np.arange(m.shape[1]))):
		ax1.text(s, c, int(m[c][s]), va='center', ha='center', rotation=90)
		ax2.text(s, c, int(M1[c][s]*mp), va='center', ha='center', rotation=90)
		ax3.text(c, s, int(Vmat[s][c]), va='center', ha='center')

for c in range (0,len(np.arange(U.shape[0]))):
	for s in range (0,len(np.arange(U.shape[1]))):
		ax4.text(s, c, int(U[c][s]*mp4), va='center', ha='center')

for c in range (0,len(np.arange(Wt.shape[0]))):
	for s in range (0,len(np.arange(Wt.shape[1]))):
		ax6.text(s, c, int(Wt[c][s]*mp6), va='center', ha='center')

plt.show()
##################################################################################

