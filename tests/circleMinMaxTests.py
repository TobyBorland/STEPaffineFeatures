import numpy as np
# from numpy import np.dot
# from math import np.sqrt

from collections import namedtuple  # PyCharm bug will flag this as type error
Point = namedtuple("Point", "x y z")
# eps machine precision constant
eps = np.finfo(float).eps

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
import mpl_toolkits.mplot3d.art3d as art3d

def rotation_matrix(d):
    """
    Calculates a rotation matrix given a vector d. The direction of d
    corresponds to the rotation axis. The length of d corresponds to
    the sin of the angle of rotation.

    Variant of: http://mail.scipy.org/pipermail/numpy-discussion/2009-March/040806.html
    """
    sin_angle = np.linalg.norm(d)

    if sin_angle == 0:
        return np.identity(3)

    d /= sin_angle

    eye = np.eye(3)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  -d[1]],
                  [-d[2],     0,  d[0]],
                  [d[1], -d[0],    0]], dtype=np.float64)

    M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
    return M

def pathpatch_2d_to_3d(pathpatch, z = 0., normal = 'z'):
    """
    Transforms a 2D Patch to a 3D patch using the given normal vector.

    The patch is projected into they XY plane, rotated about the origin
    and finally translated by z.
    """
    if type(normal) is str: #Translate strings to normal vectors
        index = "xyz".index(normal)
        normal = np.roll((1.0,0,0), index)

    normal /= np.linalg.norm(normal) #Make sure the vector is normalised

    path = pathpatch.get_path() #Get the path and the associated transform
    trans = pathpatch.get_patch_transform()

    path = trans.transform_path(path) #Apply the transform

    pathpatch.__class__ = art3d.PathPatch3D #Change the class
    pathpatch._code3d = path.codes #Copy the codes
    pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color

    verts = path.vertices #Get the vertices in 2D

    d = np.cross(normal, (0, 0, 1)) #Obtain the rotation vector
    M = rotation_matrix(d) #Get the rotation matrix

    pathpatch._segment3d = np.array([np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts])

def arrayTypeCast(xyz):
    if type(xyz) == list:
        xyz = np.array([xyz[0], xyz[1], xyz[2]])
    elif type(xyz) == Point:
        xyz = np.array([xyz.x, xyz.y, xyz.z])
    elif type(xyz) == np.ndarray:
        pass
    else:
        raise RuntimeError("input type error")
    return xyz

p = np.array([10.0, 10.0, 20.0])
circleCentre = np.array([5., 0., -5.])
#circleCentre = np.array([0., 0., 0.])
circleFeatureAxis = np.array([1.0, 0.0, 0.5])
circleLocalXaxis = np.array([0.0, 1.0, 0.])
circleRadius = 5.0

minmax = 'min'
#minmax = 'max'

# pmin = pointCircleMinDisp(p, cCentre, cLocalXaxis, cFeatureAxis, cRadius)
# def pointCircleMinDisp(p, circleCentre, circleLocalYaxis, circleLocalXaxis, circleRadius, minmax="min"):
#     # Point on circle (including interior points) closest to point "p"
#     # todo: deal with partial arcs, max disp
#
#     p = arrayTypeCast(p)
#     circleLocalXaxis = arrayTypeCast(circleLocalXaxis)
#     circleLocalYaxis = arrayTypeCast(circleLocalYaxis)
#     circleCentre = arrayTypeCast(circleCentre)

# define plane with circle points in the plane and a normal vector cNorm
# via the cross product: n = (U − CC) × (V - CC)
#cNorm = np.cross(circleFeatureAxis - circleCentre, circleLocalXaxis - circleCentre)
cNorm = np.cross(circleFeatureAxis, circleLocalXaxis)
cNorm = cNorm / np.linalg.norm(cNorm)  # unit normal vector

# distance d of point p to the plane is dot product of plane unit normal cNorm
# with the difference vector from the point in the plane P and S: d = (S − P) · u
# The distance is signed: It is positive when S is on the side of the plane where u faces and
# negative when it is on the other side. (It is zero, it S is in the plane, of course.)
# d = (p - plane[0]) * u

d = np.dot((p - circleCentre), cNorm)
#d = np.dot(p, cNorm)

# You can get the point S′, which is S projected to the plane, by subtracting d · u from S:
# S′ = S − d · u = S − ((S − P) · u) · u

pUV = p - np.dot(d, cNorm)  # point projected to circle along orthogonal
ppCentreDisp = np.linalg.norm(pUV - circleCentre)

x = circleCentre[0] + circleRadius / ppCentreDisp * (pUV[0] - circleCentre[0])
y = circleCentre[1] + circleRadius / ppCentreDisp * (pUV[1] - circleCentre[1])
z = circleCentre[2] + circleRadius / ppCentreDisp * (pUV[2] - circleCentre[2])

if minmax == "min":
    if ppCentreDisp > circleRadius:
        #return np.array([x, y, z])
        nUV = np.array([x, y, z])
    else:
        #return pp
        nUV = pUV

if minmax == "max":  # opposite side of circle
    #return np.array([circleCentre[0] - x, circleCentre[1] - y, circleCentre[2] - z])
    nUV = np.array([circleCentre[0] - x, circleCentre[1] - y, circleCentre[2] - z])

#nXY = nUV + circleCentre

ax = plt.figure().add_subplot(projection='3d')
ax.set_xlim(-20, 20)
ax.set_ylim(-20, 20)
ax.set_zlim(-20, 20)

ax.plot3D((p[0],), (p[1],), (p[2],), label="p", marker="x", markersize=5)
ax.plot3D((pUV[0], p[0]), (pUV[1], p[1]), (pUV[2], p[2]), label="n", marker="o", markersize=5, color='red')
#ax.plot3D((p[0], nXY[0]), (p[1], nXY[1]), (p[2], nXY[2]), 'gray')
#ax.plot3D((nXY[0], nnXY[0]), (nXY[1], nnXY[1]), (nXY[2], nnXY[2]), 'blue')
#ax.plot3D((nXY[0], p[0]), (nXY[1], p[1]), (nXY[2], p[2]), 'gray')
ax.plot3D((nUV[0], pUV[0]), (nUV[1], pUV[1]), (nUV[2], pUV[2]), 'gray')
ax.plot3D((circleCentre[0], 20*cNorm[0] + circleCentre[0]), (circleCentre[1], 20*cNorm[1] + circleCentre[1]), (circleCentre[2], 20*cNorm[2] + circleCentre[2]), 'green')
#ellps = Ellipse((eCentre[0], eCentre[1]), width=2*eMajorRadius, height=2*eMinorRadius, angle=0, ec='black', fc=None, fill=False)
circ = Circle((circleCentre[0], circleCentre[1]), radius=circleRadius, angle=0, ec='black', fc=None, fill=False)
ax.add_patch(circ)
#art3d.pathpatch_2d_to_3d(circ, z=circleCentre[2], zdir=(cNorm[0], cNorm[1], cNorm[2]))
pathpatch_2d_to_3d(circ, z=circleCentre[2], normal=(cNorm[0], cNorm[1], cNorm[2]))
#pathpatch_2d_to_3d(circ, z=circleCentre[2], normal=(circleCentre[0], circleCentre[1], circleCentre[2]+1))
#ellps = Ellipse((circleCentre[0], circleCentre[1]), width=20, height=12, angle=0, ec='black', fc=None, fill=False)

# circ = Circle((circleCentre[0], circleCentre[1]), radius=circleRadius, angle=0, ec='black', ec='red', fc=None, fill=False)
#ax.add_patch(ellps)
# #art3d.pathpatch_2d_to_3d(ellps, z=eCentre[2], zdir=(eNorm[0], eNorm[1], eNorm[2]))
#pathpatch_2d_to_3d(ellps, z=circleCentre[2], normal=(cNorm[0], cNorm[1], cNorm[2]))
plt.show()

pass
i=1


#2731 = CIRCLE( '', #3090, 0.00127000000000000 );
#3090 = AXIS2_PLACEMENT_3D( '', #4329, #4330, #4331 );
#4329 = CARTESIAN_POINT( '', ( 0.00812804782754277, -0.0222699757367726, 0.106053472961653 ) );
#4330 = DIRECTION( '', ( 1.83697019872103E-16, -1.00000000000000, 0.00000000000000 ) );
#4331 = DIRECTION( '', ( 1.00000000000000, 1.83697019872103E-16, 0.00000000000000 ) ); # ref_direction: the direction used to determine the direction of the local X axis.

# cCentre: CIRCLE[0] -> AXIS2_PLACEMENT_3D[0] -> CARTESIAN_POINT[:]
# cLocalYaxis: CIRCLE[0] -> AXIS2_PLACEMENT_3D[1] -> DIRECTION[:]
# cLocalXaxis: CIRCLE[0] -> AXIS2_PLACEMENT_3D[2] -> DIRECTION[:]
# cRadius: CIRCLE[1]

p = np.array([1., 1., 20.])
cCentre = np.array([0.0, 0.0, 0.0])
cFeatureAxis = np.array([1.0, 0.0, 0.2])
cLocalXaxis = np.array([0.0, 1.0, 0.2])
cRadius = 20.

# define a plane via circle points in the plane and a normal vector cNorm
# via the cross product: n = (U − CC) × (V - CC)
cNorm = np.cross(cFeatureAxis - cCentre, cLocalXaxis - cCentre)
# set to unit normal vector
cNorm = cNorm / np.linalg.norm(cNorm)

# distance d of point p to the plane is dot product of plane unit normal cNorm
# with the difference vector from the point in the plane P and S: d = (S − P) · u
# The distance is signed: It is positive when S is on the side of the plane where u faces and
# negative when it is on the other side. (It is zero, it S is in the plane, of course.)
# d = (p - plane[0]) * u

d = np.dot((p - cCentre), cNorm)

# subtracting d · u from S to find point S′ is S projected to the plane,
# S′ = S − d · u = S − ((S − P) · u) · u

#     pp = p - u * d

# pp = p - np.dot((p - cPlane), d)
pp = p - np.dot(d, cNorm)

ppCentreDisp = np.linalg.norm(pp - cCentre)

# If Q is at the centre, all points on the circle are equally close.
if np.allclose(pp, cCentre):
    pass  # ignore

if ppCentreDisp <= cRadius:
    print(pp)
else:
    x = cCentre[0] + cRadius / ppCentreDisp * (pp[0] - cCentre[0])
    y = cCentre[1] + cRadius / ppCentreDisp * (pp[1] - cCentre[1])
    z = cCentre[2] + cRadius / ppCentreDisp * (pp[2] - cCentre[2])

    print(np.array([x, y, z]))
