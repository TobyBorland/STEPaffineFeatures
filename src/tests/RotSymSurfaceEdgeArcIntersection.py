
import numpy as np

axisPoint = np.array([ 1, 0, 10])
axisDir = np.array([ 0, 0, 1])
refDir = np.array([1 ,0 ,0 ])
maxPoint = np.array([ 11, 0, 20])
edgeArcAxisPoint = np.array([ 5, 0, 20])
edgeArcAxisDir = np.array([1 ,0 ,0 ])

# define plane of maxSurfacePoint and rot-sym axisDir for cylinder, cone, toroid, etc,
# transform edge curve from global coords to local coords defined by maxSurfacePlane st maxSurfacePlane is coplanar with xy-plane
# find x,y of edge arc points where arc z=0, also where constant max/minPoint radius from rotsymAxis along z

# pel['typeName'] == 'CIRCLE':
# edgeArcAxisPoint = pel['axisPoint']

# project maxPoint, minPoint to surface rotSym axisDir
maxPointProjAxisPoint = axisPoint + (np.dot(maxPoint - axisPoint, axisDir)/np.dot(axisDir, axisDir))*axisDir

# approach here to solve in local coords
# rsNorm = np.cross(refDir, axisDir) # values from surface calc
# rsNorm = rsNorm / np.linalg.norm(rsNorm)
# rotM = np.array([refDir, rsNorm, axisDir])
#
# localEdgeArcAxisPoint = np.matmul(rotM, edgeArcAxisPoint - maxPointProjAxisPoint)

# Compare the intersect plane with the circle plane.
# If they are parallel (i.e. absolute value of dot product of normals equals 1) there are zero solutions,
# np.dot(axisDir, edgeArcAxisDir) == 0.0 #or eps

# unless they are the same plane, in which case there are infinite solutions.
#
# Otherwise, project* the plane in the circle plane, yielding a 2D-line.
# Also you may need to project the circle center in the circle plane,
# yielding a 2D-coordinate for the circle center,
# depending on where you center the coordinate system.
#
# Solve this as a 2D line-circle intersection with the 2D center and radius.
# This yields 0, 1 or 2 solutions which I projected back into 3D space (in the plane of the circle).
#
# *Project = Take two arbitrary orthogonal axes in the plane (that are orthogonal to the normal).
#pv1 = maxPointProjAxisPoint - maxPoint

lineDir = np.cross(axisDir, edgeArcAxisDir)

# direction orthogonal to plane of (circle axis & maxPoint plane axis) projecting circle centrePoint to maxPoint plane
circleProjDir = np.cross(lineDir, edgeArcAxisDir)

d = np.dot((maxPointProjAxisPoint - edgeArcAxisPoint), axisDir) / np.dot(circleProjDir, axisDir)

lineP = edgeArcAxisPoint + d * lineDir

# this point should lie on the center of a segment defining the intersection of the edge arc with a circle with maxPoint on its radius

# |maxPointProjAxisPoint - maxPoint| < |maxPointProjAxisPoint - lineP| => no intersection of edge curves and surface at maxPoint radius, something wrong
# |maxPointProjAxisPoint - maxPoint| == |maxPointProjAxisPoint - lineP| => tangent intersection of edge curves and surface at maxPoint radius, something wrong

R=5

# E is the starting point of the ray,
# L is the end point of the ray,
# C is the center of sphere you're testing against
# r is the radius of that sphere
# Compute:
# d = L - E #( Direction vector of ray, from start to end )
# f = E - C #( Vector from center sphere to ray start )


#float discriminant = b*b-4*a*c;

f = maxPointProjAxisPoint - lineP
a = np.dot(lineDir, lineDir)
b = 2 * np.dot(f, lineDir)
c = np.dot(f, f) - R*R

d = b*b-4*a*c

if d < 0:
    pass
    # no intersection
if d == 0: # eps
    # tangent
    pass

d = np.sqrt( d )

# // either solution may be on or off the ray so need to test both
# // t1 is always the smaller value, because BOTH discriminant and
# // a are nonnegative.
t1 = (-b - d)/(2*a)
t2 = (-b + d)/(2*a)

P1 = lineP + t1 * lineDir
P2 = lineP + t2 * lineDir
#This is a parametric equation:
# Px = Ex + tdx
# Py = Ey + tdy

# Q = G* np.sqrt(Q)/2*A
# B = (-B * G)
#
# P1 = D * (B + Q) + O
# P2 = D * (B - Q) + O

# Each 2D dimension is the dot product of the corresponding axis and the 3D coordinate.
# To "reverse project", multiply each dimension with the corresponding axis and add them together.
# Also to move the reverse projected plane into the circle plane
# (restoring the information that was lost when 2D projecting),
# add the dot product of the circle center and the circle normal.


# intersection(&p, &point, &direction)
#   t = dot(p._normal, point - p._point) / dot(p._normal, -direction);
#   return point + t * direction;
# }
_1=1

def rotnMat1(v1, v2):
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    if np.allclose(v1, v2):
        print("!")
        return None
    ax = np.cross(v1, v2)
    ax = ax/np.linalg.norm(ax)
    m1 = np.array([v1, ax, np.cross(ax, v1)])
    m2 = np.array([v2, ax, np.cross(ax, v2)])
    return np.matmul(m2, m1.T)

def rotnMat3(v1, v2):
   vv1 = v1/(np.sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]))
   vv2 = v2/(np.sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]))

   R=np.array([[vv2[0]*vv1[0], vv2[0]*vv1[1], vv2[0]*vv1[2]],
               [vv2[1]*vv1[0], vv2[1]*vv1[1], vv2[1]*vv1[2]],
               [vv2[2]*vv1[0], vv2[2]*vv1[1], vv2[2]*vv1[2]] ])

   return(R)

def rotation_matrix(A,B):
# a and b are in the form of numpy array

   ax = A[0]
   ay = A[1]
   az = A[2]

   bx = B[0]
   by = B[1]
   bz = B[2]

   au = A/(np.sqrt(ax*ax + ay*ay + az*az))
   bu = B/(np.sqrt(bx*bx + by*by + bz*bz))

   R=np.array([[bu[0]*au[0], bu[0]*au[1], bu[0]*au[2]], [bu[1]*au[0], bu[1]*au[1], bu[1]*au[2]], [bu[2]*au[0], bu[2]*au[1], bu[2]*au[2]] ])


   return(R)

def rotnMat2(v1, v2):
    # test for two unique vectors.
    if np.allclose(v1, v2):
        return
    mx = v1/np.linalg.norm(v1)
    mz = np.cross(v1, v2)
    mz = mz/np.linalg.norm(mz)
    my = np.cross(mz, v1)
    mz = my/np.linalg.norm(my)
    return np.array([mx, my, mz])


