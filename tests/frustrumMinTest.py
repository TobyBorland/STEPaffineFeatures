import numpy as np

# Orthogonal frustum.  Let E be the origin, D be the direction vector, U be
# the up vector, and R be the right vector.  Let u > 0 and r > 0 be the
# extents in the U and R directions, respectively.  Let n and f be the
# extents in the D direction with 0 < n < f.  The four corners of the frustum
# in the near plane are E + n*D + s0*u*U + s1*r*R where |s0| = |s1| = 1 (four
# choices).  The four corners of the frustum in the far plane are
# E + f*D + (f/n)*(s0*u*U + s1*r*R) where |s0| = |s1| = 1 (four choices).

# class Frustum3:
#     # Construction and destruction.
#     # The default constructor sets the following values:
#     # origin (E) to (0,0,0),
#     # dVector (D) to (0,0,1),
#     # uVector (U) to (0,1,0),
#     # rVector (R) to (1,0,0),
#     # dMin (n) to 1,
#     # dMax (f) to 2,
#     # uBound (u) to 1,
#     # rBound (r) to 1.
#
#     def __init__(self,
#                  inOrigin=np.zeros(3),
#                  inDVector=np.array([0., 0., 1.]),
#                  inUVector=np.array([0., 1., 0.]),
#                  inRVector=np.array([1., 0., 0.]),
#                  inDMin=1.,
#                  inDMax=2.,
#                  inUBound=1.,
#                  inRBound=1.):
#         self.origin = inOrigin
#         self.dVector = inDVector
#         self.uVector = inUVector
#         self.rVector = inRVector
#         self.dMin = inDMin
#         self.dMax = inDMax
#         self.uBound = inUBound
#         self.rBound = inRBound
#         self.mDRatio = self.dMax / self.dMin
#         self.mMTwoUF = (-2) * self.uBound * self.dMax
#         self.mMTwoRF = (-2) * self.rBound * self.dMax

        # def GetDRatio():
        #     return self.mDRatio
        #
        # def GetMTwoUF():
        #     return self.mMTwoUF
        #
        # def ComputeVertices(vertex):
        #     dScaled = self.dMin * self.dVector
        #     uScaled = self.uBound * self.uVector
        #     rScaled = self.rBound * self.rVector
        #     vertex[0] = dScaled - uScaled - rScaled
        #     vertex[1] = dScaled - uScaled + rScaled
        #     vertex[2] = dScaled + uScaled + rScaled
        #     vertex[3] = dScaled + uScaled - rScaled
        #     for i in range(4):
        #         ip = i + 4
        #         vertex[ip] = self.origin + self.mDRatio * vertex[i]
        #         vertex[i] += self.origin

        # # Comparisons to support sorted containers.
        # def __eq__(self, frustum):
        #     return (
        #         self.origin == frustum.origin
        #         and self.dVector == frustum.dVector
        #         and self.uVector == frustum.uVector
        #         and self.rVector == frustum.rVector
        #         and self.dMin == frustum.dMin
        #         and self.dMax == frustum.dMax
        #         and self.uBound == frustum.uBound
        #         and self.rBound == frustum.rBound
        #     )
        #
        # def __ne__(self, frustum):
        #     return not self.__eq__(frustum)
        #
        # def __lt__(self, frustum):
        #     if self.origin < frustum.origin:   return True
        #     if self.origin > frustum.origin:   return False
        #     if self.dVector < frustum.dVector: return True
        #     if self.dVector > frustum.dVector: return False
        #     if self.uVector < frustum.uVector: return True
        #     if self.uVector > frustum.uVector: return False
        #     if self.rVector < frustum.rVector: return True
        #     if self.rVector > frustum.rVector: return False
        #     if self.dMin < frustum.dMin:       return True
        #     if self.dMin > frustum.dMin:       return False
        #     if self.dMax < frustum.dMax:       return True
        #     if self.dMax > frustum.dMax:       return False
        #     if self.uBound < frustum.uBound:   return True
        #     if self.uBound > frustum.uBound:   return False
        #     return self.rBound < frustum.rBound
        #
        # def __le__(self, frustum):
        #     return not frustum.__lt__(self)
        #
        # def __gt__(self, frustum):
        #     return frustum.__lt__(self)
        #
        # def __ge__(self, frustum):
        #     return not self.__lt__(frustum)

# The orthogonal view frustum has origin E.
# Its coordinate axes are determined by left vector L, up vector U,
# and direction vector D. The vectors in that order form a right-handed orthonormal system.

#  The input point is stored in the member closest[0]. The frustum
#  point closest to it is stored in the member closest[1].

# class DCPQuery:
#     def __init__(self):
#         self.distance = 0
#         self.sqrDistance = 0
#         self.closest = [np.zeros(3), np.zeros(3)]

    # def __call__(self, point: np.ndarray, frustum: np.ndarray) -> Tuple[float, float, List[np.ndarray]]:
    #     # Compute coordinates of point with respect to frustum coordinate system.
    #     result = DCPQuery()


#     def frustrum(Origin=np.zeros(3),
#                  dVector=np.array([0., 0., 1.]),
#                  uVector=np.array([0., 1., 0.]),
#                  rVector=np.array([1., 0., 0.]),
#                  dMin=1.,
#                  dMax=2.,
#                  uBound=1.,
#                  rBound=1.):

f_origin = np.array([5., 5., 5.]) #np.zeros(3)
dVector = np.array([0., 0., 1.])
uVector = np.array([0., 1., 0.])
rVector = np.array([1., 0., 0.])
dMin = 1.
dMax = 2.
uBound = 1.
rBound = 1.
mDRatio = dMax / dMin
mMTwoUF = (-2) * uBound * dMax
mMTwoRF = (-2) * rBound * dMax


point = np.array([10., 0., 0.])

# diff = point - frustum[0]
# test = [
#     np.dot(diff, frustum[1]),
#     np.dot(diff, frustum[2]),
#     np.dot(diff, frustum[3])
# ]

diff = point - f_origin
test = [
    np.dot(diff, rVector),
    np.dot(diff, uVector),
    np.dot(diff, dVector)
]

#  Perform calculations in octant with non-negative R and U coordinates.
rSignChange = False
if test[0] < 0.:
    rSignChange = True
    test[0] = -test[0]

uSignChange = False
if test[1] < 0.:
    uSignChange = True
    test[1] = -test[1]

#  Frustum derived parameters.
rmin = rBound
rmax = mDRatio * rmin
umin = uBound
umax = mDRatio * umin
dmin = dMin
dmax = dMax
rminSqr = rmin * rmin
uminSqr = umin * umin
dminSqr = dmin * dmin
minRDDot = rminSqr + dminSqr
minUDDot = uminSqr + dminSqr
minRUDDot = rminSqr + minUDDot
maxRDDot = mDRatio * minRDDot
maxUDDot = mDRatio * minUDDot
maxRUDDot = mDRatio * minRUDDot

#  Algorithm computes closest point in all cases by determining
#  in which Voronoi region of the vertices, edges, and faces of
#  the frustum that the test point lives.

#  The naming conventions for the frustum components are N for near,
#  F for far, U for up, and L for left. The top face of the frustum is
#  labeled the F-face. It has two edges, the UF-edge that is in the
#  direction of L and the LF-edge that is in the direction of U.
#  It also has a vertex, the LUF-vertex at (f `/n, fµ/n, f).
#  The bottom face of the frustum is labeled the N-face. It has two edges,
#  the UN-edge that is in the direction of L and the LN-ege that is in
#  the direction of U. It also has a vertex, the LUN-vertex at (`, µ, n).
#  The remaining two faces are the L-face whose normal is (n, 0, −`) and
#  the U-face whose normal is (0, n, −µ).
#  Finally there is the LU-edge that is shared by the L-face and the U-face.

closest = np.zeros(3)

if test[2] >= dmax:
    if test[0] <= rmax:
        if test[1] <= umax:
            # F-face
            closest[0] = test[0]
            closest[1] = test[1]
            closest[2] = dmax
        else:
            # UF-edge
            closest[0] = test[0]
            closest[1] = umax
            closest[2] = dmax
    else:
        if test[1] <= umax:
            # LF-edge
            closest[0] = rmax
            closest[1] = test[1]
            closest[2] = dmax
        else:
            # LUF-vertex
            closest[0] = rmax
            closest[1] = umax
            closest[2] = dmax

elif test[2] <= dmin:
    if test[0] <= rmin:
        if test[1] <= umin:
            # N-face
            closest[0] = test[0]
            closest[1] = test[1]
            closest[2] = dmin
        else:
            udDot = umin * test[1] + dmin * test[2]
            if udDot >= maxUDDot:
                # UF-edge
                closest[0] = test[0]
                closest[1] = umax
                closest[2] = dmax
            elif udDot >= minUDDot:
                # U-face
                uDot = dmin * test[1] - umin * test[2]
                t = uDot / minUDDot
                closest[0] = test[0]
                closest[1] = test[1] - t * dmin
                closest[2] = test[2] + t * umin
            else:
                # UN-edge
                closest[0] = test[0]
                closest[1] = umin
                closest[2] = dmin

    else:
        if test[1] <= umin:
            rdDot = rmin * test[0] + dmin * test[2]
            if rdDot >= maxRDDot:
                # LF edge
                closest[0] = rmax
                closest[1] = test[1]
                closest[2] = dmax
            elif rdDot >= minRDDot:
                # L-face
                rDot = dmin * test[0] - rmin * test[2]
                t = rDot / minRDDot
                closest[0] = test[0] - t * dmin
                closest[1] = test[1]
                closest[2] = test[2] + t * rmin
            else:
                # LN-edge
                closest[0] = rmin
                closest[1] = test[1]
                closest[2] = dmin

        else:
            rudDot = rmin * test[0] + umin * test[1] + dmin * test[2]
            rEdgeDot = umin * rudDot - minRUDDot * test[1]
            if rEdgeDot >= 0:
                rdDot = rmin * test[0] + dmin * test[2]
                if rdDot >= maxRDDot:
                    # LF-edge
                    closest[0] = rmax
                    closest[1] = test[1]
                    closest[2] = dmax
                elif rdDot >= minRDDot:
                    # L-face
                    rDot = dmin * test[0] - rmin * test[2]
                    t = rDot / minRDDot
                    closest[0] = test[0] - t * dmin
                    closest[1] = test[1]
                    closest[2] = test[2] + t * rmin
                else:
                    # LN-edge
                    closest[0] = rmin
                    closest[1] = test[1]
                    closest[2] = dmin

            else:
                uEdgeDot = rmin * rudDot - minRUDDot * test[0]
                if uEdgeDot >= 0:
                    udDot = umin * test[1] + dmin * test[2]
                    if udDot >= maxUDDot:
                        # UF-edge
                        closest[0] = test[0]
                        closest[1] = umax
                        closest[2] = dmax
                    elif udDot >= minUDDot:
                        # U-face
                        uDot = dmin * test[1] - umin * test[2]
                        t = uDot / minUDDot
                        closest[0] = test[0]
                        closest[1] = test[1] - t * dmin
                        closest[2] = test[2] + t * umin
                    else:
                        # UN-edge
                        closest[0] = test[0]
                        closest[1] = umin
                        closest[2] = dmin

                else:
                    if rudDot >= maxRUDDot:
                        # LUF-vertex
                        closest[0] = rmax
                        closest[1] = umax
                        closest[2] = dmax
                    elif rudDot >= minRUDDot:
                        # LU-edge
                        t = rudDot / minRUDDot
                        closest[0] = t * rmin
                        closest[1] = t * umin
                        closest[2] = t * dmin
                    else:
                        # LUN-vertex
                        closest[0] = rmin
                        closest[1] = umin
                        closest[2] = dmin

else:
    rDot = dmin * test[0] - rmin * test[2]
    uDot = dmin * test[1] - umin * test[2]
    if rDot <= 0:
        if uDot <= 0:
            # point inside frustum
            closest = test
        else:
            udDot = umin * test[1] + dmin * test[2]
            if udDot >= maxUDDot:
                # UF-edge
                closest[0] = test[0]
                closest[1] = umax
                closest[2] = dmax
            else:
                # U-face
                t = uDot / minUDDot
                closest[0] = test[0]
                closest[1] = test[1] - t * dmin
                closest[2] = test[2] + t * umin

    else:
        if uDot <= 0:
            rdDot = rmin * test[0] + dmin * test[2]
            if rdDot >= maxRDDot:
                # LF-edge
                closest[0] = rmax
                closest[1] = test[1]
                closest[2] = dmax
            else:
                # L-face
                t = rDot / minRDDot
                closest[0] = test[0] - t * dmin
                closest[1] = test[1]
                closest[2] = test[2] + t * rmin

        else:
            rudDot = rmin * test[0] + umin * test[1] + dmin * test[2]
            rEdgeDot = umin * rudDot - minRUDDot * test[1]
            if rEdgeDot >= 0:
                rdDot = rmin * test[0] + dmin * test[2]
                if rdDot >= maxRDDot:
                    # LF-edge
                    closest[0] = rmax
                    closest[1] = test[1]
                    closest[2] = dmax
                else: # assert ( rdDot >= min RDDot )
                    # L-face
                    t = rDot / minRDDot
                    closest[0] = test[0] - t * dmin
                    closest[1] = test[1]
                    closest[2] = test[2] + t * rmin

            else:
                uEdgeDot = rmin * rudDot - minRUDDot * test[0]
                if uEdgeDot >= 0:
                    udDot = umin * test[1] + dmin * test[2]
                    if udDot >= maxUDDot:
                        # UF-edge
                        closest[0] = test[0]
                        closest[1] = umax
                        closest[2] = dmax
                    else: #  assert( udDot >= minUDDot )
                        # U-face
                        t = uDot / minUDDot
                        closest[0] = test[0]
                        closest[1] = test[1] - t * dmin
                        closest[2] = test[2] + t * umin
                else:
                    if rudDot >= maxRUDDot:
                        # LUF-vertex
                        closest[0] = rmax
                        closest[1] = umax
                        closest[2] = dmax
                    else: #  assert( rudDot >= minRUDDot )
                        # LU-edge
                        t = rudDot / minRUDDot
                        closest[0] = t * rmin
                        closest[1] = t * umin
                        closest[2] = t * dmin

diff = test - closest

# convert back to original quadrant
if rSignChange:
    closest[0] = -closest[0]
if uSignChange:
    closest[1] = -closest[1]

r_closest = [np.zeros(3), np.zeros(3)]
# convert back to original coords
r_closest[0] = point
r_closest[1] = (
    f_origin
    + closest[0] * rVector
    + closest[1] * uVector
    + closest[2] * dVector
)
r_sqrDistance = np.dot(diff, diff)
r_distance = np.sqrt(r_sqrDistance)
print (r_distance)
#return result
