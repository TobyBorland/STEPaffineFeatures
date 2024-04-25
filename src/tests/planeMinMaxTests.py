import numpy as np
import pymesh

# from numpy import np.dot
# from math import np.sqrt

from collections import namedtuple  # PyCharm bug will flag this as type error
Point = namedtuple("Point", "x y z")
# eps machine precision constant
eps = np.finfo(float).eps

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

#2766 = ELLIPSE( '', #3190, 0.000359210200000000, 0.000254000000000000 );
#3190 = AXIS2_PLACEMENT_3D( '', #4359, #4360, #4361 );
#4359 = CARTESIAN_POINT( '', ( -0.0104140000000000, -0.0222699757367726, 0.107069472961654 ) );
#4360 = DIRECTION( '', ( 0.707106781186548, -0.707106781186548, 0.00000000000000 ) );
#4361 = DIRECTION( '', ( 0.707106781186547, 0.707106781186547, 0.00000000000000 ) );

# *)
# ENTITY Ellipse
#   SUBTYPE OF (Conic);
#   first_semi_axis : positive_length_measure;
#   second_semi_axis : positive_length_measure;
# END_ENTITY;
# (*

# # Calculates the point on the ellipsoid's boundary closest to given point.
# def pointEllipsoidMinDisp(p, EllipseCentre, EllipseMajorAxis, EllipseMinorAxis, EllipseMajorDisp, EllipseMinorDisp):
#     #  Algorithm by Dr. Robert Nurnberg
#     #  http://wwwf.imperial.ac.uk/~rn/distance2ellipse.pdf
#
#     # ellipseCentre: ELLIPSE[0] -> AXIS2_PLACEMENT_3D[0] -> CARTESIAN_POINT[:]
#     # ellipseLocalYaxis: ELLIPSE[0] -> AXIS2_PLACEMENT_3D[1] -> DIRECTION[:]
#     # ellipseLocalXaxis: ELLIPSE[0] -> AXIS2_PLACEMENT_3D[2] -> DIRECTION[:]
#     # ellipseMajorRadius: ELLIPSE[1]
#     # ellipseMinorRadius: ELLIPSE[2]

# ref_direction: the direction used to determine the direction of the local X axis.

#     # /// Initializes ellipsoid instance using center point and three orthogonal vectors.
#     # /// <param name="Center">Center point.</param>
#     # /// <param name="v1">First semiaxis.</param>
#     # /// <param name="v2">Second semiaxis.</param>
#     # /// <param name="v3">Third semiaxis.</param>
#     # public Ellipsoid(Point3d Center, Vector3d v1, Vector3d v2, Vector3d v3)
#
#     #  Ellipsoid object defined by center point and three mutually orthogonal vectors.
#
#     # Length of the major semiaxis: A = v1.Norm
#     # Length of the intermediate semiaxis: B = v2.Norm
#     # Length of the minor semiaxis: C  = v3.Norm
#
#     #Coord3d local_coord = new Coord3d(this.Center, this._v1, this._v2);
#
#     # local coordinate system using origin point and two vectors.
#
#     # /// <param name="p">Origin of the coordinate system.</param>
#     # /// <param name="v1">Vector oriented along the X axis.</param>
#     # /// <param name="v2">Vector in the XY plane.</param>
#     # /// <param name="name">Name of the coordinate system.</param>
#     # public Coord3d(Point3d p, Vector3d v1, Vector3d v2, string name = "")
#     {
#
#     # if (v1.IsParallelTo(v2))
#     #    /// Check if two objects are parallel
#     # /// </summary>
#     # public bool IsParallelTo(ILinearObject obj)
#     # {
#     # Vector3d v = obj.Direction;
#     # if ((this._coord != v._coord))
#     # v = v.ConvertTo(this._coord);
#
#     return GeometRi3D.AlmostEqual(this.Normalized.Cross(v.Normalized).Norm, 0.0);
#     {
#         throw new Exception("Vectors are parallel");
#     }
#
#     v1 = v1.ConvertToGlobal().Normalized;
#     v = _coord.Axes.Transpose() * v
#     v._coord = Coord3d.GlobalCS
#
#     Vector3d v3 = v1.Cross(v2).Normalized;
#     v2 = v3.Cross(v1).Normalized;
#
#
#     local_coord = (EllipseCentre, v1, v2)
#
#     # Move origin to ellipse centre, and rotate point coordinate to local coordinates based on ellipse major & minor axis.
#     # Assuming direction_ratios of AXIS2_PLACEMENT_3D are normalised
#
#     #local_unit_x = EllipseMajorAxis/np.linalg.norm(EllipseMajorAxis)
#     local_unit_x = EllipseMajorAxis
#     local_z_vec = np.cross(EllipseMajorAxis, EllipseMinorAxis)
#     local_unit_z = local_z_vec/np.linalg.norm(local_z_vec)
#
#     local_y_vec = np.cross(local_z_vec, local_x_vec)
#     local_unit_y = local_y_vec/np.linalg.norm(local_y_vec)
#
#     # print('local_unit_x = {}'.format(local_unit_x))
#     # print('local_unit_y = {}'.format(local_unit_y))
#     # print('local_unit_z = {}'.format(local_unit_z))
#
#     rot_matrix = np.zeros((4, 4))
#     rot_matrix[3, 3] = 1.0
#     rot_matrix[0:3, 0] = local_unit_x
#     rot_matrix[0:3, 1] = local_unit_y
#     rot_matrix[0:3, 2] = local_unit_z
#     determinant = np.linalg.det(rot_matrix)
#     assert np.isclose(determinant, 1.0)
#
#     p = rot_matrix * (p - EllipseCentre)
#
#     #assert np.isclose(p.x, 0)
#     if (AlmostEqual(p.x, 0) && AlmostEqual(p.y, 0)): #(Math.Abs(a - b) <= _tolerance)
#         # Center point, choose any minor-axis
#         return Point(0, C, 0) #, local_coord)
#
#     theta = np.atan2(A * p.y, B * p.x)
#     phi = np.atan2(p.z, C * np.sqrt((p.x * p.x)/(A * A) + (p.y * p.y)/(B * B)))
#     iter = 0
#     max_iter = 100
#     n0 = p
#
#     while (iter < max_iter):
#         iter += 1
#         ct = np.cos(theta)
#         st = np.sin(theta)
#         cp = np.cos(phi)
#         sp = np.sin(phi)
#
#         F1 = (A * A - B * B) * ct * st * cp - p.x * A * st + p.y * B * ct
#         F2 = (A * A * ct * ct + B * B * st * st - C * C) * sp * cp - p.x * A * sp * ct - p.y * B * sp * st + p.z * C * cp
#
#         a11 = (A * A - B * B) * (ct * ct - st * st) * cp - p.x * A * ct - p.y * B * st
#         a12 = -(A * A - B * B) * ct * st * sp
#         a21 = -2.0 * (A * A - B * B) * ct * st * cp * sp + p.x * A * sp * st - p.y * B * sp * ct
#         a22 = (A * A * ct * ct + B * B * st * st - C * C) * (cp * cp - sp * sp) - p.x * A * cp * ct - p.y * B * cp * st - p.z * C * sp
#
#         det = a11 * a22 - a12 * a21
#         if (det == 0):
#             print("Zero determinant")
#
#         # Calc reverse matrix B[ij] = A[ij]^-1
#         b11 = a22 / det
#         b12 = -a12 / det
#         b21 = -a21 / det
#         b22 = a11 / det
#
#         theta = theta - (b11 * F1 + b12 * F2)
#         phi = phi - (b21 * F1 + b22 * F2)
#
#         n = Point(A * np.cos(phi)*np.cos(theta), B * np.cos(phi) * np.sin(theta), C * np.sin(phi), local_coord);
#
#         if (np.abs(n0 - n) < 5 * eps): # tolerance = 1E-12
#             return n
#
#         n0 = n
#
#     return n0


#     # ellipseCentre: ELLIPSE[0] -> AXIS2_PLACEMENT_3D[0] -> CARTESIAN_POINT[:]
#     # ellipseLocalYaxis: ELLIPSE[0] -> AXIS2_PLACEMENT_3D[1] -> DIRECTION[:]
#     # ellipseLocalXaxis: ELLIPSE[0] -> AXIS2_PLACEMENT_3D[2] -> DIRECTION[:]
#     # ellipseMajorRadius: ELLIPSE[1]
#     # ellipseMinorRadius: ELLIPSE[2]


#  Calculates the point on the ellipse boundary closest to given point.
#def pointEllipseMinDisp(p, EllipseCentre, EllipseMajorAxis, EllipseMinorAxis, EllipseMajorDisp, EllipseMinorDisp):


#2731 = CIRCLE( '', #3090, 0.00127000000000000 );
#3090 = AXIS2_PLACEMENT_3D( '', #4329, #4330, #4331 );
#4329 = CARTESIAN_POINT( '', ( 0.00812804782754277, -0.0222699757367726, 0.106053472961653 ) );
#4330 = DIRECTION( '', ( 1.83697019872103E-16, -1.00000000000000, 0.00000000000000 ) );
#4331 = DIRECTION( '', ( 1.00000000000000, 1.83697019872103E-16, 0.00000000000000 ) ); # ref_direction: the direction used to determine the direction of the local X axis.

#     # circleCentre: CIRCLE[0] -> AXIS2_PLACEMENT_3D[0] -> CARTESIAN_POINT[:]
#     # circleLocalYaxis: CIRCLE[0] -> AXIS2_PLACEMENT_3D[1] -> DIRECTION[:]
#     # circleLocalXaxis: CIRCLE[0] -> AXIS2_PLACEMENT_3D[2] -> DIRECTION[:]
#     # circleRadius: CIRCLE[1]






# Point on sphere's surface closest to target point "p".
# ClosestPoint(Point3d p)
# {
# if (p == this.Center)
# {
# // return any point on surface
# return this.Center.Translate(this.R * new Vector3d(1, 0, 0));
# }
# else
# {
# return this.Center.Translate(this.R * new Vector3d(this.Center, p).Normalized);
# }
# }



def planeMinPoint(p, V, minmax="min"):

    def pointTriangleMinDistance(TRI, P):
        # return distance between a point and triangle
        # SYNTAX
        #   dist = pointTriangleDistance(TRI,P)
        #   [dist,PP0] = pointTriangleDistance(TRI,P)
        #
        # DESCRIPTION
        #   Calculate the distance of a given point P from a triangle TRI.
        #   Point P is a row vector of the form 1x3. The triangle is a matrix
        #   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
        #   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
        #   to the triangle TRI.
        #   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
        #   closest point PP0 to P on the triangle TRI.
        #
        # Author: Gwolyn Fischer
        # Release: 1.0
        # Release date: 09/02/02
    
        # The algorithm is based on
        # "David Eberly, 'Distance Between Point and Triangle in 3D',
        # Geometric Tools, LLC, (1999)"
        # http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
        #
        #        ^t
        #  \     |
        #   \reg2|
        #    \   |
        #     \  |
        #      \ |
        #       \|
        #        *P2
        #        |\
        #        | \
        #  reg3  |  \ reg1
        #        |   \
        #        |reg0\
        #        |     \
        #        |      \ P1
        # -------*-------*------->s
        #        |P0      \
        #  reg4  | reg5    \ reg6
        # rewrite triangle in normal form
        B = TRI[0, :]
        E0 = TRI[1, :] - B
        # E0 = E0/np.sqrt(sum(E0.^2)); %normalize vector
        E1 = TRI[2, :] - B
        # E1 = E1/np.sqrt(sum(E1.^2)); %normalize vector
        D = B - P
        a = np.dot(E0, E0)
        b = np.dot(E0, E1)
        c = np.dot(E1, E1)
        d = np.dot(E0, D)
        e = np.dot(E1, D)
        f = np.dot(D, D)
    
        # print "{0} {1} {2} ".format(B,E1,E0)
        det = a * c - b * b
        s = b * e - c * d
        t = b * d - a * e
    
        # Teribble tree of conditionals to determine in which region of the diagram
        # shown above the projection of the point into the triangle-plane lies.
        if (s + t) <= det:
            if s < 0.0:
                if t < 0.0:
                    # region4
                    if d < 0:
                        t = 0.0
                        if -d >= a:
                            s = 1.0
                            sqrDistance = a + 2.0 * d + f
                        else:
                            s = -d / a
                            sqrDistance = d * s + f
                    else:
                        s = 0.0
                        if e >= 0.0:
                            t = 0.0
                            sqrDistance = f
                        else:
                            if -e >= c:
                                t = 1.0
                                sqrDistance = c + 2.0 * e + f
                            else:
                                t = -e / c
                                sqrDistance = e * t + f
    
                                # of region 4
                else:
                    # region 3
                    s = 0
                    if e >= 0:
                        t = 0
                        sqrDistance = f
                    else:
                        if -e >= c:
                            t = 1
                            sqrDistance = c + 2.0 * e + f
                        else:
                            t = -e / c
                            sqrDistance = e * t + f
                            # of region 3
            else:
                if t < 0:
                    # region 5
                    t = 0
                    if d >= 0:
                        s = 0
                        sqrDistance = f
                    else:
                        if -d >= a:
                            s = 1
                            sqrDistance = a + 2.0 * d + f
                            # GF 20101013 fixed typo d*s ->2*d
                        else:
                            s = -d / a
                            sqrDistance = d * s + f
                else:
                    # region 0
                    invDet = 1.0 / det
                    s = s * invDet
                    t = t * invDet
                    sqrDistance = (
                        s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f
                    )
        else:
            if s < 0.0:
                # region 2
                tmp0 = b + d
                tmp1 = c + e
                if tmp1 > tmp0:  # minimum on edge s+t=1
                    numer = tmp1 - tmp0
                    denom = a - 2.0 * b + c
                    if numer >= denom:
                        s = 1.0
                        t = 0.0
                        sqrDistance = a + 2.0 * d + f
                        # GF 20101014 fixed typo 2*b -> 2*d
                    else:
                        s = numer / denom
                        t = 1 - s
                        sqrDistance = (
                            s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f
                        )
    
                else:  # minimum on edge s=0
                    s = 0.0
                    if tmp1 <= 0.0:
                        t = 1
                        sqrDistance = c + 2.0 * e + f
                    else:
                        if e >= 0.0:
                            t = 0.0
                            sqrDistance = f
                        else:
                            t = -e / c
                            sqrDistance = e * t + f
                            # of region 2
            else:
                if t < 0.0:
                    # region6
                    tmp0 = b + e
                    tmp1 = a + d
                    if tmp1 > tmp0:
                        numer = tmp1 - tmp0
                        denom = a - 2.0 * b + c
                        if numer >= denom:
                            t = 1.0
                            s = 0
                            sqrDistance = c + 2.0 * e + f
                        else:
                            t = numer / denom
                            s = 1 - t
                            sqrDistance = (
                                s * (a * s + b * t + 2.0 * d)
                                + t * (b * s + c * t + 2.0 * e)
                                + f
                            )
    
                    else:
                        t = 0.0
                        if tmp1 <= 0.0:
                            s = 1
                            sqrDistance = a + 2.0 * d + f
                        else:
                            if d >= 0.0:
                                s = 0.0
                                sqrDistance = f
                            else:
                                s = -d / a
                                sqrDistance = d * s + f
                else:
                    # region 1
                    numer = c + e - b - d
                    if numer <= 0:
                        s = 0.0
                        t = 1.0
                        sqrDistance = c + 2.0 * e + f
                    else:
                        denom = a - 2.0 * b + c
                        if numer >= denom:
                            s = 1.0
                            t = 0.0
                            sqrDistance = a + 2.0 * d + f
                        else:
                            s = numer / denom
                            t = 1 - s
                            sqrDistance = (
                                s * (a * s + b * t + 2.0 * d)
                                + t * (b * s + c * t + 2.0 * e)
                                + f
                            )
    
        # account for numerical round-off error
        if sqrDistance < 0:
            sqrDistance = 0
    
        dist = np.sqrt(sqrDistance)
    
        PP0 = B + s * E0 + t * E1
        return dist, PP0

    disp = []
    vertices = [] #np.array([])
    p = arrayTypeCast(p)
    for v in V:
        v = arrayTypeCast(v)
        vertices.append(v)
        if minmax=="max":
            d = np.linalg.norm(v - p)
            disp.append(d)

    if minmax=="max":
        maxDisp = max(disp)
        maxPoint = vertices[disp.index(maxDisp)]
        return maxDisp, maxPoint

    tri = pymesh.triangle()
    tri.points = np.array(vertices)
    tri.split_boundary = False
    tri.verbosity = 0
    tri.run() # Execute triangle
    minPointIntersect = []
    for t in tri.faces:
        tt = np.array([vertices[t[0]], vertices[t[1]], vertices[t[2]]])
        d, pd = pointTriangleMinDistance(tt, p)
        disp.append(d)
        minPointIntersect.append(pd)

    minDisp = min(disp)
    minPoint = minPointIntersect[disp.index(minDisp)]
    return minDisp, minPoint

if __name__ == '__main__':

    vertices = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0] ]);
    tri = pymesh.triangle()
    tri.points = vertices
    # tri.max_area = 0.05;
    tri.split_boundary = False
    tri.verbosity = 0
    tri.run() # Execute triangle.

    P = np.array([0.5, -0.3, 0.5])
    disp = []
    minPointIntersect = []
    for t in tri.faces:
        tt = np.array([vertices[t[0]], vertices[t[1]], vertices[t[2]]])
        d, pd = pointTriangleMinDistance(tt, P)
        disp.append(d)
        minPointIntersect.append(pd)

    minDisp = min(disp)
    minPoint = minPointIntersect[disp.index(minDisp)]
pass
