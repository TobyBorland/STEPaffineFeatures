#2766 = ELLIPSE( '', #3190, 0.000359210200000000, 0.000254000000000000 );
#3190 = AXIS2_PLACEMENT_3D( '', #4359, #4360, #4361 );
#4359 = CARTESIAN_POINT( '', ( -0.0104140000000000, -0.0222699757367726, 0.107069472961654 ) );
#4360 = DIRECTION( '', ( 0.707106781186548, -0.707106781186548, 0.00000000000000 ) );
#4361 = DIRECTION( '', ( 0.707106781186547, 0.707106781186547, 0.00000000000000 ) );

# ENTITY Axis2_Placement_3d
# 	SUBTYPE OF (Placement);
# 	axis : OPTIONAL Direction;
# 	ref_direction : OPTIONAL Direction;
# DERIVE
# 	p : LIST [3:3] OF Direction := Build_Axes(axis,ref_direction);
# WHERE
# 	WR1 : SELF\Placement.location.dim = 3;
# 	WR2 : (NOT (EXISTS (axis))) OR (axis.dim = 3);
# 	WR3 : (NOT (EXISTS (ref_direction))) OR (ref_direction.dim = 3);
# 	WR4 : (NOT (EXISTS (axis))) OR (NOT (EXISTS (ref_direction))) OR
#           (Cross_Product(axis,ref_direction).magnitude > 0.0);
# END_ENTITY;

# *)
# ENTITY Ellipse
#   SUBTYPE OF (Conic);
#   first_semi_axis : positive_length_measure;
#   second_semi_axis : positive_length_measure;
# END_ENTITY;
# (*

# Point on ellipse (including interior points) closest to target point "p".
import numpy as np
import pymesh
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


# #Coord3d local_coord = new Coord3d(this.Center, this._v1, this._v2);
#     # /// Initializes coordinate system using origin point and two vectors.
#     # /// <param name="p">Origin of the coordinate system.</param>
#     # /// <param name="v1">Vector oriented along the X axis.</param>
#     # /// <param name="v2">Vector in the XY plane.</param>
#     public Coord3d(Point3d p, Vector3d v1, Vector3d v2, string name = "")
#     {
#         v1 = v1.ConvertToGlobal().Normalized;
#
#             v = _coord.Axes.Transpose() * v
#             v._coord = Coord3d.GlobalCS
#
#         Vector3d v3 = v1.Cross(v2).Normalized;
#         v2 = v3.Cross(v1).Normalized;
#
#         _origin = p.ConvertToGlobal();
#         _axes = new Matrix3d(v1, v2, v3);
#     }

# public Coord3d(Point3d p, double[] d1, double[] d2, string name = "")
#         {
#             Vector3d v1 = new Vector3d(d1);
#             Vector3d v2 = new Vector3d(d2);
#             if (v1.IsParallelTo(v2))
#             {
#                 throw new Exception("Vectors are parallel");
#             }
#
#             v1 = v1.Normalized;
#             Vector3d v3 = v1.Cross(v2).Normalized;
#             v2 = v3.Cross(v1).Normalized;
#
#             _origin = p.ConvertToGlobal();
#             _axes = new Matrix3d(v1, v2, v3);

#         }

#
#     # /// <summary>
#     # /// Convert vector to global coordinate system
#     # /// </summary>
#     public Vector3d ConvertToGlobal()
#     {
#         if (_coord == null || object.ReferenceEquals(_coord, Coord3d.GlobalCS))
#         {
#             return this.Copy();
#         }
#         else
#         {
#             Vector3d v = this.Copy();
#             v = _coord.Axes.Transpose() * v;
#             v._coord = Coord3d.GlobalCS;
#             return v;
#         }
#     }

#p = p.ConvertTo(local_coord);


# #def create_rotation_matrix_from_points(local_origin, local_x_dir, local_xy_plane):
#     # """
#     # Create a 4x4 homogeneous rotation matrix
#     # :param local_origin: np_array of x,y,z wrt global cs
#     # :param local_x_dir:  np_array of x,y,z point. vector from local_origin to this point defines X
#     # :param local_xy_plane: point in xy plane
#     # :return: rot_matrix: homogeneous rotation matrix
#     # """
#     local_x_vec = local_x_dir - local_origin
#     local_unit_x = unit_vector(local_x_vec) / np.linalg.norm(xyz_array)
#
#     local_xy_vec = local_xy_plane - local_origin
#     local_z_vec = np.cross(local_x_vec, local_xy_vec)
#     local_unit_z = unit_vector(local_z_vec)
#
#     local_y_vec = np.cross(local_z_vec, local_x_vec)
#     local_unit_y = unit_vector(local_y_vec)
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
#     # print('rot_matrix = \n{}'.format(rot_matrix))
#     # print('determinant = {}\n'.format(determinant))
#
#     return rot_matrix

# p = np.array([8., 8., 20.])
# eCentre = np.array([0.0, 0.0, 0.0])
# eCentre = np.array([5.0, 0.0, 0.0])
# eFeatureAxis = np.array([1.0, 0.0, 0.2])
# eLocalXaxis = np.array([0.0, 1.0, 0.2])
# eMajorRadius = 20.
# eMinorRadius = 10.

def pointEllipseMinDisp(p, eCentre, eFeatureAxis, eLocalXYaxis, eMajorRadius, eMinorRadius, minmax="min"):
    # Robert Nurnberg method, see http://wwwf.imperial.ac.uk/~rn/distance2ellipse.pdf
    # Move origin to ellipse centre, and rotate point coordinate to local coordinates
    # based on ellipse major & minor axis.
    # Assuming direction_ratios of AXIS2_PLACEMENT_3D are normalised

    eFeatureAxis = arrayTypeCast(eFeatureAxis)
    eLocalXYaxis = arrayTypeCast(eLocalXYaxis)
    eCentre = arrayTypeCast(eCentre)

    U = eFeatureAxis # - eCentre
    U = U / np.linalg.norm(U)

    eNorm = np.cross(U, eLocalXYaxis)
    eNorm = eNorm / np.linalg.norm(eNorm)

    V = np.cross(eNorm, U)
    V = V / np.linalg.norm(V)

    rotM = np.array([U, V, eNorm])
    pUV = np.matmul(rotM, p - eCentre)

    # project to ellipse plane
    pUV = np.array([pUV[0], pUV[1], 0])

    theta = np.arctan2(eMajorRadius * pUV[1], eMinorRadius * pUV[0])
    i = 0
    n0 = pUV
    radiusConst = eMajorRadius**2 - eMinorRadius**2
    pUVdisp = np.linalg.norm(pUV)

    while (i < 100): # and not strike:
        i += 1
        f = radiusConst * np.cos(theta) * np.sin(theta) - \
            pUV[0] * eMajorRadius * np.sin(theta) + \
            pUV[1] * eMinorRadius * np.cos(theta)

        f_ = radiusConst * (np.cos(theta)**2 - np.sin(theta)**2) - \
             pUV[0] * eMajorRadius * np.cos(theta) - \
             pUV[1] * eMinorRadius * np.sin(theta)

        theta = theta - f / f_

        n = np.array([eMajorRadius * np.cos(theta), eMinorRadius * np.sin(theta), 0])

        if np.allclose(n0, n): # tolerance = 1E-12
            if minmax == "min":
                if np.linalg.norm(n) > pUVdisp: # test if p < n, inside ellipse
                    return np.matmul(rotM.T, pUV) + eCentre
                else:
                    return np.matmul(rotM.T, n) + eCentre
            else:
                # furthest point on ellipse must be edge point on line that contains ellipse centre and closest point
                return np.matmul(rotM.T, np.array([-n[0], -n[1], 0])) + eCentre

        n0 = n

def pointPlaneMinMaxDisp(p, V, minmax="min"):
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
                        s * (a * s + b * t + 2.0 * d)
                        + t * (b * s + c * t + 2.0 * e)
                        + f
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
                            s * (a * s + b * t + 2 * d)
                            + t * (b * s + c * t + 2 * e)
                            + f
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
    vertices = []  # np.array([])
    p = arrayTypeCast(p)
    for v in V:
        v = arrayTypeCast(v)
        vertices.append(v)
        if minmax == "max":
            d = np.linalg.norm(v - p)
            disp.append(d)

    if minmax == "max":
        maxDisp = max(disp)
        maxPoint = vertices[disp.index(maxDisp)]
        return maxDisp, maxPoint

    tri = pymesh.triangle()
    tri.points = np.array(vertices)
    tri.split_boundary = False
    tri.verbosity = 0
    tri.run()  # Execute triangle
    minPointIntersect = []
    for t in tri.faces:
        tt = np.array([vertices[t[0]], vertices[t[1]], vertices[t[2]]])
        d, pd = pointTriangleMinDistance(tt, p)
        disp.append(d)
        minPointIntersect.append(pd)

    minDisp = min(disp)
    minPoint = minPointIntersect[disp.index(minDisp)]
    return minDisp, minPoint

def pointCircleMinMaxDisp(p, circleCentre, circleLocalXYaxis, circleFeatureAxis, circleRadius, minmax="min"):
    # Point on circle (including interior points) closest to point "p"

    p = arrayTypeCast(p)
    circleFeatureAxis = arrayTypeCast(circleFeatureAxis)
    circleLocalXYaxis = arrayTypeCast(circleLocalXYaxis)
    circleCentre = arrayTypeCast(circleCentre)

    # define plane with circle points in the plane and a normal vector cNorm
    # via the cross product: n = (U − CC) × (V - CC)
    #cNorm = np.cross(circleFeatureAxis - circleCentre, circleLocalXaxis - circleCentre)
    cNorm = np.cross(circleFeatureAxis, circleLocalXYaxis)
    cNorm = cNorm / np.linalg.norm(cNorm)  # unit normal vector

    # distance d of point p to the plane is dot product of plane unit normal cNorm
    # with the difference vector from the point in the plane P and S: d = (S − P) · u
    # distance is signed: positive when S is on the side of the plane where u faces and
    # negative on the other side. zero, S is in the plane
    # d = (p - plane[0]) * u

    d = np.dot((p - circleCentre), cNorm)
    #d = np.dot(p, cNorm)

    # point S′ is S projected to the plane, by subtracting d · u from S:
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

    return nUV

def orthoProj2line(projPoint, linePoint1, linePoint2):
    """
    project a point to an infinite line defined by 2 points, or point and vector direction
    projPoint: point to be projected
    linePoint: point on line
    linePoint2: second point on line (vector would not be offset from local origin)
    """
    projPoint = arrayTypeCast(projPoint)
    linePoint1 = arrayTypeCast(linePoint1)
    linePoint2 = arrayTypeCast(linePoint2)

    return (np.dot(linePoint2 - linePoint1, projPoint) / np.linalg.norm(projPoint)**2 ) * projPoint

def pointSegmentMinMaxDisp(p, endPoint1, endPoint2, minmax="min"):
    p = arrayTypeCast(p)
    endPoint1 = arrayTypeCast(endPoint1)
    endPoint1 = arrayTypeCast(endPoint1)
    pp = orthoProj2line(p, endPoint1, endPoint2)

    if (endPoint1[0] > endPoint2[0]) and (endPoint1[0] > pp[0]) and (pp[0] > endPoint2[0]):
        minP = pp
    if (endPoint1[0] < endPoint2[0]) and (endPoint1[0] < pp[0]) and (pp[0] < endPoint2[0]):
        minP = pp
    if (endPoint1[0] < endPoint2[0]) and (pp[0] < endPoint1[0]):
        minP = endPoint1


def pointSegmentMinDisp(p, a, b):
    # https://web.archive.org/web/20220121145748/http://geomalgorithms.com/index.html
    s = b - a
    w = p - a
    ps = np.dot(w, s)
    if ps <= 0:
        return a, np.linalg.norm(w)
    l2 = np.dot(s, s)
    if ps >= l2:
        closest = b
    else:
        closest = a + ps / l2 * s
    return closest, np.linalg.norm(p - closest)

# // dist_Point_to_Segment(): get the distance of a point to a segment
# //     Input:  a Point P and a Segment S (in any dimension)
# //     Return: the shortest distance from P to S
# float
# dist_Point_to_Segment( Point P, Segment S)
# {
#      Vector v = S.P1 - S.P0;
#      Vector w = P - S.P0;
#
#      double c1 = dot(w,v);
#      if ( c1 <= 0 )
#           return d(P, S.P0);
#
#      double c2 = dot(v,v);
#      if ( c2 <= c1 )
#           return d(P, S.P1);
#
#      double b = c1 / c2;
#      Point Pb = S.P0 + b * v;
#      return d(P, Pb);
# }

# def lineseg_dist(p, a, b):
#
#     # normalized tangent vector
#     d = np.divide(b - a, np.linalg.norm(b - a))
#
#     # signed parallel distance components
#     s = np.dot(a - p, d)
#     t = np.dot(p - b, d)
#
#     # clamped parallel distance
#     h = np.maximum.reduce([s, t, 0])
#
#     # perpendicular distance component
#     c = np.cross(p - a, d)
#
#     return np.hypot(h, np.linalg.norm(c))

# def rotation_matrix(d):
#     """
#     Calculates a rotation matrix given a vector d. The direction of d
#     corresponds to the rotation axis. The length of d corresponds to
#     the sin of the angle of rotation.
#
#     Variant of: http://mail.scipy.org/pipermail/numpy-discussion/2009-March/040806.html
#     """
#     sin_angle = np.linalg.norm(d)
#
#     if sin_angle == 0:
#         return np.identity(3)
#
#     d /= sin_angle
#
#     eye = np.eye(3)
#     ddt = np.outer(d, d)
#     skew = np.array([[    0,  d[2],  -d[1]],
#                   [-d[2],     0,  d[0]],
#                   [d[1], -d[0],    0]], dtype=np.float64)
#
#     M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
#     return M
#
# def pathpatch_2d_to_3d(pathpatch, z = 0, normal = 'z'):
#     """
#     Transforms a 2D Patch to a 3D patch using the given normal vector.
#
#     The patch is projected into they XY plane, rotated about the origin
#     and finally translated by z.
#     """
#     if type(normal) is str: #Translate strings to normal vectors
#         index = "xyz".index(normal)
#         normal = np.roll((1.0,0,0), index)
#
#     normal /= np.linalg.norm(normal)
#
#     path = pathpatch.get_path()
#     trans = pathpatch.get_patch_transform()
#
#     path = trans.transform_path(path)
#
#     pathpatch.__class__ = art3d.PathPatch3D
#     pathpatch._code3d = path.codes #Copy the codes
#     pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color
#
#     verts = path.vertices #Get the vertices in 2D
#
#     # U = np.cross(normal, np.array([0., 1., 0.]))
#     # U = U / np.linalg.norm(U)
#     #
#     # V = np.cross(normal, U)
#     # V = V / np.linalg.norm(V)
#     #
#     # rotM = np.array([U, V, normal])
#     # A = np.array([np.matmul(rotM, np.array([x, y, 0])) + np.array([0, 0, z]) for x, y in verts])
#
#     d = np.cross(normal, (0, 0, 1))
#     M = rotation_matrix(d)
#
#     pathpatch._segment3d = np.array([np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts])
#     #pathpatch._segment3d = A
#
#     #return(A)
#
# import matplotlib.pyplot as plt
# from matplotlib.patches import Ellipse
# import mpl_toolkits.mplot3d.art3d as art3d
#
# ax = plt.figure().add_subplot(projection='3d')
# ax.set_xlim(-20, 20)
# ax.set_ylim(-20, 20)
# ax.set_zlim(-20, 20)
#
# ax.plot3D((p[0],), (p[1],), (p[2],), label="p", marker="x", markersize=5)
# ax.plot3D((pUV[0], n[0]), (pUV[1], n[1]), (pUV[2], n[2]), label="n", marker="o", markersize=5, color='red')
# #ax.plot3D((p[0], pUV[0]), (p[1], pUV[1]), (p[2], pUV[2]), 'gray')
# ax.plot3D((nXY[0], nnXY[0]), (nXY[1], nnXY[1]), (nXY[2], nnXY[2]), 'blue')
# ax.plot3D((nXY[0], p[0]), (nXY[1], p[1]), (nXY[2], p[2]), 'gray')
# ax.plot3D((eCentre[0], 20*eNorm[0]), (eCentre[1], 20*eNorm[1]), (eCentre[2], 20*eNorm[2]), 'green')
#
# ellps = Ellipse((eCentre[0], eCentre[1]), width=2*eMajorRadius, height=2*eMinorRadius, angle=0, ec='black', fc=None, fill=False)
# ax.add_patch(ellps)
# #art3d.pathpatch_2d_to_3d(ellps, z=eCentre[2], zdir=(eNorm[0], eNorm[1], eNorm[2]))
# pathpatch_2d_to_3d(ellps, z=eCentre[2], normal=(eNorm[0], eNorm[1], eNorm[2]))
# #ellps2 = pathpatch_2d_to_3d(ellps, z=eCentre[2], normal=(eNorm[0], eNorm[1], eNorm[2]))
# #ax.plot(ellps2[:,0],ellps2[:,1],ellps2[:,2])
# ellps = Ellipse((eCentre[0], eCentre[1]), width=2*eMajorRadius, height=2*eMinorRadius, angle=0, ec='red', fc=None, fill=False)
# ax.add_patch(ellps)
# #art3d.pathpatch_2d_to_3d(ellps, z=eCentre[2], zdir=(eNorm[0], eNorm[1], eNorm[2]))
# pathpatch_2d_to_3d(ellps, z=eCentre[2], normal=(eCentre[0], eCentre[1], eCentre[2]+1))
# #ellps3 = pathpatch_2d_to_3d(ellps, z=eCentre[2], normal=(eCentre[0], eCentre[1], eCentre[2]+1))
# #ax.plot(ellps3[:,0],ellps3[:,1],ellps3[:,2])
# plt.show()

#======================================================================

# sphere -> requires that surface must include orthogonal point that lies on line between sphere centre and originPoint
# ditto ellipsoid?

# 201 = ADVANCED_FACE( '', ( #429 ), #430, .T. );
    # 429 = FACE_OUTER_BOUND( '', #1347, .T. );
        # 1347 = EDGE_LOOP( '', ( #2247, #2248, #2249, #2250 ) );

            # 2247 = ORIENTED_EDGE( '', *, *, #2662, .F. );
                # 2662 = EDGE_CURVE( '', #2998, #2792, #2999, .T. );
                    # 2998 = VERTEX_POINT( '', #4110 );
                        # 4110 = CARTESIAN_POINT( '', ( 0.006858, -0.046046, 0.101750 ) );
                    # 2792 = VERTEX_POINT( '', #3284 );
                        # 3284 = CARTESIAN_POINT( '', ( 0.005741, -0.047307, 0.101750 ) );
                    # 2999 = CIRCLE( '', #4111, 0.00127 );
                        # 4111 = AXIS2_PLACEMENT_3D( '', #4530, #4531, #4532 );
                            # 4530 = CARTESIAN_POINT( '', ( 0.005588, -0.046046, 0.101750 ) );
                            # 4531 = DIRECTION( '', ( 0.0, 0.0, -1.0 ) );
                            # 4532 = DIRECTION( '', ( -1.0, 0.0, 0.0 ) );

            # 2248 = ORIENTED_EDGE( '', *, *, #2663, .F. );
                # 2663 = EDGE_CURVE( '', #2952, #2998, #3000, .F. );
                    # 2952 = VERTEX_POINT( '', #3859 );
                        # 3859 = CARTESIAN_POINT( '', ( 0.006858, -0.046046, 0.102243 ) );
                    # 2998 = VERTEX_POINT( '', #4110 );
                        # 4110 = CARTESIAN_POINT( '', ( 0.006858, -0.046046, 0.101750 ) );
                    # 3000 = LINE( '', #4112, #4113 );
                        # 4112 = CARTESIAN_POINT( '', ( 0.006858, -0.046046, 0.081923 ) );
                        # 4113 = VECTOR( '', #4533, 1.0 );
                        # 4533 = DIRECTION( '', ( 0.0, 0.0, 1.0 ) );

            # 2249 = ORIENTED_EDGE( '', *, *, #2626, .F. );
                # 2626 = EDGE_CURVE( '', #2787, #2952, #2953, .T. );
                    # 2787 = VERTEX_POINT( '', #3278 );
                        # 3278 = CARTESIAN_POINT( '', ( 0.005741, -0.047307, 0.102243 ) );
                    # 2952 = VERTEX_POINT( '', #3859 );
                        # 3859 = CARTESIAN_POINT( '', ( 0.006858, -0.046046, 0.102243 ) );
                    # 2953 = CIRCLE( '', #3860, 0.001270 );
                        # 3860 = AXIS2_PLACEMENT_3D( '', #4500, #4501, #4502 );
                        # 4500 = CARTESIAN_POINT( '', ( 0.005588, -0.046046, 0.102243 ) );
                        # 4501 = DIRECTION( '', ( 0.0, 0.0, 1.0 ) );
                        # 4502 = DIRECTION( '', ( 1.0, 0.0, 0.0 ) );

            # 2250 = ORIENTED_EDGE( '', *, *, #2523, .F. );
                # 2523 = EDGE_CURVE( '', #2792, #2787, #2794, .T. );
                    # 2792 = VERTEX_POINT( '', #3284 );
                        # 3284 = CARTESIAN_POINT( '', ( 0.005741, -0.047307, 0.101750 ) );
                    # 2787 = VERTEX_POINT( '', #3278 );
                        # 3278 = CARTESIAN_POINT( '', ( 0.005741, -0.047307, 0.102243 ) );
                    # 2794 = LINE( '', #3286, #3287 );
                        # 3286 = CARTESIAN_POINT( '', ( 0.005741, -0.047307, 0.101750 ) );
                        # 3287 = VECTOR( '', #4381, 1.00 );
                        # 4381 = DIRECTION( '', ( 0.0, 0.0, 1.0 ) );

            # 430 = CYLINDRICAL_SURFACE( '', #1348, 0.001270 );
                # 1348 = AXIS2_PLACEMENT_3D( '', #2251, #2252, #2253 );
                    # 2251 = CARTESIAN_POINT( '', ( 0.005588, -0.046046, 0.101750 ) );
                    # 2252 = DIRECTION( '', ( -0.0, -0.0, -1.0 ) );
                    # 2253 = DIRECTION( '', ( 1.0, 0.0, 0.0 ) );

# cylCentre: CYLINDRICAL_SURFACE[0] -> AXIS2_PLACEMENT_3D[0] -> CARTESIAN_POINT[:]
# cylFeatureAxis: CYLINDRICAL_SURFACE[0] -> AXIS2_PLACEMENT_3D[1] -> DIRECTION[:]
# cylLocalXaxis: CYLINDRICAL_SURFACE[0] -> AXIS2_PLACEMENT_3D[2] -> DIRECTION[:]
# cylRadius: CYLINDRICAL_SURFACE[1]

# cylinder test, check for orthogonal projection of origin to cylinder axis
# if projected point is within all edges ->
# both VERTEX of any line EDGE are greater or less than projected ray thru origin and orthogonal projected point
# circle/ellipse EDGE requires furthest point check & test that vertices & furthest point are all to one side of proj ray

def ref2index(s):
    return(int(s[1:]) - 1)

originPoint = np.array([1., 1., 20.])
cylCentre = np.array([0.0, 0.0, 0.0])
cylFeatureAxis = np.array([1.0, 0.0, 0.2])
cylLocalXaxis = np.array([0.0, 1.0, 0.2])
cylRadius = 20.

def CP2point(STEPobject, pointRefString):
    '''
    return np.array of STEP AP21 CARTESIAN_POINT, DIRECTION, VECTOR
    '''
    #assert STEPobject[ref2index(pointRefString)].type_name == "CARTESIAN_POINT"
    R3 = STEPobject[ref2index(pointRefString)].params[1]
    return np.array([R3[0], R3[1], R3[2]])

# strategy for cylinders & cones (& rotationally-symmetric spheres, surfaces of revolution or ellipsoids), determine curved edges
# if rot axis is between origin and surface, maxima lies on a line segment between maxima points on curved edges
# if rot axis is beyond surface (concave), minima lies on a line segment between minima points on curved edges
# otherwise maxima/minima lie on edges.

minmax="max"
if minmax=="max":
#    if AdvancedFaceSurfaces[0]['SurfaceTypeName'] == 'CYLINDRICAL_SURFACE':
# get rot axis
cylinderRef = AdvancedFaceSurfaces[n]["SurfaceTypeRef"]
axisRef = STEP_entities[ref2index(cylinderRef)].params[1]
assert(STEP_entities[ref2index(axisRef)].type_name =='AXIS2_PLACEMENT_3D')

axisOriginRef = STEP_entities[ref2index(axisRef)].params[1]
assert(STEP_entities[ref2index(axisOriginRef)].type_name =='CARTESIAN_POINT')
# R3 = STEP_entities[ref2index(axisOriginRef)].params[1]
# p = np.array([R3[0], R3[1], R3[2]])
axisOrigin = CP2point(STEP_entities, axisOriginRef)

axisDirRef = STEP_entities[ref2index(axisRef)].params[2]
assert(STEP_entities[ref2index(axisDir1Ref)].type_name =='DIRECTION')
axisOrigin = CP2point(STEP_entities, axisOriginRef)

axisLXdirRef = STEP_entities[ref2index(axisRef)].params[3]
assert(STEP_entities[ref2index(axisLXdirRef)].type_name =='DIRECTION')
axisLXdir = CP2point(STEP_entities, axisLXdirRef)

# axis: the Direction that defines the second axis of the Axis_placement. The value of this attribute need not be specified.
# ref_direction: the direction used to determine the direction of the local X axis.
# The value of this attribute need not be specified. If axis or ref_direction is omitted, these directions are taken from the geometric coordinate system

# get the circular edges from an edge loop
# get min/max point
# from edge orientation, vertex points, determine if min, max point (or both) lie on segment.

# 1293 = EDGE_LOOP( '', ( #2202 ) );
# 2202 = ORIENTED_EDGE( '', *, *, #2631, .F. );
# 2631 = EDGE_CURVE( '', #2959, #2959, #2960, .F. );

# 2959 = VERTEX_POINT( '', #3960 );
# 3960 = CARTESIAN_POINT( '', ( -0.0101600000000000, -0.0212872087123223, 0.105578161583723 ) );

# 2959 = VERTEX_POINT( '', #3960 );
# 3960 = CARTESIAN_POINT( '', ( -0.0101600000000000, -0.0212872087123223, 0.105578161583723 ) );

# 2960 = CIRCLE( '', #3961, 0.000475311400000000 );
# 3961 = AXIS2_PLACEMENT_3D( '', #4506, #4507, #4508 );
# 4506 = CARTESIAN_POINT( '', ( -0.0101600000000000, -0.0212872087123223, 0.106053472961653 ) );
# 4507 = DIRECTION( '', ( 1.00000000000000, 0.00000000000000, 0.00000000000000 ) );
# 4508 = DIRECTION( '', ( 0.00000000000000, 0.00000000000000, -1.00000000000000 ) );

##############

# 		# 1347 = EDGE_LOOP( '', ( #2247, #2248, #2249, #2250 ) );
#
# 				#2247 = ORIENTED_EDGE( '', *, *, #2662, .F. );
# 					#2662 = EDGE_CURVE( '', #2998, #2792, #2999, .T. );
# 						#2998 = VERTEX_POINT( '', #4110 );
# 							#4110 = CARTESIAN_POINT( '', ( 0.00685804782754277, -0.0460465840350917, 0.101750792961653 ) );
#
# 						#2792 = VERTEX_POINT( '', #3284 );
# 							#3284 = CARTESIAN_POINT( '', ( 0.00574104793264291, -0.0473073341958450, 0.101750792961653 ) );
#
# 						#2999 = CIRCLE( '', #4111, 0.00127 );
# 							#4111 = AXIS2_PLACEMENT_3D( '', #4530, #4531, #4532 );
# 								#4530 = CARTESIAN_POINT( '', ( 0.00558804782754277, -0.0460465840350917, 0.101750792961653 ) );
# 								#4531 = DIRECTION( '', ( 0.00000000000000, 0.00000000000000, -1.00000000000000 ) );
# 								#4532 = DIRECTION( '', ( -1.00000000000000, 0.00000000000000, 0.00000000000000 ) );

# orientation: a BOOLEAN flag. If TRUE, the topological orientation as used coincides with the orientation,
# from start vertex to end vertex, of the edge_definition, if FALSE the vertices are reversed in order.
# edge_start: the start vertex of the Oriented_edge. This is derived from the vertices of the edge_definition after taking account of the orientation.
# edge_end: the end vertex of the Oriented_edge. This is derived from the vertices of the edge_definition after taking account of the orientation.

cylinderRef = AdvancedFaceSurfaces[-1]["EdgeLoopList"]
			
cylOrthoPoint = orthoProj2line(originPoint, cylCentre, cylFeatureAxis)
# extend intersection to cylinder surface
localCylAngle = np.arctan2(cylOrthoPoint.y, cylOrthoPoint.x)
cylExtentConvex = [cylOrthoPoint.x + cylRadius*np.cos(localCylAngle), cylOrthoPoint.y + cylRadius*np.sin(localCylAngle), originPoint.z] #Point
cylExtentConcave = [cylOrthoPoint.x - cylRadius*np.cos(localCylAngle), cylOrthoPoint.y - cylRadius*np.sin(localCylAngle), originPoint.z] #Point

edgeCrossingCount = 0
for oe in local_edge_loop: # EDGE_LOOP-ORIENTED_EDGE list of ORIENTED_EDGE
    if oe.ref == 'LINE':
        if (oe.vertex1.x > cylOrthoPoint.x and oe.vertex2.x < cylOrthoPoint.x) or \
           (oe.vertex1.x < cylOrthoPoint.x and oe.vertex2.x > cylOrthoPoint.x) or \
           (oe.vertex1.y > cylOrthoPoint.y and oe.vertex2.y > cylOrthoPoint.y) or \
           (oe.vertex1.y < cylOrthoPoint.y and oe.vertex2.y < cylOrthoPoint.y):
            edgeCrossingCount += 1

    if oe.ref == 'CIRCLE':
        circleSegMax = pointCircleMinDisp(p, circleCentre, circleLocalXYaxis, circleFeatureAxis, circleRadius, minmax="max")

    if edgeCrossingCount % 2 > 0:
        print("edge error")
    elif edgeCrossingCount > 1:
        if local_surface.ref == 'PLANE':
            # do triangles check for minimum
            # do vertex test for maxima
            if search == 'minima':



#-----------------------------------------------------------------------------

# Intersection of two circles.
# Returns 'null' (no intersection) or object of type 'Circle3d', 'Point3d' or 'Segment3d'.
# In 2D (coplanar circles) the segment will define two intersection points.

def IntersectionWith(c1, c2):
    # Relative tolerance ================================
    if (!GeometRi3D.UseAbsoluteTolerance):
        tol = GeometRi3D.Tolerance
        GeometRi3D.Tolerance = tol * Max(this.R, c.R)
        GeometRi3D.UseAbsoluteTolerance = true
        object result = this.IntersectionWith(c)
        GeometRi3D.UseAbsoluteTolerance = false
        GeometRi3D.Tolerance = tol
        return result
    #====================================================

    if (c1.Norm.IsParallelTo(c2.Norm)):
        if (this.Center.BelongsTo(new Plane3d(c2.Center, c2.Normal))):
            # Coplanar objects
            # Search 2D intersection of two circles

            # Equal circles
            if (c1.Center == c2.Center && GeometRi3D.AlmostEqual(c1.R, c2.R)):
                return c1

            d = c1.Center.DistanceTo(c2.Center)
            # Separated circles
            if (GeometRi3D.Greater(d, c1.R + c2.R)):
                return null

            # One circle inside the other
            if (d < Abs(c1.R - c2.R) - GeometRi3D.Tolerance):
                if (this.R > c.R):
                    return c2
                else:
                    return c1

            # Outer tangency
            if (GeometRi3D.AlmostEqual(d, c1.R + c2.R)):
                vec = new Vector3d(c1.Center, c2.Center)
                return c1.Center.Translate(c1.R * vec.Normalized)

            # Inner tangency
            if (Abs(Abs(c1.R - c2.R) - d) < GeometRi3D.Tolerance):
                vec = new Vector3d(c1.Center, c2.Center)
                if (c1.R > c2.R):
                    return c1.Center.Translate(c1.R * vec.Normalized)
            else:
                return c1.Center.Translate(-c1.R * vec.Normalized)

            # intersecting circles
            # Create local CS with origin in circle's center
            vec1 = Vector3d(c1.Center, c2.Center)
            vec2 = vec1.Cross(c1.Normal)
            local_cs = Coord3d(c1.Center, vec1, vec2)
            x = 0.5 * (d * d - c1.R * c1.R + c2.R * c2.R) / d
            y = 0.5 * Sqrt((-d + c2.R - c1.R) * (-d - c2.R + c1.R) * (-d + c2.R + c1.R) * (d + c2.R + c1.R)) / d
            p1 = Point3d(x, y, 0, local_cs)
            p2 = Point3d(x, -y, 0, local_cs)
            return Segment3d(p1, p2)

        else:
            # parallel objects
            return null
    else:
        # Check 3D intersection
        plane = Plane3d(c1.Center, c1.Normal)
        obj = plane.IntersectionWith(c2)
        if (obj == null):
            return null
        elif (obj.GetType() == typeof(Point3d)):
            p = (Point3d)obj
            if (p.BelongsTo(c2)):
                return p
            else:
                return null
        else:
            s = (Segment3d)obj
            return s.IntersectionWith(c1)

     def _PointLocation(p):
         if (GeometRi3D.UseAbsoluteTolerance):
             s = new Plane3d(c1.Center, c1.Normal)
             proj = p.ProjectionTo(s)
             if (GeometRi3D.AlmostEqual(p.DistanceTo(proj), 0)):
                 if ( GeometRi3D.Smaller(p.DistanceTo(c1.Center), c1.R)):
                     return 1; # Point is strictly inside
                 elif (GeometRi3D.AlmostEqual(p.DistanceTo(c1.Center), c1.R) ):
                     return 0 # Point is on boundary
                 else:
                     return -1 # Point is outside
                else:
                    return -1 # Point is outside
            else
                double tol = GeometRi3D.Tolerance;
                GeometRi3D.Tolerance = tol * c1.R;
                GeometRi3D.UseAbsoluteTolerance = true;
                int result = c1._PointLocation(p);
                GeometRi3D.UseAbsoluteTolerance = false;
                GeometRi3D.Tolerance = tol;
                return result;
            }
        }


# # /// <param name="p">Target plane</param>
# # /// <param name="point_on_circle">Closest point on circle</param>
# # /// <param name="point_on_plane">Closest point on plane</param>
# def DistanceTo(p, c):
#     if (this.IsParallelTo(p)): #this.Normal.IsOrthogonalTo(obj.Direction)
#         point_on_circle = this.Center
#         point_on_plane = point_on_circle.ProjectionTo(p)
#         return point_on_circle.DistanceTo(point_on_plane)
#
#     v1 = np.cross(cNorm, pNorm)
#     v2 = np.cross(cNorm, v1)
#     l = Line3d(cCenter, v2)
#     intersection_point = l.IntersectionWith(p)
#     Vector3d r1 = new Vector3d(l.Point);
#     Vector3d s1 = l.Direction;
#     n2 = lNorm
#
#     SetCoord(r1.Coord)
#     r1 = r1 - ((r1 * n2) + this.D) / (s1 * n2) * s1;
# return r1.ToPoint;
#
#     if (intersection_point.DistanceTo(this) <= GeometRi3D.DefaultTolerance):
#         point_on_circle = intersection_point
#         point_on_plane = intersection_point
#         return 0
#     else:
#         v1 = Vector3d(cCenter, intersection_point).Normalized
#         point_on_circle = cCenter.Translate(cRadius * v1)
#         point_on_plane = point_on_circle.ProjectionTo(p)
#         return point_on_circle.DistanceTo(point_on_plane)

# Intersection of circle with plane.
# Returns 'null' (no intersection) or object of type 'Circle3d', 'Point3d' or 'Segment3d'.
def IntersectionWith(c, s):
# Relative tolerance ================================
# if (!GeometRi3D.UseAbsoluteTolerance)
#
#     double tol = GeometRi3D.Tolerance;
# GeometRi3D.Tolerance = tol * this.R;
# GeometRi3D.UseAbsoluteTolerance = true;
# object result = this.IntersectionWith(s);
# GeometRi3D.UseAbsoluteTolerance = false;
# GeometRi3D.Tolerance = tol;
# return result;
#
# //====================================================

# def intersectPlanePlane(p1, p2):
#     #  intersection of 2-planes: a variation based on 3-plane version
#     #  see: Graphics Gems 1, pg 305
#     #  planes normal need not be unit length
#     #output args r_point r_normal
#     #  logically the 3rd plane, but we only use the normal component.
# p3_normal = np.linalg.norm(np.cross(p1, np.linalg.norm(p2)))
# det = p3_normal.length_squared()
#
# # determinant = 0, => parallel planes, no intersection
# # note: you may want to check against an epsilon value here.
# if (det != 0.0) {
#                 // calculate the final (point, normal)
# r_point = ((p3_normal.cross(p2.normal) * p1.d) +
#            (p1.normal.cross(p3_normal) * p2.d)) / det;
# r_normal = p3_normal;
# return true;
# }
# else {
# return false;
# }
# }


# convert centroid, cMinPoint, cMaxPoint edgeVertex1 & edgeVertex2 to local coordinates with axis as z

# axisDir, refDir guaranteed unit normals?


if (cNorm.IsParallelTo(s.Norm)):
    if (cCenter.BelongsTo(s)):
        # coplanar objects
        return c
    else:
        # parallel objects
        return null
else:
    # create plane on circle plane
    # get l, line along planes intersection
    #l = s.IntersectionWith(new Plane3d(circleCenter, circleNormal))

    # create transformation to local coords
    local_coord = Coord3d(cCenter, l.Direction, cNorm.Cross(l.Direction))
    p = l.Point.ConvertTo(local_coord);

    if (GeometRi3D.Greater(Abs(p.Y), cRadius)):
        return null
    elif (GeometRi3D.AlmostEqual(p.Y, this.R)):
        return Point3d(0, this.R, 0, local_coord)
    elif (GeometRi3D.AlmostEqual(p.Y, -this.R)):
        return Point3d(0, -this.R, 0, local_coord)
    else:
        d = Sqrt(Math.Pow(this.R, 2) - Math.Pow(p.Y, 2))
        p1 = Point3d(-d, p.Y, 0, local_coord)
        p2 = Point3d(d, p.Y, 0, local_coord)
        return Segment3d(p1, p2);

# coordinate system using origin point and two vectors
# <param name="p">Origin of the coordinate system.
# <param name="v1">Vector oriented along the X axis.</param>
# <param name="v2">Vector in the XY plane.</param>
def Coord3d(Point3d p, v1, v2):
if (v1.IsParallelTo(v2)):
    throw Exception("Vectors are parallel")

v1 = v1.ConvertToGlobal().Normalized
v3 = v1.Cross(v2).Normalized
v2 = v3.Cross(v1).Normalized

_origin = p.ConvertToGlobal()
_axes = Matrix3d(v1, v2, v3)

            # intersection of two planes.
            # Returns 'null' (no intersection) or object of type 'Line3d' or 'Plane3d'.
        def IntersectionWith(Plane3d s2):
        v = this.Normal.Cross(s2.Normal).ConvertToGlobal()
        if (v.Norm < GeometRi3D.Tolerance):
            # Planes are coplanar
            if (this.Point.BelongsTo(s2)):
                return this
            else:
                return null

        else:
            # Find the common point for two planes by intersecting with third plane
            # (using the 'most orthogonal' plane)
            # This part needs to be rewritten
            if (Abs(v.X) >= Abs(v.Y) && Abs(v.X) >= Abs(v.Z)):
                p = (Point3d)Coord3d.GlobalCS.YZ_plane.IntersectionWith(this, s2)
                return new Line3d(p, v)
            elif (Abs(v.Y) >= Abs(v.X) && Abs(v.Y) >= Abs(v.Z)):
                p = (Point3d)Coord3d.GlobalCS.XZ_plane.IntersectionWith(this, s2)
                return new Line3d(p, v)
            else:
                p = (Point3d)Coord3d.GlobalCS.XY_plane.IntersectionWith(this, s2)
                return new Line3d(p, v)


# define plane of maxSurfacePoint and rot-sym axisDir for cylinder, cone, toroid, etc,
# transform edge curve from global coords to local coords defined by maxSurfacePlane st maxSurfacePlane is coplanar with xy-plane
# find x,y of edge arc points where arc z=0, also where constant max/minPoint radius from rotsymAxis along z

if pel['typeName'] == 'CIRCLE':
    edgeArcAxisPoint = pel['axisPoint']

    # project maxPoint, minPoint to surface rotSym axisDir
    maxPointProjAxisPoint = axisPoint + (np.dot(maxPoint - axisPoint, axisDir)/np.dot(axisDir, axisDir))*axisDir

    rsNorm = np.cross(axisDir, refDir) # values from surface calc
    rsNorm = rsNorm / np.linalg.norm(rsNorm)
    rotM = np.array([refDir, rsNorm, axisDir])

    localEdgeArcAxisPoint = np.matmul(rotM, edgeArcAxisPoint - maxPointProjAxisPoint)

        # cMinPointUV = np.matmul(rotM, cMinPoint - axisPoint)
        # cMaxPointUV = np.matmul(rotM, cMaxPoint - axisPoint)
        # vertex1UV = np.matmul(rotM, vertex1 - axisPoint)
        # vertex2UV = np.matmul(rotM, vertex2 - axisPoint)
        #
        # # sanity check that Z-values are identical
        # zTest = [cMinPointUV[2], cMaxPointUV[2], vertex1UV[2], vertex2UV[2]]
        # if not np.isclose(max(zTest), min(zTest)):
        #     print(max(zTest) - min(zTest))

        # convert to polar angle
        cMinPointUV_angle = np.arctan2(cMinPointUV[1], cMinPointUV[0])
        cMaxPointUV_angle = np.arctan2(cMaxPointUV[1], cMaxPointUV[0])
        vertex1UV_angle = np.arctan2(vertex1UV[1], vertex1UV[0])
        vertex2UV_angle = np.arctan2(vertex2UV[1], vertex2UV[0])

        if not VOR:
            if ((cMinPointUV_angle > vertex1UV_angle) and (cMinPointUV_angle < vertex2UV_angle)) or \
                    (vertex1UV_angle == vertex2UV_angle):
                parsedEdge['minPoint'] = cMinPoint
            else:
                if np.linalg.norm(vertex1 - centroid) < np.linalg.norm(vertex2 - centroid):
                    parsedEdge['minPoint'] = vertex1
                else:
                    parsedEdge['minPoint'] = vertex2

            if ((cMaxPointUV_angle > vertex1UV_angle) and (cMaxPointUV_angle < vertex2UV_angle)) or \
                    (vertex1UV_angle == vertex2UV_angle):
                parsedEdge['maxPoint'] = cMaxPoint
            else:
                if np.linalg.norm(vertex1 - centroid) > np.linalg.norm(vertex2 - centroid):
                    parsedEdge['maxPoint'] = vertex1
                else:
                    parsedEdge['maxPoint'] = vertex2

        else: # vertex order reversed
            if (cMinPointUV_angle > vertex2UV_angle) and (cMinPointUV_angle < vertex1UV_angle):
                parsedEdge['minPoint'] = cMinPoint
            else:
                if np.linalg.norm(vertex1 - centroid) < np.linalg.norm(vertex2 - centroid):
                    parsedEdge['minPoint'] = vertex1
                else:
                    parsedEdge['minPoint'] = vertex2

            if (cMaxPointUV_angle > vertex2UV_angle) and (cMaxPointUV_angle < vertex1UV_angle):
                parsedEdge['maxPoint'] = cMaxPoint
            else:
                if np.linalg.norm(vertex1 - centroid) > np.linalg.norm(vertex2 - centroid):
                    parsedEdge['maxPoint'] = vertex1
                else:
                    parsedEdge['maxPoint'] = vertex2



