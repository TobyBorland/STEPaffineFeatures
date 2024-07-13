import numpy as np


def intersectLinePlane(l1, l2, pPoint, pNorm, tol=1e-6):
    """Computes the intersection point of a line and a plane

    Parameters
    ----------
    line : [point, point] | :class:`~compas.geometry.Line`
        Two points defining the line.
    plane : [point, vector] | :class:`~compas.geometry.Plane`
        The base point and normal defining the plane.
    tol : float, optional
        A tolerance for membership verification.

    Returns
    -------
    [float, float, float] | None
        The intersection point between the line and the plane,
        or None if the line and the plane are parallel.

    """
    # l1, l2 = line
    # pPoint, pNorm = plane

    L = l2 - l1
    cosl = np.dot(pNorm, L)

    if np.abs(cosl) <= tol:
        # if the dot product (cosine of the angle between segment and plane)
        # is close to zero the line and the normal are almost perpendicular
        # hence there is no intersection
        return None

    # based on the ratio = -dot_vectors(n, ab) / dot_vectors(n, oa)
    # there are three scenarios
    # 1) 0.0 < ratio < 1.0: the intersection is between l1 and l2
    # 2) ratio < 0.0: the intersection is on the other side of a
    # 3) ratio > 1.0: the intersection is on the other side of l2

    L = L * (-np.dot(pNorm, l1 - pPoint) / cosl)
    return l1 + L

def intersectPlanePlane(p1point, p1norm, p2point, p2norm, tol=1e-6):
    """Computes the intersection of two planes

    Parameters
    ----------
    plane1 : [point, vector] | :class:`~compas.geometry.Plane`
        The base point and normal (normalized) defining the 1st plane.
    plane2 : [point, vector] | :class:`~compas.geometry.Plane`
        The base point and normal (normalized) defining the 2nd plane.
    tol : float, optional
        A tolerance for membership verification.

    Returns
    -------
    tuple[[float, float, float], [float, float, float]] | None
        Two points defining the intersection line.
        None if the planes are parallel.

    """
    # o1, p1norm = plane1
    # p2point, p2norm = plane2

    if np.abs(np.dot(p1norm, p2norm)) >= 1 - tol:
        return None

    # direction of intersection line
    d = np.cross(p1norm, p2norm)
    # vector in plane 1 perpendicular to the direction of the intersection line
    v1 = np.cross(d, p1norm)
    # point on plane 1
    p1 = p1point + v1

    x1 = intersectLinePlane(p1point, p1, p2point, p2norm, tol=tol)
    x2 = x1 + d
    return x1, x2

def intersectSphereLine(sPoint, sRadius, l1, l2):
    # https://github.com/compas-dev/compas/blob/main/src/compas/geometry/intersections/intersections.py#L648
    """Computes the intersection of a sphere and a line.

    Parameters
    ----------
    sphere : [point, radius] | :class:`~compas.geometry.Sphere`
        A sphere defined by a point and a radius.
    line : [point, point] | :class:`~compas.geometry.Line`
        A line defined by two points.

    Returns
    -------
    tuple[[float, float, float], [float, float, float]] | [float, float, float] | None
        Two points (if the line goes through the sphere), one point (if the line is tangent to the sphere), or None (otherwise).

    Notes
    -----
    There are 3 cases of sphere-line intersection:

    1. they intersect in 2 points
    2. they intersect in 1 point (line tangent to sphere), or
    3. they do not intersect.

    Examples
    --------
    >>> from compas.geometry import allclose

    >>> sphere = (3.0, 7.0, 4.0), 10.0
    >>> line = (1.0, 0, 0.5), (2.0, 1.0, 0.5)
    >>> x1, x2 = intersection_sphere_line(sphere, line)

    >>> allclose(x1, [11.634, 10.634, 0.500], 1e-3)
    True
    >>> allclose(x2, [-0.634, -1.634, 0.50], 1e-3)
    True

    """
    # l1, l2 = line
    # sPoint, sRadius = sphere

    a = (l2[0] - l1[0]) ** 2 + (l2[1] - l1[1]) ** 2 + (l2[2] - l1[2]) ** 2
    b = 2.0 * (
        (l2[0] - l1[0]) * (l1[0] - sPoint[0]) +
        (l2[1] - l1[1]) * (l1[1] - sPoint[1]) + (
                l2[2] - l1[2]) * (l1[2] - sPoint[2])
    )

    c = (
        sPoint[0] ** 2
        + sPoint[1] ** 2
        + sPoint[2] ** 2
        + l1[0] ** 2
        + l1[1] ** 2
        + l1[2] ** 2
        - 2.0 * (sPoint[0] * l1[0] + sPoint[1] * l1[1] + sPoint[2] * l1[2])
        - sRadius**2
    )

    i = b * b - 4.0 * a * c

    if i < 0.0:  # case 3: no intersection
        return None
    elif i == 0.0:  # case 2: one intersection
        mu = -b / (2.0 * a)
        ipt = (
            l1[0] + mu * (l2[0] - l1[0]),
            l1[1] + mu * (l2[1] - l1[1]),
            l1[2] + mu * (l2[2] - l1[2]),
        )
        return ipt
    elif i > 0.0:  # case 1: two intersections
        # 1.
        mu = (-b + np.sqrt(i)) / (2.0 * a)
        ipt1 = (
            l1[0] + mu * (l2[0] - l1[0]),
            l1[1] + mu * (l2[1] - l1[1]),
            l1[2] + mu * (l2[2] - l1[2]),
        )
        # 2.
        mu = (-b - np.sqrt(i)) / (2.0 * a)
        ipt2 = (
            l1[0] + mu * (l2[0] - l1[0]),
            l1[1] + mu * (l2[1] - l1[1]),
            l1[2] + mu * (l2[2] - l1[2]),
        )
        return ipt1, ipt2

def intersectPlaneCircle(pPoint, pNorm, cPoint, cNorm, cRadius):
    """Computes the intersection of a plane and a circle.

    Parameters
    ----------
    plane : [point, vector] | :class:`~compas.geometry.Plane`
        A plane defined by a point and normal vector.
    circle : [plane, float] | :class:`~compas.geometry.Circle`
        A circle defined by a plane and radius.

    Returns
    -------
    tuple[[float, float, float], [float, float, float]] | [float, float, float] | None
        Two points (secant intersection), one point (tangent intersection), or None (otherwise).

    Notes
    -----
    There are 4 cases of plane-circle intersection:

    1. they intersect in 2 points (secant),
    2. they intersect in 1 point (tangent),
    3. they do not intersect, or
    4. they coincide (circle.plane == plane).

    Examples
    --------
    >>> plane = (0, 0, 0), (0, 0, 1)
    >>> circle = ((0, 0, 0), (0, 1, 0)), 10.0
    >>> x1, x2 = intersection_plane_circle(plane, circle)
    >>> x1
    (-10.0, 0.0, 0.0)
    >>> x2
    (10.0, 0.0, 0.0)

    """
    #cPlane, cRadius = circle
    #line = intersectPlanePlane(plane, cPlane)
    l1, l2 = intersectPlanePlane(pPoint, pNorm, cPoint, cNorm)
    if (l1 is None) or (l2 is None):
        return None
    #cPoint = cPlane[0]
    #sphere = cPoint, cRadius
    return intersectSphereLine(cPoint, cRadius, l1, l2)


def intersection_ellipse_line_xy(ellipse, line):
    """Computes the intersection of an ellipse and a line in the XY plane.

    Parameters
    ----------
    ellipse : tuple[float, float]
        The major and minor of the ellipse.
    line : [point, point] | :class:`~compas.geometry.Line`
        A line defined by two points, with at least XY coordinates.

    Returns
    -------
    tuple[[float, float, float], [float, float, float]] | [float, float, float] | None
        Two points, if the line goes through the ellipse.
        One point, if the line is tangent to the ellipse.
        None, otherwise.

    References
    ----------
    Based on [1]_.

    .. [1] C# Helper. *Calculate where a line segment and an ellipse intersect in C#*.
           Available at: http://csharphelper.com/blog/2017/08/calculate-where-a-line-segment-and-an-ellipse-intersect-in-c/

    Examples
    --------
    >>> ellipse = 6., 2.5
    >>> p1 = (4.1, 2.8, 0.)
    >>> p2 = (3.4, -3.1, 0.)
    >>> i1, i2 = intersection_ellipse_line_xy(ellipse, [p1, p2])

    """
    x1, y1 = line[0][0], line[0][1]
    x2, y2 = line[1][0], line[1][1]

    a, b = ellipse

    A = (x2 - x1) ** 2 / a**2 + (y2 - y1) ** 2 / b**2
    B = 2 * x1 * (x2 - x1) / a**2 + 2 * y1 * (y2 - y1) / b**2
    C = x1**2 / a**2 + y1**2 / b**2 - 1

    discriminant = B**2 - 4 * A * C
    if discriminant == 0:
        t = -B / (2 * A)
        return (x1 + (x2 - x1) * t, y1 + (y2 - y1) * t, 0.0)
    elif discriminant > 0:
        t1 = (-B + np.sqrt(discriminant)) / (2 * A)
        t2 = (-B - np.sqrt(discriminant)) / (2 * A)
        p1 = (x1 + (x2 - x1) * t1, y1 + (y2 - y1) * t1, 0.0)
        p2 = (x1 + (x2 - x1) * t2, y1 + (y2 - y1) * t2, 0.0)
        return p1, p2
    else:
        return None



def intersectPlaneCircle2(pPoint, pNorm, cPoint, cNorm, cRadius):
    """Computes the intersection of a plane and a circle.

    Parameters
    ----------
    plane : [point, vector] | :class:`~compas.geometry.Plane`
        A plane defined by a point and normal vector.
    circle : [plane, float] | :class:`~compas.geometry.Circle`
        A circle defined by a plane and radius.

    Returns
    -------
    tuple[[float, float, float], [float, float, float]] | [float, float, float] | None
        Two points (secant intersection), one point (tangent intersection), or None (otherwise).

    Notes
    -----
    There are 4 cases of plane-circle intersection:

    1. they intersect in 2 points (secant),
    2. they intersect in 1 point (tangent),
    3. they do not intersect, or
    4. they coincide (circle.plane == plane).

    Examples
    --------
    >>> plane = (0, 0, 0), (0, 0, 1)
    >>> circle = ((0, 0, 0), (0, 1, 0)), 10.0
    >>> x1, x2 = intersection_plane_circle(plane, circle)
    >>> x1
    (-10.0, 0.0, 0.0)
    >>> x2
    (10.0, 0.0, 0.0)

    """
    #l1, l2 = intersectPlanePlane(pPoint, pNorm, cPoint, cNorm)

    tol=1e-6
    if np.abs(np.dot(pNorm, cNorm)) >= 1 - tol:
        return None
    # direction of intersection line
    d = np.cross(pNorm, cNorm)
    # vector in plane 1 perpendicular to the direction of the intersection line
    v1 = np.cross(d, pNorm)
    # point on plane 1
    p1 = pPoint + v1

    L = p1 - pPoint
    cosl = np.dot(cNorm, L)

    if np.abs(cosl) <= tol:
        # if the dot product (cosine of the angle between segment and plane)
        # is close to zero the line and the normal are almost perpendicular
        return None

    L = L * (-np.dot(cNorm, pPoint - cPoint) / cosl)
    l1 = pPoint + L
    l2 = l1 + d

    if (l1 is None) or (l2 is None):
        return None

    a = (l2[0] - l1[0]) ** 2 + (l2[1] - l1[1]) ** 2 + (l2[2] - l1[2]) ** 2
    b = 2.0 * (
        (l2[0] - l1[0]) * (l1[0] - cPoint[0]) +
        (l2[1] - l1[1]) * (l1[1] - cPoint[1]) + (
                l2[2] - l1[2]) * (l1[2] - cPoint[2])
    )

    c = (
        cPoint[0] ** 2
        + sPoint[1] ** 2
        + sPoint[2] ** 2
        + l1[0] ** 2
        + l1[1] ** 2
        + l1[2] ** 2
        - 2.0 * (cPoint[0] * l1[0] + cPoint[1] * l1[1] + cPoint[2] * l1[2])
        - cRadius**2
    )

    i = b * b - 4.0 * a * c

    if i < 0.0:  # case 3: no intersection
        return None
    elif i == 0.0:  # case 2: one intersection
        mu = -b / (2.0 * a)
        ipt = (
            l1[0] + mu * (l2[0] - l1[0]),
            l1[1] + mu * (l2[1] - l1[1]),
            l1[2] + mu * (l2[2] - l1[2]),
        )
        return ipt
    elif i > 0.0:  # case 1: two intersections
        # 1.
        mu = (-b + np.sqrt(i)) / (2.0 * a)
        ipt1 = (
            l1[0] + mu * (l2[0] - l1[0]),
            l1[1] + mu * (l2[1] - l1[1]),
            l1[2] + mu * (l2[2] - l1[2]),
        )
        # 2.
        mu = (-b - np.sqrt(i)) / (2.0 * a)
        ipt2 = (
            l1[0] + mu * (l2[0] - l1[0]),
            l1[1] + mu * (l2[1] - l1[1]),
            l1[2] + mu * (l2[2] - l1[2]),
        )
        return ipt1, ipt2











pP = np.array([0, 0, 0])
pN = np.array([0, 0, 1])
cP = np.array([0, 0, 0])
cN = np.array([0, 1, 0])
cR = 10.0
x1, x2 = intersectPlaneCircle(pP, pN, cP, cN, cR)

# Intersection of 2-planes: a variation based on the 3-plane version.
# see: Graphics Gems 1 pg 305
# Note that the 'normal' components of the planes need not be unit length
# bool isect_plane_plane_to_normal_ray(
#         const Plane& p1, const Plane& p2,
#         // output args
#         Vector3f& r_point, Vector3f& r_normal)
# {
#     // logically the 3rd plane, but we only use the normal component.
#   p3_normal = p1.normal.cross(p2.normal);
#     const float det = p3_normal.length_squared();
#
#     // If the determinant is 0, that means parallel planes, no intersection.
#     // note: you may want to check against an epsilon value here.
#     if (det != 0.0) {
#         // calculate the final (point, normal)
#         r_point = ((p3_normal.cross(p2.normal) * p1.d) +
#                    (p1.normal.cross(p3_normal) * p2.d)) / det;
#         r_normal = p3_normal;
#         return true;
#     }
#     else {
#         return false;
#     }
# }
#
# def intersectPlanePlane(p1, p2):
#     #  intersection of 2-planes: a variation based on 3-plane version
#     #  see: Graphics Gems 1, pg 305
#     #  planes normal need not be unit length
#     #output args r_point r_normal
#     #  logically the 3rd plane, but we only use the normal component.
#     p1norm = p2/np.linalg.norm(p1)
#     p2norm = p2/np.linalg.norm(p2)
#     p3norm = np.cross(p1, p2)
#     p3norm = p3 / np.linalg.norm(p3)
#     det = np.linalg.norm(p3norm)
#     det = det*det
#
#
#     # determinant = 0, => parallel planes, no intersection
#     # note: you may want to check against an epsilon value here.
#     if (det > 1e-16):
#         # calculate the final (point, normal)
#         rPoint = np.cross(p3norm, p2norm) * p1.d +  (p1.normal.cross(p3_normal) * p2.d)) / det;
# r_normal = p3_normal;
# return true;
# }
# else:
#     return False
#
#
#
# # Circle intersection
#
# c1 = np.array([])
#
# cDiff = (c1 - c2)
# cDiffDisp = np.linalg.norm(cDiff)
#
# if cDiffDisp > (r1 + r2):
#     print(" no intersection")
# if cDiffDisp == (r1 + r2):
#     print(" tangent intersection")
#
# p = (c1*r2 + c2*r1) / (r1 + r2)
#
# if cDiffDisp < (r2 - r1):
#     print(" no intersection")
# if cDiffDisp == (r2 - r1):
#     print(" tangent intersection")
#
# p = (c1 - c2) * r2 /(r2 - r1)
#
# # otherwise two intersection points
# cDiffNorm = cDiff/cDiffDisp
# cDiffOrthog = cDiffNorm * n1.Normalized
# q = cDiffDisp*cDiffDisp + r2*r2 - r1*r1
# # dx = 1/2 * q / cDiffDisp
# # dy = 1/2 * Sqrt(4 * cDiffDisp^2 * r2^2 - q^2) / cDiffDisp
# # p1,2 = c1 + cDiffNorm * dx +/- cDiffOrthog * dy