
# *****SURFACE*****

# Plane - position - Axis_placement_3d
# (detect if plane normal passes through origin)

# Conical_surface - position - Axis_placement_3d
#                 - radius - Positive_length_measure
#                 - semi_angle - plane_angle_measure
# ( is cone axis parallel to ray through origin?,
#  are cone surfaces orthogonal to ray through origin?)

# Cylindrical_surface - position - Axis_placement_3d
#                     - radius - Positive_length_measure
# ( is cylindrical surface perpendicular to ray through origin?
# ( if not, extract point from bounding curve)

# Spherical_surface - position - Axis_placement_3d
#                   - radius - Positive_length_measure
# ( point on surface cotangent to centre and origin point)

# Swept_surface - swept_curve - Curve

# Swept_surface - Extruded_surface - extrusion_axis - Direction

# Swept_surface - Surface_of_revolution - axis_direction - Direction
#                                       - axis_point - Cartesian_point

# Toroidal_surface - position - Axis_placement_3d
#                  - radius - Positive_length_measure
#                  - minor_radius - Positive_length_measure
# ( if toroid is complete at outer/inner radius, increment rut-symm maxima/minima )

# Bounded_surface - B_spline_surface - U_degree - integer
#                                    - V_degree - integer
#                                    - control_points - Cartesian_point
#                                    - U_closed - boolean
#                                    - V_closed - boolean
# (part of complex_entity listing, need to find footpoint for maxima/minima)

# Bounded_surface - B_spline_surface - Rational_B_spline_surface - real

# Bounded_surface - B_spline_surface - Surface_with_explicit_knots - U_knot_multiplicities - integer
#                                                                  - U_knot_values - Parameter_value
#                                                                  - V_knot_multiplicities - integer

# Bounded_surface - B_spline_surface - Surface_with_implicit_knots - Knot_specification

# Bounded_surface - Curve_bounded_surface - basis_surface - Surface
#                                         - boundaries - Boundary_curve
#                                         - implicit_outer - boolean

# Bounded_surface - Rectangular_composite_surface - segments - Surface_patch
#                                                 - using_surfaces - Surface_patch
# Bounded_surface - parent_surface - Surface_patch - u_sense - boolean
#                                                  - v_sense - boolean
#                                                  - u_transition - Surface_transition_code
#                                                  - v_transition - Surface_transition_code

# Bounded_surface - Trimmed_surface - u1 - Parameter_value - real
#                                   - u2 - Parameter_value - real
#                                   - v1 - Parameter_value - real
#                                   - v2 - Parameter_value - real
#                                   - usense - boolean
#                                   - vsense - boolean

# *****CURVE*****

# Conic - position - Axis_placement

# Conic - Circle - Closed_curve (optional) - Closed_composite_curve
#                - radius - Positive_length_measure
# ( does circle axis pass thru origin? )

# Conic - Ellipse - first_semi_axis - Positive_length_measure
#                 - second_semi_axis - Positive_length_measure

# Conic - Hyperbola - semi_axis - Positive_length_measure
#                   - imaginary_semi_axis - Positive_length_measure

# Conic - Parabola - focal_distance - Length_measure

# Line - point - Cartesian_point
#      - line_direction - Direction

# Bounded_curve - segment_curve - Composite_curve_segment

# Bounded_curve - Composite_curve - curve_segment - Composite_curve_segment
# Bounded_curve - Composite_curve - using_curves - Composite_curve_segment

# Bounded_curve - Closed_composite_curve - Boundary_curve

# Bounded_curve - B_spline_curve - degree - integer
#                                - control_points - Cartesian_point
#                                - closed - boolean

# Bounded_curve - B_spline_curve - Rational_b_spline_curve - weight_values - real

# Bounded_curve - B_spline_curve - Curve_with_explicit_knots - knot_values - Parameter_value
#                                                            - knot_multiplicities - integer

# Bounded_curve - B_spline_curve - Curve_with_implicit_knots - knot_type - Knot_specification

# Trimmed_curve - start_point - Cartesian_point
#               - end_point - Cartesian_point


from stepcode.Part21 import Parser
import os, glob, re
#import verb
from geomdl import BSpline
from geomdl import NURBS
#from geomdl.visualization import VisMPL
from geomdl import exchange
from geomdl import operations
from geomdl.visualization import VisMPL
from geomdl.elements import Vertex, Triangle
from geomdl import tessellate
from geomdl import operations as op
from geomdl import utilities
import copy

from collections import namedtuple  # PyCharm bug will flag this as type error

Point = namedtuple("Point", "x y z")

from scipy import spatial # numpy version < 1.24
import numpy as np
import pymesh

# eps machine precision constant
eps = np.finfo(float).eps

eps_STEP_AP21 = 1e-6  # STEP precision seems to end here

# need a precision factor that reflects increasing accuracy arising from iterative calculations of centroid --------------------------------------<<<<<<<<<<<<<<<<<
# also scale factor affecting numerical precision
# also STEP seems to be E-6 decimal places at best
iterativeScaledPrecision = 1e-4  # for testing

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


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)

# does one find a centroid through axis-intersection?
# a number of axes that pass through a common point?
# equidistant corners?
def medianPoint(pointArray):
    """
    Returns centroid determined as the median of cartesian values in an input array of points

    :param pointArray: array of Cartesian points
    :return: point
    """

    # from statistics import median # apparently slow/accurate
    # from numpy import median
    if len(pointArray) < 2:
        raise RuntimeError("medianPoint() malformed input")

    xRange = [xyz[0] for xyz in pointArray]
    xCentroid = max(xRange) - ((max(xRange) - min(xRange)) / 2)
    yRange = [xyz[1] for xyz in pointArray]
    yCentroid = max(yRange) - ((max(yRange) - min(yRange)) / 2)
    zRange = [xyz[2] for xyz in pointArray]
    zCentroid = max(zRange) - ((max(zRange) - min(zRange)) / 2)
    return np.array([xCentroid, yCentroid, zCentroid])


def pointDisp(p1, p2):
    """
    Return displacement(s) between input points or list of points

    :param p1: single point
    :param p2: single or list of points
    :return: disp value or list of values
    """

    if not isinstance(p2, tuple):
        raise RuntimeError("pointDisp second variable must be single point")
    if isinstance(p1, tuple):
        # noinspection PyUnresolvedReferences
        return sqrt((p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2)
    elif isinstance(p1, list):
        return [
            np.sqrt((pp.x - p2.x) ** 2 + (pp.y - p2.y) ** 2 + (pp.z - p2.z) ** 2)
            for pp in p1
        ]
    else:
        raise RuntimeError("pointDisp error")


def pointInArc(p, v1, v2, _refDir, _auxDir, _normDir, aCentreP, rotSym=True):
    # (note: if VOR, reverse v1, v2 entries)
    # test if points are on arc between vertex1 & vertex2
    # convert aCentreP, p, v1 & v2 to local coordinates with arcAxisDir as z
    # p may be point or list, returns boolean
    # if not rotSym, no checks for equal radius (e.g. ellipse)

    rotM = np.array([_refDir, _auxDir, _normDir])

    v1UV = np.matmul(rotM, v1 - aCentreP)
    v1UV_angle = np.arctan2(v1UV[1], v1UV[0])

    v2UV = np.matmul(rotM, v2 - aCentreP)
    v2UV_angle = np.arctan2(v2UV[1], v2UV[0])

    if isinstance(p, list):
        pUV_angle = []
        for pp in p:
            ppUV = np.matmul(rotM, pp - aCentreP)
            ppUV_angle = np.arctan2(ppUV[1], ppUV[0])
            pUV_angle.append(ppUV_angle)
    elif isinstance(p, np.ndarray):
        pUV = np.matmul(rotM, p - aCentreP)
        pUV_angle = np.arctan2(pUV[1], pUV[0])

    # sanity check that Z-values are identical
    if isinstance(p, list):
        zTest = [pp[2] for pp in pUV]
        zTest = zTest + [v1UV[2], v2UV[2]]
    else:
        zTest = [pUV[2], v1UV[2], v2UV[2]]

    # if not np.isclose(max(zTest), min(zTest)):
    planeTestFail = 0
    if max(zTest) - min(zTest) > eps_STEP_AP21:
        planeTestFail += 1
        #print("pointInArc() transform error: " + str(max(zTest) - min(zTest)))
        #raise RuntimeError("error")

    # check if vertex order reversed before calling function
    if isinstance(p, list):
        returnList = []
        for ppUV_angle in pUV_angle:
            if (((ppUV_angle > v1UV_angle) and (ppUV_angle < v2UV_angle)) or
                    (v1UV_angle == v2UV_angle)):
                returnList.append(True)
            else:
                returnList.append(False)

        if any(returnList):
            print("pointInArc() transform error: " + str(max(zTest) - min(zTest)))

        return returnList

    else:
        returnBool = (((pUV_angle >= v1UV_angle) and (pUV_angle <= v2UV_angle)) or
                (v1UV_angle == v2UV_angle))

        if returnBool and not planeTestFail:
            print("pointInArc() transform deviation: " + str(max(zTest) - min(zTest)))

        return returnBool

# def checkFilePath(filePath):
#     if os.path.isfile(filePath) and os.access(filePath, os.R_OK):
#         return filePath
#     else:
#         raise Exception("{0} is not a readable filepath".format(filePath))

primitivesDir = os.path.normpath(
    r"/media/foobert/Dell2HDD/PY_DELL2/SWKS_Rhino_compare_circle/primitives"
)
primitivesDir = os.path.normpath(
    r"/home/foobert/STEP_test_files/primitives"
)
testDir = os.path.normpath(
    r"/home/foobert/STEP_test_files"
)

#filepath = primitivesDir + "/Cube/unit_cube_inc8.0_blend.06.stp"
#filepath = primitivesDir + "/Cylinder/unit_cyl.stp"
#filepath = primitivesDir + "/Primitive_Cone-PartCone.step"

#filepath = r"/home/foobert/Downloads/RR_STEP_test_1.step"
#filepath = r"/home/foobert/Downloads/RR_STEP_test_N.step"
#filepath = r"/home/foobert/Downloads/00000001_1ffb81a71e5b402e966b9341_step_000.step"
#filepath = r"/home/foobert/Downloads/00000010_b4b99d35e04b4277931f9a9c_step_000.step"
# filepath = r"/home/foobert/Downloads/LEA-M8F(AP203).STEP"

filepath = primitivesDir + "/Cube/unit_cube.stp"
filepath = testDir + "/TiltedCylinder2.step"
#filepath = testDir + "/TiltedConeAP203.step"
#filepath = testDir + "/OffsetCone_AP214_noPlacement_noParametric.step"
filepath = testDir + "/OffsetCone-PartCone.step"
#filepath = testDir + "/TiltedCylinder4_AP214_PC.step"
#filepath = testDir + "/Cylinder5_AP214_PC.step"
filepath = testDir + "/Drexel_blended_cylinder.step"
#filepath = testDir + "/Synth_ellipse_plain.step"


try:
    # filepath sanity check
    checkPath = glob.glob(filepath)
    if len(checkPath) != 1:
        raise Exception("pathname failure")
    else:
        with open(filepath) as f:
            STEPdata = f.read()

except:
    raise Exception("file failure")

STEPdata = STEPdata.replace(";\n", ";;")
STEPdata = STEPdata.replace("\n", "")
STEPdata = STEPdata.replace(";;", ";\n")

P = Parser()
# model = P.parse(/home/foobert/.config/JetBrains/PyCharmCE2023.1/scratches/simple_schema.exp)
model = P.parse(STEPdata)

STEP_entities = model.sections[0].entities

ClosedShells = []
FaceInstances = []
for cs in STEP_entities:
    if hasattr(cs, "type_name"):  # stepcode.Part21.ComplexEntity type missing attribute
        if cs.type_name == "CLOSED_SHELL":
            ClosedShells.append(cs)

def ref2index(s):
    return int(s[1:]) - 1

# encountering step files with non-contiguous line and #ref indices
if ref2index(STEP_entities[-1].ref) is not len(STEP_entities):
    STEP_entities.sort(key=lambda x: int(x.ref[1:]))
    padded_STEP_entities = [0] * int(STEP_entities[-1].ref[1:])
    for s in STEP_entities:
        padded_STEP_entities[ref2index(s.ref)] = s
    STEP_entities = padded_STEP_entities

# for cs in ClosedShells:
#     FaceRefs = [ref2index(r) for r in cs.params[-1]]  # get list of face refs
#     FaceRefs = [
#         fr
#         for fr in FaceRefs
#         if STEP_entities[fr].type_name in ["ADVANCED_FACE", "FACE_SURFACE"]
#     ]
#     # 'CARTESIAN_POINT', 'DIRECTION' # 'MANIFOLD_SOLID_BREP'
#
#     # recursively extract face edge points
#     # NURBS, conics, circles, holes need alternative object recording centres, directions, knots, etc
#     # for fr in FaceRefs:

# first pass to find origin point
# find all CARTESIAN_POINT values child to VERTEX_POINT


def pointProjectAxis(p, axP, axN):
    '''
    project point p to axis defined by point axP and normal axN
    p: point to project
    axP: axis point
    axN: axis normal vector
    '''
    vpp = p - axP
    vpn = axP - axN
    return vpn * np.dot(vpp, vpn) / np.dot(vpn, vpn)


def pointCircleMinMaxDisp(p, centrePoint, normAxis, radius, interior=False):
    # Point on circle (including interior points) furthest/closest to point "p"

    # normAxis = auxDir / np.linalg.norm(normAxis)  # unit normal vector

    # distance d of point p to the plane is dot product of plane unit normal, normAxis
    # with the difference vector from the point in the plane P and S: d = (S − P) · u
    # distance is signed: positive when S is on the side of the plane where u faces and
    # negative on the other side. zero, S is in the plane
    # d = (p - plane[0]) * u

    d = np.dot((p - centrePoint), normAxis)

    # point S′ is S projected to the plane, by subtracting d · u from S:
    # S′ = S − d · u = S − ((S − P) · u) · u

    pUV = p - np.dot(d, normAxis)  # point projected to circle along orthogonal
    ppCentreDisp = np.linalg.norm(pUV - centrePoint)

    # if point pUV, projected to circle plane, falls outside circle radius, replace with edge point.

    #centreNormal = False
    if np.isclose(ppCentreDisp, 0):
        # p lies on normal through circle centrePoint,
        return centrePoint, None

    pUVdir = (pUV - centrePoint)/ppCentreDisp

    nUVmin = centrePoint + radius * pUVdir
    nUVmax = centrePoint - radius * pUVdir

    #nUVmin = centrePoint + radius / (ppCentreDisp * (pUV - centrePoint))
    #nUVmax = centrePoint - radius / (ppCentreDisp * (pUV - centrePoint))

    if ppCentreDisp >= radius and not interior:
        return nUVmin, nUVmax
    else:
        # if ppCentreDisp < eps_STEP_AP21:
        #     return centrePoint, centrePoint
        # else:
        centreOffset = pUV - centrePoint
        centreOffset = centreOffset/np.linalg.norm(centreOffset)
        nUVmax = centrePoint - radius * centreOffset
        #
        nUVmin = pUV

    return nUVmin, nUVmax,


def pointEllipseMinMaxDisp(
                            p,
                            eCentre,
                            eLocalXaxis,
                            eLocalYaxis,
                            eNormalAxis,
                            eMajorRadius,
                            eMinorRadius,
                            interior=False
):
    # Robert Nurnberg method, see http://wwwf.imperial.ac.uk/~rn/distance2ellipse.pdf
    # Move origin to ellipse centre, and rotate point coordinate to local coordinates
    # based on ellipse major & minor axis.
    # Assuming direction_ratios of AXIS2_PLACEMENT_3D are normalised

    # eFeatureAxis = arrayTypeCast(eFeatureAxis)
    # eLocalXYaxis = arrayTypeCast(eLocalXYaxis)
    # eCentre = arrayTypeCast(eCentre)

    # cLocalYaxis = np.cross(eLocalXaxis - eCentre, eNormalAxis - eCentre) # STEP AP21 gives axes wrt local origin
    # eLocalYaxis = np.cross(eLocalXaxis, eNormalAxis)
    # eLocalYaxis = eLocalYaxis / np.linalg.norm(eLocalYaxis)  # unit normal vector

    rotM = np.array([eLocalXaxis, eLocalYaxis, eNormalAxis])
    pUV = np.matmul(rotM, p - eCentre)

    # project to ellipse plane
    pUV = np.array([pUV[0], pUV[1], 0])

    theta = np.arctan2(eMajorRadius * pUV[1], eMinorRadius * pUV[0])
    i = 0
    n0 = pUV
    radiusConst = eMajorRadius**2 - eMinorRadius**2
    pUVdisp = np.linalg.norm(pUV)

    while i < 100:
        i += 1
        f = (
            radiusConst * np.cos(theta) * np.sin(theta)
            - pUV[0] * eMajorRadius * np.sin(theta)
            + pUV[1] * eMinorRadius * np.cos(theta)
        )

        f_ = (
            radiusConst * (np.cos(theta) ** 2 - np.sin(theta) ** 2)
            - pUV[0] * eMajorRadius * np.cos(theta)
            - pUV[1] * eMinorRadius * np.sin(theta)
        )

        theta = theta - f / f_

        n = np.array([eMajorRadius * np.cos(theta), eMinorRadius * np.sin(theta), 0])

        if np.allclose(n0, n):  # tolerance = 1E-12
            if (
                np.linalg.norm(n) > pUVdisp
            ) and not interior:  # test if p < n, inside ellipse
                nUVmin = np.matmul(rotM.T, pUV) + eCentre
            else:
                nUVmin = np.matmul(rotM.T, n) + eCentre

            # furthest point on ellipse must be edge point on line that contains ellipse centre and closest point
            nUVmax = np.matmul(rotM.T, np.array([-n[0], -n[1], 0])) + eCentre

        n0 = n
    return nUVmin, nUVmax



def splineCurveMinMaxPointDisp(curve, p, maxSearch=True): #TODO: find multiple maxima, verb_eps_constant
    # based on Nurbs Book, Piegl & Tiller p.230
    # same idea to find orthogonal tangent, but at maximum/minimum displacement from a centroid point
    # revised to test minima segment candidates using 3-point defined circle parameters
    def radiusCentre3points2(p1, p2, p3):
        # radius + centre from 3 point circle

        t = p2 - p1
        u = p3 - p1
        v = p3 - p2

        w = np.cross(t, u)  # triangle normal
        wsl = np.linalg.norm(w)
        if (wsl < 10e-14): return False  # triangle area too small (additionally check points for colinearity)

        wsl = np.dot(w, w)
        iwsl2 = 1. / (2. * wsl)
        tt = np.dot(t, t)
        uu = np.dot(u, u)

        circCenter = p1 + (u * tt * (np.dot(u, v)) - t * uu * (np.dot(t, v))) * iwsl2
        circRadius = np.sqrt(tt * uu * (np.dot(v, v)) * iwsl2 * 0.5)
        # circAxis   = w / np.sqrt(wsl)

        return circRadius, circCenter

    d_min = np.inf
    u = 0.0
    curveSampleFactor = 2

    numSamples = curve.ctrlpts_size * curve.degree * curveSampleFactor
    span = (curve.knotvector[-1] - curve.knotvector[0]) / (numSamples - 1)
    kvs = np.array([curve.knotvector[0] + (span * i) for i in range(0, numSamples)])
    pts = curve.evaluate_list(kvs)

    # curvature approach fails without C2 continuity at minima/maxima?
    # not much of an issue with STEP spline limitations
    for i in range(1, len(pts) - 2):
        u1 = kvs[i]
        p0 = np.array(pts[i - 1])
        p1 = np.array(pts[i])
        p2 = np.array(pts[i + 1])

        dcc, cc = radiusCentre3points2(p0, p1, p2)

        if (cc == np.inf):
            # curvature method failure, use displacement minima , dcc_p = np.inf
            print("spline discrete point 3-point circle curvature failure ")
            du_p = [np.linalg.norm(p - i) for i in pts]
            if maxSearch:
                #u = pts[du_p.index(max(du_p))][0]
                u = kvs[np.where(du_p == max(du_p))[0][0]]
            else:
                #u = pts[du_p.index(min(du_p))][0]
                u = kvs[np.where(du_p == min(du_p))[0][0]]
            break
            # du_p = [du0_p, du1_p, du2_p]
            # if (max(du_p) < d_min):
            #     d_min = max(du_p)
            #     offset = du_p.index(min(du_p)) - 1
            #     u = pts[i+offset][0]

        dcc_p = np.linalg.norm(p - cc) # distance circle centre

        if maxSearch:
            dua = max([np.linalg.norm(p - i) for i in [p0, p1, p2]])
            # centre of 3 point circle is nearer or further than max/min sample point from p
            #orthoSign = dua > (dcc_p + (2*dcc)) # circle centre closer to p than circle points => maximum curve
            if (dcc < d_min) and (dua > dcc_p):
                d_min = dcc # find the smallest radius
                u = u1
        else: #minSearch
            dua = min([np.linalg.norm(p - i) for i in [p0, p1, p2]])
            # centre of 3 point circle is nearer or further than max/min sample point from p
            #orthoSign = dua > (dcc_p - (2*dcc)) # circle centre further from p than circle points => minimum curve
            if (dcc < d_min) and (dua > dcc_p):
                d_min = dcc
                u = u1

        #print("dua: " + repr(dua) + "    dcc: " + repr(dcc))

    eps1 = 0.0001   # Euclidean distance measure
    eps2 = 0.0005   # zero cosine measure
    verb_EPSILON = 1E-10
    minu = curve.knotvector[0]
    maxu = curve.knotvector[-1]

    curveClosed = np.linalg.norm(np.array(curve.ctrlpts[0]) - np.array(curve.ctrlpts[-1]))
    closed = False
    if (curveClosed < verb_EPSILON):
        if (curveClosed*curveClosed < verb_EPSILON):
            closed = True
    cu = u

    # def n(u2,e1,d):
    #     #   Newton's method: 	 u* = u - f / f'
    #     #   use product rule to form derivative, f':   f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
    #     #   d:  ( C(u) - p )
    #
    #     #f1 = verb.verb_core_Vec.dot(e1[1], d)       # C'(u) * ( C(u) - p )
    #     #s0 = verb.verb_core_Vec.dot(e1[2], d)       # C"(u) * ( C(u) - p )
    #     #s1 = verb.verb_core_Vec.dot(e1[1], e1[1])   # C'(u) * C'(u)
    #     f1 = np.dot(e1[1], d)       # C'(u) * ( C(u) - p )
    #     s0 = np.dot(e1[2], d)       # C"(u) * ( C(u) - p )
    #     s1 = np.dot(e1[1], e1[1])   # C'(u) * C'(u)
    #     df = s0 + s1
    #     return u2 - f1 / df
    #     #return (u2 - np.dot(e1[1], d)) / (np.dot(e1[2], d) + np.dot(e1[1], e1[1]))

    def n(u2, e1, d):
        #   Newton's method: 	 u* = u - f / f'
        #   use product rule to form derivative, f':   f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
        #   d:  ( C(u) - p )
        return u2 - (np.dot(e1[1], d) / (np.dot(e1[2], d) + np.dot(e1[1], e1[1])))

    maxits = 5
    i = 0
    while(i < maxits):
        # Newton iteration here is to find orthogonal tangent and is max/min agnostic
        e = np.array(curve.derivatives(cu, order=2) )   #f(cu)
        dif = e[0] - p                                  #   C(u) - p
        c1v = np.linalg.norm(dif)                       #   |C(u) - p|
        c2n = np.dot(e[1], dif)                         #   C'(u) * (C(u) - P)
        c2d = np.linalg.norm(e[1]) * c1v                #   |C'(u)||C(u) - P|
        c2v = c2n / c2d                                 #   |C'(u) * (C(u) - P)| / |C'(u)||C(u) - P|
        #c2v = np.dot(e[1], dif) / (np.linalg.norm(e[1]) * np.linalg.norm(dif))

        if (c1v < eps1) and (np.abs(c2v) < eps2):
            return cu
        #ct = n(cu, e, dif)                           #   u* = u - f / f'
        ct = cu - (np.dot(e[1], dif) / (np.dot(e[2], dif) + np.dot(e[1], e[1])))

        if ct < minu:
            if closed:
                #ct = maxu - (ct - minu) # NURBS book
                ct = maxu - (minu - ct)
            else:
                ct = minu
        elif ct > maxu:
            if closed:
                ct = minu + (ct - maxu)
            else:
                ct = maxu

        c3v = np.linalg.norm(np.multiply(ct - cu, e[1]))
        if c3v < eps1:
            return cu

        cu = ct
        i += 1

    return cu

# def pointNURBSmaxDisp(p, curve):
#     # based on verb & Nurbs Book, Piegl & Tiller p.230
#     # Orthogonal tangents work for both maxima and minima =>
#     # Newton-Raphson will work if correct subdivision of curve selected.
#     # same idea to find orthogonal tangent, but at maximum displacement from a centroid point
#     # revised to test minima segment candidates using 3-point defined circle parameters
#
#     def radiusCentre3points(p1, p2, p3):
#         # radius + centre from 3 point circle
#
#         p1 = np.array([p1[0], p1[1], p1[2]])
#         p2 = np.array([p2[0], p2[1], p2[2]])
#         p3 = np.array([p3[0], p3[1], p3[2]])
#
#         t = p2 - p1
#         u = p3 - p1
#         v = p3 - p2
#
#         w = np.cross(t, u)  # triangle normal
#         wsl = np.linalg.norm(w)
#
#         # triangle area too small (additionally check points for colinearity)
#         if wsl < 10e-14:
#             return np.inf, np.inf
#
#         wsl = np.dot(w, w)
#         iwsl2 = 1.0 / (2.0 * wsl)
#         tt = np.dot(t, t)
#         uu = np.dot(u, u)
#
#         circCenter = p1 + (u * tt * (np.dot(u, v)) - t * uu * (np.dot(t, v))) * iwsl2
#         circRadius = np.sqrt(tt * uu * (np.dot(v, v)) * iwsl2 * 0.5)
#         # circAxis   = w / np.sqrt(wsl)
#
#         #print(circCenter)
#
#         return circRadius, [circCenter[0], circCenter[1], circCenter[2]]
#
#     d_min = np.inf  # d_max = 0. for finding minima
#     u = 0.0
#     pts = verb.verb_eval_Tess.rationalCurveRegularSample(
#         curve, (len(curve.controlPoints) * curve.degree * 2), True
#     )
#
#     for i in range(1, len(pts) - 2):
#         u1 = pts[i][0]
#         p0 = pts[i - 1][1:]
#         p1 = pts[i][1:]
#         p2 = pts[i + 1][1:]
#         dcc, cc = radiusCentre3points(p0, p1, p2)
#         du0_p = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, p0))
#         du1_p = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, p1))
#         du2_p = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, p2))
#         if (cc == np.inf):
#             # curvature method failure, use displacement minima , dcc_p = np.inf
#             du_p = [verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, i[1:])) for i in pts]
#             u = pts[du_p.index(max(du_p))][0]
#             print("spline discrete point 3-point circle curvature failure ")
#             break
#             # du_p = [du0_p, du1_p, du2_p]
#             # if (max(du_p) < d_min):
#             #     d_min = max(du_p)
#             #     offset = du_p.index(min(du_p)) - 1
#             #     u = pts[i+offset][0]
#         else:
#             dcc_p = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, cc))
#             dua = (du0_p + du1_p + du2_p) / 3
#             # centre of 3 point circle is nearer or further than sample point barycentre from p
#             orthoSign = dua > dcc_p  # closer means maximum curve
#             # print("dua: "+repr(dua)+"    dcc: "+repr(dcc))
#             if (dcc < d_min) and orthoSign:  # dcc > d_max, initialised as 0.0 for minima
#                 d_min = dcc
#                 u = u1
#
#     #   solve:  C'(u) * ( C(u) - P ) = 0 = f(u)    C(u) is the curve, p is the point, * is a dot product
#     #   Newton-Raphson method:  u* = u - f / f'
#     #   derivative via product rule, f':   	f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
#     #   Piegl & Tiller suggested halting criteria
#     #
#     #    |C(u) - p| < e1                Euclidean distance measure
#     #
#     #    |C'(u)*(C(u) - P)|
#     #    ------------------  < e2       zero cosine measure
#     #    |C'(u)| |C(u) - P|
#     #
#     #     1) check 2 & 3
#     #     2) if at least one of these is not, compute new value, otherwise halt
#     #     3) ensure the parameter stays within range
#     #    			* if not closed, don't allow outside of range a-b
#     #    			* if closed (e.g. circle), allow to move back to beginning
#     #     4)  if |(u* - u)C'(u)| < e1, halt
#
#     eps1 = 0.0001  # Euclidean distance measure
#     eps2 = 0.0005  # zero cosine measure
#     minu = curve.knots[0]
#     maxu = verb.verb_core_ArrayExtensions.last(curve.knots)
#     closed = (
#         verb.verb_core_Vec.normSquared(
#             verb.verb_core_Vec.sub(
#                 curve.controlPoints[0],
#                 verb.verb_core_ArrayExtensions.last(curve.controlPoints),
#             )
#         )
#         < verb.verb_core_Constants.EPSILON
#     )
#     cu = u
#
#     def f(u1):
#         return verb.verb_eval_Eval.rationalCurveDerivatives(curve, u1, 2)
#
#     def n(u2, e1, d):
#         #   Newton's method: 	 u* = u - f / f'
#         #   use product rule to form derivative, f':   f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
#         #   d:  ( C(u) - p )
#
#         f1 = verb.verb_core_Vec.dot(e1[1], d)  # C'(u) * ( C(u) - p )
#         s0 = verb.verb_core_Vec.dot(e1[2], d)  # C"(u) * ( C(u) - p )
#         s1 = verb.verb_core_Vec.dot(e1[1], e1[1])  # C'(u) * C'(u)
#         df = s0 + s1
#         return u2 - f1 / df
#
#     maxits = 5
#     i = 0
#     while i < maxits:
#         e = f(cu)
#         dif = verb.verb_core_Vec.sub(e[0], p)  #   C(u) - p
#         c1v = verb.verb_core_Vec.norm(dif)  #   |C(u) - p|
#         c2n = verb.verb_core_Vec.dot(e[1], dif)  #   C'(u) * (C(u) - P)
#         c2d = verb.verb_core_Vec.norm(e[1]) * c1v  #   |C'(u)||C(u) - P|
#         c2v = c2n / c2d  #   |C'(u) * (C(u) - P)| / |C'(u)||C(u) - P|
#         c1 = c1v < eps1
#         c2 = np.abs(c2v) < eps2
#         if c1 and c2:
#             return cu
#         ct = n(cu, e, dif)  #   u* = u - f / f'
#
#         if ct < minu:
#             if closed:
#                 # ct = maxu - (ct - minu) # NURBS book
#                 ct = maxu - (minu - ct)
#             else:
#                 ct = minu
#         elif ct > maxu:
#             if closed:
#                 ct = minu + (ct - maxu)
#             else:
#                 ct = maxu
#
#         c3v = verb.verb_core_Vec.norm(verb.verb_core_Vec.mul(ct - cu, e[1]))
#         if c3v < eps1:
#             return cu
#         # print(ct)
#         cu = ct
#         i += 1
#
#     return cu  # np.array([cu[0], cu[1, cu[2])


def intersectSegmentPlane(
    vertex0, vertex1, planeNorm, planePoint, precision=eps_STEP_AP21
):
    """
    Return intersection point of segment and plane or None if no intersection.

    :param vertex0: segment end point
    :param vertex1: second segment end point
    :param planePoint: point on plane
    :param planeNorm: normal vector defining plane direction (does not need to be normalized)
    :param precision: constant for STEP numerical float comparison
    :return: point of segment plane intersection or none
    """
    s = vertex1 - vertex0
    SprojP = np.dot(planeNorm, s)
    if abs(SprojP) > precision:
        F = -np.dot(planeNorm, (vertex0 - planePoint)) / SprojP
        if 0.0 <= F <= 1.0:
            return vertex0 + np.dot(s, F)  # point intersects segment
        else:
            return None  # intersection beyond segment
    return None  #   segment is parallel to plane

def intersectPlaneCircle(pPoint, pNorm, cPoint, cNorm, cRadius, tol=eps_STEP_AP21):
    """Cintersection of a plane and a circle.
    plane point
    plane norm
    circle point
    circle norm
    circle radius

    There are 4 cases of plane-circle intersection:

    1. intersect in 2 points (secant),
    2. intersect in 1 point (tangent),
    3. do not intersect, or
    4. coincide (circle.plane == plane).
    """
    if np.abs(np.dot(pNorm, cNorm)) >= 1 - tol:
        return None

    d = np.cross(pNorm, cNorm)  # direction of intersection line
    # vector in plane 1 perpendicular to the direction of the intersection line
    v1 = np.cross(d, pNorm)
    p1 = pPoint + v1    # point on plane 1
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
        return [ipt1, ipt2]


def intersectArcPlane( planeNormal,
                        planePoint,
                        arcCentrePoint,
                        arcRefDir,
                        arcAuxDir,
                        arcNormal,
                        arcRadius,
                        arcVertex1,
                        arcVertex2 ):

    # planePoint = centroidProjAxisPoint for most calls

    # points within the edge
    validatedPoints = []

    # determine if circle and plane intersect, if so return intersection line/segment (cosine projection?)
    circleIntersectPoints = intersectPlaneCircle(planePoint, planeNormal, arcCentrePoint, arcNormal, arcRadius)
    if circleIntersectPoints is None:
        return []
    if len(circleIntersectPoints) == 0:
        return [] # no intersection
    if len(circleIntersectPoints) == 1:
        print("tangent exception")
        pass
    # check projected segment is between arc vertices,
    for cip in circleIntersectPoints:

        if pointInArc(cip, arcVertex1, arcVertex2, arcRefDir, arcAuxDir, arcNormal, arcCentrePoint, rotSym=True):
            validatedPoints.append(cip)

    return validatedPoints


def pointSegmentMinDisp(p, s0, s1):
    #   https://web.archive.org/web/20220121145748/http://geomalgorithms.com/index.html
    # final two boolean values indicate whether point is at or beyond endpoints of segments
    v = s1 - s0
    w = p - s0

    c1 = np.dot(w, v)
    if ( c1 <= 0 ):
      return s0, np.linalg.norm(w), True, False

    c2 = np.dot(v, v)
    if ( c2 <= c1 ):
      return s1, np.linalg.norm(p - s1), False, True

    b = c1 / c2
    pb = s0 + b * v
    return pb, np.linalg.norm(p - pb), False, False


# def pointSegmentMinDisp(p, a, b):
#     # https://web.archive.org/web/20220121145748/http://geomalgorithms.com/index.html
#     s = b - a
#     w = p - a
#     ps = np.dot(w, s)
#     if ps <= 0:
#         return a, np.linalg.norm(w), True
#     l2 = np.dot(s, s)
#     if ps >= l2:
#         closest = b
#         vertexEnd = True
#     else:
#         closest = a + ps / l2 * s
#         vertexEnd = False
#     return closest, np.linalg.norm(p - closest), vertexEnd


def planeMinMaxPoint(p, V, minmax="min"):
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
    vertices = []  # np.array([])
    #p = arrayTypeCast(p)
    for v in V:
        #v = arrayTypeCast(v)
        vertices.append(v)
        if minmax == "max":
            d = np.linalg.norm(v - p)
            disp.append(d)

    if minmax == "max":
        maxDisp = max(disp)
        maxPoint = vertices[disp.index(maxDisp)]
        return maxPoint, maxDisp

    tri = pymesh.triangle()
    tri.points = vertices  # np.array(vertices)
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
    return minPoint, minDisp


def cleanSubRefs(refStr):
    """
    flatten parameter sets and extract #d string references
    """
    flattenList = (
        lambda irregularL: [e for i in irregularL for e in flattenList(i)]
        if type(irregularL) is list
        else [irregularL]
    )
    sp = flattenList(refStr)  # flatten sublists
    sp = [i for i in sp if isinstance(i, str)]
    sp = [re.match("#[0-9]+", i) for i in sp]
    sp = [i.string for i in sp if i is not None]
    return sp


def curveEnclosingRectangle(v1, v2, centreP, normV, radius):
    # in order to establish union of partial ellipse/circle, create box encompassing circle/ellipse
    # in the case of ellipse, use eMajorRadius
    # vertex midpoint
    midP = (v2 + v1)/2
    midPv = np.cross(v2 - v1, normV)
    # midpoint - centroid vector, falls apart if coincident with centre point
    #midPv = midP - centreP
    midPv = midPv / np.linalg.norm(midPv)
    farP = centreP + midPv * radius
    orthoPv = np.cross(normV, midPv)
    orthoPv = orthoPv / np.linalg.norm(orthoPv)
    p1 = midP + radius * orthoPv
    p2 = midP - radius * orthoPv
    p3 = farP - radius * orthoPv
    p4 = farP + radius * orthoPv
    return [p1, p2, p3, p4]


def CP2point(STEPobject, pointRefString):
    """
    return np.array of STEP AP21 CARTESIAN_POINT, DIRECTION, VECTOR
    """
    # assert STEPobject[ref2index(pointRefString)].type_name == "CARTESIAN_POINT"
    #R3 = STEPobject[ref2index(pointRefString)].params[1]
    R3 = STEPobject[ref2index(pointRefString)].params[-1]
    return np.array([R3[0], R3[1], R3[2]])


def axis2Placement3D(strRef, STEPobj):
    # process STEP AXIS2_PLACEMENT_3D entity
    assert STEPobj[ref2index(strRef)].type_name == "AXIS2_PLACEMENT_3D"
    # D = {}

    if '#' in STEPobj[ref2index(strRef)].params[0]:  # check for name string
        # not a robust approach as certain values are optional
        offset = 0
    else:
        offset = 1

    axisPoint = STEPobj[ref2index(strRef)].params[offset]
    assert STEPobj[ref2index(axisPoint)].type_name == "CARTESIAN_POINT"
    axisPoint = CP2point(STEPobj, axisPoint)
    # D["axisPoint"] = axisPoint

    auxDir = STEP_entities[ref2index(strRef)].params[offset+1]
    assert STEP_entities[ref2index(auxDir)].type_name == "DIRECTION"
    auxDir = CP2point(STEP_entities, auxDir)
    # D["auxDir"] = auxDir

    refDir = STEP_entities[ref2index(strRef)].params[offset+2]
    assert STEP_entities[ref2index(refDir)].type_name == "DIRECTION"
    refDir = CP2point(STEP_entities, refDir)
    # D["refDir"] = refDir

    return axisPoint, auxDir, refDir


def axis2Placement3D_2(strRef, STEPobj):
    # process STEP AXIS2_PLACEMENT_3D entity
    # axis: the Direction that defines the second axis of the Axis_placement. (Y or V, not normal)
    # The value of this attribute need not be specified.
    # ref_direction: the direction used to determine the direction of the local X axis. (or U)
    # The value of this attribute need not be specified.
    # If axis or ref_direction is omitted, these directions are taken from the geometric coordinate system
    # If both axis and ref_direction are provided then the vector product of axis and ref_direction shall not be a null vector.

    if STEPobj[ref2index(strRef)].type_name != "AXIS2_PLACEMENT_3D":
        print("axis2Placement3D_2 assignment failure: " + STEPobj[ref2index(strRef)].type_name)
        return None

    axisPoint = None
    normDir = None
    refDir = None
    subRefList = cleanSubRefs(STEPobj[ref2index(strRef)].params)

    if (STEPobj[ref2index(subRefList[0])].type_name == "CARTESIAN_POINT"):
        axisPoint = CP2point(STEPobj, subRefList[0])

        if len(subRefList) > 1:
            if (STEPobj[ref2index(subRefList[1])].type_name == "DIRECTION"):
                normDir = CP2point(STEPobj, subRefList[1])

                # relying on order specified
                # https://www.steptools.com/stds/smrl/data/modules/elemental_geometric_shape/sys/4_info_reqs.htm#elemental_geometric_shape_arm.axis_placement
                if len(subRefList) > 2:
                    if (STEPobj[ref2index(subRefList[2])].type_name == "DIRECTION"):
                        refDir = CP2point(STEPobj, subRefList[2])

    return axisPoint, normDir, refDir

# # test to determine all 'EDGE_CURVE' are children of 'ADVANCED_FACE'
# for r in STEP_entities:
#     if hasattr(r, "type_name"):
#         if r.type_name == "EDGE_LOOP":
#             rp = cleanSubRefs(r.params)  # extract all '#d'
#             for p in rp:
#                 if hasattr(STEP_entities[ref2index(p)], "type_name"):
#                     if STEP_entities[ref2index(p)].type_name != "ORIENTED_EDGE":
#                         ref = STEP_entities[ref2index(p)].type_name
#                         print(
#                             "ORIENTED_EDGE child of EDGE_LOOP assumption failure: "
#                             + ref
#                         )

# 189 = ADVANCED_FACE( '', ( #409, #410 ), #411, .T. )
# 409 = FACE_OUTER_BOUND( '', #1293, .T. );
# 1293 = EDGE_LOOP( '', ( #2202 ) );
# 2202 = ORIENTED_EDGE( '', *, *, #2631, .F. );
# 2631 = EDGE_CURVE( '', #2959, #2959, #2960, .F. );

# The geometry used in the definition of the face shall be restricted.
# The face geometry shall be an elementary surface, swept surface, or b-spline surface;
# The geometry of all bounding edges of the face shall be fully defined as edge curves;
# The types of curve used to define the geometry of edges shall be restricted to lines,
# conics, polylines, surface curves, or b spline curves;
# All vertices used in the face definition shall be of type vertex point with geometry defined by a cartesian point;
# The use of oriented paths in the definition of the edge loops of the advanced face is prohibited;
# If the face geometry is a swept surface then the swept curve used in the definition shall be a line,
# conic, polyline, or a b-spline curve;
# For any vertex loop used to bound the face, the loop vertex shall be a
# vertex whose geometry shall be defined by a cartesian point;

# The face bounds shall be defined by either edge loops or vertex loops;

# If a surface curve is used as part of a face bound then the associated geometry attribute shall
# reference pcurves not surfaces;
# If a polyline is used either to define a swept surface or as part of a face bound,
# it shall contain at least three points;
# Any instance of advanced face which has the geometry of a complete spherical surface shall be bounded by a vertex loop
# located at the point where the z axis of the placement coordinate system leaves the surface of the sphere,
# (that is, at the North Pole of the sphere).

def insideOutsideSurfaceTest(localCentroid, mPoint, AFSdata):
    # project ray from localCentroid toward mPoint in order to determine whether this ray intersects the surface defined in AFS
    # the angle of edge vertices of an intercepted face are calculated with respect to this point and summed
    # return true if point within surface edges
    sNormDen = np.linalg.norm(mPoint - localCentroid)
    if sNormDen > eps_STEP_AP21:
        sNorm = (mPoint - localCentroid) / sNormDen
    else:
        # point is at localCentroid, find normal from surface
        print("insideOutsideSurfaceTest() fail, median localCentroid at surface, not solved")
        _1=1
    edgeTotalAngle = 0
    vertexLoopCount = 0
    for edge in AFSdata["FaceOuterBoundEdgeLoopList"]:
        v1 = edge["vertex1"]
        v2 = edge["vertex2"]

        #if not np.isclose(v1, v2).all():
        if edge["vertex1ref"] == edge["vertex2ref"]: # loop, discard
            if edge['typeName'] == 'CIRCLE':
                vertexLoopCount +=1
        else:
            # project vertex points to a plane defined by mPoint and sNorm
            m_v1 = v1 - mPoint
            disp = np.dot(m_v1, sNorm) # scalar distance from point to plane along the normal
            p_v1 = v1 - disp * sNorm

            m_v2 = v2 - mPoint
            disp = np.dot(m_v2, sNorm)
            p_v2 = v2 - disp * sNorm

            # cosine of the angle
            v1v2cos = np.dot(p_v1, p_v2) / np.linalg.norm(p_v1 - mPoint) / np.linalg.norm(p_v2 - mPoint)
            edgeTotalAngle += np.arccos(np.clip(v1v2cos, -1, 1))
            # zero value implies mPoint lies on edge

    # # find matching inverted edges
    # matchingEdgeCount = 0
    # for edgeSource in AFSdata['FaceOuterBoundEdgeLoopList']:
    #     v1ref = edgeSource['vertex1ref']
    #     v2ref = edgeSource['vertex2ref']
    #     if (v1ref != v2ref):
    #         for edgeTarget in AFSdata['FaceOuterBoundEdgeLoopList']:
    #                 if (v2ref == edgeTarget['vertex1ref']) and (v1ref == edgeTarget['vertex2ref']) or \
    #                         (v2ref == edgeTarget['vertex2ref']) and (v1ref == edgeTarget['vertex1ref']):
    #                     if edgeSource["typeName"] == edgeTarget["typeName"]:
    #                         if edgeSource.get("subTypeName") is not None:
    #                             if edgeTarget.get("subTypeName") is not None:
    #                                 if edgeSource["subTypeName"] == edgeTarget["subTypeName"]:
    #                                     matchingEdgeCount +=1
    #                         else:
    #                             matchingEdgeCount +=1

    # surfaces that curve > 90deg can give false negative, assume inside/outside test is to detect nearest surface

    # exception for cylinder, conic, sphere
    if AFS['SurfaceTypeName'] in ["CYLINDRICAL_SURFACE",]:
        # edges angle sum to zero from 2x loops and seam edge ---------------------needs testing
        if (np.abs(edgeTotalAngle) < eps_STEP_AP21) and (vertexLoopCount == 2): # (len(AFSdata["FaceOuterBoundEdgeLoopList"]) < 4):
            # vertex angles that sum to zero will still indicate a crossing, for either CW, or ACW direction configuration
            # say like point is on the LHS of a clockwise loop?, or whatever
            # direction of edge closest to surface point?

            # if method is uniquely used for closest surface point on a plane encompassing projection and cylinder axis
            # can then ignore far-side intersection
            return True
        else:
            return False
    if AFS['SurfaceTypeName'] in ["CONICAL_SURFACE",]:
        # edges angle sum to zero from 2x loops and seam edge ---------------------needs testing
        if (np.abs(edgeTotalAngle) < eps_STEP_AP21) and (vertexLoopCount == 2): # (len(AFSdata["FaceOuterBoundEdgeLoopList"]) < 4):
            return True
        else:
            return False
    if AFS['SurfaceTypeName'] in ["PLANE", ]: # full circle
        if (np.abs(edgeTotalAngle) < eps_STEP_AP21) and (len(AFSdata["FaceOuterBoundEdgeLoopList"]) < 2):
            return True
        else:
            return False

    # may have to take VOR into account depending on orientation of surface,
    # e.g. concave or convex spherical surface
    # assume that minima point is always closest to centroid?
    #if not edge["VOR"]:
    if not AFSdata["SurfaceNormalOutwards"]:
        if (edgeTotalAngle - (2 * np.pi)) < eps_STEP_AP21:
            return True
        else:
            return False
    else:
        if (edgeTotalAngle + (2 * np.pi)) < eps_STEP_AP21:
            return True
        else:
            return False


def lineParse(localCentroid, parsedEdge_, STEP_entities_):
    # ENTITY Line
    #   SUBTYPE OF (Curve);
    #   point : Cartesian_point;
    #   line_direction : Direction;

    # if parsedEdge_.get("subEdgeType") is not None:
    #     edgeType_ = parsedEdge_["subEdgeType"]
    # else:
    edgeType_ = parsedEdge_["edgeType"]

    cleanEdgeType = cleanSubRefs(edgeType_)
    for cet in cleanEdgeType:
        if (STEP_entities_[ref2index(cet)].type_name == "CARTESIAN_POINT"):
            linePoint = CP2point(STEP_entities_, cet)
            parsedEdge_["linePoint"] = linePoint

        if (STEP_entities_[ref2index(cet)].type_name == "VECTOR"):
            lineVectorDisp = STEP_entities_[ref2index(cet)].params[-1]
            parsedEdge_["lineVectorDisp"] = lineVectorDisp
            lineVectorDir = STEP_entities_[ref2index(cet)].params[-2]
            if (STEP_entities_[ref2index(lineVectorDir)].type_name == "DIRECTION"):
                lineVectorDir = CP2point(STEP_entities_, lineVectorDir)
                parsedEdge_["lineVectorDir"] = lineVectorDir

    minPoint, minPointCentroidDisp, v1ext, v2ext = pointSegmentMinDisp(localCentroid, parsedEdge_["vertex1"], parsedEdge_["vertex2"])
    # todo, if there are multiple minima, return list? store minpoint as list?
    parsedEdge_["minPoint"] = minPoint
    parsedEdge_["minPointCentroidDisp"] = minPointCentroidDisp
    parsedEdge_['vertex1extremaMin'] = v1ext
    parsedEdge_['vertex2extremaMin'] = v2ext

    if not v1ext and v2ext:
        parsedEdge_['vertex1extremaMax'] = True
        parsedEdge_['vertex2extremaMax'] = False
        parsedEdge_['maxPoint'] = parsedEdge_['vertex1']
        parsedEdge_['maxPointCentroidDisp'] = parsedEdge_['vertex1centroidDisp']
    if v1ext and not v2ext:
        parsedEdge_['vertex1extremaMax'] = False
        parsedEdge_['vertex2extremaMax'] = True
        parsedEdge_['maxPoint'] = parsedEdge_['vertex2']
        parsedEdge_['maxPointCentroidDisp'] = parsedEdge_['vertex2centroidDisp']
    if not v1ext and not v2ext:
        minPointv1disp = np.linalg.norm(minPoint - parsedEdge_['vertex1'])
        minPointv2disp = np.linalg.norm(minPoint - parsedEdge_['vertex2'])
        parsedEdge_['vertex1extremaMin'] = False
        parsedEdge_['vertex2extremaMin'] = False
        parsedEdge_['vertex1extremaMax'] = False
        parsedEdge_['vertex2extremaMax'] = False
        if minPointv1disp > minPointv2disp:
            parsedEdge_["maxPoint"] = parsedEdge_['vertex2']
            parsedEdge_['maxPointCentroidDisp'] = np.linalg.norm(parsedEdge_['vertex2'] - localCentroid)
            parsedEdge_['minPointCentroidDisp'] = parsedEdge_['vertex2centroidDisp']
        elif minPointv1disp < minPointv2disp:
            parsedEdge_['maxPoint'] = parsedEdge_['vertex1']
            parsedEdge_['maxPointCentroidDisp'] = np.linalg.norm(parsedEdge_['vertex1'] - localCentroid)
            parsedEdge_['minPointCentroidDisp'] = parsedEdge_['vertex1centroidDisp']
        elif (minPointv1disp - minPointv2disp) < eps_STEP_AP21: # midpoint
            # todo maxPoint = [parsedEdge_["vertex1"], parsedEdge_["vertex2"]]
            parsedEdge_['maxPoint'] = parsedEdge_['vertex1']
            parsedEdge_['maxPointCentroidDisp'] = parsedEdge_['vertex1centroidDisp']
            parsedEdge_['vertex1extremaMax'] = True
            parsedEdge_['vertex2extremaMax'] = True

    # if parsedEdge_.get('vertex1extremaMax') is None:
    #     _1=1

    # if parsedEdge_["vertex2centroidDisp"] > parsedEdge_["vertex1centroidDisp"]:  # multiple maxima for equidistant vector?
    #     maxPoint = parsedEdge_["vertex2"]
    #     maxPointCentroidDisp = parsedEdge_["vertex2centroidDisp"]
    # parsedEdge_["maxPoint"] = maxPoint
    # parsedEdge_["maxPointCentroidDisp"] = maxPointCentroidDisp
    # parsedEdge_['vertexExtremaMax'] = True


def circleParse(localCentroid, parsedEdge, STEP_entities_):

    # if parsedEdge_.get("subEdgeType") is not None:
    #     edgeType_ = parsedEdge_["subEdgeType"]
    # else:
    edgeType_ = parsedEdge["edgeType"]

    radius = edgeType_[-1]
    parsedEdge["radius"] = radius

    axisPoint, normDir, refDir = axis2Placement3D_2(edgeType_[-2], STEP_entities_)
    parsedEdge["axisPoint"] = axisPoint
    parsedEdge["normDir"] = normDir
    parsedEdge["refDir"] = refDir
    # axisDir, refDir guaranteed unit normals?
    auxDir = np.cross(normDir, refDir)
    auxDir = auxDir / np.linalg.norm(auxDir)
    parsedEdge["auxDir"] = auxDir

    # test if cylinder is rotationally symmetrical, if localCentroid is close to axisDir through axisPoint
    rotSymDisp = np.linalg.norm(
        (localCentroid - axisPoint) - np.dot((localCentroid - axisPoint), normDir) * normDir
    )

    # # project localCentroid to disc plane
    # # distance d of point p to the plane is dot product of plane unit normal, normAxis
    # # with the difference vector from the point in the plane P and S: d = (S − P) · u
    # # distance is signed: positive when S is on the side of the plane where u faces and
    # # negative on the other side. zero, S is in the plane
    #
    # Cdisp = np.dot((localCentroid - axisPoint), normDir)
    # Cproj = localCentroid - np.dot(Cdisp, normDir)  # point projected to plane along orthogonal
    # ppCentreDisp = np.linalg.norm(Cproj - axisPoint)


    # if axisPoint[0]==20 and axisPoint[1]==0 and axisPoint[2]==10:
    #     _1=1


    # min/max point on arc wrt centroid
    minPoint, maxPoint = pointCircleMinMaxDisp(
        localCentroid,
        axisPoint,
        normDir,
        radius
    )

    v1 = parsedEdge["vertex1"]
    v2 = parsedEdge["vertex2"] # v1 always equals v2 for full circle
    parsedEdge['minPoint'] = minPoint
    parsedEdge['vertex1extremaMax'] = False
    parsedEdge['vertex2extremaMax'] = False
    parsedEdge['vertex1extremaMin'] = False
    parsedEdge['vertex2extremaMin'] = False

    # if maxPoint[0]<0:
    #     _1=1

    # centroid orthogonal to centrePoint
    if (maxPoint is None) or (rotSymDisp < eps_STEP_AP21):
        # case of minPoint equidistant from circle centre - edge case
        # assign rotSymCentre feature point at centre
        parsedEdge['minPoint'] = None
        parsedEdge["rotSymCentre"] = axisPoint
        parsedEdge["rotSymRadius"] = radius
        parsedEdge["minPointCentroidDisp"] = np.linalg.norm(centroid - axisPoint)
        parsedEdge['vertex1extremaMax'] = True
        parsedEdge['vertex2extremaMax'] = True
        parsedEdge["maxPointCentroidDisp"] = np.sqrt(radius**2 + np.linalg.norm(axisPoint - localCentroid)**2)
        parsedEdge['maxPoint'] = None # recalculation possibility
    elif parsedEdge['vertex1ref'] == parsedEdge['vertex2ref']: # full circle
        parsedEdge['maxPoint'] = maxPoint
        if np.linalg.norm(minPoint - v1) < eps_STEP_AP21:
            parsedEdge["minPointCentroidDisp"] = parsedEdge['vertex1centroidDisp']
            parsedEdge["maxPointCentroidDisp"] = np.linalg.norm(maxPoint - localCentroid)
            parsedEdge['vertex1extremaMin'] = True
            parsedEdge['vertex2extremaMin'] = True
        elif np.linalg.norm(maxPoint - v1) < eps_STEP_AP21:
            parsedEdge["maxPointCentroidDisp"] = parsedEdge['vertex1centroidDisp']
            parsedEdge["minPointCentroidDisp"] = np.linalg.norm(minPoint - localCentroid)
            parsedEdge['vertex1extremaMax'] = True
            parsedEdge['vertex2extremaMax'] = True
        else:
            parsedEdge["minPointCentroidDisp"] = np.linalg.norm(minPoint - localCentroid)
            parsedEdge["maxPointCentroidDisp"] = np.linalg.norm(maxPoint - localCentroid)
    else: # partial circle arc
        if pointInArc(minPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
            parsedEdge["minPoint"] = minPoint
            parsedEdge["minPointCentroidDisp"] = np.linalg.norm(minPoint - localCentroid)
            if np.linalg.norm(minPoint - parsedEdge["vertex1"]) < eps_STEP_AP21:
                parsedEdge['vertex1extremaMin'] = True
            if np.linalg.norm(minPoint - parsedEdge["vertex2"]) < eps_STEP_AP21:
                parsedEdge['vertex2extremaMin'] = True
        if pointInArc(maxPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
            parsedEdge["maxPoint"] = maxPoint
            parsedEdge["maxPointCentroidDisp"] = np.linalg.norm(maxPoint - localCentroid)
            if np.linalg.norm(maxPoint - v1) < eps_STEP_AP21:
                parsedEdge['vertex1extremaMax'] = True
            if np.linalg.norm(maxPoint - v2) < eps_STEP_AP21:
                parsedEdge['vertex2extremaMax'] = True


def ellipseParse(localCentroid, parsedEdge_, STEP_entities_):
    # first_semi_axis: half the length of the first diameter of the Ellipse.
    # second_semi_axis: half the length of the second diameter of the Ellipse.

    # if parsedEdge_.get("subEdgeType") is not None:
    #     edgeType_ = parsedEdge_["subEdgeType"]
    # else:
    edgeType_ = parsedEdge_["edgeType"]

    cleanEdgeType = cleanSubRefs(edgeType_)
    majorRadius = cleanEdgeType[-2]
    parsedEdge_["majorRadius"] = majorRadius

    minorRadius = cleanEdgeType[-1]
    parsedEdge_["minorRadius"] = minorRadius

    axisPoint, normDir, refDir = axis2Placement3D_2(edgeType_[-2], STEP_entities_)
    parsedEdge_["axisPoint"] = axisPoint
    parsedEdge_["normDir"] = normDir
    parsedEdge_["refDir"] = refDir
    # axisDir, refDir guaranteed unit normals?
    auxDir = np.cross(normDir, refDir)
    auxDir = auxDir / np.linalg.norm(auxDir)
    parsedEdge_["auxDir"] = auxDir

    minPoint, maxPoint = pointEllipseMinMaxDisp(
        localCentroid,
        axisPoint,
        refDir,
        auxDir,
        normDir,
        majorRadius,
        minorRadius
    )

    # test if minPoint, maxPoint is on segment between vertex1 & vertex2
    v1 = pel["vertex1"]
    v2 = pel["vertex2"]

    # if pointInArc(minPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=False):
    #     parsedEdge_["minPoint"] = minPoint
    #     parsedEdge_['vertexExtremaMin'] = False
    # else:  # one of the vertices is a minima
    #     parsedEdge_['vertexExtremaMin'] = True
    #     if parsedEdge_['vertex1centroidDisp'] < parsedEdge_['vertex2centroidDisp']:
    #         parsedEdge_["minPoint"] = parsedEdge_['vertex1']
    #         parsedEdge_["minPointCentroidDisp"] = parsedEdge_['vertex1centroidDisp']
    #     else:
    #         parsedEdge_['minPoint'] = parsedEdge_['vertex2']
    #         parsedEdge_["minPointCentroidDisp"] = parsedEdge_['vertex2centroidDisp']
    #
    # if pointInArc(maxPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
    #     parsedEdge_["maxPoint"] = maxPoint
    #     parsedEdge_['vertexExtremaMax'] = False
    # else:  # one of the vertices is a maxima
    #     parsedEdge_['vertexExtremaMax'] = True
    #     if parsedEdge_['vertex1centroidDisp'] > parsedEdge_['vertex2centroidDisp']:
    #         parsedEdge_["maxPoint"] = parsedEdge_['vertex1']
    #         parsedEdge_["maxPointCentroidDisp"] = parsedEdge_['vertex1centroidDisp']
    #     else:
    #         parsedEdge_['maxPoint'] = parsedEdge_['vertex2']
    #         parsedEdge_["maxPointCentroidDisp"] = parsedEdge_['vertex2centroidDisp']

    if parsedEdge_['vertex1ref'] == parsedEdge_['vertex2ref']:  # full ellipse? untested
        if np.linalg.norm(minPoint - v1) < eps_STEP_AP21:
            parsedEdge_['vertex1extremaMin'] = True
            parsedEdge_['vertex2extremaMin'] = True
            parsedEdge_['vertex1extremaMax'] = False
            parsedEdge_['vertex2extremaMax'] = False
        if np.linalg.norm(maxPoint - v1) < eps_STEP_AP21:
            parsedEdge_['vertex1extremaMin'] = False
            parsedEdge_['vertex2extremaMin'] = False
            parsedEdge_['vertex1extremaMax'] = True
            parsedEdge_['vertex2extremaMax'] = True

    else:  # partial ellipse arc
        if pointInArc(minPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
            parsedEdge_["minPoint"] = minPoint
            parsedEdge_["minPointCentroidDisp"] = np.linalg.norm(minPoint - localCentroid)
            if np.linalg.norm(minPoint - v1) < eps_STEP_AP21:
                parsedEdge_['vertex1extremaMin'] = True
                parsedEdge_['vertex2extremaMin'] = False
            if np.linalg.norm(minPoint - v2) < eps_STEP_AP21:
                parsedEdge_['vertex1extremaMin'] = False
                parsedEdge_['vertex2extremaMin'] = True

        if pointInArc(maxPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
            parsedEdge_["maxPoint"] = maxPoint
            parsedEdge_["maxPointCentroidDisp"] = np.linalg.norm(maxPoint - localCentroid)
            if np.linalg.norm(maxPoint - v1) < eps_STEP_AP21:
                parsedEdge_['vertex1extremaMax'] = True
                parsedEdge_['vertex2extremaMax'] = False
            if np.linalg.norm(maxPoint - v2) < eps_STEP_AP21:
                parsedEdge_['vertex1extremaMax'] = False
                parsedEdge_['vertex2extremaMax'] = True


        # if pointInArc(minPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
        #     parsedEdge_["minPoint"] = minPoint
        #     parsedEdge_["minPointCentroidDisp"] = np.linalg.norm(minPoint - centroid)
        #     parsedEdge_['vertex1extremaMin'] = False
        #     parsedEdge_['vertex2extremaMin'] = False
        # else:  # one of the vertices is a minima (only valid for straight line)
        #     if parsedEdge_['vertex1centroidDisp'] < parsedEdge_['vertex2centroidDisp']:
        #         parsedEdge_["minPoint"] = parsedEdge_['vertex1']
        #         parsedEdge_["minPointCentroidDisp"] = parsedEdge_['vertex1centroidDisp']
        #         parsedEdge_['vertex1extremaMin'] = True
        #         parsedEdge_['vertex2extremaMin'] = False
        #     else:
        #         parsedEdge_['minPoint'] = parsedEdge_['vertex2']
        #         parsedEdge_["minPointCentroidDisp"] = parsedEdge_['vertex2centroidDisp']
        #         parsedEdge_['vertex1extremaMin'] = False
        #         parsedEdge_['vertex2extremaMin'] = True
        #
        # if pointInArc(maxPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
        #     parsedEdge_["maxPoint"] = maxPoint
        #     parsedEdge_["maxPointCentroidDisp"] = np.linalg.norm(maxPoint - centroid)
        #     parsedEdge_['vertex1extremaMax'] = True
        #     parsedEdge_['vertex2extremaMax'] = True
        # else:  # one of the vertices is a maxima (only valid for straight line)
        #     if parsedEdge_['vertex1centroidDisp'] > parsedEdge_['vertex2centroidDisp']:
        #         parsedEdge_["maxPoint"] = parsedEdge_['vertex1']
        #         parsedEdge_["maxPointCentroidDisp"] = parsedEdge_['vertex1centroidDisp']
        #         parsedEdge_['vertex1extremaMax'] = True
        #         parsedEdge_['vertex2extremaMax'] = False
        #     else:
        #         parsedEdge_['maxPoint'] = parsedEdge_['vertex2']
        #         parsedEdge_["maxPointCentroidDisp"] = parsedEdge_['vertex2centroidDisp']
        #         parsedEdge_['vertex1extremaMax'] = False
        #         parsedEdge_['vertex2extremaMax'] = True


def BsplineCurveWithKnotsParse(localCentroid, parsedEdge_, STEP_entities_):
    # ENTITY B_Spline_Curve
    # 	SUPERTYPE OF ((ONEOF (Uniform_Curve, B_Spline_Curve_With_Knots, Quasi_Uniform_Curve, Bezier_Curve) ANDOR Rational_B_Spline_Curve))
    # 	SUBTYPE OF (Bounded_Curve);
    # 	degree : INTEGER;
    # 	control_points_list : LIST [2:?] OF Cartesian_Point;
    # 	curve_form : B_Spline_Curve_Form;
    # 	closed_curve : LOGICAL;
    # 	self_intersect : LOGICAL;
    # ENTITY B_Spline_Curve_With_Knots
    # 	SUBTYPE OF (B_Spline_Curve);
    # 	knot_multiplicities : LIST [2:?] OF INTEGER;
    # 	knots : LIST [2:?] OF Parameter_Value;
    # 	knot_spec : Knot_Type;

    # if parsedEdge_.get("subEdgeType") is not None:
    #     edgeType_ = parsedEdge_["subEdgeType"]
    # else:
    edgeType_ = parsedEdge_["edgeType"]

    if type(edgeType_[0]) == str:  # check for name string
        offset = 1
    else:
        offset = 0
    curveDegree = edgeType_[offset]
    controlPointsListRefs = edgeType_[offset + 1]
    curveForm = edgeType_[offset + 2]
    closedCurve = "T" in edgeType_[offset + 3]
    selfIntersect = "T" in edgeType_[offset + 4]
    knotMultiplicities = edgeType_[offset + 5]
    knots = edgeType_[offset + 6]
    knotSpec = edgeType_[offset + 7]

    parsedEdge_["curveDegree"] = curveDegree
    parsedEdge_["controlPointsListRefs"] = controlPointsListRefs
    parsedEdge_["curveForm"] = curveForm
    parsedEdge_["closedCurve"] = closedCurve
    parsedEdge_["selfIntersect"] = selfIntersect
    parsedEdge_["knotMultiplicities"] = knotMultiplicities
    parsedEdge_["knots"] = knots
    parsedEdge_["knotSpec"] = knotSpec

    # extract control points
    controlPointsList = []
    for cpl in controlPointsListRefs:
        if (STEP_entities_[ref2index(cpl)].type_name == "CARTESIAN_POINT"):
            controlPoint = CP2point(STEP_entities_, cleanSubRefs(cpl)[0])
            controlPointsList.append(controlPoint)
    parsedEdge_["controlPointsList"] = controlPointsList

    #controlPointsList = controlPointsList[1:-1] #----------------------------TEST-TEST-TEST!!!
    #knotMultiplicities = knotMultiplicities[1:-1] #----------------------------TEST-TEST-TEST!!!

    _knotvector = []
    for kml in range(len(knotMultiplicities)):
        for i in range(knotMultiplicities[kml]):
            _knotvector.append(knots[kml])
    parsedEdge_["knotvector"] = _knotvector

    if (len(_knotvector) - len(controlPointsList) - curveDegree) != 1:
        print("maldefined B-spline!")

    #  rule for closed periodic curves is : number of input control points + curve degree = the total number of control points required.
    # (Reflect.field(Math,"fabs")(((vec[i] if i >= 0 and i < len(vec) else None) - rep)) > verb_core_Constants.EPSILON)
    # NURBS curve is rational with non-unity weights

    # Creates a rational B_spline curve on the basis <Knots, Multiplicities> of degree <Degree>.
    # Raises ConstructionError subject to the following conditions 0 < Degree <= MaxDegree.
    # Weights.Length() == Poles.Length()
    # Knots.Length() == Mults.Length() >= 2
    # Knots(i) < Knots(i+1) (Knots are increasing)
    # 1 <= Mults(i) <= Degree
    #
    # On a non periodic curve the first and last multiplicities may be Degree+1
    # (this is even recommended if you want the curve to start and finish on the first and last pole).
    # On a periodic curve the first and the last multicities must be the same.
    # on non-periodic curves
    # Poles.Length() == Sum(Mults(i)) - Degree - 1 >= 2
    #
    # on periodic curves
    # Poles.Length() == Sum(Mults(i)) except the first or last

    # The knot vector must be non-decreasing and of length (degree + 1) * 2 or greater
    # [ (degree + 1 copies of the first knot), internal non-decreasing knots, (degree + 1 copies of the last knot) ]

    # knot multiplicity must be explicit for verb library # byKnotsControlPointsWeights
    # BsplineKnotCurve = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(
    #     curveDegree, knotvector, [cp.tolist() for cp in controlPointsList]
    # ) #altered line 6281 verb.py

    BsplineKnotCurve = BSpline.Curve(normalize_kv=False)
    BsplineKnotCurve.degree = 3
    BsplineKnotCurve.ctrlpts = controlPointsList
    BsplineKnotCurve.knotvector = _knotvector

    # BsplineKnotCurve = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(
    #     curveDegree, knotvector, [cp.tolist() for cp in controlPointsList]
    # ) #altered line 6281 verb.py

    minPoint = BsplineKnotCurve.closestPoint(centroid)
    minPoint = np.array([minPoint[0], minPoint[1], minPoint[2]])

    # OpenCascade
    # Creates a rational B_spline curve on the basis <Knots, Multiplicities> of degree <Degree>.
    # Raises ConstructionError subject to the following conditions:
    # 0 < Degree <= MaxDegree.
    # Weights.Length() == Poles.Length()
    # Knots.Length() == Mults.Length() >= 2
    # Knots(i) < Knots(i+1) (Knots are increasing)
    # 1 <= Mults(i) <= Degree
    #
    # On a non periodic curve the first and last multiplicities may be Degree+1
    # (this is even recommanded if you want the curve to start and finish on the first and last pole).
    # On a periodic curve the first and the last multicities must be the same.
    # on non-periodic curves
    # Poles.Length() == Sum(Mults(i)) - Degree - 1 >= 2
    # on periodic curves
    # Poles.Length() == Sum(Mults(i)) except the first or last

    # K, the knots sequence of this BSpline curve.
    # In this sequence, knots with a multiplicity greater than 1 are repeated.
    # In the case of a non-periodic curve the length of the sequence must be equal to the sum of the NbKnots multiplicities of the knots of the curve (where NbKnots is the number of knots of this BSpline curve). This sum is also equal to : NbPoles + Degree + 1 where NbPoles is the number of poles and Degree the degree of this BSpline curve. In the case of a periodic curve, if there are k periodic knots, the period is Knot(k+1) - Knot(1). The initial sequence is built by writing knots 1 to k+1, which are repeated according to their corresponding multiplicities. If Degree is the degree of the curve, the degree of continuity of the curve at the knot of index 1 (or k+1) is equal to c = Degree + 1 - Mult(1). c knots are then inserted at the beginning and end of the initial sequence:
    #
    # the c values of knots preceding the first item Knot(k+1) in the initial sequence are inserted at the beginning; the period is subtracted from these c values;
    # the c values of knots following the last item Knot(1) in the initial sequence are inserted at the end; the period is added to these c values. The length of the sequence must therefore be equal to: NbPoles + 2*Degree - Mult(1) + 2. Example For a non-periodic BSpline curve of degree 2 where:
    # the array of knots is: { k1 k2 k3 k4 },
    # with associated multiplicities: { 3 1 2 3 }, the knot sequence is: K = { k1 k1 k1 k2 k3 k3 k4 k4 k4 } For a periodic BSpline curve of degree 4 , which is "C1" continuous at the first knot, and where :
    # the periodic knots are: { k1 k2 k3 (k4) } (3 periodic knots: the points of parameter k1 and k4 are identical, the period is p = k4 - k1),
    # with associated multiplicities: { 3 1 2 (3) }, the degree of continuity at knots k1 and k4 is: Degree + 1 - Mult(i) = 2. 2 supplementary knots are added at the beginning and end of the sequence:
    # at the beginning: the 2 knots preceding k4 minus the period; in this example, this is k3 - p both times;
    # at the end: the 2 knots following k1 plus the period; in this example, this is k2 + p and k3 + p. The knot sequence is therefore: K = { k3-p k3-p k1 k1 k1 k2 k3 k3 k4 k4 k4 k2+p k3+p } Exceptions Standard_DimensionError if the array K is not of the appropriate length.Returns the knots sequence.

    # for each spline you have to satisfy the equation: m = n + p + 1, where
    # m + 1 is the number of the knots (each knot counted separately, no multiplicity)
    # n + 1 is the number of control points
    # p is the curve degree
    #
    # OCC documentation for Geom_BSplieCurve:
    # Knot(i+k) = Knot(i) + period
    # Pole(i+p) = Pole(i)
    # exclude the repeating control point from description of the spline and the 'redundant' knots.

    # Going from B-Splines to NURBS means the addition of weights and non-uniform knots.
    # Weights are straightforward since they correspond directly to points.
    # As for knots, an unclosed NURBS spline requires n + degree + 1 knots, where n is the number of points.
    # A closed NURBS curve requires only n + 1 knots. The periodicity is defined by equivalence of the first and last knot,
    # i.e. knot[0] == knot[-1]

    # an unclosed NURBS spline requires n + degree + 1 knots, where n is the number of points.
    # A closed NURBS curve requires only n + 1 knots.

    # STEP files mostly use NURBS when defining curves.
    # But the "b_spline_curve_with_knots" command allows for uniform b splines too.
    # It does not discern between the two.

    # STEP files almost only use "open b splines", which means the curve starts and ends at the first and last control points.
    # So the task was to find the control points in between, which there are n-1 of. But there are only n-3 interpolating points,
    # so I input the last two equations as the tangents at the start and end.
    # This is now a linear system of n-1 vector equations, which is easily solved.

    # constructing closed B-spline curves by wrapping knots.
    # to construct a closed B-spline curve C(u) of degree p defined by n+1 control points P0, P1, ..., Pn.
    # The following is the construction procedure:
    # Add a new control point Pn+1 = P0. Therefore, the number of control points is n+2.
    # Find an appropriate knot sequence of n+1 knots u0, u1, ..., un. These knots are not necessarily uniform, an advantage over the method discussed above.
    # Add p+2 knots and wrap around the first p+2 knots: un+1 = u0, un+2 = u1, ..., un+p = up-1, un+p+1 = up, un+p+2 = up+1 as shown in the following diagram. In this way, we have n+p+2 = (n+1) + p + 1 knots

    # https://github.com/FreeCAD/FreeCAD/blob/69097667df47b2cc86d8688d2dbb545319e33e68/src/Mod/Sketcher/App/Sketch.cpp#L1384

    #     // If periodic, startpoint and endpoint do not play a role in the solver, this can remove
    #     // unnecessary DoF of determining where in the curve the start and the stop should be. However,
    #     // since start and end points are placed above knots, removing them leads to that knot being
    #     // unusable.

    #     // If periodic, startpoint and endpoint do not play a role in the solver, this can remove
    #     // unnecessary DoF of determining where in the curve the start and the stop should be. However,
    #     // since start and end points are placed above knots, removing them leads to that knot being
    #     // unusable.

    #     // WARNING: This is only valid where the multiplicity of the endpoints conforms with a BSpline
    #     // only then the startpoint is the first control point and the endpoint is the last control
    #     // point accordingly, it is never the case for a periodic BSpline. NOTE: For an external
    #     // B-spline (i.e. fixed=true) we must not set the coincident constraints as the points are not
    #     // movable anyway. See #issue 0003176: Sketcher: always over-constrained when referencing
    #     // external B-Spline

# https://git.dev.opencascade.org/gitweb/?p=occt.git;a=blob;f=src/BSplCLib/BSplCLib.cxx;h=1414f5008b3b2a79583a0af67e0a05fcd52569cb;hb=HEAD
    # 1917   // index is the first pole of the current curve for insertion schema
    # 1918
    # 1919   if (Periodic) index = -Mults(Mults.Lower());
    # 1920   else          index = -Degree-1;

    # 1922   // on Periodic curves the first knot and the last knot are inserted later
    # 1923   // (they are the same knot)
    # 1924   firstmult = 0;  // multiplicity of the first-last knot for periodic


    # Newton root finding over a bunch of points probably makes more sense for mapping local maxima/minima
    # TODO: pointNURBSmaxDisp() may be more accurate for minima

    #maxPointU = pointNURBSmaxDisp(localCentroid, BsplineKnotCurve._data)  # TODO: only finds one maxima
    #maxPoint = verb.verb_eval_Eval.dehomogenize(verb.verb_eval_Eval.curvePoint(BsplineKnotCurve._data, maxPointU))
    #maxPoint = np.array([maxPoint[0], maxPoint[1], maxPoint[2]])
    maxPointUV = splineCurveMinMaxPointDisp(BsplineKnotCurve, localCentroid, maxSearch=True)
    maxPoint = BsplineKnotCurve.derivatives(maxPointUV[0], maxPointUV[1], order=0)[0][0]

    minPointCentroidDisp = np.linalg.norm(minPoint - localCentroid)
    parsedEdge_["minPointCentroidDisp"] = minPointCentroidDisp
    parsedEdge_["minPoint"] = minPoint
    maxPointCentroidDisp = np.linalg.norm(maxPoint - localCentroid)
    parsedEdge_["maxPointCentroidDisp"] = maxPointCentroidDisp
    parsedEdge_["maxPoint"] = maxPoint

    v1 = parsedEdge_["vertex1"]
    v2 = parsedEdge_["vertex2"]

    # independent check for vertices,
    parsedEdge_['vertex1extremaMin'] = False
    parsedEdge_['vertex2extremaMin'] = False
    parsedEdge_['vertex1extremaMax'] = False
    parsedEdge_['vertex2extremaMax'] = False

    if parsedEdge_['vertex1centroidDisp'] < minPointCentroidDisp:
        parsedEdge_["minPoint"] = parsedEdge_['vertex1']
        parsedEdge_["minPointCentroidDisp"] = parsedEdge_['vertex1centroidDisp']
        parsedEdge_['vertex1extremaMin'] = True

    if parsedEdge_['vertex2centroidDisp'] < minPointCentroidDisp:
        parsedEdge_["minPoint"] = parsedEdge_['vertex2']
        parsedEdge_["minPointCentroidDisp"] = parsedEdge_['vertex2centroidDisp']
        parsedEdge_['vertex2extremaMin'] = True

    if parsedEdge_['vertex1centroidDisp'] > maxPointCentroidDisp:
        parsedEdge_["maxPoint"] = parsedEdge_['vertex1']
        parsedEdge_["maxPointCentroidDisp"] = parsedEdge_['vertex1centroidDisp']
        parsedEdge_['vertex1extremaMax'] = True

    if parsedEdge_['vertex2centroidDisp'] > maxPointCentroidDisp:
        parsedEdge_["maxPoint"] = parsedEdge_['vertex2']
        parsedEdge_["maxPointCentroidDisp"] = parsedEdge_['vertex2centroidDisp']
        parsedEdge_['vertex2extremaMax'] = True


# def subEdgeParse(edgeInstance_, parsedEdge_, STEP_entities_):
#     # SEAM_EDGE,
#     parsedEdge_["subTypeName"] = STEP_entities_[ref2index(edgeInstance_)].type_name


def commonEdgeParse(edgeInstance_, STEP_entities_):
    parsedEdge = {}
    if (STEP_entities_[ref2index(edgeInstance_)].type_name == "ORIENTED_EDGE"):
        # if not hasattr(STEP_entities_[ref2index(edgeLoopInstance)].type_name, 'ORIENTED_EDGE'): ----------------------
        VOR = "F" in STEP_entities_[ref2index(edgeInstance_)].params[-1]  # vertexOrderReversed
        parsedEdge["VOR"] = VOR

    edgeCurveRef = STEP_entities_[ref2index(edgeInstance_)].params[-2]
    assert STEP_entities_[ref2index(edgeCurveRef)].type_name == "EDGE_CURVE"
    # edge_geometry: the Curve defining the geometric shape of the edge.
    # same_sense: a BOOLEAN variable giving the relationship between the topological sense
    #             of the edge and the parametric sense of the curve.
    edgeCurveParams = STEP_entities_[ref2index(edgeCurveRef)].params
    # edgeSameSense = 'T' in edgeCurveParams[-1]; parsedEdge['sameSense'] = edgeSameSense
    edgeCurveParams = cleanSubRefs(edgeCurveParams)

    if (STEP_entities_[ref2index(edgeCurveParams[0])].type_name == "VERTEX_POINT"):
        vertex1 = STEP_entities_[ref2index(edgeCurveParams[0])].params
        vertex1 = CP2point(STEP_entities_, cleanSubRefs(vertex1)[0])
        if not VOR:
            parsedEdge["vertex1"] = vertex1
            parsedEdge["vertex1ref"] = STEP_entities_[ref2index(edgeCurveParams[0])].ref
            vertex1centroidDisp = np.linalg.norm(centroid - vertex1)
            parsedEdge["vertex1centroidDisp"] = vertex1centroidDisp
        else:
            parsedEdge["vertex2"] = vertex1
            parsedEdge["vertex2ref"] = STEP_entities_[ref2index(edgeCurveParams[0])].ref
            vertex2centroidDisp = np.linalg.norm(centroid - vertex1)
            parsedEdge["vertex2centroidDisp"] = vertex2centroidDisp

    if (STEP_entities_[ref2index(edgeCurveParams[1])].type_name == "VERTEX_POINT"):
        vertex2 = STEP_entities_[ref2index(edgeCurveParams[1])].params
        vertex2 = CP2point(STEP_entities_, cleanSubRefs(vertex2)[0])
        if not VOR:
            parsedEdge["vertex2"] = vertex2
            parsedEdge["vertex2ref"] = STEP_entities_[ref2index(edgeCurveParams[1])].ref
            vertex2centroidDisp = np.linalg.norm(centroid - vertex2)
            parsedEdge["vertex2centroidDisp"] = vertex2centroidDisp
        else:
            parsedEdge["vertex1"] = vertex2
            parsedEdge["vertex1ref"] = STEP_entities_[ref2index(edgeCurveParams[0])].ref
            vertex2centroidDisp = np.linalg.norm(centroid - vertex2)
            parsedEdge["vertex1centroidDisp"] = vertex2centroidDisp

    edgeType = STEP_entities_[ref2index(edgeCurveParams[2])].params
    edgeTypeName = STEP_entities_[ref2index(edgeCurveParams[2])].type_name
    edgeRef = STEP_entities_[ref2index(edgeCurveParams[2])].ref
    parsedEdge["typeName"] = edgeTypeName
    parsedEdge["edgeRef"] = edgeRef
    parsedEdge["edgeType"] = edgeType

    return parsedEdge


def simpleEdgeParse(c, parsedEdge_, STEP_entities_):
    #print("simpleEdgeParse edge type: " + parsedEdge_["typeName"])

    if (parsedEdge_["typeName"] == "LINE"):
        lineParse(c, parsedEdge_, STEP_entities_)

    elif (parsedEdge_["typeName"] == "CIRCLE"):
        circleParse(c, parsedEdge_, STEP_entities_)

    elif (parsedEdge_["typeName"] == "ELLIPSE"):
        ellipseParse(c, parsedEdge_, STEP_entities_)

    elif (parsedEdge_["typeName"] == "B_SPLINE_CURVE_WITH_KNOTS"):
        BsplineCurveWithKnotsParse(c, parsedEdge_, STEP_entities_)

    else: print("simpleEdgeParse unknown edge type: " + parsedEdge_["typeName"])


def sphereSurfaceParse(AFS, localCentroid):
    # ENTITY Spherical_surface
    #   SUBTYPE OF (Surface);
    #   position : Axis_placement_3d;
    #   radius : positive_length_measure;

    # Attribute definitions:
    # position: an Axis_placement_3d that defines the location and orientation of the surface for the purposes of parameterization.
    # The centre of the Spherical_surface is at the location
    # radius: the radius of the Spherical_surface.

    #  The axis is the placement Z axis direction and the ref_direction is an approximation to
    # the placement X axis direction.= normDir

    ParsedSurface = {}
    if (AFS.get("maxPoint") is not None) and (AFS.get("minPoint") is not None): # recalculation of max/min
        centrePoint = AFS["axisPoint"]
        normDir = AFS["normDir"] # Z
        refDir = AFS["refDir"] # X
        auxDir = AFS["auxDir"]  # Y
    else:
        centrePoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        radius = AFS["SurfaceParams"][-1]

    ParsedSurface["radius"] = radius
    ParsedSurface["axisPoint"] = centrePoint
    ParsedSurface["normDir"] = normDir  # Z
    ParsedSurface["refDir"] = refDir  # X
    ParsedSurface["auxDir"] = auxDir

    # sphere surface maxima & minima exist on a vector between centroid and sphere centre point
    # minima/maxima = centre point -/+ radius

    # "The local z-axis corresponds to the normal of planar and spherical surfaces and to the axis of cylindrical,
    # conical and toroidal surfaces."
    # Jaider Oussama et al Int. Journal of Engineering Research and Applications www.ijera.com
    # ISSN : 2248-9622, Vol. 4, Issue 5( Version 6), May 2014, pp.14-25
    # refDir, auxDir doesn't mean much in stock implementaions.

    sphereCentreCentroidDir = centrePoint - localCentroid
    sphereCentreCentroidDir = sphereCentreCentroidDir / np.linalg.norm(sphereCentreCentroidDir)

    # determine if maxPoint/minPoint are inside or outside an edge_loop defined boundary
    # CCW directionality of edge vertices determines outward surface normal
    # surface minPoint/maxPoint always orthonormal to localCentroid

    maxPoint = centrePoint + sphereCentreCentroidDir * radius
    maxPointCentroidDisp = maxPoint - localCentroid

    if not insideOutsideSurfaceTest(localCentroid, maxPoint, AFS):
        # find maxima among edges
        maxPointCentroidDisp = [edge['maxPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
        maxPoint = AFS["FaceOuterBoundEdgeLoopList"][maxPointCentroidDisp.index(max(maxPointCentroidDisp))]
        maxPointCentroidDisp = max(maxPointCentroidDisp)
    else:
        ParsedSurface["maxPoint"] = maxPoint
        ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp

    minPoint = centrePoint - sphereCentreCentroidDir * radius
    minPointCentroidDisp = minPoint - localCentroid

    if not insideOutsideSurfaceTest(localCentroid, minPoint, AFS):
        # find maxima among edges
        minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
        minPoint = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]
        minPointCentroidDisp = min(minPointCentroidDisp)
    else:
        ParsedSurface["minPoint"] = minPoint
        ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp

    AFS['ParsedSurface'] = ParsedSurface
    #return ParsedSurface


def toroidSurfaceParse(AFS, localCentroid):
    # ENTITY Toroidal_surface
    #   SUBTYPE OF (Surface);
    #   position : Axis_placement_3d;
    #   radius : positive_length_measure;
    #   minor_radius : positive_length_measure;

    # Attribute definitions:
    # position: an Axis_placement_3d that defines the location and orientation of the surface.
    # The centre of the Toroidal_surface is at the location. The plane of the ref_direction and the axis is a
    # central plane of symmetry of the Toroidal_surface. The normal to this plane is the central axis of the surface.
    # radius: the distance from the central axis of the Toroidal_surface to the centre of one of the circles
    # produced by intersecting the Toroidal_surface with a plane containing the central axis.
    # The surface can be produced by sweeping a circle of radius minor_radius along a circle of radius radius.
    # minor_radius: the radius of one of the circles produced by intersecting the Toroidal_surface with a plane
    # containing the central axis.

    ParsedSurface = {}

    if (AFS.get("maxPoint") is not None) and (AFS.get("minPoint") is not None): # recalculation of max/min
        axisPoint = AFS["axisPoint"]
        normDir = AFS["normDir"] # Z
        refDir = AFS["refDir"] # X
        auxDir = AFS["auxDir"]  # Y
    else:
        axisPoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        minorRadius = AFS["SurfaceParams"][-1]
        majorRadius = AFS["SurfaceParams"][-2]

    ParsedSurface["axisPoint"] = axisPoint
    ParsedSurface["normDir"] = normDir  # Z
    ParsedSurface["refDir"] = refDir  # X
    ParsedSurface["auxDir"] = auxDir

    ParsedSurface["minorRadius"] = minorRadius
    ParsedSurface["majorRadius"] = majorRadius

    # axisPoint, auxDir, refDir = axis2Placement3D(AFS["SurfaceParams"][1], STEP_entities)
    # ParsedSurface["axisPoint"] = axisPoint
    # ParsedSurface["auxDir"] = auxDir
    # ParsedSurface["refDir"] = refDir
    # normDir = np.cross(refDir, auxDir)
    # normDir = normDir / np.linalg.norm(normDir)
    # ParsedSurface["normDir"] = normDir

    # localCentroid = np.array([1, 1, 0])
    # axisPoint = np.array([0, 0, 10])
    # axisDir = np.array([0, 0, 1])
    # refDir = np.array([1, 0, 0])
    # tMajorRadius = 10
    # tMinorRadius = 1

    # test if toroid is rotationally symmetrical, if localCentroid is close to normDir through axisPoint

    rotSymDisp = np.linalg.norm(
        (centroid - axisPoint) - np.dot((centroid - axisPoint), normDir) * normDir
    )

    # need a precision factor that reflects increasing accuracy arising from iterative calculations of localCentroid --------------------------------------<<<<<<<<<<<<<<<<<
    # does one find a localCentroid through axis-intersection?
    # also scale factor affecting numerical precision
    # also STEP seems to be E-6 decimal places at best
    # iterativeScaledPrecision = 1e-4  # for testing

    if rotSymDisp < eps_STEP_AP21:  # iterativeScaledPrecision:
        #ParsedSurface["rotSymMax"] = axisPoint
        ParsedSurface["rotSymCentre"] = axisPoint
    else:
        # create axes for circles based at points of majorMinPoint, majorMaxPoint
        # minima case, create new axis at majorMinPoint
        # get min, max of major radius
        majorMinPoint, majorMaxPoint = pointCircleMinMaxDisp(
            localCentroid,
            axisPoint,
            normDir,
            majorRadius
        )

        # vector from toroid centre to minima point on majorRadius
        majorMinRadialDir = (majorMinPoint - axisPoint) / np.linalg.norm(
            majorMinPoint - axisPoint
        )
        # vector tangential to major radius
        majorMinTangentDir = np.cross(majorMinRadialDir, normDir)
        majorMinTangentDir = majorMinTangentDir / np.linalg.norm(
            majorMinTangentDir
        )

        # point on minor circle cross-section at minimum disp from localCentroid
        minorMinPoint, _ = pointCircleMinMaxDisp(
            localCentroid,
            majorMinPoint,
            majorMinTangentDir,
            minorRadius,
        )

        #ParsedSurface["minPoint"] = minorMinPoint

        # same for maxima
        majorMaxRadialDir = (majorMaxPoint - axisPoint) / np.linalg.norm(
            majorMaxPoint - axisPoint
        )
        majorMaxTangentDir = np.cross(majorMaxRadialDir, auxDir)
        majorMaxTangentDir = majorMaxTangentDir / np.linalg.norm(
            majorMaxTangentDir
        )

        minorMaxPoint, _ = pointCircleMinMaxDisp(
            localCentroid,
            majorMaxPoint,
            majorMaxTangentDir,
            minorRadius,
        )

        # np.linalg.norm(centroid - tMinorMinPoint) < np.linalg.norm(centroid - tMinorMaxPoint)
        #ParsedSurface["maxPoint"] = minorMaxPoint

        # est whether maxima/minima are within defined edges
        # call edge loop tests after surface max/min points extraction, line 1384

        # determine if maxPoint/minPoint are inside or outside an edge_loop defined boundary
        # CCW directionality of edge vertices determines outward surface normal
        # surface minPoint/maxPoint always orthonormal to localCentroid

        if not insideOutsideSurfaceTest(localCentroid, minorMaxPoint, AFS):
            # find maxima among edges
            maxPointCentroidDisp = [edge['maxPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
            minorMaxPoint = AFS["FaceOuterBoundEdgeLoopList"][maxPointCentroidDisp.index(max(maxPointCentroidDisp))]
            maxPointCentroidDisp = max(maxPointCentroidDisp)
        else:
            ParsedSurface["maxPoint"] = minorMaxPoint
            ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp

        if not insideOutsideSurfaceTest(localCentroid, minorMinPoint, AFS):
            # find minima among edges
            minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
            minorMinPoint = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]
            minPointCentroidDisp = min(minPointCentroidDisp)
        else:
            ParsedSurface["minPoint"] = minorMinPoint
            ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp

    AFS['ParsedSurface'] = ParsedSurface
    #return ParsedSurface


def cylinderSurfaceParse(AFS, localCentroid):
    # ENTITY Cylindrical_surface
    #   SUBTYPE OF (Surface);
    #   position : Axis_placement_3d;
    #   radius : positive_length_measure;

    # position: an Axis_placement_3d that defines the location and orientation of the surface.
    # The axis of the Cylindrical_surface passes through the location and is normal to the plane of ref_direction and axis.
    # radius: the radius of the circular curve of intersection between the surface of the cylinder and a plane
    # perpendicular to the axis of the Cylindrical_surface.

    ParsedSurface = {}

    if (AFS.get("maxPoint") is not None) and (AFS.get("minPoint") is not None): # recalculation of max/min
        axisPoint = AFS["axisPoint"]
        normDir = AFS["normDir"] # Z
        refDir = AFS["refDir"] # X
        auxDir = AFS["auxDir"]  # Y
    else:
        axisPoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        radius = AFS["SurfaceParams"][-1]

    ParsedSurface["radius"] = radius
    ParsedSurface["axisPoint"] = axisPoint
    ParsedSurface["normDir"] = normDir  # Z
    ParsedSurface["refDir"] = refDir  # X
    ParsedSurface["auxDir"] = auxDir

    # test if cylinder is rotationally symmetrical, if localCentroid is close to axisDir through axisPoint
    rotSymDisp = np.linalg.norm(
        (localCentroid - axisPoint) - np.dot((localCentroid - axisPoint), normDir) * normDir
    )

    if rotSymDisp < eps_STEP_AP21:
        # ParsedSurface["rotSymMax"] = axisPoint at circle
        minSurfaceDispPoint = localCentroid  # for creating a minima plane intersection
        #centroidProjAxisPoint = centroid
        # maxima is circular edge centre, determined from vertices
        # minima is orthogonal to centroid
        # shouldn't there also be a radius measure?? check: not in doctorate schema

        # todo, update concept where the region between a ridge and groove is a sphere
        # todo: truth determined via wandering hillclimb point path, or grad

        # test if any edge crosses cylinder midsection
        # edges should all have been calculated prior to surface

        # circle edge cannot be declared groove or ridge prior to discovery of adjacent surfaces
        # rotationally-symmetric edge must be assigned indeterminate status with centroid & radius


        for edge in AFS["FaceOuterBoundEdgeLoopList"]:
            if np.isclose(edge["minPointCentroidDisp"], radius, eps): # eps_STEP_AP21 ?? todo could this fail for small values?
                if ParsedSurface.get("rotSymRadius") is None:
                    #ParsedSurface["rotSymCentre"] = np.array([,, edge['minPoint'][3]])
                    ParsedSurface["rotSymRadius"] = radius
                    ParsedSurface["minPoint"] = None # a band rather than a point

                    # project centroid to symmetric axis to find rotSymCentre
                    #centroidProjAxisPoint = normDir * np.dot((axisPoint - localCentroid), normDir) / np.dot(normDir, normDir)
                    centroidProjAxisPoint = pointProjectAxis(localCentroid, axisPoint, normDir)
                    ParsedSurface["rotSymCentre"] = centroidProjAxisPoint
                    # minPointCentroidDisp/maxPointCentroidDisp exist where centroid is on rotSym axis
                    ParsedSurface['minPointCentroidDisp'] = np.sqrt(np.linalg.norm(centroidProjAxisPoint - localCentroid)**2 + radius**2)
                    # if np.isclose(edge["minPointCentroidDisp"], radius, eps):
                    # ignore maxPoints at cylinder extrema and rely on circle edges?---------------------------------------------
    else:
        # cylinder/conic maxima is just a search of greatest disp of vertices from localCentroid
        # if localCentroid is outside full cylindrical surface,
        #   cylinder minima -> test for intersection with surface on orthogonal projection from axis to localCentroid
        # if localCentroid within surface, |centroid - projAxisDir| < cRadius

        # test:
        # axisPoint = np.array([1, 0, 0])
        # localCentroid = np.array([2, 2, 0])
        # axisDir = np.array([1.5, 0, 0])
        # cRadius = 1

        # orthogonal projection of localCentroid to cylinder axis
        # vCA = localCentroid - axisPoint
        # vNA = axisPoint - normDir
        # centroidProjAxisPoint = vNA*np.dot(vCA, vNA)/np.dot(vNA, vNA)

        centroidProjAxisPoint = pointProjectAxis(localCentroid, axisPoint, normDir)



        # centroidProjAxisPoint = axisPoint + (
        #             np.dot((centroid - axisPoint), normDir) / np.dot(normDir, normDir)) * normDir
        centroidAxisDir = centroidProjAxisPoint - localCentroid
        centroidAxisDisp = np.linalg.norm(centroidAxisDir)
        centroidAxisDir = centroidAxisDir / centroidAxisDisp

        if centroidAxisDisp > radius:
            surfaceMinPoint = localCentroid + (centroidAxisDisp - radius) * centroidAxisDir
        elif centroidAxisDisp < radius:
            surfaceMinPoint = localCentroid - (radius - centroidAxisDisp) * centroidAxisDir
        elif centroidAxisDisp - radius < eps:
            surfaceMinPoint = localCentroid

        # insideOutsideSurfaceTest() method to determine if surfaceMinPoint or surfaceMaxPoint within defined surface edges
        # note the max point will be on one of the edges

        if not insideOutsideSurfaceTest(localCentroid, surfaceMinPoint, AFS):
            # find minima among edges
            minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
            surfaceMinPoint = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]["minPoint"]
            minPointCentroidDisp = min(minPointCentroidDisp)
            ParsedSurface["minPoint"] = surfaceMinPoint
            ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp
        else:
            ParsedSurface["minPoint"] = surfaceMinPoint
            ParsedSurface["minPointCentroidDisp"] = np.linalg.norm(surfaceMinPoint - localCentroid)

    # find maxima among edges
    maxPointCentroidDisp = [edge['maxPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"] if edge.get('maxPointCentroidDisp') is not None]
    ParsedSurface["maxPoint"] = AFS["FaceOuterBoundEdgeLoopList"][maxPointCentroidDisp.index(max(maxPointCentroidDisp))]["maxPoint"]
    ParsedSurface["maxPointCentroidDisp"] = max(maxPointCentroidDisp)
    # todo --what happens with two equal maxpoints!

    # # shift below to function?=========================================================???????????????
    #
    # intersectedEdgeSet = []
    # intersectPoint = []
    # #rotSymFeature = []
    # rotSymFeatureEdgeSet = []
    # ParsedSurface["rotSymEdge"] = []
    #
    #
    # # determine whether any edge in edge loop contains an intersection with a radius circle or plane containing this point
    # for pel in parsedEdgeLoop:
    #
    #     v1 = pel["vertex1"]
    #     v2 = pel["vertex2"]
    #
    #     if pel["typeName"] == "LINE":
    #         # test whether a plane containing minPoint thru axisDir intersects an edge segment
    #         edgeIntersect = intersectSegmentPlane(v1, v2,
    #                                               normDir,
    #                                               minSurfaceDispPoint)
    #         if edgeIntersect is not None:
    #             intersectedEdgeSet.append(pel)
    #             intersectPoint.append(edgeIntersect)
    #
    #     if pel["typeName"] == "CIRCLE":
    #
    #         edgeArcIntersect = intersectArcPlane(normDir,
    #                                              centroidProjAxisPoint,
    #                                              pel["axisPoint"],
    #                                              pel["refDir"],
    #                                              pel["auxDir"],
    #                                              pel["normDir"],
    #                                              pel["radius"],
    #                                              v1,
    #                                              v2)
    #
    #         if len(edgeArcIntersect) > 0:
    #             intersectPoint = intersectPoint + edgeArcIntersect
    #             intersectedEdgeSet.append(pel)
    #
    #         # create a rotationally symmetric edge centre point
    #         # perhaps easier to determine ridges/grooves following minima maxima chase?
    #         # requires rot-sym surface axis to pass through centroid
    #
    #         # if rotSymDisp < eps_STEP_AP21:
    #         #     # this is a circle on a cylinder with axis colinear with centroid
    #         #     pass
    #
    #         if (np.linalg.norm(pel["normDir"] - ParsedSurface["normDir"]) < eps_STEP_AP21 or
    #                 np.linalg.norm(pel["normDir"] + ParsedSurface["normDir"]) < eps_STEP_AP21):
    #             # arc circle norm is parallel to surface cylinder norm, (default case for simple cylinder)
    #             # test whether centroid projects to arc centre
    #
    #             arcDisp = np.dot((centroid - pel["axisPoint"]), pel["normDir"])
    #
    #             # centroid projected to arc plane
    #             projCentroid = centroid - np.dot(arcDisp, pel["normDir"])
    #             if (np.linalg.norm(projCentroid - pel["axisPoint"]) < eps_STEP_AP21 or
    #                     np.linalg.norm(projCentroid + pel["axisPoint"]) < eps_STEP_AP21):
    #                 # arc on a cylinder cannot be determined as groove or ridge without surrounding surface data
    #                 ParsedSurface["rotSymEdge"].append(pel["axisPoint"])
    #                 rotSymFeatureEdgeSet.append(pel)
    #                 # rotSymFeature.append(pel["axisPoint"])
    #
    #             _1 = 1
    #
    #         # # circleYaxis = np.cross(pel["axisDir"], pel["refDir"])
    #         # # circleYaxis = circleYaxis/np.linalg.norm(circleYaxis)
    #         #
    #         # circleEdge = verb.verb_geom_Circle( pel["axisPoint"].tolist(),
    #         #                                     pel["refDir"].tolist(),
    #         #                                     pel["auxDir"].tolist(),
    #         #                                     pel["cRadius"] )
    #         #
    #         # circleMinSurfaceRadius = verb.verb_geom_Circle( axisPoint.tolist(),
    #         #                                                 refDir.tolist(),
    #         #                                                 auxDir.tolist(),
    #         #                                                 centroidAxisDisp )
    #         #
    #         # edgeArcIntersect = verb.verb_geom_Intersect.curves(circleEdge, circleMinSurfaceRadius)
    #         #
    #         # # general intersections of edge arc with plane of minSurfaceDisp point
    #         # # edgeArcIntersect1, edgeArcIntersect2 = intersectArcRotSymSurface( surfaceAxisDir,
    #         # #                                                                   surfaceAxisPoint,
    #         # #                                                                   minSurfaceDisp,
    #         # #                                                                   pel["axisDir"],
    #         # #                                                                   pel["axisPoint"],
    #         # #                                                                   pel["cRadius"] )
    #         #
    #         # # test if circle intersections lies in vertex defined arc
    #         # if edgeArcIntersect1 is not None:
    #         #     if pointInArc( edgeArcIntersect1,
    #         #                    v1, v2,
    #         #                    pel["axisDir"],
    #         #                    pel["refDir"],
    #         #                    pel["axisPoint"]):
    #         #         intersectedEdgeSet.append(pel)
    #         #         intersectPoint.append(edgeIntersect)
    #         #
    #         # if edgeArcIntersect2 is not None:
    #         #     if pointInArc( edgeArcIntersect2,
    #         #                    v1, v2,
    #         #                    pel["axisDir"],
    #         #                    pel["refDir"],
    #         #                    pel["axisPoint"]):
    #         #         if pel not in intersectedEdgeSet:
    #         #             intersectedEdgeSet.append(pel)
    #         #         intersectPoint.append(edgeIntersect)
    #
    #     if pel["typeName"] == "ELLIPSE":
    #         # general intersections of edge arc with plane of minSurfaceDisp point
    #         # edgeArcIntersect1, edgeArcIntersect2 = intersectArcRotSymSurface( surfaceAxisDir,
    #         #                                                                   surfaceAxisPoint,
    #         #                                                                   minSurfaceDisp,
    #         #                                                                   pel["axisDir"],
    #         #                                                                   pel["axisPoint"],
    #         #                                                                   pel["eMajorRadius"],
    #         #                                                                   pel["eMinorRadius"])
    #
    #         # class EllipseArc extends NurbsCurve {
    #         #     //* Length 3 array representing the center of the arc
    #         #     //* Length 3 array representing the xaxis
    #         #     //* Length 3 array representing the perpendicular yaxis
    #         #     //* Minimum angle of the EllipseArc
    #         #     //* Maximum angle of the EllipseArc
    #         #
    #         #     public function new(   center : Point,
    #         #                            xaxis : Vector,
    #         #                            yaxis : Vector,
    #         #                            minAngle : Float,
    #         #                            maxAngle : Float )
    #
    #         # class Ellipse extends EllipseArc {
    #         #
    #         #     //Create an ellipse
    #         #     //* Length 3 array representing the center of the circle
    #         #     //* Length 3 array representing the xaxis
    #         #     //* Length 3 array representing the perpendicular yaxis
    #         #
    #         #     public function new(   center : Point,
    #         #                             xaxis : Vector,
    #         #                             yaxis : Vector )
    #
    #         # normal direction is inferred from xaxis X yaxis rather than STEP axisDir convention
    #         # major/minor radius is defined as xaxis/yaxis
    #         # to convert from STEP, find major/minor axis vectors and scale to minorRadius majorRadius values
    #         # have to assume AP21 majorRadius is aligned with refDir until proven otherwise
    #
    #         ellipseYaxis = pel["auxDir"]
    #         ellipseYaxis = ellipseYaxis / np.linalg.norm(ellipseYaxis)
    #         ellipseYaxis = ellipseYaxis * pel["minorRadius"]
    #
    #         ellipseXaxis = pel["refDir"]
    #         ellipseXaxis = ellipseXaxis / np.linalg.norm(ellipseXaxis)
    #         ellipseXaxis = ellipseXaxis * pel["majorRadius"]
    #
    #         ellipseEdge = verb.verb_geom_Ellipse(pel["axisPoint"].tolist(),
    #                                              ellipseXaxis.tolist(),
    #                                              ellipseYaxis.tolist())
    #
    #         # surfaceYaxisDir = np.cross(surfaceAxisDir, surfaceRefDir)
    #         # surfaceYaxisDir = surfaceYaxisDir/np.linalg.norm(surfaceYaxisDir)
    #
    #         # circleMinSurfaceRadius = verb.verb_geom_Circle( axisPoint.tolist(),
    #         #                                                 refDir.tolist(),
    #         #                                                 auxDir.tolist(),
    #         #                                                 centroidAxisDisp )
    #
    #         # verb core Plane
    #
    #         auxDirP = np.cross(centroidAxisDir, auxDir)  # surface axisDir
    #         p1 = centroidProjAxisPoint + centroidAxisDir * centroidAxisDisp
    #         p3 = centroidProjAxisPoint - centroidAxisDir * centroidAxisDisp
    #         p2 = centroidProjAxisPoint + auxDirP * centroidAxisDisp
    #         p4 = centroidProjAxisPoint - auxDirP * centroidAxisDisp
    #
    #         intersectPlane = verb.verb_eval_Make.fourPointSurface(p1, p2, p3, p4)
    #
    #         # class Intersect {
    #         #
    #         #     //Determine the intersection of two curves
    #         #     //* ICurve object
    #         #     //* ICurve object
    #         #     //* tolerance for the intersection
    #         #     //
    #         #     //**returns**
    #         #     //* a possibly empty array of CurveCurveIntersection objects
    #         #
    #         #     public static function curves( first : ICurve, second : ICurve, tol : Float = 1e-3  )
    #
    #         # edgeArcIntersect = verb.verb_geom_Intersect.curveAndSurface(ellipseEdge, intersectPlane)
    #         # edgeArcIntersect = verb.verb_geom_Intersect.curves(ellipseEdge, circleMinSurfaceRadius)
    #         edgeArcIntersect = verb.verb_eval_Intersect.curveAndSurface(ellipseEdge.asNurbs(), intersectPlane, tol=1e-3)
    #
    #         if len(edgeArcIntersect) > 0:
    #             intersectPoint = intersectPoint + edgeArcIntersect
    #             intersectedEdgeSet.append(pel)
    #
    #         # test if circle intersections lies in vertex defined arc
    #         # if edgeArcIntersect1 is not None:
    #         #     if pointInArc( edgeArcIntersect1,
    #         #                    v1, v2,
    #         #                    pel["axisDir"],
    #         #                    pel["refDir"],
    #         #                    pel["axisPoint"]):
    #         #         intersectedEdgeSet.append(pel)
    #         #         intersectPoint.append(edgeIntersect)
    #         #
    #         # if edgeArcIntersect2 is not None:
    #         #     if pointInArc( edgeArcIntersect2,
    #         #                    v1, v2,
    #         #                    pel["axisDir"],
    #         #                    pel["refDir"],
    #         #                    pel["axisPoint"]):
    #         #         if pel not in intersectedEdgeSet:
    #         #             intersectedEdgeSet.append(pel)
    #         #         intersectPoint.append(edgeIntersect)
    #
    #     if pel["typeName"] == "B_SPLINE_CURVE_WITH_KNOTS":
    #         # try intersecting spline curve with plane containing minPoint
    #         # find closestPoint() < eps
    #
    #         curveDegree = pel["curveDegree"]
    #         controlPointsList = pel["controlPointsList"]
    #         knotvector = pel["knotvector"]
    #
    #         NURBS_edge = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(curveDegree,
    #                                                                            knotvector,
    #                                                                            [cp.tolist() for cp in
    #                                                                             controlPointsList])
    #
    #         refDirP = np.cross(centroidAxisDir, axisDir)  # surface axisDir
    #         p1 = centroidProjAxisPoint + centroidAxisDir * minSurfaceDisp
    #         p3 = centroidProjAxisPoint - centroidAxisDir * minSurfaceDisp
    #         p2 = centroidProjAxisPoint + refDirP * minSurfaceDisp
    #         p4 = centroidProjAxisPoint - refDirP * minSurfaceDisp
    #
    #         #   Generate the control points, weights, and knots of a surface defined by 4 points
    #         #   first point in counter-clockwise form
    #         #   second point in counter-clockwise form
    #         #   third point in counter-clockwise form
    #         #   forth point in counter-clockwise form
    #
    #         #   returns: NurbsSurfaceData object
    #         #   fourPointSurface( p1 : Point, p2 : Point, p3 : Point, p4 : Point, degree : Int = 3 )
    #
    #         intersectPlane = verb.verb_eval_Make.fourPointSurface(p1, p2, p3, p4)
    #         edgeArcIntersect = verb_eval_Intersect.curveAndSurface(NURBS_edge, intersectPlane, tol=1e-3)
    #         #     //Get the intersection of a NURBS curve and a NURBS surface without an estimate
    #         #     //
    #         #     //**params**
    #         #     //
    #         #     //* NurbsCurveData
    #         #     //* NurbsSurfaceData
    #         #     //* tolerance for the curve intersection
    #         #     //
    #         #     //**returns**
    #         #     //
    #         #     //* array of CurveSurfaceIntersection objects
    #         #
    #         #     public static function curveAndSurface( curve : NurbsCurveData,
    #         #                                               surface : NurbsSurfaceData,
    #         #                                               tol : Float = 1e-3,
    #         #                                               crvBbTree : IBoundingBoxTree<NurbsCurveData> = null,
    #         #                                               srfBbTree : IBoundingBoxTree<NurbsSurfaceData> = null ) : Array<CurveSurfaceIntersection>  {
    #         #
    #         #         crvBbTree = crvBbTree != null ? crvBbTree : new LazyCurveBoundingBoxTree( curve );
    #         #         srfBbTree = srfBbTree != null ? srfBbTree : new LazySurfaceBoundingBoxTree( surface );
    #
    #         if len(edgeArcIntersect) > 0:
    #             intersectPoint = intersectPoint + edgeArcIntersect
    #             intersectedEdgeSet.append(pel)
    #
    # if len(intersectPoint) == 2 and len(intersectedEdgeSet) == 2:
    #     if pointInArc(centroidProjAxisPoint,
    #                   intersectPoint[0],
    #                   intersectPoint[1],
    #                   ParsedSurface["refDir"],
    #                   ParsedSurface["auxDir"],
    #                   ParsedSurface["normDir"],
    #                   ParsedSurface["axisPoint"],
    #                   rotSym=True):
    #         if rotSymDisp < eps_STEP_AP21:
    #             ParsedSurface["rotSymEdge"].append(centroidProjAxisPoint)  # -----------------------------------------
    #             # minima or maxima
    #         else:
    #             ParsedSurface["minPoint"] = minSurfaceDispPoint
    #         ParsedSurface["minPointCentroidDisp"] = np.linalg.norm(
    #             centroid - minSurfaceDispPoint)  # minPointCentroidDisp
    #
    # # if len(intersectPoint) > 2:
    #
    # # # can't determine in advance whether rotSymMax or rotSymCentre => requires subsequent comparison to adjacent nodes
    # # # same as maxPoint minPoint comparison

    AFS['ParsedSurface'] = ParsedSurface


def coneSurfaceParse(AFS, localCentroid):
    # ENTITY Conical_surface
    #   SUBTYPE OF (Surface);
    #   position : Axis_placement_3d;
    #   radius : positive_length_measure;
    #   semi_angle : plane_angle_measure;

    # position: an Axis_placement_3d that defines the location and orientation of the surface.
    # The axis of the Conical_surface passes through the location and is normal to the plane of ref_direction and axis.
    # radius: the radius of the circular curve of intersection between the cone and a plane perpendicular to the
    # axis of the cone passing through the location. This will have value 0.0 if the location is at the apex of the cone.
    # semi_angle: the cone semi-angle, this is the angle between the axis of the Conical_surface and the generating line.

    # note that normal is defined in apex->location direction
    # https://quaoar.su/files/standards/Standard%20-%202003%20-%20ISO%2010303-42.pdf p88

    ParsedSurface = {}

    if (AFS.get("maxPoint") is not None) and (AFS.get("minPoint") is not None): # recalculation of max/min
        axisPoint = AFS["axisPoint"]
        normDir = AFS["normDir"] # Z
        refDir = AFS["refDir"] # X
        auxDir = AFS["auxDir"]  # Y
    else:
        axisPoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        semiAngle = AFS["SurfaceParams"][-1]
        radius = AFS["SurfaceParams"][-2]

    ParsedSurface["radius"] = radius
    ParsedSurface["axisPoint"] = axisPoint
    ParsedSurface["normDir"] = normDir  # Z
    ParsedSurface["refDir"] = refDir  # X
    ParsedSurface["auxDir"] = auxDir

    ParsedSurface["semiAngle"] = semiAngle

    # test if cone is rotationally symmetrical, if localCentroid is close to axisDir through axisPoint
    rotSymDisp = np.linalg.norm(
        (localCentroid - axisPoint) - np.dot((localCentroid - axisPoint), normDir) * normDir
    )

    if rotSymDisp < eps_STEP_AP21:
        # ParsedSurface["rotSymMax"] = axisPoint at circle
        #minSurfaceDispPoint = localCentroid  #
        #centroidProjAxisPoint = centroid
        # maxima is circular edge centre, determined from vertices

        # edges should all have been calculated prior to surface
        # circle edge cannot be declared groove or ridge prior to discovery of adjacent surfaces
        # rotationally-symmetric edge must be assigned indeterminate status with centroid & radius

        edgeRadiusDisp = [edge['radius'] for edge in AFS["FaceOuterBoundEdgeLoopList"] if edge.get('radius') is not None]
        ParsedSurface["rotSymCentre"] = AFS["FaceOuterBoundEdgeLoopList"][edgeRadiusDisp.index(min(edgeRadiusDisp))]["axisPoint"]
        ParsedSurface['rotSymRadius'] = min(edgeRadiusDisp)
        ParsedSurface["minPoint"] = None  # a band rather than a point
        ParsedSurface["maxPoint"] = None
        # only rotSymCentre where there is no explicit edge
        # todo worthwhile finding maxima?

        #     for edge in AFS["FaceOuterBoundEdgeLoopList"]:
        #         if np.isclose(edge["minPointCentroidDisp"], radius, eps): # eps_STEP_AP21 ?? todo could this fail for small values?
        #             if ParsedSurface.get("rotSymRadius") is None:
        #                 #ParsedSurface["rotSymCentre"] = np.array([,, edge['minPoint'][3]])
        #                 ParsedSurface["rotSymRadius"] = radius
        #                 ParsedSurface["minPoint"] = None # a band rather than a point
        #
        #                 # project centroid to symmetric axis to find rotSymCentre
        #                 #centroidProjAxisPoint = normDir * np.dot((axisPoint - localCentroid), normDir) / np.dot(normDir, normDir)
        #                 centroidProjAxisPoint = pointProjectAxis(localCentroid, axisPoint, normDir)
        #                 ParsedSurface["rotSymCentre"] = centroidProjAxisPoint
        #                 # minPointCentroidDisp/maxPointCentroidDisp exist where centroid is on rotSym axis
        #                 ParsedSurface['minPointCentroidDisp'] = np.sqrt(np.linalg.norm(centroidProjAxisPoint - localCentroid)**2 + radius**2)


    else:
        # cylinder/conic maxima is a search of greatest disp of vertices from localCentroid
        # if localCentroid is outside full cylindrical surface,
        #   cylinder minima -> test for intersection with surface on orthogonal projection from axis to localCentroid
        # if localCentroid within surface, |centroid - projAxisDir| < cRadius

        # test:
        # axisPoint = np.array([1, 0, 0])
        # localCentroid = np.array([2, 2, 0])
        # axisDir = np.array([1.5, 0, 0])
        # cRadius = 1

        # # orthogonal projection of localCentroid to vector between conic apex and radius at axisPoint
        # # vCA = localCentroid - axisPoint
        # # vNA = axisPoint - normDir
        # # centroidProjAxisPoint = vNA*np.dot(vCA, vNA)/np.dot(vNA, vNA)
        #
        # centroidProjAxisPoint = pointProjectAxis(localCentroid, axisPoint, normDir)
        #
        # # centroidProjAxisPoint = axisPoint + (
        # #             np.dot((centroid - axisPoint), normDir) / np.dot(normDir, normDir)) * normDir
        # centroidAxisDir = centroidProjAxisPoint - localCentroid
        # centroidAxisDisp = np.linalg.norm(centroidAxisDir)
        # centroidAxisDir = centroidAxisDir / centroidAxisDisp
        #
        # if centroidAxisDisp > radius:
        #     surfaceMinPoint = localCentroid + (centroidAxisDisp - radius) * centroidAxisDir
        # elif centroidAxisDisp < radius:
        #     surfaceMinPoint = localCentroid - (radius - centroidAxisDisp) * centroidAxisDir
        # elif centroidAxisDisp - radius < eps:
        #     surfaceMinPoint = localCentroid
        #
        # # insideOutsideSurfaceTest() method to determine if surfaceMinPoint or surfaceMaxPoint within defined surface edges
        # # note the max point will be on one of the edges
        #
        # if not insideOutsideSurfaceTest(localCentroid, surfaceMinPoint, AFS):
        #     # find minima among edges
        #     minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
        #     ParsedSurface["minPoint"] = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]["minPoint"]
        #     ParsedSurface["minPointCentroidDisp"] = min(minPointCentroidDisp)
        # else:
        #     ParsedSurface["minPoint"] = surfaceMinPoint
        #     ParsedSurface["minPointCentroidDisp"] = np.linalg.norm(surfaceMinPoint - localCentroid)


        # get adjacent disp between axisPoint and cone/triangle apex
        axisPoint_apexPoint_disp = radius / np.tan(semiAngle)
        # find point at cone/triangle apex
        apexPoint = axisPoint + axisPoint_apexPoint_disp * (-normDir)

        # create a vector in the base plane that defines radius,
        # and in the plane of axisPoint, apex & localCentroid

        #centroidAxisPointDir = axisPoint - localCentroid
        #centroidAxisPointDir = centroidAxisPointDir / np.linalg.norm(centroidAxisPointDir)
        #coneEdgeDir = np.cross(centroidAxisPointDir, normDir)
        #centroidConeEdgeDir = np.cross(coneEdgeDir, normDir)
        # determine if centroid is closer to one side of cone
        oopDir = np.cross(localCentroid - axisPoint, -normDir)
        coneEdgeDir = np.cross(oopDir, -normDir)
        coneEdgeDir /= np.linalg.norm(coneEdgeDir)
        minConeEdgePoint = axisPoint + (radius * coneEdgeDir)

        # maxConeEdgePoint = axisPoint + radius * localCentroidConeEdgeDir
        # find minimum disp between localCentroid and segment defined between apexPoint, minConeEdgePoint

        minPoint, minPointCentroidDisp, _1, _2 = pointSegmentMinDisp(localCentroid, apexPoint, minConeEdgePoint)

        minNorm = (minPoint - localCentroid) / np.linalg.norm(minPoint - localCentroid) # same as maxNorm tho..

        if not insideOutsideSurfaceTest(localCentroid, minPoint, AFS):
            # find minima among edges
            minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
            ParsedSurface["minPoint"] = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]
            ParsedSurface["minPointCentroidDisp"] = min(minPointCentroidDisp)
        else:
            ParsedSurface["minPoint"] = minPoint
            ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp

        # find maxima among (arc/circle) edges
        maxPointCentroidDisp = [edge['maxPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"] if edge.get('maxPointCentroidDisp') is not None]
        print(maxPointCentroidDisp)
        ParsedSurface["maxPoint"] = AFS["FaceOuterBoundEdgeLoopList"][maxPointCentroidDisp.index(max(maxPointCentroidDisp))]["maxPoint"]
        ParsedSurface["maxPointCentroidDisp"] = max(maxPointCentroidDisp)

    AFS['ParsedSurface'] = ParsedSurface


def planeSurfaceParse(AFS, localCentroid):
    # ENTITY Plane
    #   SUBTYPE OF (Surface);
    #   position : Axis_placement_3d;

    ParsedSurface = {}

    if (AFS.get("maxPoint") is not None) and (AFS.get("minPoint") is not None): # recalculation of max/min
        axisPoint = AFS["axisPoint"]
        normDir = AFS["normDir"] # Z
        refDir = AFS["refDir"] # X
        #auxDir = AFS["auxDir"]  # Y
    else:
        axisPoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)

    ParsedSurface["axisPoint"] = axisPoint
    ParsedSurface["normDir"] = normDir  # Z
    ParsedSurface["refDir"] = refDir  # X

    # rotSymFeatureEdgeSet = []
    # rotSymFeature = []

    # circle segments, ellipse segments have to be evaluated independently
    # not sure how B-splines might work, max point box?
    # circleIndex = [i for i, fobell in enumerate(AFS['FaceOuterBoundEdgeLoopList']) if
    #                'CIRCLE' in fobell['typeName'] or 'CIRCLE' in fobell['subTypeName']]
    circleIndex = [i for i, fobell in enumerate(AFS['FaceOuterBoundEdgeLoopList']) if
                   'CIRCLE' in fobell['typeName']]
    # ellipseIndex = [i for i, fobell in enumerate(AFS['FaceOuterBoundEdgeLoopList']) if
    #                'ELLIPSE' in fobell['typeName'] or 'ELLIPSE' in fobell['subTypeName']]
    ellipseIndex = [i for i, fobell in enumerate(AFS['FaceOuterBoundEdgeLoopList']) if
                   'ELLIPSE' in fobell['typeName']]

    # get list of vertices wrt surface
    #surfaceTypes = []
    surfaceVertices = [fobell['vertex1'] for fobell in AFS['FaceOuterBoundEdgeLoopList']]
        #surfaceTypes.append(fobell['typeName'])

    if len(surfaceVertices) > 2:
        # strategy for planar polygonal surfaces is to divide into non-overlapping triangles
        # where corners are polygon vertices and test for nearest/furthest point on triangle
        minPoint, minPointCentroidDisp = planeMinMaxPoint(centroid, surfaceVertices, minmax="min")
        maxPoint, maxPointCentroidDisp = planeMinMaxPoint(centroid, surfaceVertices, minmax="max")

        # not required to test within edges

    if len(circleIndex) > 0:
        localCircleMinP = []
        for ci in circleIndex:
            fobell = AFS['FaceOuterBoundEdgeLoopList'][ci]
            cMinPoint, cMaxPoint = pointCircleMinMaxDisp(centroid,
                                                 fobell['axisPoint'],
                                                 fobell['normDir'],
                                                 fobell['radius'],
                                                 interior=True)
            if fobell["vertex1ref"] == fobell["vertex2ref"]: # full circle
                localCircleMinP.append(cMinPoint)
            else:
                boxPoints = curveEnclosingRectangle(fobell['vertex1'],
                                                    fobell['vertex2'],
                                                    fobell['axisPoint'],
                                                    fobell['normDir'],
                                                    fobell['radius'])
                rMinPoint, _ = planeMinMaxPoint(centroid, boxPoints, minmax="min")

                # check circle-centre, circle maxima, box maxima are collinear, and BM-CC > CM-CC
                # cross product of (BM-CC) and (CM-CC) = 0 => points BM, CM and CC collinear.
                # CC < CM < BM => dot product of (BM-CC) and (CM-CC) is positive
                # and is less than the square of the distance between a and b.
                CM_CC = cMinPoint - fobell['axisPoint']
                BM_CC = rMinPoint - fobell['axisPoint']
                if (np.linalg.norm(BM_CC) < eps_STEP_AP21) and (np.linalg.norm(CM_CC) < eps_STEP_AP21):
                    # coincident case
                    localCircleMinP.append(cMinPoint)
                elif np.linalg.norm(np.cross(BM_CC, CM_CC)) < eps_STEP_AP21: #=tol
                    if 0 < np.dot(BM_CC, CM_CC) < (BM_CC*BM_CC):
                        localCircleMinP.append(cMinPoint)

        if len(ellipseIndex) > 0:
            localEllipseMinP = []
            for ei in ellipseIndex:
                fobell = AFS['FaceOuterBoundEdgeLoopList'][ei]
                eMinPoint, eMaxPoint = pointEllipseMinMaxDisp(centroid,
                                                              fobell['refDir'],
                                                              fobell['auxDir'],
                                                              fobell['normDir'],
                                                              fobell['majorRadius'],
                                                              fobell['minorRadius'],
                                                              interior=True)

                if fobell["vertex1ref"] == fobell["vertex2ref"]:  # full circle
                    localEllipseMinP.append(cMinPoint)
                else:
                    boxPoints = curveEnclosingRectangle(fobell['vertex1'],
                                                        fobell['vertex2'],
                                                        fobell['axisPoint'],
                                                        fobell['normDir'],
                                                        fobell['radius'])
                    rMinPoint, rMinPointDisp = planeMinMaxPoint(centroid, boxPoints, minmax="min")

                    EM_CC = eMinPoint - fobell['axisPoint']
                    BM_CC = rMinPoint - fobell['axisPoint']
                    if np.cross(BM_CC, EM_CC) < eps_STEP_AP21:  # =tol
                        if 0 < np.dot(BM_CC, EM_CC) < (BM_CC * BM_CC):
                            localEllipseMinP.append(eMinPoint)

        if len(circleIndex) > 0: # extract closest point to localCentroid from contenders
            localCircleMinDisp = [np.linalg.norm(lcmp - localCentroid) for lcmp in localCircleMinP]
            localCircleMinP = localCircleMinP[localCircleMinDisp.index(min(localCircleMinDisp))]
            localCircleMinDisp = min(localCircleMinDisp)
            if len(surfaceVertices) > 2:
                if localCircleMinDisp < minPointDisp:
                    minPoint = localCircleMinP
                    minPointCentroidDisp = localCircleMinDisp
            else:
                minPoint = localCircleMinP
                minPointCentroidDisp = localCircleMinDisp

        if len(ellipseIndex) > 0:
            localEllipseMinDisp = [np.linalg.norm(lemp - localCentroid) for lemp in localEllipseMinP]
            localEllipseMinP = localEllipseMinP[localEllipseMinDisp.index(min(localEllipseMinDisp))]
            localEllipseMinDisp = min(localEllipseMinDisp)
            if len(surfaceVertices) > 2:
                if localEllipseMinDisp < minPointDisp:
                    minPoint = localEllipseMinP
                    minPointDisp = localEllipseMinDisp
            else:
                minPoint = localEllipseMinP
                minPointDisp = localEllipseMinDisp

    # note local minPoint of planar surface is generally a global minPoint
    ParsedSurface["minPoint"] = minPoint
    ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp

    # maxPoint is determined from vertex and edge maxima, exception for rotSymCentre
    maxPoint = [fobell['maxPoint'] for fobell in AFS['FaceOuterBoundEdgeLoopList'] if 'maxPoint' in fobell.keys()]
    maxPointCentroidDisp = [fobell['maxPointCentroidDisp'] for fobell in AFS['FaceOuterBoundEdgeLoopList']if 'maxPointCentroidDisp' in fobell.keys()]
    if len(maxPoint) > 0 and len(maxPointCentroidDisp) > 0:
        maxPoint = maxPoint[maxPointCentroidDisp.index(max(maxPointCentroidDisp))]
        maxPointCentroidDisp = max(maxPointCentroidDisp)
        ParsedSurface["maxPoint"] = maxPoint
        ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp

    rotSymEdges = [fobell['rotSymCentre'] for fobell in AFS['FaceOuterBoundEdgeLoopList'] if 'rotSymCentre' in fobell.keys()]
    if len(rotSymEdges) == 1:
        ParsedSurface["rotSymCentre"] = rotSymEdges[0]
    elif len(rotSymEdges) > 1:
        rsm_x = [rsm[0] for rsm in rotSymEdges]
        rsm_y = [rsm[1] for rsm in rotSymEdges]
        rsm_z = [rsm[2] for rsm in rotSymEdges]
        if (all(x == rsm_x[0] for x in rsm_x) and
                all(y == rsm_y[0] for y in rsm_y) and
                all(z == rsm_z[0] for z in rsm_z)):
            ParsedSurface["rotSymCentre"] = rotSymEdges[0]
        else:
            print("indeterminate rotSym")

    # TEMPORARY TEST--------------------------------------------------------------------------------
    # if maxPoint is None: # circle centre orthogonal to centroid, maxPoint not defined
    #     _1=1

    ParsedSurface["maxPoint"] = maxPoint
    ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp
    if maxPoint is not None: # circle centre orthogonal to centroid, maxPoint not defined
        if not insideOutsideSurfaceTest(localCentroid, maxPoint, AFS):
            # find maxima among edges
            maxPointCentroidDisp = [edge['maxPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
            maxPoint = AFS["FaceOuterBoundEdgeLoopList"][maxPointCentroidDisp.index(max(maxPointCentroidDisp))]
            maxPointCentroidDisp = max(maxPointCentroidDisp)
        # else:
        #     ParsedSurface["maxPoint"] = maxPoint
        #     ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp

    ParsedSurface["minPoint"] = minPoint
    ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp
    if not insideOutsideSurfaceTest(localCentroid, minPoint, AFS):
        # find maxima among edges
        minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
        minPoint = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]
        minPointCentroidDisp = min(minPointCentroidDisp)

    AFS['ParsedSurface'] = ParsedSurface


def primitiveSurfaceParse(AFS, c):
    if AFS["SurfaceTypeName"] not in [
        "SPHERICAL_SURFACE",
        "TOROIDAL_SURFACE",
        "CYLINDRICAL_SURFACE",
        "CONICAL_SURFACE",
        "PLANE"
    ]:
        print("Untreated surface type: " + AFS["SurfaceTypeName"])
    else:
        if AFS["SurfaceTypeName"] == "SPHERICAL_SURFACE":
            sphereSurfaceParse(AFS, c) #ParsedSurface

        if AFS["SurfaceTypeName"] == "TOROIDAL_SURFACE":
            toroidSurfaceParse(AFS, c)

        if (AFS["SurfaceTypeName"] == "CYLINDRICAL_SURFACE"):
            cylinderSurfaceParse(AFS, c)

        if (AFS["SurfaceTypeName"] == "CONICAL_SURFACE"): # still must test with simple STEP file
            coneSurfaceParse(AFS, c)

        if AFS["SurfaceTypeName"] == "PLANE":
            planeSurfaceParse(AFS, c)



        # if vp.type_name == "VERTEX_POINT":
        #     cpRef = ref2index(vp.params[-1])
        #     if STEP_entities[cpRef].type_name == "CARTESIAN_POINT":
        #         cp = Point(
        #             STEP_entities[cpRef].params[-1][0],
        #             STEP_entities[cpRef].params[-1][1],
        #             STEP_entities[cpRef].params[-1][1],
        #         )
        #        vertexPoints.append(cp)


AdvancedFaceSurfaces = []
for se in STEP_entities:
    if hasattr(se, "type_name"):  # stepcode.Part21.ComplexEntity type missing attribute
        if (se.type_name == "ADVANCED_FACE") or (se.type_name == "FACE_SURFACE"):
            SurfaceNormalOutwards = 'T' in se.params[-1]
            afRefs = cleanSubRefs(se.params)
            SurfaceClass = {
                "SurfaceNormalOutwards": SurfaceNormalOutwards,
                "SurfaceTypeName": None,
                'SurfaceRef': None,
                "SurfaceParams": None,
                "EdgeLoopList": [],
            }
            for se2 in afRefs:
                if hasattr(STEP_entities[ref2index(se2)], "type_name"): # should prolly assume [0] is FACE_*-----------------<<<<<<
                    afTypeName = STEP_entities[ref2index(se2)].type_name
                    if (afTypeName == "FACE_OUTER_BOUND") or (afTypeName == "FACE_BOUND"):  # no risk of complex entities from here on in
                        # todo assume 'FACE_INNER_BOUND' or 'HOLE' etc can run the same edge loop from here ----------------------<<<<<<
                        se3 = STEP_entities[ref2index(se2)].params
                        se3 = cleanSubRefs(se3)[0]
                        if STEP_entities[ref2index(se3)].type_name == "EDGE_LOOP":
                            SurfaceClass["EdgeLoopList"] = STEP_entities[ref2index(se3)].params[-1]
                    else:
                        # print(afTypeName)
                        SurfaceClass["SurfaceTypeName"] = afTypeName
                        SurfaceClass["SurfaceRef"] = se2
                        SurfaceClass["SurfaceParams"] = STEP_entities[
                            ref2index(se2)
                        ].params

                elif not hasattr(STEP_entities[ref2index(se2)], "type_name"):
                    print("complex entity @todo")
            AdvancedFaceSurfaces.append(SurfaceClass)

# axis: the Direction that defines the second axis of the Axis_placement. The value of this attribute need not be specified.
# ref_direction: the direction used to determine the direction of the local X axis.
# The value of this attribute need not be specified. If axis or ref_direction is omitted, these directions are taken from the geometric coordinate system

# orientation: a BOOLEAN flag. If TRUE, the topological orientation as used coincides with the orientation,
# from start vertex to end vertex, of the edge_definition, if FALSE the vertices are reversed in order.
# edge_start: the start vertex of the Oriented_edge. This is derived from the vertices of the edge_definition after taking account of the orientation.
# edge_end: the end vertex of the Oriented_edge. This is derived from the vertices of the edge_definition after taking account of the orientation.

# extract edge data and convert into useful format

#
#10 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#15),#113);
#11 = AXIS2_PLACEMENT_3D('',#12,#13,#14);
#12 = CARTESIAN_POINT('',(0.,0.,0.));
#13 = DIRECTION('',(0.,0.,1.));
#14 = DIRECTION('',(1.,0.,-0.));

outermostPoints = []
# need to exclude Advanced_brep_shape_representation which seems to contain a [0,0,0] origin point
ABSRpointRef = None
for ABSR in STEP_entities:
    if hasattr(ABSR, "type_name"):  # stepcode.Part21.ComplexEntity type missing attribute
        if ABSR.type_name == 'ADVANCED_BREP_SHAPE_REPRESENTATION':
            # get
            cleanABSR = cleanSubRefs(ABSR.params)
            if STEP_entities[ref2index(cleanABSR[0])].type_name == 'AXIS2_PLACEMENT_3D':
                ABSRpointRef = cleanSubRefs(STEP_entities[ref2index(cleanABSR[0])].params)[0]
                #_1=1


    for vp in STEP_entities:
        if hasattr(vp, "type_name"):  # stepcode.Part21.ComplexEntity type missing attribute
            if vp.type_name == "CARTESIAN_POINT":
                if (len(vp.params[-1]) == 3):
                    if (vp.ref == ABSRpointRef):
                        pass
                    else:
                        outermostPoints.append(np.array([vp.params[-1][0], vp.params[-1][1], vp.params[-1][2]]))



            #
#
centroid = medianPoint(outermostPoints)
# continue to ignore local maxima at both radial extents and NURB surface maxima?
# for all elements, surfaces and bounding edges, iterate maxima discovery relative to centroid


for AFS in AdvancedFaceSurfaces:
    EdgeLoops = AFS["EdgeLoopList"]
    # testing a specific edge loop instance 44CONE 0SPLINE
    # EdgeLoops = AdvancedFaceSurfaces[0]["EdgeLoopList"]

    parsedEdgeLoop = []
    for edgeInstance in EdgeLoops:

        parsedEdge = commonEdgeParse(edgeInstance, STEP_entities)
        simpleEdgeParse(centroid, parsedEdge, STEP_entities)
        parsedEdgeLoop.append(parsedEdge)

        if parsedEdge["typeName"] in ["SEAM_CURVE", "SURFACE_CURVE"]:
            #   SEAM_CURVE (not sure if this is worthwhile retaining)
            #   associated geometry in following list does not seem to be defined elsewhere

            #   Attribute	                Type	                                        Defined By
            #   name	                    label (STRING)	                                representation_item
            #   curve_3d	                curve (ENTITY)	                                surface_curve
            #   associated_geometry	        LIST OF pcurve_or_surface (SELECT)	            surface_curve
            #   master_representation	    preferred_surface_curve_representation (ENUM)	surface_curve

            cleanEdgeType = cleanSubRefs(parsedEdge["edgeType"])
            edgeTypeName = STEP_entities[ref2index(cleanEdgeType[0])].type_name
            # if subEdgeTypeName in [
            #                             "LINE",
            #                             "CIRCLE",
            #                             "ELLIPSE",
            #                             "B_SPLINE_CURVE_WITH_KNOTS",
            #                         ]:
            parsedEdge["superTypeName"] = parsedEdge["typeName"]

            edgeRef = STEP_entities[ref2index(cleanEdgeType[0])].ref
            parsedEdge["superEdgeRef"] = parsedEdge["edgeRef"]
            parsedEdge["edgeRef"] = edgeRef
            parsedEdge["superEdgeType"] = parsedEdge["edgeType"]
            parsedEdge["edgeType"] = STEP_entities[ref2index(cleanEdgeType[0])].params # todo: change edgeParams
            parsedEdge["superTypeName"] = parsedEdge["typeName"]
            parsedEdge["typeName"] = STEP_entities[ref2index(edgeRef)].type_name
            simpleEdgeParse(centroid, parsedEdge, STEP_entities)

            if parsedEdge["edgeRef"] not in [pel["edgeRef"] for pel in parsedEdgeLoop]:
                parsedEdgeLoop.append(parsedEdge)

            for cet in cleanEdgeType[1:]: # get associated surfaces, LIST OF pcurve_or_surface
                if STEP_entities[ref2index(cet)].type_name == "PCURVE":
                    # extract surface data from second parameter term
                    surfaceRef = STEP_entities[ref2index(cet)].params[1]

                    # todo:  shouldn't recalculate on every pass?
                    # todo:  create function collateField(surface?) to .. dict data

                    # surfaceRefCollection = []
                    # for AFS2 in AdvancedFaceSurfaces:
                    #     surfaceRefCollection.append(AFS2["SurfaceRef"])
                    # if surfaceRef not in surfaceRefCollection:
                    if surfaceRef not in [AFS2["SurfaceRef"] for AFS2 in AdvancedFaceSurfaces]:
                        _1=1
                        print("undetected surface! PCURVE")
                else:
                    # todo: code for other surface option missing
                    print("undetermined surfacetype failure line2824")
                    _1 = 1

    AFS["FaceOuterBoundEdgeLoopList"] = parsedEdgeLoop
    primitiveSurfaceParse(AFS, centroid)

    # The standard CSG primitives are the cone, eccentric_cone, cylinder, sphere, torus,
    # block, right_angular_wedge, ellipsoid, tetrahedron and pyramid.

    #        b_spline_surface,
    #        b_spline_surface_with_knots,
    #        bezier_surface,
    #        conical_surface,
    #        curve_bounded_surface,
    #        cylindrical_surface,
    #        degenerate_toroidal_surface,
    #        offset_surface,
    #        quasi_uniform_surface,
    #        rational_b_spline_surface,
    #        rectangular_composite_surface,
    #        rectangular_trimmed_surface,
    #        spherical_surface,
    #        surface,
    #        surface_of_linear_extrusion,
    #        surface_of_revolution,
    #        surface_replica,
    #        swept_surface,
    #        toroidal_surface,
    #        uniform_surface,

    # #35797 =( BOUNDED_SURFACE ( )  B_SPLINE_SURFACE ( 3, 1, (
    # ( #6498, #18391 ),
    #     ( #10581, #38916 ),
    #         ( #14677, #43009 ),
    #             ( #18762, #47098 ) ),
    #                 .UNSPECIFIED., .F., .F., .F. )
    #             B_SPLINE_SURFACE_WITH_KNOTS ( ( 4, 4 ),
    #         ( 2, 2 ),
    #         ( 0.0000000000000000000, 1.000000000000000000 ),
    #         ( 0.0000000000000000000, 1.000000000000000000 ),
    #         .UNSPECIFIED. )
    #     GEOMETRIC_REPRESENTATION_ITEM ( )  RATIONAL_B_SPLINE_SURFACE ( (
    #     ( 1.000000000000000000, 1.000000000000000000),
    #     ( 0.3333333333333333700, 0.3333333333333333700),
    #     ( 0.3333333333333333700, 0.3333333333333333700),
    #     ( 1.000000000000000000, 1.000000000000000000) ) )
    # REPRESENTATION_ITEM ( '' )  SURFACE ( )  );

# iteratively calculate max/min points and median centroid point
# extract max values for all surfaces to add to median vertex calculation for centroid
outermostPoints = []
lastCentroid = centroid
for AFS in AdvancedFaceSurfaces:
    if AFS["ParsedSurface"].get('maxPoint') is not None:
        outermostPoints.append(AFS["ParsedSurface"]["maxPoint"])
        # print("AFS[ParsedSurface][maxPoint]")
        # print(AFS["ParsedSurface"]["maxPoint"])
    if AFS["ParsedSurface"].get('minPoint') is not None:
        outermostPoints.append(AFS["ParsedSurface"]["minPoint"])
        # print("AFS[ParsedSurface][minPoint]")
        # print(AFS["ParsedSurface"]["minPoint"])
    if AFS["ParsedSurface"].get('axisPoint') is not None:
        outermostPoints.append(AFS["ParsedSurface"]["axisPoint"])
        # print("AFS[ParsedSurface][axisPoint]")
        # print(AFS["ParsedSurface"]["axisPoint"])
    for edge in AFS['FaceOuterBoundEdgeLoopList']:
        if edge.get('maxPoint') is not None:
            outermostPoints.append(edge["maxPoint"])
            # print("maxPoint")
            # print(edge["maxPoint"])
        if edge.get('minPoint') is not None:
            outermostPoints.append(edge["minPoint"])
            # print("minPoint")
            # print(edge["minPoint"])
        if edge.get('axisPoint') is not None:
            outermostPoints.append(edge["axisPoint"])
            # print("axisPoint")
            # print(edge["axisPoint"])

centroid = medianPoint(outermostPoints)

# iterate centroid-features relationship
while (np.linalg.norm(lastCentroid - centroid) > eps):
    # recalculate surfaces & edges max/min points according to new centroid
    for AFS in AdvancedFaceSurfaces:
        # update surface after edges
        for edge in AFS['FaceOuterBoundEdgeLoopList']:
            simpleEdgeParse(centroid, edge, STEP_entities)
        primitiveSurfaceParse(AFS, centroid)

    lastCentroid = centroid
    for AFS in AdvancedFaceSurfaces:
        if AFS["ParsedSurface"]["maxPoint"] is not None:
            outermostPoints.append(AFS["ParsedSurface"]["maxPoint"])
            for edge in AFS['FaceOuterBoundEdgeLoopList']:
                outermostPoints.append(edge["maxPoint"])

    centroid = medianPoint(outermostPoints)
    print("centroid update: " + str(np.linalg.norm(lastCentroid - centroid)))

# todo: weed out surfaces that only touch endpoints? May be already true

# find adjoining surfaces to every surface edge
for AFS in AdvancedFaceSurfaces:
    for edgeSource in AFS['FaceOuterBoundEdgeLoopList']:
        edgeSource['edgeAdjSurfacesIndex'] = {}
        for AFS2index, AFS2 in enumerate(AdvancedFaceSurfaces):
            for edgeTarget in AFS2['FaceOuterBoundEdgeLoopList']:
                if edgeSource['edgeRef'] == edgeTarget['edgeRef']:
                    if AFS2['SurfaceRef'] not in edgeSource['edgeAdjSurfacesIndex'].keys():
                        edgeSource['edgeAdjSurfacesIndex'][AFS2['SurfaceRef']] = AFS2index

# find adjoining surfaces to every vertex 1 surface point
for AFS in AdvancedFaceSurfaces:
    for edgeSource in AFS['FaceOuterBoundEdgeLoopList']:
        v1ref = edgeSource['vertex1ref']
        edgeSource['vertex1SurfacesIndex'] = {}
        for AFS2index, AFS2 in enumerate(AdvancedFaceSurfaces):
            for edgeTarget in AFS2['FaceOuterBoundEdgeLoopList']:
                if (v1ref == edgeTarget['vertex1ref']) or (v1ref == edgeTarget['vertex2ref']):
                    if AFS2index not in edgeSource['vertex1SurfacesIndex'].keys():
                        edgeSource['vertex1SurfacesIndex'][AFS2['SurfaceRef']] = AFS2index

# find adjoining surfaces to every vertex 2 surface point
for AFS in AdvancedFaceSurfaces:
    for edgeSource in AFS['FaceOuterBoundEdgeLoopList']:
        v2ref = edgeSource['vertex2ref']
        edgeSource['vertex2SurfacesIndex'] = {}
        for AFS2index, AFS2 in enumerate(AdvancedFaceSurfaces):
            for edgeTarget in AFS2['FaceOuterBoundEdgeLoopList']:
                if (v1ref == edgeTarget['vertex1ref']) or (v1ref == edgeTarget['vertex2ref']):
                    if AFS2index not in edgeSource['vertex2SurfacesIndex'].keys():
                        edgeSource['vertex2SurfacesIndex'][AFS2['SurfaceRef']] = AFS2index

# find adjoining edges to every vertex 1 surface point
for AFS in AdvancedFaceSurfaces:
    for edgeSource in AFS['FaceOuterBoundEdgeLoopList']:
        v1ref = edgeSource['vertex1ref']
        edgeSource['vertex1edgesIndex'] = {}
        for AFS2index, AFS2 in enumerate(AdvancedFaceSurfaces):
            for edgeTargetIndex, edgeTarget in enumerate(AFS2['FaceOuterBoundEdgeLoopList']):
                if (v1ref == edgeTarget['vertex1ref']) or (v1ref == edgeTarget['vertex2ref']):
                    if edgeTarget['edgeRef'] not in edgeSource['vertex1edgesIndex'].keys():
                        edgeSource['vertex1edgesIndex'][edgeTarget['edgeRef']] = (AFS2index, edgeTargetIndex)

# find adjoining edges to every vertex 2 surface point
for AFS in AdvancedFaceSurfaces:
    for edgeSource in AFS['FaceOuterBoundEdgeLoopList']:
        v2ref = edgeSource['vertex2ref']
        edgeSource['vertex2edgesIndex'] = {}
        for AFS2index, AFS2 in enumerate(AdvancedFaceSurfaces):
            for edgeTargetIndex, edgeTarget in enumerate(AFS2['FaceOuterBoundEdgeLoopList']):
                if (v2ref == edgeTarget['vertex1ref']) or (v2ref == edgeTarget['vertex2ref']):
                    if edgeTarget['edgeRef'] not in edgeSource['vertex2edgesIndex'].keys():
                        edgeSource['vertex2edgesIndex'][edgeTarget['edgeRef']] = (AFS2index, edgeTargetIndex)

# rank rotationally symmetric edges as grooves/spherical/ridges based on adjoining surfaces
# namely for the 2 surfaces adjoining an arc/circle, measure the displacement of the surface relative to this axis of rotation
# project each point on adjoining surfaces to rotational axis to determine average displacement
# how to determine the extent of relevant points and edges on each side? Limit to adjoining surface.
# (not check edges completely enclose rotational axis)
# average the displacement to this axis of vertices & edges
for AFSindex, AFS in enumerate(AdvancedFaceSurfaces):
    # commence with searching local minima around edges associated with AFS surface
    # no real reason to separate vertices into points closer to localCentroid than arc/circle, and those on the far side.
    for edge in AFS['FaceOuterBoundEdgeLoopList']:
        #
        if (edge.get('rotSymCentre') is not None) and (edge.get('rotSymRadius') is not None):
            adjSurfaceSet = edge['edgeAdjSurfacesIndex'].values()
            # presuming it is not required to identify edges

            # normalise axis normal
            rotSymCentroidDisp = np.linalg.norm(edge['rotSymCentre'] - centroid)
            if rotSymCentroidDisp > eps:
                rotAxis = (edge['rotSymCentre'] - centroid) / rotSymCentroidDisp
            else:
                # if centre to rotSym feature coincides with centroid, get rotAxis from other parsed data
                if AFS['ParsedSurface'].get('normDir') is not None:
                    rotAxis = AFS['ParsedSurface']['normDir']
                else:
                    print("no surface axisDir value, cross()?")

            # adjacent surfaces minima,  nearPointAxisDisp > rotSymCentroidDisp
            nearPointAxisDisp = []
            farPointAxisDisp = []
            #ddump=[]

            for adjSurfaceIndex in adjSurfaceSet:
                adjSurface = AdvancedFaceSurfaces[adjSurfaceIndex]
                for adjSurfaceEdge in adjSurface['FaceOuterBoundEdgeLoopList']:
                    if edge['edgeRef'] != adjSurfaceEdge['edgeRef']:

                        # project all points, vertex1, vertex2, maxPoint, minPoint to rotational axis
                        adjSurfaceVertices = [adjSurfaceEdge['vertex1'], adjSurfaceEdge['vertex2']]
                        if not (adjSurfaceEdge['vertex1extremaMax'] and adjSurfaceEdge['vertex2extremaMax']):  # also add maxPoint
                            if adjSurfaceEdge['maxPoint'] is not None:
                                adjSurfaceVertices.append(adjSurfaceEdge['maxPoint'])
                        if not (adjSurfaceEdge['vertex1extremaMin'] and adjSurfaceEdge['vertex2extremaMin']):  # also add minPoint
                            if adjSurfaceEdge['minPoint'] is not None:
                                adjSurfaceVertices.append(adjSurfaceEdge['minPoint'])

                        for asv in adjSurfaceVertices:
                            #pAxisProj = pointProjectAxis(asv, centroid, rotAxis)
                            #pAxisProj = rotAxis * np.dot((asv - centroid), rotAxis) / np.dot(rotAxis, rotAxis)

                            if rotSymCentroidDisp > eps:
                                pAxisProj = rotAxis * np.dot((asv - centroid), rotAxis) / np.dot(rotAxis, rotAxis)
                            else:
                                pAxisProj = rotAxis * np.dot(asv, rotAxis) / np.dot(rotAxis, rotAxis)

                            pAxisDisp = np.linalg.norm(asv - pAxisProj)
                            #ddump.append(np.dot((asv - edge['rotSymCentre']), rotAxis))

                            rotSymEdgeDisp = np.dot((asv - edge['rotSymCentre']), rotAxis)
                            #print(rotSymEdgeDisp)
                            if rotSymEdgeDisp > eps: # eps_
                                #if np.linalg.norm(pAxisProj - centroid) > rotSymCentroidDisp:
                                farPointAxisDisp.append(pAxisDisp)
                            elif rotSymEdgeDisp < -eps:
                                nearPointAxisDisp.append(pAxisDisp)

            if len(nearPointAxisDisp)==0 or len(farPointAxisDisp)==0:
                # disc at end of object, defined as rotSymMax
                edge['rotSymMax'] = True
            else:
                NPADA = np.average(nearPointAxisDisp)
                FPADA = np.average(farPointAxisDisp)
                # if (NPADA >= edge['radius']) and (FPADA >= edge['radius']):
                #     # groove
                #     edge['rotSymMin'] = True
                #
                # if (NPADA < edge['radius']) and (FPADA < edge['radius']):
                #     # ridge
                #     edge['rotSymMax'] = True

                # conical midpoints
                avgRadius = np.abs(NPADA - FPADA)/2 + min([NPADA, FPADA])
                if avgRadius > edge['radius']:
                    # ridge
                    edge['rotSymMax'] = True
                else:
                    edge['rotSymMin'] = True
        # else:
        #     edge['rotSymCentre'] = None
        #     edge['rotSymRadius'] = None

    # the other class of minima integral to surfaces, e.g. ring minima on cylinders
    if (AFS['ParsedSurface'].get('rotSymCentre') is not None) and (AFS['ParsedSurface'].get('rotSymRadius') is not None):

        # normalise axis normal
        rotSymCentroidDisp = np.linalg.norm(AFS['ParsedSurface']['rotSymCentre'] - centroid)
        if rotSymCentroidDisp > eps:
            rotAxis = (AFS['ParsedSurface']['rotSymCentre'] - centroid) / rotSymCentroidDisp
        else:
            # if centre to rotSym feature coincides with centroid, get rotAxis from other parsed data
            if AFS['ParsedSurface'].get('normDir') is not None:
                rotAxis = AFS['ParsedSurface']['normDir']
            else:
                print("no surface axisDir value, cross()?")

        # todo what happens with several rotationally symmetric features within a surface?
        # todo list check?

        # adjacent surfaces minima,  nearPointAxisDisp > rotSymCentroidDisp
        nearPointAxisDisp = []
        farPointAxisDisp = []

        # use edges in surface FaceOuterBoundEdgeLoopList
        for adjSurfaceEdge in AFS['FaceOuterBoundEdgeLoopList']:
            # project all points, vertex1, vertex2, maxPoint, minPoint to rotational axis
            adjSurfaceVertices.append(adjSurfaceEdge['vertex1'])
            adjSurfaceVertices.append(adjSurfaceEdge['vertex2'])
            if not (adjSurfaceEdge['vertex1extremaMax'] and adjSurfaceEdge[
                'vertex2extremaMax']):  # also add maxPoint
                if adjSurfaceEdge['maxPoint'] is not None:
                    adjSurfaceVertices.append(adjSurfaceEdge['maxPoint'])
            if not (adjSurfaceEdge['vertex1extremaMin'] and adjSurfaceEdge[
                'vertex2extremaMin']):  # also add minPoint
                if adjSurfaceEdge['minPoint'] is not None:
                    adjSurfaceVertices.append(adjSurfaceEdge['minPoint'])

            for asv in adjSurfaceVertices:
                #pAxisProj = pointProjectAxis(asv, centroid, rotAxis)

                #pAxisProj = rotAxis * np.dot((asv - centroid), rotAxis) / np.dot(rotAxis, rotAxis)
                #pAxisProj = rotAxis * np.dot(asv, rotAxis) / np.dot(rotAxis, rotAxis)

                if rotSymCentroidDisp > eps:
                    pAxisProj = rotAxis * np.dot((asv - centroid), rotAxis) / np.dot(rotAxis, rotAxis)
                else:
                    pAxisProj = rotAxis * np.dot(asv, rotAxis) / np.dot(rotAxis, rotAxis)

                pAxisDisp = np.linalg.norm(asv - pAxisProj)
                rotSymEdgeDisp = np.dot((asv - AFS['ParsedSurface']['rotSymCentre']), rotAxis)
                if rotSymEdgeDisp > eps:  # eps_
                    # if np.linalg.norm(pAxisProj - centroid) > rotSymCentroidDisp:
                    farPointAxisDisp.append(pAxisDisp)
                elif rotSymEdgeDisp < -eps:
                    nearPointAxisDisp.append(pAxisDisp)

            if len(nearPointAxisDisp) == 0 or len(farPointAxisDisp) == 0:
                # disc at end of object, defined as rotSymMax
                AFS['ParsedSurface']['rotSymMax'] = True
            else:
                NPADA = np.average(nearPointAxisDisp)
                FPADA = np.average(farPointAxisDisp)
                # if (NPADA >= edge['radius']) and (FPADA >= edge['radius']):
                #     # groove
                #     edge['rotSymMin'] = True
                #
                # if (NPADA < edge['radius']) and (FPADA < edge['radius']):
                #     # ridge
                #     edge['rotSymMax'] = True

                # conical midpoints
                avgRadius = np.abs(NPADA - FPADA) / 2 + min([NPADA, FPADA])
                if avgRadius > AFS['ParsedSurface']['rotSymRadius']:
                    # groove
                    AFS['ParsedSurface']['rotSymMax'] = True
                else:
                    AFS['ParsedSurface']['rotSymMin'] = True

# for every local minima or maxima feature point, find the surrounding local minima points by way of geometrical features
# (i.e. opposing minima on a thin plate don't count; adjacency is not calculated by cartesian distance but surface topography)

# create similar dict structure listing local minima/maxima tuple (index1, index2) key and index list of surrounding local maxima/minima
# ignore overall minima within a surface and edge boundaries

# surface maxima/minima may be associated with related edge maxima/minima (note B-spline surface may have several max/min)
# find the minima/maxima associated with every edge

def array3x1Match(a, A):
    T = [np.isclose(a[0], b[0], atol=eps) and np.isclose(a[1], b[1], atol=eps) and np.isclose(a[2], b[2], atol=eps) for b in A]
    return any(T)

surfaceMinPoint = []
surfaceMinFeatureCentroidDisp = []
surfaceMinRotSym = []

surfaceMaxPoint = []
surfaceMaxFeatureCentroidDisp = []
surfaceMaxRotSym = []

for AFSindex, AFS in enumerate(AdvancedFaceSurfaces):
    if AFS["SurfaceTypeName"] == "SPHERICAL_SURFACE":
        # edge and surface point only equidistant from origin if centre coincides with origin
        _1=1
        print("write code for SPHERICAL_SURFACE")
    else:
        # commence with searching local minima around edges associated with AFS surface
        for edge in AFS['FaceOuterBoundEdgeLoopList']:
            if edge['superTypeName'] != 'SEAM_CURVE': # exclude SEAM_CURVE data
                edgeMinima = {}
                edgeMinimaCentroidDisp = {}

                # in the instance that a vertex is the edge minima, compare minima associated with adjoining edges & surfaces
                # adjEdgeSet, adjSurfaceSet contain the indices of adjacent surface & edges

                if edge['vertex1extremaMin']:
                    adjSurfaceSet = edge['vertex1SurfacesIndex'].values()
                    adjEdgeSet = edge['vertex1edgesIndex'].values()
                elif edge['vertex2extremaMin']:
                    adjSurfaceSet = edge['vertex2SurfacesIndex'].values()
                    adjEdgeSet = edge['vertex2edgesIndex'].values()
                else:  # minPoint somewhere between vertices on edge todo midpoint=> surfaces adjoining endpoints
                    adjSurfaceSet = edge['edgeAdjSurfacesIndex'].values()
                    # adjEdgeSet is union of set of edges at each vertex & remove identity edge
                    adjEdgeSet = ([edge['vertex1edgesIndex'][i] for i in edge['vertex1edgesIndex'].keys() if i is not edge['edgeRef']] +
                                  [edge['vertex2edgesIndex'][i] for i in edge['vertex2edgesIndex'].keys() if i is not edge['edgeRef']])

                # vertexX, adjacent surfaces minima
                for adjSurfaceIndex in adjSurfaceSet:
                    adjSurface = AdvancedFaceSurfaces[adjSurfaceIndex]
                    if adjSurface['ParsedSurface'].get('minPoint') is not None:
                        edgeMinima[adjSurface['SurfaceRef']] = adjSurface['ParsedSurface']['minPoint'] # multiple minima problem ?-----------------
                        edgeMinimaCentroidDisp[adjSurface['SurfaceRef']] = adjSurface['ParsedSurface']['minPointCentroidDisp']
                    # 'rotSymCentre' is centre of rotationally symmetrical surface feature (e.g. cylinder minimum)
                    if adjSurface['ParsedSurface'].get('rotSymCentre') is not None:
                        if adjSurface['ParsedSurface'].get('rotSymMin') is not None:
                            if adjSurface['ParsedSurface'].get('minPointCentroidDisp') is not None:
                                edgeMinima[adjSurface['SurfaceRef']] = adjSurface['ParsedSurface']['rotSymCentre']
                                edgeMinimaCentroidDisp[adjSurface['SurfaceRef']] = adjSurface['ParsedSurface']['minPointCentroidDisp']
                            else:
                                print("unexpected!")

                # vertexX, adjacent edges minima, get all relevant minima adjacent to origin vertex/edge/surface
                for adjSurfaceEdgeIndex in adjEdgeSet:
                    adjEdge = AdvancedFaceSurfaces[adjSurfaceEdgeIndex[0]]['FaceOuterBoundEdgeLoopList'][adjSurfaceEdgeIndex[1]]
                    if adjEdge.get('minPoint') is not None:
                        edgeMinima[adjEdge['edgeRef']] = adjEdge['minPoint']
                        edgeMinimaCentroidDisp[adjEdge['edgeRef']] = adjEdge['minPointCentroidDisp']
                    if adjEdge.get('rotSymMin') is not None:
                        edgeMinima[adjEdge['edgeRef']] = adjEdge['rotSymCentre']
                        edgeMinimaCentroidDisp[adjEdge['edgeRef']] = adjEdge['rotSymRadius']

                # test if local minimum relative to surrounding minima
                if len(edgeMinimaCentroidDisp) < 1:
                    deltaMin = 0.0
                else:
                    deltaMin = edge['minPointCentroidDisp'] - min(edgeMinimaCentroidDisp.values())
                if (np.abs(deltaMin) * 2 < eps_STEP_AP21) or (deltaMin <= 0):
                    if edge.get('rotSymCentre') is not None: # rotSym feature
                        if not array3x1Match(edge['rotSymCentre'], surfaceMinRotSym):  # check already exists
                            surfaceMinRotSym.append(edge['rotSymCentre'])
                            surfaceMinFeatureCentroidDisp .append(edge['minPointCentroidDisp'])
                    if edge['minPoint'] is not None:
                        if not array3x1Match(edge['minPoint'], surfaceMinPoint):
                            surfaceMinPoint.append(edge['minPoint'])
                            surfaceMinFeatureCentroidDisp .append(edge['minPointCentroidDisp'])
                        edge['localMinima'] = True
                    else:
                        edge['localMinima'] = False

        # test minima of surface patch (AFS) against those of immediate surface edges
        # note that multiple local surface minima can have an overall minima value
        # todo this approach breaks where there are multiple equal minima within a single surface: prepare for list representation
        edgeMinima = {}
        edgeMinimaCentroidDisp = {}
        for edge in AFS['FaceOuterBoundEdgeLoopList']: # todo 'faceInnerBoundEdgeLoopList'
            if edge['minPoint'] is not None:
                edgeMinima[edge['edgeRef']] = edge['minPoint']
                edgeMinimaCentroidDisp[edge['edgeRef']] = edge['minPointCentroidDisp']

        # test minima of surface patch against those of surface edges
        # groove representation here? - 'minPointCentroidDisp' will work, 'minArcPoint' represents arc centre? minFeatureCentroidDisp
        if (AFS['ParsedSurface'].get('minPointCentroidDisp') is not None):
            if len(edgeMinimaCentroidDisp) < 1:
                deltaMin = 0.0
            else:
                deltaMin = AFS['ParsedSurface']['minPointCentroidDisp'] - min(edgeMinimaCentroidDisp.values())

            if (np.abs(deltaMin) < eps_STEP_AP21) or (deltaMin < 0):
                if AFS['ParsedSurface']['minPoint'] is None: # rotSym feature
                    if len(surfaceMinPoint) > 0:
                        if len(surfaceMinRotSym) > 0:
                            if not array3x1Match(AFS['ParsedSurface']['rotSymCentre'], surfaceMinRotSym):  # check already exists
                                surfaceMinRotSym.append(AFS['ParsedSurface']['rotSymCentre'])
                                surfaceMinFeatureCentroidDisp.append(AFS['ParsedSurface']['minPointCentroidDisp'])
                        else: # first value
                            surfaceMinRotSym.append(AFS['ParsedSurface']['rotSymCentre'])
                            surfaceMinFeatureCentroidDisp.append(AFS['ParsedSurface']['minPointCentroidDisp'])
                else:
                    if surfaceMinPoint is not None:
                        if not array3x1Match(AFS['ParsedSurface']['minPoint'], surfaceMinPoint): # check already exists
                            surfaceMinPoint.append(AFS['ParsedSurface']['minPoint'])
                            surfaceMinFeatureCentroidDisp.append(AFS['ParsedSurface']['minPointCentroidDisp'])
                    else:
                        surfaceMinPoint.append(AFS['ParsedSurface']['minPoint'])
                        surfaceMinFeatureCentroidDisp.append(AFS['ParsedSurface']['minPointCentroidDisp'])
                    AFS['ParsedSurface']['localMinima'] = True
                    for edge in AFS['FaceOuterBoundEdgeLoopList']:
                        edge['localMinima'] = False
            # else:
            #     AFS['ParsedSurface']['localMinima'] = False

# print(surfaceMinPoint)
# print(surfaceMinFeatureCentroidDisp )

for AFSindex, AFS in enumerate(AdvancedFaceSurfaces):
    if AFS["SurfaceTypeName"] == "SPHERICAL_SURFACE":
        _1=1
        print("write code for SPHERICAL_SURFACE")
        # edge and surface point only equidistant from origin if centre coincides with origin
    else:
        # commence with searching local minima around edges associated with AFS surface
        for edge in AFS['FaceOuterBoundEdgeLoopList']:
            if edge['superTypeName'] != 'SEAM_CURVE':  # exclude SEAM_CURVE data
                edgeMaxima = {}
                edgeMaximaCentroidDisp = {}


            # # in the instance that a vertex is the edge maxima, compare maxima associated with adjoining edges & surfaces
            # # adjEdgeSet, adjSurfaceSet contain the indices of adjacent surface & edges
            # determine whether maxPoint is vertex or edge; multiple equal maxima, see hill-climbing categorisation
            # multiple local maxima suggest saddle points and/or surrounding local minima

                if edge['vertex1extremaMax']:
                    adjSurfaceSet = edge['vertex1SurfacesIndex'].values()
                    adjEdgeSet = edge['vertex1edgesIndex'].values()
                elif edge['vertex2extremaMax']:
                    adjSurfaceSet = edge['vertex2SurfacesIndex'].values()
                    adjEdgeSet = edge['vertex2edgesIndex'].values()
                else:  # maxPoint somewhere between vertices on edge
                    adjSurfaceSet = edge['edgeAdjSurfacesIndex'].values()
                    # adjEdgeSet is union of set of edges at each vertex
                    adjEdgeSet = ([edge['vertex1edgesIndex'][i] for i in edge['vertex1edgesIndex'].keys() if i is not edge['edgeRef']] +
                                  [edge['vertex2edgesIndex'][i] for i in edge['vertex2edgesIndex'].keys() if i is not edge['edgeRef']])

                # vertexX, adjacent surfaces maxima
                for adjSurfaceIndex in adjSurfaceSet:
                    adjSurface = AdvancedFaceSurfaces[adjSurfaceIndex]
                    if adjSurface['ParsedSurface'].get('maxPoint') is not None:
                        edgeMaxima[adjSurface['SurfaceRef']] = adjSurface['ParsedSurface']['maxPoint']
                        edgeMaximaCentroidDisp[adjSurface['SurfaceRef']] = adjSurface['ParsedSurface']['maxPointCentroidDisp'] # todo rename maxFeatureCentroidDisp?
                    # 'rotSymCentre' is centre to rotationally-symmetric feature, meaning that it's a placeholder for no-point-minima condition
                    if adjSurface['ParsedSurface'].get('rotSymCentre') is not None:
                        if adjSurface['ParsedSurface'].get('rotSymMax') is not None:
                            if adjSurface['ParsedSurface'].get('maxPointCentroidDisp') is not None:
                                edgeMinimaCentroidDisp[adjSurface['SurfaceRef']] = adjSurface['ParsedSurface']['maxPointCentroidDisp']
                            else:
                                print("unexpected! - frustum instance")

                # vertexX, adjacent edges maxima
                for adjSurfaceEdgeIndex in adjEdgeSet:
                    adjEdge = AdvancedFaceSurfaces[adjSurfaceEdgeIndex[0]]['FaceOuterBoundEdgeLoopList'][
                        adjSurfaceEdgeIndex[1]]
                    if adjEdge.get('maxPoint') is not None:
                        edgeMaxima[adjEdge['edgeRef']] = adjEdge['maxPoint']
                        edgeMaximaCentroidDisp[adjEdge['edgeRef']] = adjEdge['maxPointCentroidDisp']
                    if adjEdge.get('rotSymMax') is not None:
                        edgeMinima[adjEdge['edgeRef']] = adjEdge['rotSymCentre']
                        edgeMinimaCentroidDisp[adjEdge['edgeRef']] = adjEdge['rotSymRadius']


                # test if local maximum relative to surrounding maxima
                if len(edgeMaximaCentroidDisp) < 1:
                    deltaMax = 0.0
                else:
                    deltaMax = edge['maxPointCentroidDisp'] - max(edgeMaximaCentroidDisp.values())
                if (np.abs(deltaMax) < eps_STEP_AP21) or (deltaMax > 0): # EPS greater than or equal to
                    if edge.get('rotSymCentre') is not None: # rotSym feature
                        if not array3x1Match(edge['rotSymCentre'], surfaceMaxRotSym):  # check already exists
                            surfaceMaxRotSym.append(edge['rotSymCentre'])
                            surfaceMaxFeatureCentroidDisp .append(edge['maxPointCentroidDisp'])
                    if edge['maxPoint'] is not None:
                        if not array3x1Match(edge['maxPoint'], surfaceMaxPoint):
                            surfaceMaxPoint.append(edge['maxPoint'])
                            surfaceMaxFeatureCentroidDisp.append(edge['maxPointCentroidDisp'])
                        edge['localMaxima'] = True
                    else:
                        edge['localMaxima'] = False

        # test minima of surface patch (AFS) against those of immediate surface edges
        edgeMaxima = {}
        edgeMaximaCentroidDisp = {}
        for edge in AFS['FaceOuterBoundEdgeLoopList']:
            if edge['maxPoint'] is not None:
                edgeMaxima[edge['edgeRef']] = edge['maxPoint']
                edgeMaximaCentroidDisp[edge['edgeRef']] = edge['maxPointCentroidDisp']

        # test maxima of surface patch against those of surface edges
        #if (AFS['ParsedSurface'].get('maxPointCentroidDisp') is not None) and (len(edgeMaximaCentroidDisp) > 0):
        if (AFS['ParsedSurface'].get('maxPointCentroidDisp') is not None):
            if len(edgeMaximaCentroidDisp) < 1:
                deltaMax = 0.0
            else:
                deltaMax = AFS['ParsedSurface']['maxPointCentroidDisp'] - min(edgeMaximaCentroidDisp.values())

            if (np.abs(deltaMax) < eps_STEP_AP21) or (deltaMax > 0):
                if AFS['ParsedSurface']['maxPoint'] is None: # rotSym feature
                    if len(surfaceMaxPoint) > 0:
                        if len(surfaceMaxRotSym) > 0:
                            if not array3x1Match(AFS['ParsedSurface']['rotSymCentre'], surfaceMaxRotSym):  # check already exists
                                surfaceMaxRotSym.append(AFS['ParsedSurface']['rotSymCentre'])
                                surfaceMaxFeatureCentroidDisp.append(AFS['ParsedSurface']['maxPointCentroidDisp'])
                        else: # first value
                            surfaceMaxRotSym.append(AFS['ParsedSurface']['rotSymCentre'])
                            surfaceMaxFeatureCentroidDisp.append(AFS['ParsedSurface']['maxPointCentroidDisp'])
                    else:
                        if not array3x1Match(AFS['ParsedSurface']['maxPoint'], surfaceMaxPoint):
                            surfaceMaxPoint.append(AFS['ParsedSurface']['maxPoint'])
                            surfaceMaxFeatureCentroidDisp.append(AFS['ParsedSurface']['maxPointCentroidDisp'])
                        AFS['ParsedSurface']['localMaxima'] = True
                        for edge in AFS['FaceOuterBoundEdgeLoopList']:
                            edge['localMaxima'] = False


# maxima/minima on edges, but not at vertices => traverse adjoining surfaces and edges until first min/max
# maxima/minima at vertices => traverse adjoining surfaces and edges until first min/max


# if maxPoint/minPoint at an edge are at an edge vertex, this requires that all surfaces adjacent to this vertex are identified & tested

# # determine local maxima/minima for each surface
# for AFS in AdvancedFaceSurfaces:
#     mostMinPoint = AFS['FaceOuterBoundEdgeLoopList'][0]['minPoint']
#     mostMinPointCentroidDisp = AFS['FaceOuterBoundEdgeLoopList'][0]['minPointCentroidDisp']
#     mostMinPointEdgeIndex = 0
#     for ind, fobell in enumerate(AFS['FaceOuterBoundEdgeLoopList']):
#         if 'minPointCentroidDisp' in fobell:
#             fobell['localMinima'] = False
#             if fobell['minPointCentroidDisp'] < mostMinPointCentroidDisp:
#                 mostMinPointCentroidDisp = fobell['minPointCentroidDisp']
#                 mostMinPoint = fobell['minPoint']
#                 mostMinPointEdgeIndex = ind
#
#     if AFS['ParsedSurface']['minPointCentroidDisp'] < mostMinPointCentroidDisp:  # halt at local minima
#         AFS['ParsedSurface']['localMinima'] = True
#     else:
#         AFS['ParsedSurface']['localMinima'] = False
#         AFS['FaceOuterBoundEdgeLoopList']['mostMinPointEdgeIndex']['localMinima'] = True
#
#     mostMaxPoint = AFS['FaceOuterBoundEdgeLoopList'][0]['maxPoint']
#     mostMaxPointCentroidDisp = AFS['FaceOuterBoundEdgeLoopList'][0]['maxPointCentroidDisp']
#     mostMaxPointEdgeIndex = 0
#     for ind, fobell in enumerate(AFS['FaceOuterBoundEdgeLoopList']):
#         if 'maxPointCentroidDisp' in fobell:
#             fobell['localMaxima'] = False
#             if fobell['maxPointCentroidDisp'] > mostMaxPointCentroidDisp:
#                 mostMaxPointCentroidDisp = fobell['maxPointCentroidDisp']
#                 mostMaxPoint = fobell['maxPoint']
#                 mostMaxPointEdgeIndex = ind
#
#     if AFS['ParsedSurface']['maxPointCentroidDisp'] > mostMaxPointCentroidDisp:
#         AFS['ParsedSurface']['localMaxima'] = True
#     else:
#         AFS['ParsedSurface']['localMaxima'] = False
#         AFS['FaceOuterBoundEdgeLoopList'][mostMinPointEdgeIndex]['localMaxima'] = True

# traverse data and construct a set of minima points, including their nearest minima points
# repeat for maxima points
# philosophy is to find local minima within surrounding ring of nearest



# def searchMinima(surfaceInstance, surfaceSet):
#     if surfaceInstance['ParsedSurface']['localMinima']: # halt search
#         return surfaceInstance['ParsedSurface']['minPoint'], surfaceInstance['ParsedSurface']['minPointCentroidDisp']
#     else:
#         for edge in surfaceInstance['FaceOuterBoundEdgeLoopList']:
#             if edge['vertexExtrema']:
#                 if edge['localMinima']: # a vertex may also be a boundary local minima, but requires checking adjacent surfaces
#                     if np.allclose(edge['minPoint'], edge['vertex1'], eps_STEP_AP21):
#                         for adjSurfaceInd in edge['vertex1Surfaces']:
#                             adjSurface = surfaceSet[adjSurfaceInd]
#                             if not adjSurface['ParsedSurface']['localMinima']:
#                                 searchMinima(surfaceInstance, surfaceSet)
#
#                     elif np.allclose(edge['minPoint'], edge['vertex2'], eps_STEP_AP21):
#                         for adjSurfaceInd in edge['vertex2Surfaces']:
#                             adjSurface = surfaceSet[adjSurfaceInd]
#                             searchMinima(adjSurface, surfaceSet)
#                             # if not adjSurface['ParsedSurface']['localMinima']:
#                             #     searchMinima(adjSurface, surfaceSet)
#
#                     else: # minPont/maxPoint is not at vertices
#                         adjSurfaceInd = edge['AdjoiningEdgeSurface']
#                         adjSurface = surfaceSet[adjSurfaceInd]
#                         searchMinima(adjSurface, surfaceSet)
#
# minimaFeatures = []
# for AFS in AdvancedFaceSurfaces:
#     minimaFeatures.append(searchMinima(AFS, AdvancedFaceSurfaces))
#
# def searchMaxima(surfaceInstance, surfaceSet):
#     if surfaceInstance['ParsedSurface']['localMaxima']: # halt search, usually spherical, toroidal, NURBS
#         return surfaceInstance['ParsedSurface']['maxPoint'], surfaceInstance['ParsedSurface']['maxPointCentroidDisp']
#     else:
#         for edge in surfaceInstance['FaceOuterBoundEdgeLoopList']:
#             if edge['vertexExtrema']:
#                 if edge['localMaxima']: # a vertex may also be a boundary local minima, but requires checking adjacent surfaces
#                     if np.allclose(edge['maxPoint'], edge['vertex1'], eps_STEP_AP21):
#                         for adjSurfaceInd in edge['vertex1Surfaces']:
#                             adjSurface = surfaceSet[adjSurfaceInd]
#                             if not adjSurface['ParsedSurface']['localMaxima']:
#                                 searchMaxima(surfaceInstance, surfaceSet)
#
#                     elif np.allclose(edge['maxPoint'], edge['vertex2'], eps_STEP_AP21):
#                         for adjSurfaceInd in edge['vertex1Surfaces']:
#                             adjSurface = surfaceSet[adjSurfaceInd]
#                             searchMaxima(adjSurface, surfaceSet)
#                             # if not adjSurface['ParsedSurface']['localMaxima']:
#                             #     searchMaxima(adjSurface, surfaceSet)
#
#                     else: # minPont/maxPoint is not at vertices
#                         adjSurfaceInd = edge['AdjoiningEdgeSurface']
#                         adjSurface = surfaceSet[adjSurfaceInd]
#                         searchMaxima(adjSurface, surfaceSet)
#
# maximaFeatures = []
# for AFS in AdvancedFaceSurfaces:
#     maximaFeatures.append(searchMaxima(AFS, AdvancedFaceSurfaces))
#
#     # for maxima/minima of each surface, if minima/maxima is on edge,
#     # then query adjoining surface(s) for greater/lesser maxima/minima
#
#     # from AdvancedFaceSurfaces[], recalculate centroid from discovered vertices?

# identify minima/maxima points in hierarchy from local minima maxima points.
# for every point in all points, find nearest neighbours and determine maximal disp from centroid
# calculate rotSym minima maxima from surrounding vertices

# for finding nearest vertices, can use KD-tree, e.g.
# from scipy import spatial # numpy < 1.24
# A = np.random.random((10,3))*100
# pt = [6, 30]  # <-- the point to find
# A[spatial.KDTree(A).query(pt)[1]] # <-- the nearest point
# distance,index = spatial.KDTree(A).query(pt)

# # test ensuring all children of 'ORIENTED_EDGE' is 'EDGE_CURVE
# OrientedEdgeChildren = []
# for s in STEP_entities:
#     if hasattr(s, 'type_name'):  # stepcode.Part21.ComplexEntity type missing attribute
#         if s.type_name == 'ORIENTED_EDGE':
#             sp = cleanSubRefs(s.params)
#             for ec in sp:
#                 if hasattr(STEP_entities[ref2index(ec)], 'type_name'):
#                     #if STEP_entities[ref2index(ec)].type_name != 'VERTEX_POINT':
#                     OrientedEdgeChildren.append(STEP_entities[ref2index(ec)].type_name)

# # find all children of 'EDGE_CURVE'
# EdgeCurveChildren = []
# for s in STEP_entities:
#     if hasattr(s, 'type_name'):  # stepcode.Part21.ComplexEntity type missing attribute
#         if s.type_name == "EDGE_CURVE":
#             sp = cleanSubRefs(s.params)
#             for ec in sp:
#                 if hasattr(STEP_entities[ref2index(ec)], 'type_name'):
#                     if STEP_entities[ref2index(ec)].type_name != 'VERTEX_POINT':
#                         EdgeCurveChildren.append(STEP_entities[ref2index(ec)].type_name)

# extract all rotationally symmetric features,
# test against adjoining surface for groove or ridge distinction

#STOP HERE

# replacement of verb NURB functions with geomdl or SINTEF/Splipy equivalent
# verb.verb_eval_Tess.rationalCurveRegularSample()

# verb.verb_eval_Eval.rationalCurveDerivatives()

# verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights()
# verb.verb_eval_Eval.dehomogenize()
# verb.verb_eval_Eval.curvePoint()

# numpy analogues exist
# verb.verb_core_Vec.norm()
# verb.verb_core_Vec.sub()
# verb.verb_core_Vec.dot()
# verb.verb_core_Vec.mul()
# verb.verb_core_ArrayExtensions.last()
# verb.verb_core_Constants.EPSILON
# verb.verb_core_Vec.normSquared()




def STEPgetChildSet(STEPobject, parentName, childName):
    """
    given a parent name string, return all sets of child sub-ref numbers that match child name string
    represent as DFS traversal of N-ary tree?
    """

    def cleanSubRefs(refStr):
        """
        flatten parameter sets and extract #d string references
        """
        flattenList = (
            lambda irregularL: [e for i in irregularL for e in flattenList(i)]
            if type(irregularL) is list
            else [irregularL]
        )
        sp = flattenList(refStr)  # flatten sublists
        sp = [i for i in sp if isinstance(i, str)]
        sp = [re.match("#[0-9]+", i) for i in sp]
        sp = [i.string for i in sp if i is not None]
        return sp

    # def findChild(STEPobject, root, child):
    #
    #
    # childSet = []
    # for s in STEPobject:
    #     if hasattr(s, 'type_name'):  # stepcode.Part21.ComplexEntity type missing attribute
    #         if s.type_name == parentName:
    #             childChildSet = []
    #             # extract all '#d'
    #             sp = cleanSubRefs(s.params)
    #             for ec in sp:
    #                 if hasattr(STEP_entities[ref2index(ec)], 'type_name'):
    #                     if STEP_entities[ref2index(ec)].type_name != childRef:
    #                         EdgeCurveChildren.append(STEP_entities[ref2index(ec)].type_name)


# create list of EDGE_LOOP that define individual surfaces
# child of ADVANCED_FACE & FACE_OUTER_BOUND

# ???
# P = np.array([1,0,1])
# center = np.array([1,0,0])
# radius = 1
# n2 = np.array([0,0,1])
#
# Delta = P-center
# dist = np.sqrt(np.dot(n2, Delta)**2 + (np.linalg.norm(np.cross(n2, Delta))- radius)**2)
# ???

# strategy for planar polygonal surfaces is to divide into non-overlapping triangles where corners are polygon vertices
# and test for nearest/furthest point on triangle

# import and build pymesh locally, see https://github.com/PyMesh/PyMesh/blob/main/README.md


# class STEPface(object):
#     """
#     container class for STEP AP 203 surface
#     """
#
#     def __init__(self, baseRef=""):
#         self.baseRef = baseRef
#         self.typeName = '' # 'CYLINDRICAL_SURFACE',
#
#     def addShape(
#         self,
#         shapeObject,
#         shapeName
#     ):
#         shapeViewObject = self.modelDoc.addObject("Part::FeaturePython", shapeName)


# STEP Part21 decomposition:
# https://www.mbx-if.org/documents/AP203e2_html/AP203e2.htm

# root node, CLOSED_SHELL (OPEN_SHELL), Connected_Face_Set SUPERTYPE OF (ONEOF (Closed_Shell, Open_Shell))

# rethink
# deconstruct STEP model to find all faces
# categorise circles, cylinders and other forms with axes.
# extract points from Face_outer_bound (holes in face for minima?)

# first pass, find all edge cartesian_points to generate origin based on median
# second pass, extract surface maxima from faces,

# get all VERTEX, part of EDGE_LOOP
# get median barycentre/origin

# for all EDGE_LOOP, get VERTEX, find furthest vertex
# for all surfaces bound by EDGE_LOOP, calculate if they contain a point normal to an origin-ray


# note that for planar regions, maxima/minima is at a polygon point
# unless polygon plane is perfectly normal to ray =>

# random start, find points closest to initial Desernosphere rays (rethink: more relevant to NURBS surfaces)
# get polygon point sets attached to point, find polygon point with max/min disp, repeat
# cylinder surface ?? different process, look for entity type
# efficient graph node search between interconnecting points?
# surface normal - Newell's method for surface normal, cross product == zero, origin on a line parallel to surface normal

# for any planar polygon, if 1 point is further from/closer to origin than the others, it is polygon maxima/minima
# if 2 points share distance, lacal max/min is at edge.

# for any polygon point or edge at min/max find adjacent polygon to point/edge and check for associated points/edges of lesser/greater origin displacement.
# (3 points at equivalent disp on valid for symmetric polygon centered and normal to axis thru origin)

# for any cylindrical surface, get axis position
# find shortest distance from origin cylinder axis, returned point is orthogonal
# determine if all cylinder edges are above or below this point
# if not, max disp is point + radius

# torus, find max/min point on major radius circle, repeat for a minor radius circle at max/min point

# ellipsoid, try closest point to an ellipse parallel to major axis?

# 1. get the median origin from all points
# 2. find maximum displacement of points from this point, construct ray
# 3. for every ray find closest cosine angle points
# 3a. identify closest quadrant through point signs (required to translate model to origin?)
# segmentAndPlane( p0 : Point, p1 : Point, v0 : Point, n : Point )
# 3b. distances from points to ray
# 4. identify surfaces associated with these (edge) points
# 5. for each surface, determine if ray is between bounding edge points


    #   (axis1_placement,
    #        axis2_placement_2d,
    #        axis2_placement_3d,
    #        b_spline_curve,
    #        b_spline_curve_with_knots,
    #        b_spline_surface,
    #        b_spline_surface_with_knots,
    #        bezier_curve,
    #        bezier_surface,
    #        boundary_curve,
    #        cartesian_point,
    #        cartesian_transformation_operator_3d,
    #        circle,
    #        composite_curve,
    #        composite_curve_on_surface,
    #        composite_curve_segment,
    #        conic,
    #        conical_surface,
    #        curve,
    #        curve_bounded_surface,
    #        curve_replica,
    #        cylindrical_surface,
    #        degenerate_pcurve,
    #        degenerate_toroidal_surface,
    #        dimension_count,
    #        dimension_of,
    #        direction,
    #        ellipse,
    #        evaluated_degenerate_pcurve,
    #        geometric_representation_context,
    #        geometric_representation_item,
    #        hyperbola,
    #        intersection_curve,
    #        line,
    #        offset_curve_3d,
    #        offset_surface,
    #        outer_boundary_curve,
    #        parabola,
    #        pcurve,
    #        plane,
    #        point,
    #        point_on_curve,
    #        point_on_surface,
    #        point_replica,
    #        polyline,
    #        quasi_uniform_curve,
    #        quasi_uniform_surface,
    #        rational_b_spline_curve,
    #        rational_b_spline_surface,
    #        rectangular_composite_surface,
    #        rectangular_trimmed_surface,
    #        reparametrised_composite_curve_segment,
    #        seam_curve,
    #        spherical_surface,
    #        surface,
    #        surface_curve,
    #        surface_of_linear_extrusion,
    #        surface_of_revolution,
    #        surface_replica,
    #        swept_surface,
    #        toroidal_surface,
    #        trimmed_curve,
    #        uniform_curve,
    #        uniform_surface,


    # # test each surface defined in AdvancedFaceSurfaces for maxima/minima from centroid
    # if AFS["SurfaceTypeName"] == "SPHERICAL_SURFACE":
    #     # # ENTITY Spherical_surface
    #     # #   SUBTYPE OF (Surface);
    #     # #   position : Axis_placement_3d;
    #     # #   radius : positive_length_measure;
    #     #
    #     # # Attribute definitions:
    #     # # position: an Axis_placement_3d that defines the location and orientation of the surface for the purposes of parameterization.
    #     # # The centre of the Spherical_surface is at the location
    #     # # radius: the radius of the Spherical_surface.
    #     #
    #     # #  The axis is the placement Z axis direction and the ref_direction is an approximation to
    #     # # the placement X axis direction.= normDir
    #     #
    #     # radius = AFS["SurfaceParams"][-1]
    #     # ParsedSurface["radius"] = radius
    #     # centrePoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)
    #     # ParsedSurface["axisPoint"] = centrePoint
    #     # ParsedSurface["normDir"] = normDir # Z
    #     # ParsedSurface["refDir"] = refDir # X
    #     # auxDir = np.cross(normDir, refDir)
    #     # auxDir = auxDir / np.linalg.norm(auxDir) # Y
    #     # ParsedSurface["auxDir"] = auxDir
    #     #
    #     # # sphere surface maxima & minima exist on a vector between centroid and sphere centre point
    #     # # minima/maxima = centre point -/+ radius
    #     #
    #     # # "The local z-axis corresponds to the normal of planar and spherical surfaces and to the axis of cylindrical,
    #     # # conical and toroidal surfaces."
    #     # # Jaider Oussama et al Int. Journal of Engineering Research and Applications www.ijera.com
    #     # # ISSN : 2248-9622, Vol. 4, Issue 5( Version 6), May 2014, pp.14-25
    #     # # refDir, auxDir doesn't mean much in stock implementaions.
    #     #
    #     # sphereCentreCentroidDir = centrePoint - centroid
    #     # sphereCentreCentroidDir = sphereCentreCentroidDir / np.linalg.norm(sphereCentreCentroidDir)
    #     #
    #     # # determine if maxPoint/minPoint are inside or outside an edge_loop defined boundary
    #     # # CCW directionality of edge vertices determines outward surface normal
    #     # # surface minPoint/maxPoint always orthonormal to centroid
    #     #
    #     # maxPoint = centrePoint + sphereCentreCentroidDir * radius
    #     # maxPointCentroidDisp = maxPoint - centroid
    #     # maxNorm = (maxPoint - centroid) / np.linalg.norm(maxPoint - centroid)
    #     #
    #     # if not insideOutsideSurfaceTest(maxNorm, maxPoint, AFS):
    #     #     # find maxima among edges
    #     #     maxPointCentroidDisp = [edge['maxPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
    #     #     maxPoint = AFS["FaceOuterBoundEdgeLoopList"][maxPointCentroidDisp.index(max(maxPointCentroidDisp))]
    #     #     maxPointCentroidDisp = max(maxPointCentroidDisp)
    #     # else:
    #     #     ParsedSurface["maxPoint"] = maxPoint
    #     #     ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp
    #     #
    #     # minPoint = centrePoint - sphereCentreCentroidDir * radius
    #     # minPointCentroidDisp = minPoint - centroid
    #     # minNorm = (minPoint - centroid) / np.linalg.norm(minPoint - centroid)
    #     #
    #     # if not insideOutsideSurfaceTest(minNorm, minPoint, AFS):
    #     #     # find maxima among edges
    #     #     minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
    #     #     minPoint = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]
    #     #     minPointCentroidDisp = min(minPointCentroidDisp)
    #     # else:
    #     #     ParsedSurface["minPoint"] = minPoint
    #     #     ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp
    #     #
    #     #AFS['ParsedSurface'] = sphereSurfaceParse(AFS)
    #     sphereSurfaceParse(AFS) #ParsedSurface
    #
    # if AFS["SurfaceTypeName"] == "TOROIDAL_SURFACE":
    #     # # ENTITY Toroidal_surface
    #     # #   SUBTYPE OF (Surface);
    #     # #   position : Axis_placement_3d;
    #     # #   radius : positive_length_measure;
    #     # #   minor_radius : positive_length_measure;
    #     #
    #     # # Attribute definitions:
    #     # # position: an Axis_placement_3d that defines the location and orientation of the surface.
    #     # # The centre of the Toroidal_surface is at the location. The plane of the ref_direction and the axis is a
    #     # # central plane of symmetry of the Toroidal_surface. The normal to this plane is the central axis of the surface.
    #     # # radius: the distance from the central axis of the Toroidal_surface to the centre of one of the circles
    #     # # produced by intersecting the Toroidal_surface with a plane containing the central axis.
    #     # # The surface can be produced by sweeping a circle of radius minor_radius along a circle of radius radius.
    #     # # minor_radius: the radius of one of the circles produced by intersecting the Toroidal_surface with a plane
    #     # # containing the central axis.
    #     #
    #     # minorRadius = AFS["SurfaceParams"][-1]
    #     # ParsedSurface["minorRadius"] = minorRadius
    #     # majorRadius = AFS["SurfaceParams"][-2]
    #     # ParsedSurface["majorRadius"] = majorRadius
    #     #
    #     # # axisPoint, auxDir, refDir = axis2Placement3D(AFS["SurfaceParams"][1], STEP_entities)
    #     # # ParsedSurface["axisPoint"] = axisPoint
    #     # # ParsedSurface["auxDir"] = auxDir
    #     # # ParsedSurface["refDir"] = refDir
    #     # # normDir = np.cross(refDir, auxDir)
    #     # # normDir = normDir / np.linalg.norm(normDir)
    #     # # ParsedSurface["normDir"] = normDir
    #     #
    #     # # centroid = np.array([1, 1, 0])
    #     # # axisPoint = np.array([0, 0, 10])
    #     # # axisDir = np.array([0, 0, 1])
    #     # # refDir = np.array([1, 0, 0])
    #     # # tMajorRadius = 10
    #     # # tMinorRadius = 1
    #     #
    #     # axisPoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)
    #     # parsedEdge["axisPoint"] = axisPoint
    #     # parsedEdge["normDir"] = normDir
    #     # parsedEdge["refDir"] = refDir # axisDir, refDir guaranteed unit normals?
    #     # auxDir = np.cross(normDir, refDir)
    #     # auxDir = auxDir / np.linalg.norm(auxDir)
    #     # parsedEdge["auxDir"] = auxDir
    #     #
    #     # # test if toroid is rotationally symmetrical, if centroid is close to normDir through axisPoint
    #     #
    #     # rotSymDisp = np.linalg.norm(
    #     #     (centroid - axisPoint) - np.dot((centroid - axisPoint), normDir) * normDir
    #     # )
    #     #
    #     # # need a precision factor that reflects increasing accuracy arising from iterative calculations of centroid --------------------------------------<<<<<<<<<<<<<<<<<
    #     # # does one find a centroid through axis-intersection?
    #     # # also scale factor affecting numerical precision
    #     # # also STEP seems to be E-6 decimal places at best
    #     # # iterativeScaledPrecision = 1e-4  # for testing
    #     #
    #     # if rotSymDisp < eps_STEP_AP21:  # iterativeScaledPrecision:
    #     #     ParsedSurface["rotSymMax"] = axisPoint
    #     #     ParsedSurface["rotSymMin"] = axisPoint
    #     # else:
    #     #     # create axes for circles based at points of majorMinPoint, majorMaxPoint
    #     #     # minima case, create new axis at majorMinPoint
    #     #     # get min, max of major radius
    #     #     majorMinPoint, majorMaxPoint = pointCircleMinMaxDisp(
    #     #         centroid,
    #     #         axisPoint,
    #     #         normDir,
    #     #         majorRadius
    #     #     )
    #     #
    #     #     # vector from toroid centre to minima point on majorRadius
    #     #     majorMinRadialDir = (majorMinPoint - axisPoint) / np.linalg.norm(
    #     #         majorMinPoint - axisPoint
    #     #     )
    #     #     # vector tangential to major radius
    #     #     majorMinTangentDir = np.cross(majorMinRadialDir, normDir)
    #     #     majorMinTangentDir = majorMinTangentDir / np.linalg.norm(
    #     #         majorMinTangentDir
    #     #     )
    #     #
    #     #     # point on minor circle cross-section at minimum disp from centroid
    #     #     minorMinPoint, _ = pointCircleMinMaxDisp(
    #     #         centroid,
    #     #         majorMinPoint,
    #     #         majorMinTangentDir,
    #     #         minorRadius,
    #     #     )
    #     #
    #     #     ParsedSurface["minPoint"] = minorMinPoint
    #     #
    #     #     # same for maxima
    #     #     majorMaxRadialDir = (majorMaxPoint - axisPoint) / np.linalg.norm(
    #     #         majorMaxPoint - axisPoint
    #     #     )
    #     #     majorMaxTangentDir = np.cross(majorMaxRadialDir, auxDir)
    #     #     majorMaxTangentDir = majorMaxTangentDir / np.linalg.norm(
    #     #         majorMaxTangentDir
    #     #     )
    #     #
    #     #     minorMaxPoint, _ = pointCircleMinMaxDisp(
    #     #         centroid,
    #     #         majorMaxPoint,
    #     #         majorMaxTangentDir,
    #     #         minorRadius,
    #     #     )
    #     #
    #     #     # np.linalg.norm(centroid - tMinorMinPoint) < np.linalg.norm(centroid - tMinorMaxPoint)
    #     #     ParsedSurface["maxPoint"] = minorMaxPoint
    #     #
    #     #     # need to test whether maxima/minima are within defined edges
    #     #     # call edge loop tests after surface max/min points extraction, line 1384
    #     #
    #     #     # determine if maxPoint/minPoint are inside or outside an edge_loop defined boundary
    #     #     # CCW directionality of edge vertices determines outward surface normal
    #     #     # surface minPoint/maxPoint always orthonormal to centroid
    #     #
    #     #     maxNorm = (maxPoint - centroid) / np.linalg.norm(maxPoint - centroid)
    #     #     if not insideOutsideSurfaceTest(maxNorm, maxPoint, AFS):
    #     #         # find maxima among edges
    #     #         maxPointCentroidDisp = [edge['maxPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
    #     #         maxPoint = AFS["FaceOuterBoundEdgeLoopList"][maxPointCentroidDisp.index(max(maxPointCentroidDisp))]
    #     #         maxPointCentroidDisp = max(maxPointCentroidDisp)
    #     #     else:
    #     #         ParsedSurface["maxPoint"] = maxPoint
    #     #         ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp
    #     #
    #     #     minNorm = (minPoint - centroid) / np.linalg.norm(minPoint - centroid)
    #     #     if not insideOutsideSurfaceTest(minNorm, minPoint, AFS):
    #     #         # find maxima among edges
    #     #         minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
    #     #         minPoint = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]
    #     #         minPointCentroidDisp = min(minPointCentroidDisp)
    #     #     else:
    #     #         ParsedSurface["minPoint"] = minPoint
    #     #         ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp
    #     #
    #     #     AFS['ParsedSurface'] = ParsedSurface
    #     toroidSurfaceParse(AFS)
    #
    # if (AFS["SurfaceTypeName"] == "CYLINDRICAL_SURFACE"):
    #     # # ENTITY Cylindrical_surface
    #     # #   SUBTYPE OF (Surface);
    #     # #   position : Axis_placement_3d;
    #     # #   radius : positive_length_measure;
    #     #
    #     # # position: an Axis_placement_3d that defines the location and orientation of the surface.
    #     # # The axis of the Cylindrical_surface passes through the location and is normal to the plane of ref_direction and axis.
    #     # # radius: the radius of the circular curve of intersection between the surface of the cylinder and a plane
    #     # # perpendicular to the axis of the Cylindrical_surface.
    #     #
    #     # radius = AFS["SurfaceParams"][-1]
    #     # ParsedSurface["radius"] = radius
    #     # axisPoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)
    #     # ParsedSurface["axisPoint"] = axisPoint
    #     # ParsedSurface["normDir"] = normDir
    #     # ParsedSurface["refDir"] = refDir # axisDir, refDir guaranteed unit normals?
    #     # auxDir = np.cross(normDir, refDir)
    #     # auxDir = auxDir / np.linalg.norm(auxDir)
    #     # ParsedSurface["auxDir"] = auxDir
    #     #
    #     # # ===================================================================================
    #     # # TESTING MAX/MIN POINTS ON ROTATIONALLY SYMMETRIC SURFACES
    #     # # CHECKING FOR INTERSECTION ON EDGES
    #     # # ===================================================================================
    #     #
    #     # # test if cylinder is rotationally symmetrical, if centroid is close to axisDir through axisPoint
    #     # rotSymDisp = np.linalg.norm(
    #     #     (centroid - axisPoint) - np.dot((centroid - axisPoint), normDir) * normDir
    #     # )
    #     #
    #     # if rotSymDisp < eps_STEP_AP21:
    #     #     #ParsedSurface["rotSymMax"] = axisPoint #-------------------at circle
    #     #     minSurfaceDispPoint = centroid # for creating a minima plane intersection
    #     #     centroidProjAxisPoint = centroid
    #     #     # maxima is circular edge centre, determined from vertices
    #     #     # minima is orthogonal to centroid
    #     #     # shouldn't there also be a radius measure?? check-nope
    #     # else:
    #     #     # cylinder/conic maxima is just a search of greatest disp of vertices from centroid
    #     #     # if centroid is outside full cylindrical surface,
    #     #     #   cylinder minima -> test for intersection with surface on orthogonal projection from axis to centroid
    #     #     # if centroid within surface, |centroid - projAxisDir| < cRadius
    #     #
    #     #     # test:
    #     #     # axisPoint = np.array([1, 0, 0])
    #     #     # centroid = np.array([2, 2, 0])
    #     #     # axisDir = np.array([1.5, 0, 0])
    #     #     # cRadius = 1
    #     #
    #     #     centroidProjAxisPoint = axisPoint + (np.dot((centroid - axisPoint), normDir)/np.dot(normDir, normDir)) * normDir
    #     #     centroidAxisDir = centroidProjAxisPoint - centroid
    #     #     centroidAxisDisp = np.linalg.norm(centroidAxisDir)
    #     #     centroidAxisDir = centroidAxisDir/centroidAxisDisp
    #     #
    #     #     if centroidAxisDisp > radius:
    #     #         minSurfaceDispPoint = centroid + (np.linalg.norm(centroidAxisDir) - radius) * centroidAxisDir
    #     #     else:
    #     #         minSurfaceDispPoint = centroid - (radius - np.linalg.norm(centroidAxisDir)) * centroidAxisDir
    #     #
    #     # # shift below to function?=========================================================???????????????
    #     #
    #     # intersectedEdgeSet = []
    #     # intersectPoint = []
    #     # rotSymFeature = []
    #     # rotSymFeatureEdgeSet = []
    #     # ParsedSurface["rotSymEdge"] = []
    #     # # determine whether any edge in edge loop contains an intersection with a radius circle or plane containing this point
    #     # for pel in parsedEdgeLoop:
    #     #
    #     #     v1 = pel["vertex1"]
    #     #     v2 = pel["vertex2"]
    #     #
    #     #     if pel["typeName"] == "LINE":
    #     #         # test whether a plane containing minPoint thru axisDir intersects an edge segment
    #     #         edgeIntersect = intersectSegmentPlane( v1, v2,
    #     #                                                normDir,
    #     #                                                minSurfaceDispPoint )
    #     #         if edgeIntersect is not None:
    #     #             intersectedEdgeSet.append(pel)
    #     #             intersectPoint.append(edgeIntersect)
    #     #
    #     #     if pel["typeName"] == "CIRCLE":
    #     #
    #     #         edgeArcIntersect = intersectArcPlane(normDir,
    #     #                               centroidProjAxisPoint,
    #     #                               pel["axisPoint"],
    #     #                               pel["refDir"],
    #     #                               pel["auxDir"],
    #     #                               pel["normDir"],
    #     #                               pel["radius"],
    #     #                               v1,
    #     #                               v2)
    #     #
    #     #         if len(edgeArcIntersect) > 0:
    #     #             intersectPoint = intersectPoint + edgeArcIntersect
    #     #             intersectedEdgeSet.append(pel)
    #     #
    #     #         # create a rotationally symmetric edge centre point
    #     #         # perhaps easier to determine ridges/grooves following minima maxima chase?
    #     #         # requires rot-sym surface axis to pass through centroid
    #     #
    #     #         # if rotSymDisp < eps_STEP_AP21:
    #     #         #     # this is a circle on a cylinder with axis colinear with centroid
    #     #         #     pass
    #     #
    #     #         if (np.linalg.norm(pel["normDir"] - ParsedSurface["normDir"]) < eps_STEP_AP21 or
    #     #                 np.linalg.norm(pel["normDir"] + ParsedSurface["normDir"]) < eps_STEP_AP21):
    #     #             # arc circle norm is parallel to surface cylinder norm, (default case for simple cylinder)
    #     #             # test whether centroid projects to arc centre
    #     #
    #     #             arcDisp = np.dot((centroid - pel["axisPoint"]), pel["normDir"])
    #     #
    #     #             # centroid projected to arc plane
    #     #             projCentroid = centroid - np.dot(arcDisp, pel["normDir"])
    #     #             if (np.linalg.norm(projCentroid - pel["axisPoint"]) < eps_STEP_AP21 or
    #     #                     np.linalg.norm(projCentroid + pel["axisPoint"]) < eps_STEP_AP21):
    #     #                 # arc on a cylinder cannot be determined as groove or ridge without surrounding surface data
    #     #                 ParsedSurface["rotSymEdge"].append(pel["axisPoint"])
    #     #                 rotSymFeatureEdgeSet.append(pel)
    #     #                 #rotSymFeature.append(pel["axisPoint"])
    #     #
    #     #             _1=1
    #     #
    #     #         # # circleYaxis = np.cross(pel["axisDir"], pel["refDir"])
    #     #         # # circleYaxis = circleYaxis/np.linalg.norm(circleYaxis)
    #     #         #
    #     #         # circleEdge = verb.verb_geom_Circle( pel["axisPoint"].tolist(),
    #     #         #                                     pel["refDir"].tolist(),
    #     #         #                                     pel["auxDir"].tolist(),
    #     #         #                                     pel["cRadius"] )
    #     #         #
    #     #         # circleMinSurfaceRadius = verb.verb_geom_Circle( axisPoint.tolist(),
    #     #         #                                                 refDir.tolist(),
    #     #         #                                                 auxDir.tolist(),
    #     #         #                                                 centroidAxisDisp )
    #     #         #
    #     #         # edgeArcIntersect = verb.verb_geom_Intersect.curves(circleEdge, circleMinSurfaceRadius)
    #     #         #
    #     #         # # general intersections of edge arc with plane of minSurfaceDisp point
    #     #         # # edgeArcIntersect1, edgeArcIntersect2 = intersectArcRotSymSurface( surfaceAxisDir,
    #     #         # #                                                                   surfaceAxisPoint,
    #     #         # #                                                                   minSurfaceDisp,
    #     #         # #                                                                   pel["axisDir"],
    #     #         # #                                                                   pel["axisPoint"],
    #     #         # #                                                                   pel["cRadius"] )
    #     #         #
    #     #         # # test if circle intersections lies in vertex defined arc
    #     #         # if edgeArcIntersect1 is not None:
    #     #         #     if pointInArc( edgeArcIntersect1,
    #     #         #                    v1, v2,
    #     #         #                    pel["axisDir"],
    #     #         #                    pel["refDir"],
    #     #         #                    pel["axisPoint"]):
    #     #         #         intersectedEdgeSet.append(pel)
    #     #         #         intersectPoint.append(edgeIntersect)
    #     #         #
    #     #         # if edgeArcIntersect2 is not None:
    #     #         #     if pointInArc( edgeArcIntersect2,
    #     #         #                    v1, v2,
    #     #         #                    pel["axisDir"],
    #     #         #                    pel["refDir"],
    #     #         #                    pel["axisPoint"]):
    #     #         #         if pel not in intersectedEdgeSet:
    #     #         #             intersectedEdgeSet.append(pel)
    #     #         #         intersectPoint.append(edgeIntersect)
    #     #
    #     #     if pel["typeName"] == "ELLIPSE":
    #     #         # general intersections of edge arc with plane of minSurfaceDisp point
    #     #         # edgeArcIntersect1, edgeArcIntersect2 = intersectArcRotSymSurface( surfaceAxisDir,
    #     #         #                                                                   surfaceAxisPoint,
    #     #         #                                                                   minSurfaceDisp,
    #     #         #                                                                   pel["axisDir"],
    #     #         #                                                                   pel["axisPoint"],
    #     #         #                                                                   pel["eMajorRadius"],
    #     #         #                                                                   pel["eMinorRadius"])
    #     #
    #     #         # class EllipseArc extends NurbsCurve {
    #     #         #     //* Length 3 array representing the center of the arc
    #     #         #     //* Length 3 array representing the xaxis
    #     #         #     //* Length 3 array representing the perpendicular yaxis
    #     #         #     //* Minimum angle of the EllipseArc
    #     #         #     //* Maximum angle of the EllipseArc
    #     #         #
    #     #         #     public function new(   center : Point,
    #     #         #                            xaxis : Vector,
    #     #         #                            yaxis : Vector,
    #     #         #                            minAngle : Float,
    #     #         #                            maxAngle : Float )
    #     #
    #     #         # class Ellipse extends EllipseArc {
    #     #         #
    #     #         #     //Create an ellipse
    #     #         #     //* Length 3 array representing the center of the circle
    #     #         #     //* Length 3 array representing the xaxis
    #     #         #     //* Length 3 array representing the perpendicular yaxis
    #     #         #
    #     #         #     public function new(   center : Point,
    #     #         #                             xaxis : Vector,
    #     #         #                             yaxis : Vector )
    #     #
    #     #         # normal direction is inferred from xaxis X yaxis rather than STEP axisDir convention
    #     #         # major/minor radius is defined as xaxis/yaxis
    #     #         # to convert from STEP, find major/minor axis vectors and scale to minorRadius majorRadius values
    #     #         # have to assume AP21 majorRadius is aligned with refDir until proven otherwise
    #     #
    #     #         ellipseYaxis = pel["auxDir"]
    #     #         ellipseYaxis = ellipseYaxis/np.linalg.norm(ellipseYaxis)
    #     #         ellipseYaxis = ellipseYaxis * pel["minorRadius"]
    #     #
    #     #         ellipseXaxis = pel["refDir"]
    #     #         ellipseXaxis = ellipseXaxis / np.linalg.norm(ellipseXaxis)
    #     #         ellipseXaxis = ellipseXaxis * pel["majorRadius"]
    #     #
    #     #         ellipseEdge = verb.verb_geom_Ellipse( pel["axisPoint"].tolist(),
    #     #                                               ellipseXaxis.tolist(),
    #     #                                               ellipseYaxis.tolist() )
    #     #
    #     #         # surfaceYaxisDir = np.cross(surfaceAxisDir, surfaceRefDir)
    #     #         # surfaceYaxisDir = surfaceYaxisDir/np.linalg.norm(surfaceYaxisDir)
    #     #
    #     #         # circleMinSurfaceRadius = verb.verb_geom_Circle( axisPoint.tolist(),
    #     #         #                                                 refDir.tolist(),
    #     #         #                                                 auxDir.tolist(),
    #     #         #                                                 centroidAxisDisp )
    #     #
    #     #         # verb core Plane
    #     #
    #     #         auxDirP = np.cross(centroidAxisDir, auxDir)  #  surface axisDir
    #     #         p1 = centroidProjAxisPoint + centroidAxisDir * centroidAxisDisp
    #     #         p3 = centroidProjAxisPoint - centroidAxisDir * centroidAxisDisp
    #     #         p2 = centroidProjAxisPoint + auxDirP * centroidAxisDisp
    #     #         p4 = centroidProjAxisPoint - auxDirP * centroidAxisDisp
    #     #
    #     #         intersectPlane = verb.verb_eval_Make.fourPointSurface(p1, p2, p3, p4)
    #     #
    #     #         # class Intersect {
    #     #         #
    #     #         #     //Determine the intersection of two curves
    #     #         #     //* ICurve object
    #     #         #     //* ICurve object
    #     #         #     //* tolerance for the intersection
    #     #         #     //
    #     #         #     //**returns**
    #     #         #     //* a possibly empty array of CurveCurveIntersection objects
    #     #         #
    #     #         #     public static function curves( first : ICurve, second : ICurve, tol : Float = 1e-3  )
    #     #
    #     #         #edgeArcIntersect = verb.verb_geom_Intersect.curveAndSurface(ellipseEdge, intersectPlane)
    #     #         #edgeArcIntersect = verb.verb_geom_Intersect.curves(ellipseEdge, circleMinSurfaceRadius)
    #     #         edgeArcIntersect = verb.verb_eval_Intersect.curveAndSurface(ellipseEdge.asNurbs(), intersectPlane, tol=1e-3)
    #     #
    #     #         if len(edgeArcIntersect) > 0:
    #     #             intersectPoint = intersectPoint + edgeArcIntersect
    #     #             intersectedEdgeSet.append(pel)
    #     #
    #     #         # test if circle intersections lies in vertex defined arc
    #     #         # if edgeArcIntersect1 is not None:
    #     #         #     if pointInArc( edgeArcIntersect1,
    #     #         #                    v1, v2,
    #     #         #                    pel["axisDir"],
    #     #         #                    pel["refDir"],
    #     #         #                    pel["axisPoint"]):
    #     #         #         intersectedEdgeSet.append(pel)
    #     #         #         intersectPoint.append(edgeIntersect)
    #     #         #
    #     #         # if edgeArcIntersect2 is not None:
    #     #         #     if pointInArc( edgeArcIntersect2,
    #     #         #                    v1, v2,
    #     #         #                    pel["axisDir"],
    #     #         #                    pel["refDir"],
    #     #         #                    pel["axisPoint"]):
    #     #         #         if pel not in intersectedEdgeSet:
    #     #         #             intersectedEdgeSet.append(pel)
    #     #         #         intersectPoint.append(edgeIntersect)
    #     #
    #     #     if pel["typeName"] == "B_SPLINE_CURVE_WITH_KNOTS":
    #     #         # try intersecting spline curve with plane containing minPoint
    #     #         # find closestPoint() < eps
    #     #
    #     #         curveDegree = pel["curveDegree"]
    #     #         controlPointsList = pel["controlPointsList"]
    #     #         knotvector = pel["knotvector"]
    #     #
    #     #         NURBS_edge = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights( curveDegree,
    #     #                                                                             knotvector,
    #     #                                                                             [cp.tolist() for cp in controlPointsList] )
    #     #
    #     #         refDirP = np.cross(centroidAxisDir, axisDir)  #  surface axisDir
    #     #         p1 = centroidProjAxisPoint + centroidAxisDir * minSurfaceDisp
    #     #         p3 = centroidProjAxisPoint - centroidAxisDir * minSurfaceDisp
    #     #         p2 = centroidProjAxisPoint + refDirP * minSurfaceDisp
    #     #         p4 = centroidProjAxisPoint - refDirP * minSurfaceDisp
    #     #
    #     #         #   Generate the control points, weights, and knots of a surface defined by 4 points
    #     #         #   first point in counter-clockwise form
    #     #         #   second point in counter-clockwise form
    #     #         #   third point in counter-clockwise form
    #     #         #   forth point in counter-clockwise form
    #     #
    #     #         #   returns: NurbsSurfaceData object
    #     #         #   fourPointSurface( p1 : Point, p2 : Point, p3 : Point, p4 : Point, degree : Int = 3 )
    #     #
    #     #         intersectPlane = verb.verb_eval_Make.fourPointSurface(p1, p2, p3, p4)
    #     #         edgeArcIntersect = verb_eval_Intersect.curveAndSurface(NURBS_edge, intersectPlane, tol=1e-3)
    #     #         #     //Get the intersection of a NURBS curve and a NURBS surface without an estimate
    #     #         #     //
    #     #         #     //**params**
    #     #         #     //
    #     #         #     //* NurbsCurveData
    #     #         #     //* NurbsSurfaceData
    #     #         #     //* tolerance for the curve intersection
    #     #         #     //
    #     #         #     //**returns**
    #     #         #     //
    #     #         #     //* array of CurveSurfaceIntersection objects
    #     #         #
    #     #         #     public static function curveAndSurface( curve : NurbsCurveData,
    #     #         #                                               surface : NurbsSurfaceData,
    #     #         #                                               tol : Float = 1e-3,
    #     #         #                                               crvBbTree : IBoundingBoxTree<NurbsCurveData> = null,
    #     #         #                                               srfBbTree : IBoundingBoxTree<NurbsSurfaceData> = null ) : Array<CurveSurfaceIntersection>  {
    #     #         #
    #     #         #         crvBbTree = crvBbTree != null ? crvBbTree : new LazyCurveBoundingBoxTree( curve );
    #     #         #         srfBbTree = srfBbTree != null ? srfBbTree : new LazySurfaceBoundingBoxTree( surface );
    #     #
    #     #
    #     #         if len(edgeArcIntersect) > 0:
    #     #             intersectPoint = intersectPoint + edgeArcIntersect
    #     #             intersectedEdgeSet.append(pel)
    #     #
    #     # # sanity check that points in intersectPoint are of greater centroid disp values than points in edges?
    #     #
    #     # # determine winding number of edges intersected for validity of maxima/minima
    #     # # count edge direction?
    #     # # with a central axis, easier to just check is point is in arc
    #     #
    #     # # if AFS['SurfaceTypeName'] == 'CYLINDRICAL_SURFACE': # debug
    #     # #     _1=1
    #     #
    #     # if len(intersectPoint) == 2 and len(intersectedEdgeSet) == 2:
    #     #     if pointInArc(centroidProjAxisPoint,
    #     #                   intersectPoint[0],
    #     #                   intersectPoint[1],
    #     #                   ParsedSurface["refDir"],
    #     #                   ParsedSurface["auxDir"],
    #     #                   ParsedSurface["normDir"],
    #     #                   ParsedSurface["axisPoint"],
    #     #                   rotSym=True):
    #     #         if rotSymDisp < eps_STEP_AP21:
    #     #             ParsedSurface["rotSymEdge"].append(centroidProjAxisPoint)  #-----------------------------------------
    #     #             # minima or maxima
    #     #         else:
    #     #             ParsedSurface["minPoint"] = minSurfaceDispPoint
    #     #         ParsedSurface["minPointCentroidDisp"] = np.linalg.norm(centroid - minSurfaceDispPoint)  #minPointCentroidDisp
    #     #
    #     # #if len(intersectPoint) > 2:
    #     #
    #     # # # can't determine in advance whether rotSymMax or rotSymMin => requires subsequent comparison to adjacent nodes
    #     # # # same as maxPoint minPoint comparison
    #     cylinderSurfaceParse(AFS)
    #
    # if (AFS["SurfaceTypeName"] == "CONICAL_SURFACE"): # still must test with simple STEP file
    #     # # ENTITY Conical_surface
    #     # #   SUBTYPE OF (Surface);
    #     # #   position : Axis_placement_3d;
    #     # #   radius : positive_length_measure;
    #     # #   semi_angle : plane_angle_measure;
    #     #
    #     # # position: an Axis_placement_3d that defines the location and orientation of the surface.
    #     # # The axis of the Conical_surface passes through the location and is normal to the plane of ref_direction and axis.
    #     # # radius: the radius of the circular curve of intersection between the cone and a plane perpendicular to the
    #     # # axis of the cone passing through the location. This will have value 0.0 if the location is at the apex of the cone.
    #     # # semi_angle: the cone semi-angle, this is the angle between the axis of the Conical_surface and the generating line.
    #     #
    #     # semiAngle = AFS["SurfaceParams"][-1]
    #     # ParsedSurface["semiAngle"] = semiAngle
    #     # radius = AFS["SurfaceParams"][-2]
    #     # ParsedSurface["radius"] = radius
    #     #
    #     # # axisPoint, axisDir, refDir = axis2Placement3D_2(
    #     # #     AFS["SurfaceParams"][1], STEP_entities
    #     # # )
    #     # # ParsedSurface["axisPoint"] = axisPoint
    #     # # ParsedSurface["axisDir"] = axisDir
    #     # # ParsedSurface["refDir"] = refDir
    #     #
    #     # axisPoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][1], STEP_entities)
    #     # ParsedSurface["axisPoint"] = axisPoint
    #     # ParsedSurface["normDir"] = normDir
    #     # ParsedSurface["refDir"] = refDir # axisDir, refDir guaranteed unit normals?
    #     # # auxDir = np.cross(normDir, refDir)
    #     # # auxDir = auxDir / np.linalg.norm(auxDir)
    #     # # ParsedSurface["auxDir"] = auxDir
    #     #
    #     # # get adjacent disp between axisPoint and cone/triangle apex
    #     # axisPoint_apexPoint_disp = radius/np.tan(semiAngle)
    #     # # find point at cone/triangle apex
    #     # apexPoint = axisPoint + axisPoint_apexPoint_disp * normDir
    #     # # create a vector in the base plane that defines radius,
    #     # # and in the plane of axisPoint, apex & centroid
    #     # centroidAxisPointDir = axisPoint - centroid
    #     # centroidAxisPointDir = centroidAxisPointDir/np.linalg.norm(centroidAxisPointDir)
    #     # coneEdgeDir = np.cross(centroidAxisPointDir, normDir)
    #     # centroidConeEdgeDir = np.cross(coneEdgeDir, normDir)
    #     # minConeEdgePoint = axisPoint - radius * centroidConeEdgeDir
    #     # # maxConeEdgePoint = axisPoint + radius * centroidConeEdgeDir
    #     # # find minimum disp between centroid and segment defined between apexPoint, minConeEdgePoint
    #     #
    #     # surfaceMinPoint, surfaceMinPointDisp = pointSegmentMinDisp(centroid, apexPoint, minConeEdgePoint)
    #     #
    #     # # test for intersection of edges with a plane through surfaceMinPoint with normDir normal
    #     #
    #     # # 1x intersection => point tangential to plane thru edge
    #     # # 2x edge intersections => simple angle test on intersection plane
    #     # # 2x intersection on one edge => same
    #     # # Nx intersection => winding number solution?
    #     #
    #     # auxDirP = np.cross(centroidAxisDir, auxDir)  # surface axisDir
    #     # p1 = centroidProjAxisPoint + centroidAxisDir * centroidAxisDisp
    #     # p3 = centroidProjAxisPoint - centroidAxisDir * centroidAxisDisp
    #     # p2 = centroidProjAxisPoint + auxDirP * centroidAxisDisp
    #     # p4 = centroidProjAxisPoint - auxDirP * centroidAxisDisp
    #     #
    #     # intersectPlanePatch = verb.verb_eval_Make.fourPointSurface(p1, p2, p3, p4)
    #     # # intersectPlane := centroidProjAxisPoint, auxDirP
    #     # #------------------------------------------------------------------------------------<<<<<
    #     # AFS['ParsedSurface'] = ParsedSurface
    #     coneSurfaceParse(AFS)

    # if AFS["SurfaceTypeName"] == "PLANE":
    #     # ENTITY Plane
    #     #   SUBTYPE OF (Surface);
    #     #   position : Axis_placement_3d;
    #
    #     axisPoint, normDir, refDir = axis2Placement3D_2(AFS["SurfaceParams"][-1], STEP_entities)
    #     ParsedSurface["axisPoint"] = axisPoint
    #     ParsedSurface["normDir"] = normDir
    #     ParsedSurface["refDir"] = refDir # axisDir, refDir guaranteed unit normals?
    #
    #     rotSymFeatureEdgeSet = []
    #     rotSymFeature = []
    #
    #     # circle segments, ellipse segments have to be evaluated independentely
    #     # not sure how B-splines might work, max point box?
    #     circleIndex = [i for i, fobell in enumerate(AFS['FaceOuterBoundEdgeLoopList']) if
    #                    'CIRCLE' in fobell['typeName']]
    #     ellipseIndex = [i for i, fobell in enumerate(AFS['FaceOuterBoundEdgeLoopList']) if
    #                    'ELLIPSE' in fobell['typeName']]
    #
    #     # get list of vertices wrt surface
    #     surfaceVertices = []
    #     #surfaceTypes = []
    #     for fobell in AFS['FaceOuterBoundEdgeLoopList']:
    #         surfaceVertices.append(fobell['vertex1'])
    #         #surfaceTypes.append(fobell['typeName'])
    #
    #     if len(surfaceVertices) > 2:
    #         # strategy for planar polygonal surfaces is to divide into non-overlapping triangles
    #         # where corners are polygon vertices and test for nearest/furthest point on triangle
    #         minPoint, minPointDisp = planeMinMaxPoint(centroid, surfaceVertices, minmax="min")
    #         maxPoint, maxPointDisp = planeMinMaxPoint(centroid, surfaceVertices, minmax="max")
    #
    #         # not required to test within edges
    #
    #     if len(circleIndex) > 0:
    #         localCircleMinP = []
    #         for ci in circleIndex:
    #             fobell = AFS['FaceOuterBoundEdgeLoopList'][ci]
    #             cMinPoint, cMaxPoint = pointCircleMinMaxDisp(centroid,
    #                                                  fobell['axisPoint'],
    #                                                  fobell['normDir'],
    #                                                  fobell['radius'],
    #                                                  interior=True)
    #             boxPoints = curveEnclosingRectangle(fobell['vertex1'],
    #                                                 fobell['vertex2'],
    #                                                 fobell['axisPoint'],
    #                                                 fobell['normDir'],
    #                                                 fobell['radius'])
    #             rMinPoint, _ = planeMinMaxPoint(centroid, boxPoints, minmax="min")
    #
    #             # check circle-centre, circle maxima, box maxima are collinear, and BM-CC > CM-CC
    #             # cross product of (BM-CC) and (CM-CC) = 0 => points BM, CM and CC collinear.
    #             # CC < CM < BM => dot product of (BM-CC) and (CM-CC) is positive
    #             # and is less than the square of the distance between a and b.
    #             CM_CC = cMinPoint - fobell['axisPoint']
    #             BM_CC = rMinPoint - fobell['axisPoint']
    #             if (np.linalg.norm(BM_CC) < eps_STEP_AP21) and (np.linalg.norm(CM_CC) < eps_STEP_AP21):
    #                 # coincident case
    #                 localCircleMinP.append(cMinPoint)
    #             elif np.linalg.norm(np.cross(BM_CC, CM_CC)) < eps_STEP_AP21: #=tol
    #                 if 0 < np.dot(BM_CC, CM_CC) < (BM_CC*BM_CC):
    #                     localCircleMinP.append(cMinPoint)
    #
    #         if len(ellipseIndex) > 0:
    #             localEllipseMinP = []
    #             for ei in ellipseIndex:
    #                 fobell = AFS['FaceOuterBoundEdgeLoopList'][ei]
    #                 eMinPoint, eMaxPoint = pointEllipseMinMaxDisp(centroid,
    #                                                               fobell['refDir'],
    #                                                               fobell['auxDir'],
    #                                                               fobell['normDir'],
    #                                                               fobell['majorRadius'],
    #                                                               fobell['minorRadius'],
    #                                                               interior=True)
    #
    #                 boxPoints = curveEnclosingRectangle(fobell['vertex1'],
    #                                                     fobell['vertex2'],
    #                                                     fobell['axisPoint'],
    #                                                     fobell['normDir'],
    #                                                     fobell['radius'])
    #                 rMinPoint, rMinPointDisp = planeMinMaxPoint(centroid, boxPoints, minmax="min")
    #
    #                 EM_CC = eMinPoint - fobell['axisPoint']
    #                 BM_CC = rMinPoint - fobell['axisPoint']
    #                 if np.cross(BM_CC, EM_CC) < eps_STEP_AP21:  # =tol
    #                     if 0 < np.dot(BM_CC, EM_CC) < (BM_CC * BM_CC):
    #                         localEllipseMinP.append(eMinPoint)
    #
    #         if len(circleIndex) > 0:
    #             localCircleMinDisp = [np.linalg.norm(lcmp - centroid) for lcmp in localCircleMinP]
    #             localCircleMinP = localCircleMinP[localCircleMinDisp.index(min(localCircleMinDisp))]
    #             localCircleMinDisp = min(localCircleMinDisp)
    #             if len(surfaceVertices) > 2:
    #                 if localCircleMinDisp < minPointDisp:
    #                     minPoint = localCircleMinP
    #                     minPointDisp = localCircleMinDisp
    #             else:
    #                 minPoint = localCircleMinP
    #                 minPointDisp = localCircleMinDisp
    #
    #         if len(ellipseIndex) > 0:
    #             localEllipseMinDisp = [np.linalg.norm(lemp - centroid) for lemp in localEllipseMinP]
    #             localEllipseMinP = localEllipseMinP[localEllipseMinDisp.index(min(localEllipseMinDisp))]
    #             localEllipseMinDisp = min(localEllipseMinDisp)
    #             if len(surfaceVertices) > 2:
    #                 if localEllipseMinDisp < minPointDisp:
    #                     minPoint = localEllipseMinP
    #                     minPointDisp = localEllipseMinDisp
    #             else:
    #                 minPoint = localEllipseMinP
    #                 minPointDisp = localEllipseMinDisp
    #
    #     # note local minPoint of planar surface is generally a global minPoint
    #     ParsedSurface["minPoint"] = minPoint
    #     ParsedSurface["minPointCentroidDisp"] = minPointDisp
    #
    #     # maxPoint is determined from vertex and edge maxima, exception for rotSymMax/rotSymMin
    #     maxPoint = [fobell['maxPoint'] for fobell in AFS['FaceOuterBoundEdgeLoopList'] if 'maxPoint' in fobell.keys()]
    #     maxPointCentroidDisp = [fobell['maxPointCentroidDisp'] for fobell in AFS['FaceOuterBoundEdgeLoopList']if 'maxPointCentroidDisp' in fobell.keys()]
    #     if len(maxPoint) > 0 and len(maxPointCentroidDisp) > 0:
    #         maxPoint = maxPoint[maxPointCentroidDisp.index(max(maxPointCentroidDisp))]
    #         maxPointCentroidDisp = max(maxPointCentroidDisp)
    #         ParsedSurface["maxPoint"] = maxPoint
    #         ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp
    #
    #     # can't determine in advance whether rotSymMax or rotSymMin => requires subsequent comparison to adjacent nodes
    #     # same as maxPoint minPoint comparison
    #     rotSymEdges= [fobell['rotSymEdge'] for fobell in AFS['FaceOuterBoundEdgeLoopList'] if 'rotSymEdge' in fobell.keys()]
    #     if len(rotSymEdges) == 1:
    #         ParsedSurface["rotSymEdge"] = rotSymEdges[0]
    #     elif len(rotSymEdges) > 1:
    #         rsm_x = [rsm[0] for rsm in rotSymEdges]
    #         rsm_y = [rsm[1] for rsm in rotSymEdges]
    #         rsm_z = [rsm[2] for rsm in rotSymEdges]
    #         if (all(x == rsm_x[0] for x in rsm_x) and
    #                 all(y == rsm_y[0] for y in rsm_y) and
    #                 all(z == rsm_z[0] for z in rsm_z)):
    #             ParsedSurface["rotSymEdge"] = rotSymEdges[0]
    #         else:
    #             print("indeterminate rotSym")
    #
    #     # TEMPORARY TEST--------------------------------------------------------------------------------
    #
    #     maxNorm = (maxPoint - centroid) / np.linalg.norm(maxPoint - centroid)
    #     if not insideOutsideSurfaceTest(maxNorm, maxPoint, AFS):
    #         # find maxima among edges
    #         maxPointCentroidDisp = [edge['maxPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
    #         maxPoint = AFS["FaceOuterBoundEdgeLoopList"][maxPointCentroidDisp.index(max(maxPointCentroidDisp))]
    #         maxPointCentroidDisp = max(maxPointCentroidDisp)
    #     else:
    #         ParsedSurface["maxPoint"] = maxPoint
    #         ParsedSurface["maxPointCentroidDisp"] = maxPointCentroidDisp
    #
    #     minNorm = (minPoint - centroid) / np.linalg.norm(minPoint - centroid)
    #     if not insideOutsideSurfaceTest(minNorm, minPoint, AFS):
    #         # find maxima among edges
    #         minPointCentroidDisp = [edge['minPointCentroidDisp'] for edge in AFS["FaceOuterBoundEdgeLoopList"]]
    #         minPoint = AFS["FaceOuterBoundEdgeLoopList"][minPointCentroidDisp.index(min(minPointCentroidDisp))]
    #         minPointCentroidDisp = min(minPointCentroidDisp)
    #     else:
    #         ParsedSurface["minPoint"] = minPoint
    #         ParsedSurface["minPointCentroidDisp"] = minPointCentroidDisp
    #
    #     AFS['ParsedSurface'] = ParsedSurface

    # strategy for cylinders & cones (& rotationally-symmetric spheres, surfaces of revolution or ellipsoids), determine curved edges
    # if rot axis is between origin and surface, maxima lies on a line segment between maxima points on curved edges
    # if rot axis is beyond surface (concave), minima lies on a line segment between minima points on curved edges
    # otherwise maxima/minima lie on edges.

        # extract edge entities from a surface and parse
        # parsedEdge = {}
        # if (STEP_entities[ref2index(el)].type_name == "ORIENTED_EDGE"):
        #     # if not hasattr(STEP_entities[ref2index(el)].type_name, 'ORIENTED_EDGE'): ----------------------------------
        #     VOR = "F" in STEP_entities[ref2index(el)].params[-1]  #  vertexOrderReversed
        #     parsedEdge["VOR"] = VOR
        #
        # edgeCurveRef = STEP_entities[ref2index(el)].params[-2]
        # assert STEP_entities[ref2index(edgeCurveRef)].type_name == "EDGE_CURVE"
        # # edge_geometry: the Curve defining the geometric shape of the edge.
        # # same_sense: a BOOLEAN variable giving the relationship between the topological sense
        # #             of the edge and the parametric sense of the curve.
        # edgeCurveParams = STEP_entities[ref2index(edgeCurveRef)].params
        # # edgeSameSense = 'T' in edgeCurveParams[-1]; parsedEdge['sameSense'] = edgeSameSense
        # edgeCurveParams = cleanSubRefs(edgeCurveParams)
        #
        # if (STEP_entities[ref2index(edgeCurveParams[0])].type_name == "VERTEX_POINT"):
        #     vertex1 = STEP_entities[ref2index(edgeCurveParams[0])].params
        #     vertex1 = CP2point(STEP_entities, cleanSubRefs(vertex1)[0])
        #     if not VOR:
        #         parsedEdge["vertex1"] = vertex1
        #         parsedEdge["vertex1ref"] = STEP_entities[ref2index(edgeCurveParams[0])].ref
        #         vertex1centroidDisp = np.linalg.norm(centroid - vertex1)
        #         parsedEdge["vertex1centroidDisp"] = vertex1centroidDisp
        #     else:
        #         parsedEdge["vertex2"] = vertex1
        #         parsedEdge["vertex2ref"] = STEP_entities[ref2index(edgeCurveParams[0])].ref
        #         vertex2centroidDisp = np.linalg.norm(centroid - vertex2)
        #         parsedEdge["vertex2centroidDisp"] = vertex2centroidDisp
        #
        # if (STEP_entities[ref2index(edgeCurveParams[1])].type_name == "VERTEX_POINT"):
        #     vertex2 = STEP_entities[ref2index(edgeCurveParams[1])].params
        #     vertex2 = CP2point(STEP_entities, cleanSubRefs(vertex2)[0])
        #     if not VOR:
        #         parsedEdge["vertex2"] = vertex2
        #         parsedEdge["vertex2ref"] = STEP_entities[ref2index(edgeCurveParams[1])].ref
        #         vertex2centroidDisp = np.linalg.norm(centroid - vertex2)
        #         parsedEdge["vertex2centroidDisp"] = vertex2centroidDisp
        #     else:
        #         parsedEdge["vertex1"] = vertex2
        #         parsedEdge["vertex1ref"] = STEP_entities[ref2index(edgeCurveParams[0])].ref
        #         vertex2centroidDisp = np.linalg.norm(centroid - vertex2)
        #         parsedEdge["vertex1centroidDisp"] = vertex2centroidDisp
        #
        # edgeType = STEP_entities[ref2index(edgeCurveParams[2])].params
        # edgeTypeName = STEP_entities[ref2index(edgeCurveParams[2])].type_name
        # edgeRef = STEP_entities[ref2index(edgeCurveParams[2])].ref
        # parsedEdge["typeName"] = edgeTypeName
        # parsedEdge["edgeRef"] = edgeRef
        #
        # if edgeTypeName not in [
        #     "LINE",
        #     "CIRCLE",
        #     "ELLIPSE",
        #     "B_SPLINE_CURVE_WITH_KNOTS",
        # ]:
        #     print("Untreated edge type: " + edgeTypeName)

        # if edgeTypeName == "LINE":
        # # # ENTITY Line
        # # #   SUBTYPE OF (Curve);
        # # #   point : Cartesian_point;
        # # #   line_direction : Direction;
        # #     cleanEdgeType = cleanSubRefs(edgeType)
        # #     for cet in cleanEdgeType:
        # #         if (STEP_entities[ref2index(cet)].type_name == "CARTESIAN_POINT"):
        # #             linePoint = CP2point(STEP_entities, cet)
        # #             parsedEdge["linePoint"] = linePoint
        # #
        # #         if (STEP_entities[ref2index(cet)].type_name == "VECTOR"):
        # #             lineVectorDisp = STEP_entities[ref2index(cet)].params[-1]
        # #             parsedEdge["lineVectorDisp"] = lineVectorDisp
        # #             lineVectorDir = STEP_entities[ref2index(cet)].params[-2]
        # #             if (STEP_entities[ref2index(lineVectorDir)].type_name == "DIRECTION"):
        # #                 lineVectorDir = CP2point(STEP_entities, lineVectorDir)
        # #                 parsedEdge["lineVectorDir"] = lineVectorDir
        # #
        # #     minPoint, minPointCentroidDisp, vertexExtrema = pointSegmentMinDisp(centroid, vertex1, vertex2)
        # #     parsedEdge["minPoint"] = minPoint
        # #     parsedEdge["minPointCentroidDisp"] = minPointCentroidDisp
        # #     parsedEdge['vertexExtremaMin'] = vertexExtrema
        # #
        # #     maxPoint = vertex1
        # #     maxPointCentroidDisp = vertex1centroidDisp
        # #     if vertex2centroidDisp > vertex1centroidDisp: # multiple maxima for equidistant vector?
        # #         maxPoint = vertex2
        # #         maxPointCentroidDisp = vertex2centroidDisp
        # #     parsedEdge["maxPoint"] = maxPoint
        # #     parsedEdge["maxPointCentroidDisp"] = maxPointCentroidDisp
        # #     parsedEdge['vertexExtremaMax'] = True
        #
        #     parsedEdge.update(lineParse(edgeType, STEP_entities)) #parsedEdge_ = lineParse(edgeType, STEP_entities)
        #
        # if edgeTypeName == "CIRCLE":
        #     # radius = edgeType[-1]
        #     # parsedEdge["radius"] = radius
        #     #
        #     # axisPoint, normDir, refDir = axis2Placement3D_2(edgeType[-2], STEP_entities)
        #     # parsedEdge["axisPoint"] = axisPoint
        #     # parsedEdge["normDir"] = normDir
        #     # parsedEdge["refDir"] = refDir
        #     # # axisDir, refDir guaranteed unit normals?
        #     # auxDir = np.cross(normDir, refDir)
        #     # auxDir = auxDir / np.linalg.norm(auxDir)
        #     # parsedEdge["auxDir"] = auxDir
        #     #
        #     # # min/max point on arc wrt centroid
        #     # minPoint, maxPoint = pointCircleMinMaxDisp(
        #     #     centroid,
        #     #     axisPoint,
        #     #     normDir,
        #     #     radius
        #     # )
        #     #
        #     # v1 = parsedEdge["vertex1"]
        #     # v2 = parsedEdge["vertex2"]
        #     #
        #     # if np.linalg.norm(minPoint - axisPoint) < eps_STEP_AP21:
        #     #     # case of minPoint equidistant from circle centre - edge case
        #     #     # assign rotSymMin feature point at centre
        #     #     parsedEdge["rotSymEdge"] = axisPoint
        #     #
        #     # else:
        #     #     if pointInArc(minPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
        #     #         parsedEdge["minPoint"] = minPoint
        #     #         parsedEdge["minPointCentroidDisp"] = np.linalg.norm(minPoint - centroid)
        #     #         parsedEdge['vertexExtremaMin'] = False
        #     #     else: # one of the vertices is a minima
        #     #         parsedEdge['vertexExtremaMin'] = True
        #     #         if parsedEdge['vertex1centroidDisp'] < parsedEdge['vertex2centroidDisp']:
        #     #             parsedEdge["minPoint"] = parsedEdge['vertex1']
        #     #             parsedEdge["minPointCentroidDisp"] = parsedEdge['vertex1centroidDisp']
        #     #         else:
        #     #             parsedEdge['minPoint'] = parsedEdge['vertex2']
        #     #             parsedEdge["minPointCentroidDisp"] = parsedEdge['vertex2centroidDisp']
        #     #
        #     # if np.linalg.norm(maxPoint - axisPoint) < eps_STEP_AP21:
        #     #     # case of maxPoint equidistant from circle centre - edge case
        #     #     # assign rotSymMax feature point at centre
        #     #     parsedEdge["rotSymEdge"] = axisPoint
        #     #
        #     # else:
        #     #     if pointInArc(maxPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
        #     #         parsedEdge["maxPoint"] = maxPoint
        #     #         parsedEdge["maxPointCentroidDisp"] = np.linalg.norm(maxPoint - centroid)
        #     #         parsedEdge['vertexExtremaMax'] = False
        #     #     else: # one of the vertices is a maxima
        #     #         parsedEdge['vertexExtremaMax'] = True
        #     #         if parsedEdge['vertex1centroidDisp'] > parsedEdge['vertex2centroidDisp']:
        #     #             parsedEdge["maxPoint"] = parsedEdge['vertex1']
        #     #             parsedEdge["maxPointCentroidDisp"] = parsedEdge['vertex1centroidDisp']
        #     #         else:
        #     #             parsedEdge['maxPoint'] = parsedEdge['vertex2']
        #     #             parsedEdge["maxPointCentroidDisp"] = parsedEdge['vertex2centroidDisp']
        #
        #     parsedEdge.update(circleParse(edgeType, STEP_entities))
        #
        # if edgeTypeName == "ELLIPSE":
        #     # # first_semi_axis: half the length of the first diameter of the Ellipse.
        #     # # second_semi_axis: half the length of the second diameter of the Ellipse.
        #     #
        #     # cleanEdgeType = cleanSubRefs(edgeType)
        #     #
        #     # majorRadius = cleanEdgeType[-2]
        #     # parsedEdge["majorRadius"] = majorRadius
        #     #
        #     # minorRadius = cleanEdgeType[-1]
        #     # parsedEdge["minorRadius"] = minorRadius
        #     #
        #     # # axisPoint, auxDir, refDir = axis2Placement3D_2(cleanEdgeType[-3], STEP_entities)
        #     # # parsedEdge["axisPoint"] = axisPoint
        #     # # parsedEdge["auxDir"] = auxDir
        #     # # parsedEdge["refDir"] = refDir
        #     # # normDir = np.cross(auxDir, refDir)
        #     # # normDir = normDir / np.linalg.norm(normDir)
        #     # # parsedEdge["normDir"] = normDir ------------------------------------------------------------------------------------------?????
        #     #
        #     # axisPoint, normDir, refDir = axis2Placement3D_2(edgeType[-2], STEP_entities)
        #     # parsedEdge["axisPoint"] = axisPoint
        #     # parsedEdge["normDir"] = normDir
        #     # parsedEdge["refDir"] = refDir
        #     # # axisDir, refDir guaranteed unit normals?
        #     # auxDir = np.cross(normDir, refDir)
        #     # auxDir = auxDir / np.linalg.norm(auxDir)
        #     # parsedEdge["auxDir"] = auxDir
        #     #
        #     # minPoint, maxPoint = pointEllipseMinMaxDisp(
        #     #     centroid,
        #     #     axisPoint,
        #     #     refDir,
        #     #     auxDir,
        #     #     normDir,
        #     #     majorRadius,
        #     #     minorRadius
        #     # )
        #     #
        #     # # test if minPoint, maxPoint is on segment between vertex1 & vertex2
        #     # v1 = pel["vertex1"]
        #     # v2 = pel["vertex2"]
        #     #
        #     # if pointInArc(minPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=False):
        #     #     parsedEdge["minPoint"] = minPoint
        #     #     parsedEdge['vertexExtremaMin'] = False
        #     # else: # one of the vertices is a minima
        #     #     parsedEdge['vertexExtremaMin'] = True
        #     #     if parsedEdge['vertex1centroidDisp'] < parsedEdge['vertex2centroidDisp']:
        #     #         parsedEdge["minPoint"] = parsedEdge['vertex1']
        #     #         parsedEdge["minPointCentroidDisp"] = parsedEdge['vertex1centroidDisp']
        #     #     else:
        #     #         parsedEdge['minPoint'] = parsedEdge['vertex2']
        #     #         parsedEdge["minPointCentroidDisp"] = parsedEdge['vertex2centroidDisp']
        #     #
        #     # if pointInArc(maxPoint, v1, v2, refDir, auxDir, normDir, axisPoint, rotSym=True):
        #     #     parsedEdge["maxPoint"] = maxPoint
        #     #     parsedEdge['vertexExtremaMax'] = False
        #     # else: # one of the vertices is a maxima
        #     #     parsedEdge['vertexExtremaMax'] = True
        #     #     if parsedEdge['vertex1centroidDisp'] > parsedEdge['vertex2centroidDisp']:
        #     #         parsedEdge["maxPoint"] = parsedEdge['vertex1']
        #     #         parsedEdge["maxPointCentroidDisp"] = parsedEdge['vertex1centroidDisp']
        #     #     else:
        #     #         parsedEdge['maxPoint'] = parsedEdge['vertex2']
        #     #         parsedEdge["maxPointCentroidDisp"] = parsedEdge['vertex2centroidDisp']
        #     #
        #     # # # convert centroid, cMinPoint, cMaxPoint vertex1 & vertex2 to local coordinates with axis as z
        #     # # # convert to polar coords
        #     # #
        #     # # # axisDir, refDir guaranteed unit normals?
        #     # # eNorm = np.cross(auxDir, refDir)
        #     # # eNorm = eNorm / np.linalg.norm(eNorm)
        #     # # rotM = np.array([refDir, eNorm, auxDir])
        #     # #
        #     # # eMinPointUV = np.matmul(rotM, eMinPoint - axisPoint)
        #     # # eMaxPointUV = np.matmul(rotM, eMaxPoint - axisPoint)
        #     # # vertex1UV = np.matmul(rotM, vertex1 - axisPoint)
        #     # # vertex2UV = np.matmul(rotM, vertex2 - axisPoint)
        #     # #
        #     # # # sanity check that Z-values are identical
        #     # # zTest = [eMinPointUV[2], eMaxPointUV[2], vertex1UV[2], vertex2UV[2]]
        #     # # if not np.isclose(max(zTest), min(zTest)):
        #     # #     print(max(zTest) - min(zTest))
        #     # #
        #     # # eMinPointUV_angle = np.arctan2(eMinPointUV[1], eMinPointUV[0])
        #     # # eMaxPointUV_angle = np.arctan2(eMaxPointUV[1], eMaxPointUV[0])
        #     # # vertex1UV_angle = np.arctan2(vertex1UV[1], vertex1UV[0])
        #     # # vertex2UV_angle = np.arctan2(vertex2UV[1], vertex2UV[0])
        #     # #
        #     # # if not VOR:
        #     # #     if (
        #     # #         (eMinPointUV_angle > vertex1UV_angle)
        #     # #         and (eMinPointUV_angle < vertex2UV_angle)
        #     # #     ) or (vertex1UV_angle == vertex2UV_angle):
        #     # #         parsedEdge["minPoint"] = eMinPoint
        #     # #     else:
        #     # #         if np.linalg.norm(vertex1 - centroid) < np.linalg.norm(
        #     # #             vertex2 - centroid
        #     # #         ):
        #     # #             parsedEdge["minPoint"] = vertex1
        #     # #         else:
        #     # #             parsedEdge["minPoint"] = vertex2
        #     # #
        #     # #     if (
        #     # #         (eMaxPointUV_angle > vertex1UV_angle)
        #     # #         and (eMaxPointUV_angle < vertex2UV_angle)
        #     # #     ) or (vertex1UV_angle == vertex2UV_angle):
        #     # #         parsedEdge["maxPoint"] = eMaxPoint
        #     # #     else:
        #     # #         if np.linalg.norm(vertex1 - centroid) > np.linalg.norm(
        #     # #             vertex2 - centroid
        #     # #         ):
        #     # #             parsedEdge["maxPoint"] = vertex1
        #     # #         else:
        #     # #             parsedEdge["maxPoint"] = vertex2
        #     # #
        #     # # else:  # vertex order reversed
        #     # #     if (eMinPointUV_angle > vertex2UV_angle) and (
        #     # #         eMinPointUV_angle < vertex1UV_angle
        #     # #     ):
        #     # #         parsedEdge["minPoint"] = cMinPoint
        #     # #     else:
        #     # #         if np.linalg.norm(vertex1 - centroid) < np.linalg.norm(
        #     # #             vertex2 - centroid
        #     # #         ):
        #     # #             parsedEdge["minPoint"] = vertex1
        #     # #         else:
        #     # #             parsedEdge["minPoint"] = vertex2
        #     # #
        #     # #     if (eMaxPointUV_angle > vertex2UV_angle) and (
        #     # #         eMaxPointUV_angle < vertex1UV_angle
        #     # #     ):
        #     # #         parsedEdge["maxPoint"] = cMaxPoint
        #     # #     else:
        #     # #         if np.linalg.norm(vertex1 - centroid) > np.linalg.norm(
        #     # #             vertex2 - centroid
        #     # #         ):
        #     # #             parsedEdge["maxPoint"] = vertex1
        #     # #         else:
        #     # #             parsedEdge["maxPoint"] = vertex2
        #
        #     parsedEdge.update(ellipseParse(edgeType, STEP_entities))
        #
        # # note that a NURBS edge may have more than one local maxima/minima
        # if edgeTypeName == "B_SPLINE_CURVE_WITH_KNOTS":
        #     # # ENTITY B_Spline_Curve
        #     # # 	SUPERTYPE OF ((ONEOF (Uniform_Curve, B_Spline_Curve_With_Knots, Quasi_Uniform_Curve, Bezier_Curve) ANDOR Rational_B_Spline_Curve))
        #     # # 	SUBTYPE OF (Bounded_Curve);
        #     # # 	degree : INTEGER;
        #     # # 	control_points_list : LIST [2:?] OF Cartesian_Point;
        #     # # 	curve_form : B_Spline_Curve_Form;
        #     # # 	closed_curve : LOGICAL;
        #     # # 	self_intersect : LOGICAL;
        #     # # ENTITY B_Spline_Curve_With_Knots
        #     # # 	SUBTYPE OF (B_Spline_Curve);
        #     # # 	knot_multiplicities : LIST [2:?] OF INTEGER;
        #     # # 	knots : LIST [2:?] OF Parameter_Value;
        #     # # 	knot_spec : Knot_Type;
        #     #
        #     # if type(edgeType[0]) == str: # check for name string
        #     #     offset = 1
        #     # else:
        #     #     offset = 0
        #     # curveDegree = edgeType[offset]
        #     # controlPointsListRefs = edgeType[offset+1]
        #     # curveForm = edgeType[offset+2]
        #     # closedCurve = "T" in edgeType[offset+3]
        #     # selfIntersect = "T" in edgeType[offset+4]
        #     # knotMultiplicities = edgeType[offset+5]
        #     # knots = edgeType[offset+6]
        #     # knotSpec = edgeType[offset+7]
        #     #
        #     # parsedEdge["curveDegree"] = curveDegree
        #     # parsedEdge["controlPointsListRefs"] = controlPointsListRefs
        #     # parsedEdge["curveForm"] = curveForm
        #     # parsedEdge["closedCurve"] = closedCurve
        #     # parsedEdge["selfIntersect"] = selfIntersect
        #     # parsedEdge["knotMultiplicities"] = knotMultiplicities
        #     # parsedEdge["knots"] = knots
        #     # parsedEdge["knotSpec"] = knotSpec
        #     #
        #     # # extract control points
        #     # controlPointsList = []
        #     # for cpl in controlPointsListRefs:
        #     #     if (STEP_entities[ref2index(cpl)].type_name == "CARTESIAN_POINT"):
        #     #         controlPoint = CP2point(STEP_entities, cleanSubRefs(cpl)[0])
        #     #         controlPointsList.append(controlPoint)
        #     # parsedEdge["controlPointsList"] = controlPointsList
        #     #
        #     # knotvector = []
        #     # for kml in range(len(knotMultiplicities)):
        #     #     for i in range(knotMultiplicities[kml]):
        #     #         knotvector.append(knots[kml])
        #     # parsedEdge["knotvector"] = knotvector
        #     #
        #     # if (len(knotvector) - len(controlPointsList) - curveDegree) != 1:
        #     #     print("maldefined B-spline!")
        #     #
        #     # # knot multiplicity must be explicit for verb library
        #     # BsplineKnotCurve = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(
        #     #     curveDegree, knotvector, [cp.tolist() for cp in controlPointsList]
        #     # )
        #     # minPoint = BsplineKnotCurve.closestPoint(centroid)
        #     # minPoint = np.array([minPoint[0], minPoint[1], minPoint[2]])
        #     #
        #     # # Newton root finding over a bunch of points probably makes more sense for mapping local maxima/minima
        #     # # TODO: pointNURBSmaxDisp() may be more accurate for minima
        #     # maxPointU = pointNURBSmaxDisp(
        #     #     centroid, BsplineKnotCurve._data
        #     # )  # TODO: only finds one maxima
        #     # maxPoint = verb.verb_eval_Eval.dehomogenize(
        #     #     verb.verb_eval_Eval.curvePoint(BsplineKnotCurve._data, maxPointU)
        #     # )
        #     # maxPoint = np.array([maxPoint[0], maxPoint[1], maxPoint[2]])
        #     #
        #     # minPointCentroidDisp = np.linalg.norm(minPoint - centroid)
        #     # parsedEdge["minPointCentroidDisp"] = minPointCentroidDisp
        #     # parsedEdge["minPoint"] = minPoint
        #     # parsedEdge['vertexExtremaMin'] = False
        #     # # if vertex1centroidDisp < centroidMinPointDisp:
        #     # #     minPoint = vertex1
        #     # # if vertex2centroidDisp < centroidMinPointDisp:
        #     # #     minPoint = vertex2
        #     # if parsedEdge['vertex1centroidDisp'] < minPointCentroidDisp:
        #     #     parsedEdge["minPoint"] = parsedEdge['vertex1']
        #     #     parsedEdge["minPointCentroidDisp"] = parsedEdge['vertex1centroidDisp']
        #     #     parsedEdge['vertexExtremaMin'] = True
        #     # if parsedEdge['vertex2centroidDisp'] < minPointCentroidDisp:
        #     #     parsedEdge["minPoint"] = parsedEdge['vertex2']
        #     #     parsedEdge["minPointCentroidDisp"] = parsedEdge['vertex2centroidDisp']
        #     #     parsedEdge['vertexExtremaMin'] = True
        #     #
        #     # maxPointCentroidDisp = np.linalg.norm(maxPoint - centroid)
        #     # parsedEdge["maxPointCentroidDisp"] = maxPointCentroidDisp
        #     # parsedEdge["maxPoint"] = maxPoint
        #     # parsedEdge['vertexExtremaMax'] = False
        #     # # if vertex1centroidDisp > centroidMaxPointDisp:
        #     # #     maxPoint = vertex1
        #     # # if vertex2centroidDisp > centroidMaxPointDisp:
        #     # #     maxPoint = vertex2
        #     # # parsedEdge["maxPoint"] = maxPoint
        #     # if parsedEdge['vertex1centroidDisp'] > maxPointCentroidDisp:
        #     #     parsedEdge["maxPoint"] = parsedEdge['vertex1']
        #     #     parsedEdge["maxPointCentroidDisp"] = parsedEdge['vertex1centroidDisp']
        #     #     parsedEdge['vertexExtremaMax'] = True
        #     # if parsedEdge['vertex2centroidDisp'] > maxPointCentroidDisp:
        #     #     parsedEdge["maxPoint"] = parsedEdge['vertex2']
        #     #     parsedEdge["maxPointCentroidDisp"] = parsedEdge['vertex2centroidDisp']
        #     #     parsedEdge['vertexExtremaMax'] = True
        #
        #     parsedEdge.update(BsplineCurveWithKnotsParse(edgeType, STEP_entities))
