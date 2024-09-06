# continuation of WIP_2jun24.py, addressing complex edges with multiple maxima/minima
# introduction of dict structure, featureExtrema, to maintain local extrema (x, y) value, centroid disp value, u spline value, v spline surface optional value
# remove pre-calculation of ParsedEdge['vertex1extremaMax'], ParsedEdge['vertex2extremaMax'] flags

# Python STEP parser, NIST STEPcode using NIST BSD-3 license, STEPcode v0.8 -- github.com/stepcode/stepcode
# https://github.com/stepcode/stepcode/blob/develop/src/exp2python/python/README.md
from stepcode.Part21 import Parser
import os, glob, re
from pprint import pprint # only debug printout

# NURBS-Python (geomdl) is a self-contained, object-oriented pure Python B-Spline and NURBS library
# geomdl v5.3.1
# Please refer to the following DOI link to access the article:
# https://doi.org/10.1016/j.softx.2018.12.005

# @article{bingol2019geomdl,
#   title={{NURBS-Python}: An open-source object-oriented {NURBS} modeling framework in {Python}},
#   author={Bingol, Onur Rauf and Krishnamurthy, Adarsh},
#   journal={{SoftwareX}},
#   volume={9},
#   pages={85--94},
#   year={2019},
#   publisher={Elsevier}
# }

# from geomdl import BSpline
from geomdl import NURBS
from geomdl import compatibility
from geomdl import utilities
#from geomdl import evaluators
# import copy

# from collections import namedtuple  # PyCharm bug will flag this as type error
# Point = namedtuple("Point", "x y z")

import numpy as np

# eps machine precision constant
eps = np.finfo(float).eps
eps_STEP_AP21 = 1e-6  # STEP precision seems to end here

# eps_bspline = 1E-10 # precision used in bspline/NURB surface subroutines
# machine_EPS = np.finfo(float).eps # float comparison


def FreeCADpointSyntax(A, color=(0., 0., 0.)):
    '''print the instructions to create a FreeCAD model of a list of points'''
    print('import FreeCAD as App')
    print('import Draft')
    for i, xyz in enumerate(A):
        print('P{:} = Draft.make_point( {:1.8f}, {:1.8f}, {:1.8f}, color=({:1.8f}, {:1.8f}, {:1.8f}))'.format(i, xyz[0],
                                                                                                              xyz[1],
                                                                                                              xyz[2],
                                                                                                              color[0],
                                                                                                              color[1],
                                                                                                              color[2]))
    print('App.ActiveDocument.recompute()')
    # print('setview()')


# todo: dict data arrays like ParsedEdge should become classes, once all required fields & methods identified

# ==========================general geometric functions below, shunt into an object I guess

def angleAxis2RotMat(rotAxis, theta):
    """
    Affine transformation matrix from 3x3 rotation matrix and displacement vector.
    3x3 rotation matrix derived from CCW angle theta around axis
    From: http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle & transform3d

    :param rotAxis: axis of rotation vector axis emanating from origin
    :param theta: scalar radian angle of rotation
    :return: 4x4 matrix representing affine transformation
    """

    # normalise to unit vector
    rotNorm = np.linalg.norm(rotAxis)
    if rotNorm > eps:
        rotAxis = rotAxis / rotNorm
        rotAxis = rotAxis.tolist()

    # if rotAxis.__class__.__name__ == "list":
    #     rotAxis = Point(rotAxis[0], rotAxis[1], rotAxis[2])

    s = np.sin(theta)
    c = np.cos(theta)
    C = 1 - c
    xs = rotAxis.x * s
    ys = rotAxis.y * s
    zs = rotAxis.z * s
    xC = rotAxis.x * C
    yC = rotAxis.y * C
    zC = rotAxis.z * C
    xyC = rotAxis.x * yC
    yzC = rotAxis.y * zC
    zxC = rotAxis.z * xC
    xxC = rotAxis.x * xC
    yyC = rotAxis.y * yC
    zzC = rotAxis.z * zC

    return np.array(
        [
            [xxC + c, xyC - zs, zxC + ys],
            [xyC + zs, yyC + c, yzC - xs],
            [zxC - ys, yzC + xs, zzC + c],
        ]
    )


def rotationMatrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about the given axis by theta radians.

    :param axis: axis defined as a vector point relative to cartesian origin
    :param theta: rotation angle in radians
    :return: rotation matrix
    """

    # Euler-Rodrigues formula
    # http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    # from numpy import asarray, array, dot
    # from math import sqrt, cos, sin

    # axis = np.asarray(axis)
    # theta = np.asarray(theta)
    SDAA = np.sqrt(np.dot(axis, axis))
    if SDAA == 0.0:
        raise RuntimeError("rotationMatrix() zero div error (0,0,0)?")

    axis = axis / SDAA
    a = np.cos(theta / 2)
    b, c, d = -axis * np.sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


def arrayTypeCast(xyz):
    if type(xyz) == list:
        xyz = np.array([xyz[0], xyz[1], xyz[2]])
    # elif type(xyz) == Point:
    #     xyz = np.array([xyz.x, xyz.y, xyz.z])
    elif type(xyz) == np.ndarray:
        pass
    else:
        raise RuntimeError("input type error")
    return xyz


def cart2pol(x, y):
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def array3x1match(a, A, atol=eps):
    '''test whether a numpy array of 3 elements (a) exists within a list of 3 element arrays'''
    T = [np.isclose(a[0], b[0], atol=eps) and np.isclose(a[1], b[1], atol=eps) and np.isclose(a[2], b[2], atol) for
         b in A]
    return any(T)

#debugA=[]
def medianOfFeaturePoints(Surfaces):
    """get barycentre of all AFS feature points"""

    #global debugA
    outermostPoints = []

    for afs in Surfaces:
        if afs['ParsedSurface'].get('pointFeature') is not None:
            if afs['ParsedSurface']['pointFeature'].get('xyz') is not None:
                [outermostPoints.append(mp) for mp in afs['ParsedSurface']['pointFeature']['xyz']]
        # if afs['ParsedSurface'].get('axisPoint') is not None:
        #     outermostPoints.append(afs['ParsedSurface']['axisPoint'])
        # print('afs[ParsedSurface][axisPoint]')
        # print(afs['ParsedSurface']['axisPoint'])
        for edge in afs['outerBoundEdgeLoop']:
            if edge.get('pointFeature') is not None:
                if edge['pointFeature'].get('xyz') is not None:
                    [outermostPoints.append(mp) for mp in edge['pointFeature']['xyz']]
            if edge.get('axisPoint') is not None:
                outermostPoints.append(edge['axisPoint']) # ----------------------------------------????
            # print('axisPoint')
            # print(edge['axisPoint'])

    #debugA.append(outermostPoints)

    return medianPoint(outermostPoints)


def medianPoint(pointArray):
    """
    Returns centroid determined as the median of cartesian values in an input array of points

    :param pointArray: array of Cartesian points
    :return: point
    """

    # from statistics import median # apparently slow/accurate
    # from numpy import median
    if len(pointArray) < 2:
        #raise RuntimeError("medianPoint() malformed input")
        print("medianPoint() inadequate input")
        return np.array([])

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


def radiusCentre3points(p1, p2, p3):
    # radius + centre from 3 point circle

    t = p2 - p1
    u = p3 - p1
    v = p3 - p2

    w = np.cross(t, u)  # triangle normal
    wsl = np.linalg.norm(w)
    if (wsl < 10e-14): return False, False  # triangle area too small (additionally check points for colinearity)

    wsl = np.dot(w, w)
    iwsl2 = 1. / (2. * wsl)
    tt = np.dot(t, t)
    uu = np.dot(u, u)

    circCenter = p1 + (u * tt * (np.dot(u, v)) - t * uu * (np.dot(t, v))) * iwsl2
    circRadius = np.sqrt(tt * uu * (np.dot(v, v)) * iwsl2 * 0.5)
    # circAxis   = w / np.sqrt(wsl)

    return circRadius, circCenter


def splineToPolyline(c, curveSampleFactor=2):
    '''
    Get set of sampled points from geomdl spline curve object
    :param c: geomdl curve object
    :param curveSampleFactor:
    :return cartesian points, parametric values
    '''
    numSamples = c.ctrlpts_size * c.degree * curveSampleFactor
    span = (c.knotvector[-1] - c.knotvector[0]) / (numSamples - 1)
    kvs = np.array([c.knotvector[0] + (span * i) for i in range(0, numSamples)])
    pts = c.evaluate_list(kvs)
    return (pts, kvs)


def rationalSurfaceExtremaParam_5(S, p,
                                  maxSearch=True,
                                  localExtrema=False,
                                  curvatureTest=False,
                                  uv_xyz=True,
                                  eps1=0.0001,
                                  eps2=0.0005,
                                  deltaFactor=4,
                                  eps_bspline=1E-10):
    '''
    Use either Newton-Raphson, or a hillclimbing neighbour search to find minimal or maximal distance on a surface, S, from a supplied point, p
    - not tested with surfaces closed along U & V
    maxSearch: search for most distant point from p if true, else nearest if false
    localExtrema: find any local extrema which differs from surrounding sample points, tolerance not sufficiently defined
    curvatureTest: rather than displacement of sampled points to provide seed locations for N-R/hillclimb search,
    use a combination of the centres of curvature at any point across U-axis and V-axis to identify C2 inflection point.
    uv_xyz: return u,v normalised values or cartesian x,y,z values
    eps1: point coincidence limit, halt N-R on reaching this minima
    eps2: cosine angle limit, halt N-R on reaching this minima
    delta: sample interval dividing surface into points
    '''

    # easy to find a minimum, but a maximum found through the original hillclimbing algorithm will settle on local maxima
    # for surfaces, minimize the following:
    #
    # f = Su(u,v) * r = 0
    # g = Sv(u,v) * r = 0
    #
    #   where r = S(u,v) - P
    #
    # Requires Newton-Raphson iteration, objective function is vector valued
    #
    #     J d = k
    #
    #     d =   [ u* - u, v* - v ] (alternatively u[i+1] -  u[i], v[i+1]- v[i])
    #     k = - [ f(u,v), g(u,v) ]
    #     J =     |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
    #              Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
    #
    # halting conditions:
    # point coincidence
    #
    #         |S(u,v) - p| < e1
    #
    # cosine
    #
    #    |Su(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Su(u,v)| |S(u,v) - P|
    #
    #    |Sv(u,v)*(S(u,v) - P)|
    #    ----------------------  < e2
    #    |Sv(u,v)| |S(u,v) - P|
    #
    # 1) first check 2 & 3
    # 2) if at least one of these is not, compute new value, otherwise halt
    # 3) ensure the parameter stays within range
    #     * if not closed, don't allow outside of range a-b
    #     * if closed (e.g. circle), allow to move back to beginning
    # 4)  if |(u* - u)C'(u)| < e1, halt

    # todo: take value from STEP fields
    # check if surface closed along U and V

    def closedV():
        return np.isclose([np.array(s[0]) - np.array(s[-1]) for s in S.ctrlpts2d], eps_bspline).all()

    def closedU():
        return np.isclose(np.array(S.ctrlpts2d[0]) - np.array(S.ctrlpts2d[-1]), eps_bspline).all()  # eps_STEP_AP21

    # S.delta = delta  # set evaluation delta
    # tradeoff between sample density and feature identification

    # deltaFactor = 5
    # deltaFactor = 4

    # S.delta_u = 1 / (S.ctrlpts_size_u * deltaFactor)
    # S.delta_v = 1 / (S.ctrlpts_size_v * deltaFactor)

    S.delta_u = 1 / (len(S.knotvector_u) * deltaFactor)
    S.delta_v = 1 / (len(S.knotvector_v) * deltaFactor)

    # if S.evalpts == None:  # evaluate surface points---------appears redundant
    #     S.evaluate()

    def neighbouringUV(U, V, isClosedU, isClosedV):

        for u in range(0, U - 2 + (2 * isClosedU)):
            for v in range(0, V - 2 + (2 * isClosedV)):

                if v == 0:
                    NW = (((u + 2) % U) * V)
                    N = (((u + 2) % U) * V) + 1
                    NE = (((u + 2) % U) * V) + 2

                    W = (((u + 1) % U) * V)
                    C = (((u + 1) % U) * V) + 1
                    E = (((u + 1) % U) * V) + 2

                    SW = ((u % U) * V)
                    S = ((u % U) * V) + 1
                    SE = ((u % U) * V) + 2

                if v >= 1:
                    NW = N
                    W = C  # central instance
                    SW = S

                    N = NE
                    C = E
                    S = SE

                    NE = (((u + 2) % U) * V) + ((v + 2) % V)
                    E = (((u + 1) % U) * V) + ((v + 2) % V)
                    SE = ((u % U) * V) + ((v + 2) % V)

                # print([NW, N, NE])
                # print([W, C, E])
                # print([SW, S, SE])
                # print("==========")
                # print([u, v, ((u * V) + v) ])
                # print("==========")

                # if v==V - 2:
                #     _1=1

                # if C==3:
                #     _1=1

                yield ([NW, N, NE, W, C, E, SW, S, SE], (u, v))

    # def _neighbouringUV(U, V, isClosedU, isClosedV):
    #     # generator function to yield neighbouring indices over a matrix of size U x V
    #
    #     openU = not isClosedU
    #     openV = not isClosedV
    #
    #     for u in range(openU, U - openU):
    #         for v in range(openV, V - openV):
    #
    #             if v == 0:  # closed in v direction
    #                 if u == U - 1:
    #                     NW = V - 1
    #                 else:
    #                     NW = ((u + 2) * V) + v - 1
    #
    #                 W = (V * (u + 1)) - 1
    #
    #                 if u == 0:
    #                     SW = (U * V) - 1  # V - 1
    #                 else:
    #                     SW = (u * V) - 1
    #
    #                 if u == U - 1:  # closed in u direction
    #                     N = v
    #                 else:
    #                     N = V * (u + 1)
    #
    #                 C = u * V
    #
    #                 if u == 0:
    #                     S = V * (U - 1)
    #                 else:
    #                     S = V * (u - 1)
    #
    #             if v == 1:
    #                 if u == U - 1:  # closed in u direction
    #                     NW = v - 1
    #                 else:
    #                     NW = ((u + 1) * V) + v - 1
    #
    #                 W = (u * V) + v - 1
    #
    #                 if u == 0:
    #                     SW = ((U - 1) * V) + v - 1
    #                 else:
    #                     SW = ((u - 1) * V) + v - 1
    #
    #                 if u == U - 1:  # closed in v direction
    #                     N = v
    #                 else:
    #                     N = ((u + 1) * V) + v
    #
    #                 C = (u * V) + v
    #
    #                 if u == 0:
    #                     S = ((U - 1) * V) + 1
    #                 else:
    #                     S = ((u - 1) * V) + v
    #
    #             if v > 1:
    #                 NW = N
    #                 W = C  # central instance
    #                 SW = S
    #
    #                 N = NE
    #                 C = E
    #                 S = SE
    #
    #             if v == V - 1:  # closed in v direction
    #                 if u < U - 1:
    #                     NE = (u + 1) * V
    #                 elif u == U - 1:  # closed in u direction
    #                     NE = 0
    #
    #                 if u < U - 1:
    #                     E = u * V
    #                 elif u == U - 1:
    #                     E = ((u - 1) * V) + v + 1
    #
    #                 if u == 0:
    #                     SE = ((U - 1) * V)
    #                 elif u < U - 1:
    #                     SE = (u - 1) * V
    #                 elif u == U - 1:
    #                     SE = ((u - 2) * V) + v + 1
    #
    #             elif v < V - 1:
    #                 if u < U - 1:  # only have to calculate this once per loop
    #                     NE = ((u + 1) * V) + v + 1
    #                 elif u == U - 1:  # closed in u direction
    #                     NE = v + 1
    #
    #                 if u < U - 1:  # only have to calculate this once per loop
    #                     E = (u * V) + v + 1
    #                 elif u == U - 1:
    #                     E = (u * V) + v + 1
    #
    #                 if u == 0:
    #                     SE = ((U - 1) * V) + v + 1
    #                 elif u < U - 1:  # only have to calculate this once per loop
    #                     SE = ((u - 1) * V) + v + 1
    #                 elif u == U - 1:
    #                     SE = ((u - 1) * V) + v + 1
    #
    #             # print([NW, N, NE])
    #             # print([W, C, E])
    #             # print([SW, S, SE])
    #             # print("==========")
    #
    #             yield ([NW, N, NE, W, C, E, SW, S, SE], (u, v))

    def localNeighbourSearch(S):
        extremaUV = []
        extremaP = []
        maxBound = []

        # not S.evalpts
        uSamples = np.linspace(S.knotvector_u[0],
                               S.knotvector_u[-1],
                               num=len(S.knotvector_u) * deltaFactor,
                               endpoint=~closedU())
        vSamples = np.linspace(S.knotvector_v[0],
                               S.knotvector_v[-1],
                               num=len(S.knotvector_v) * deltaFactor,
                               endpoint=~closedV())
        # uSamples = np.linspace(S.knotvector_u[0], S.knotvector_u[-1], num=len(S.knotvector_u) * deltaFactor)
        # vSamples = np.linspace(S.knotvector_v[0], S.knotvector_v[-1], num=len(S.knotvector_v) * deltaFactor)
        uvArray = np.zeros((len(uSamples) * len(vSamples), 2))

        for vi, vS in enumerate(vSamples):
            for ui, uS in enumerate(uSamples):
                uvArray[(ui * len(vSamples)) + vi] = [uS, vS]

        pArray = np.array(S.evaluate_list(uvArray))

        # test1 = S.evaluate_single(uvArray[0]) #(0.0, 0.0))
        # test2 = S.evaluate_single(uvArray[207]) #(0.0, 1.0))
        # minSample=[]

        for CC, uv in neighbouringUV(len(uSamples), len(vSamples), closedU(), closedV()):

            pNW = pArray[CC[0]]
            pN = pArray[CC[1]]
            pNE = pArray[CC[2]]
            pE = pArray[CC[3]]
            pC = pArray[CC[4]]
            pW = pArray[CC[5]]
            pSW = pArray[CC[6]]
            pS = pArray[CC[7]]
            pSE = pArray[CC[8]]

            dpuv_NW = np.linalg.norm(p - pNW)
            dpuv_N = np.linalg.norm(p - pN)
            dpuv_NE = np.linalg.norm(p - pNE)
            dpuv_E = np.linalg.norm(p - pE)
            dpuv = np.linalg.norm(p - pC)
            dpuv_W = np.linalg.norm(p - pW)
            dpuv_SW = np.linalg.norm(p - pSW)
            dpuv_S = np.linalg.norm(p - pS)
            dpuv_SE = np.linalg.norm(p - pSE)

            # testSample = [pNW, pN, pNE, pW, pC, pE, pSW, pS, pSE]
            # if ((dpuv_S <= dpuv_NW) and (dpuv_S <= dpuv_N) and (dpuv_S <= dpuv_NE) and
            #         (dpuv_S <= dpuv_W) and (dpuv_S <= dpuv_E) and
            #         (dpuv_S <= dpuv_SW) and (dpuv_S <= dpuv) and (dpuv_S <= dpuv_SE)):
            #     minSample.append(pS)

            if maxSearch:
                if ((dpuv >= dpuv_NW) and (dpuv >= dpuv_N) and (dpuv >= dpuv_NE) and
                        (dpuv >= dpuv_W) and (dpuv >= dpuv_E) and
                        (dpuv >= dpuv_SW) and (dpuv >= dpuv_S) and (dpuv >= dpuv_SE)):
                    # b_NW = np.linalg.norm(pC - pNW)
                    b_N = np.linalg.norm(pC - pN)
                    # b_NE = np.linalg.norm(pC - pNE)
                    b_E = np.linalg.norm(pC - pE)

                    b_W = np.linalg.norm(pC - pW)
                    # b_SW = np.linalg.norm(pC - pSW)
                    b_S = np.linalg.norm(pC - pS)
                    # b_SE = np.linalg.norm(pC - pSE)

                    maxBound.append(max([b_N, b_E, b_W, b_S]))
                    extremaUV.append(uvArray[CC[4]])  # uv[0]*S.sample_size_v + uv[1]
                    extremaP.append(pArray[CC[4]])

            else:  # minS
                if ((dpuv <= dpuv_NW) and (dpuv <= dpuv_N) and (dpuv <= dpuv_NE) and
                        (dpuv <= dpuv_W) and (dpuv <= dpuv_E) and
                        (dpuv <= dpuv_SW) and (dpuv <= dpuv_S) and (dpuv <= dpuv_SE)):
                    # where a point is orthogonal to a planar surface -> sphere radius test, or accept minima

                    # b_NW = np.linalg.norm(pC - pNW)
                    b_N = np.linalg.norm(pC - pN)
                    # b_NE = np.linalg.norm(pC - pNE)
                    b_E = np.linalg.norm(pC - pE)

                    b_W = np.linalg.norm(pC - pW)
                    # b_SW = np.linalg.norm(pC - pSW)
                    b_S = np.linalg.norm(pC - pS)
                    # b_SE = np.linalg.norm(pC - pSE)

                    maxBound.append(max([b_N, b_E, b_W, b_S]))
                    extremaUV.append(uvArray[CC[4]])
                    extremaP.append(pArray[CC[4]])

        return extremaUV, extremaP, maxBound

    def curvatureSearch():

        # not S.evalpts
        uSamples = np.linspace(S.knotvector_u[0],
                               S.knotvector_u[-1],
                               num=len(S.knotvector_u) * deltaFactor,
                               endpoint=~closedU())
        vSamples = np.linspace(S.knotvector_v[0],
                               S.knotvector_v[-1],
                               num=len(S.knotvector_v) * deltaFactor,
                               endpoint=~closedV())
        uvArray = np.zeros((len(uSamples) * len(vSamples), 2))

        for vi, vS in enumerate(vSamples):
            for ui, uS in enumerate(uSamples):
                uvArray[(ui * len(vSamples)) + vi] = [uS, vS]

        pArray = np.array(S.evaluate_list(uvArray))

        discreteK = [np.inf] * len(uvArray)
        for CC, uv in neighbouringUV(len(uSamples), len(vSamples), closedU(), closedV()):

            puv = np.array(pArray[CC[4]])
            pur, pucc = radiusCentre3points(pArray[CC[3]], puv, pArray[CC[5]])
            pvr, pvcc = radiusCentre3points(pArray[CC[7]], puv, pArray[CC[1]])

            # find media point of U curvature centre and V curvature, kuv
            if pur is np.inf and pvr is not np.inf:
                kuv = pvcc - puv
            elif pur is not np.inf and pvr is np.inf:
                kuv = pucc - puv
            elif pur is np.inf and pvr is np.inf:
                discreteK[(S.sample_size_u * uv[0]) + uv[1]] = np.inf
                break
            else:
                kuv = ((pucc - puv) + (pvcc - puv)) / 2

            dkuv = np.linalg.norm(kuv - puv)

            if np.linalg.norm(p - puv) > np.linalg.norm(kuv - p):
                # p is on the same side of the surface as the median curvature
                discreteK[(S.sample_size_u * uv[0]) + uv[1]] = dkuv

            if np.linalg.norm(p - puv) < np.linalg.norm(kuv - p):
                # smallest radius, with curvature centre on far side of surface from p
                discreteK[(S.sample_size_u * uv[0]) + uv[1]] = -dkuv

            # todo: examine issue with non-colinear curve centre and projecting point
            # e.g. proj(p)
            #                 if np.dot(p - p1, pucc - p1)/np.linalg.norm(pucc - p1) > 0:
            #                     discreteK[i] = pur
            #                 else:
            #                     discreteK[i] = -pur

        return discreteK

    def localCurvatureExtrema(discreteK):
        leu = []
        # localExtremaP = []
        # maxBound = []

        for CC, uv in neighbouringUV(len(S.knotvector_u) * deltaFactor,
                                     len(S.knotvector_v) * deltaFactor,
                                     closedU(), closedV()):
            surroundNodes = [CC[1], CC[3], CC[5], CC[7]]

            if not any(discreteK[ssn] == np.inf for ssn in surroundNodes):
                if maxSearch:
                    if (all((np.abs(discreteK[ssn]) >= discreteK[(S.sample_size_u * uv[0]) + uv[1]]) for ssn in
                            surroundNodes) and
                            (discreteK[(S.sample_size_u * uv[0]) + uv[1]] > 0)):
                        # this version excludes negative curves from surrounding curvature values
                        # if all(np.abs(discreteK[ssn]) > discreteK[S.sample_size_u * u + v] or
                        # (discreteK[ssn]<0) for ssn in surroundNodes) and
                        # (discreteK[S.sample_size_u * u + v] > 0):

                        pC = np.array(S.evaluate_single(uv))
                        localExtremaP(pC)
                        maxBound.append(
                            max([np.linalg.norm(pC - np.array(S.evaluate_single(uvk))) for uvk in surroundNodes]))
                        leu.append(uv)
                else:
                    if all((np.abs(discreteK[ssn]) >= discreteK[(S.sample_size_u * uv[0]) + uv[1]]) for ssn in
                           surroundNodes) and (discreteK[(S.sample_size_u * uv[0]) + uv[1]] < 0):
                        pC = np.array(S.evaluate_single(uv))
                        localExtremaP(pC)
                        maxBound.append(
                            max([np.linalg.norm(pC - np.array(S.evaluate_single(uvk))) for uvk in surroundNodes]))
                        leu.append(uv)

        return leu

    def subSearchUV(mpU, mpV, tol):
        # create simple subdivisions around a minimum point as simple hillclimb search
        tolFlag = 0
        divU = S.delta_u
        divV = S.delta_v
        dpuvlocal = np.linalg.norm(p - np.array(S.evaluate_single((mpU, mpV,))))
        # dpuvlocal = np.inner(p - np.array(S.evaluate_single((mpU, mpV,))))

        while not tolFlag:
            divU /= 2
            divV /= 2

            subDivs = [[mpU - divU, mpV + divV],
                       [mpU, mpV + divV],
                       [mpU + divU, mpV + divV],
                       [mpU - divU, mpV],
                       [mpU, mpV],
                       [mpU + divU, mpV],
                       [mpU - divU, mpV - divV],
                       [mpU, mpV - divV],
                       [mpU + divU, mpV - divV]]

            for sdv in subDivs:
                if sdv[0] > 1.0: sdv[0] = 1.0
                if sdv[0] < 0.0: sdv[0] = 0.0
                if sdv[1] > 1.0: sdv[1] = 1.0
                if sdv[1] < 0.0: sdv[1] = 0.0

            pSubDivs = S.evaluate_list(subDivs)

            # test np.inner
            # dispUV = [np.inner(p - psd) for psd in pSubDivs]

            dpuv_NW = np.linalg.norm(p - pSubDivs[0])
            dpuv_N = np.linalg.norm(p - pSubDivs[1])
            dpuv_NE = np.linalg.norm(p - pSubDivs[2])

            dpuv_W = np.linalg.norm(p - pSubDivs[3])
            # dpuv = np.linalg.norm(p - pSubDivs[4])
            dpuv_E = np.linalg.norm(p - pSubDivs[5])

            dpuv_SW = np.linalg.norm(p - pSubDivs[6])
            dpuv_S = np.linalg.norm(p - pSubDivs[7])
            dpuv_SE = np.linalg.norm(p - pSubDivs[8])

            dispUV = [dpuv_NW,
                      dpuv_N,
                      dpuv_NE,
                      dpuv_W,
                      # dpuvlocal,
                      dpuv_E,
                      dpuv_SW,
                      dpuv_S,
                      dpuv_SE]

            if maxSearch:
                if np.abs(max(dispUV) - dpuvlocal) >= tol:
                    # np.abs(max(dispUV) - dpuvlocal) >= tol:
                    maxUVindex = dispUV.index(max(dispUV))
                    mpU = subDivs[maxUVindex][0]
                    mpV = subDivs[maxUVindex][1]
                    dpuvlocal = max(dispUV)
                else:
                    tolFlag += 1
            else:
                if np.abs(min(dispUV) - dpuvlocal) >= tol:
                    minUVindex = dispUV.index(min(dispUV))
                    mpU = subDivs[minUVindex][0]
                    mpV = subDivs[minUVindex][1]
                    dpuvlocal = min(dispUV)
                else:
                    tolFlag += 1
        return (mpU, mpV,)

    def NewtonRaphson(cuv, maxits=5):

        # def f(uv):
        #     return surface.derivatives(uv[0], uv[1], 2)

        # e[0][0], the surface point itself
        # e[0][1], the 1st derivative w.r.t. v
        # e[2][1], the 2nd derivative w.r.t. u and 1st derivative w.r.t. v

        def n(uv, e, r):  # np.array inputs
            #   f = Su(u,v) * r = 0
            #   g = Sv(u,v) * r = 0

            Su = e[1][0]
            Sv = e[0][1]

            Suu = e[2][0]
            Svv = e[0][2]

            Suv = e[1][1]
            Svu = e[1][1]

            f = np.dot(Su, r)
            g = np.dot(Sv, r)
            k = [-f, -g]

            J00 = np.dot(Su, Su) + np.dot(Suu, r)
            J01 = np.dot(Su, Sv) + np.dot(Suv, r)
            J10 = np.dot(Su, Sv) + np.dot(Svu, r)
            J11 = np.dot(Sv, Sv) + np.dot(Svv, r)

            # d =   [ u* - u, v* - v ]
            # k = - [ f(u,v), g(u,v) ]
            # J =   |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
            #        Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r

            J = [[J00, J01], [J10, J11]]
            d = np.linalg.solve(J, k)
            return d + uv

        i = 0
        # e = None
        while i < maxits:
            if (cuv[0] < 0) or (cuv[0] > 1) or (cuv[1] < 0) or (cuv[1] > 1):
                return (np.inf, np.inf)  # derivatives sometimes > |2|

            # test = S.derivatives(1.0, 0.0, 0)

            e = np.array(S.derivatives(cuv[0], cuv[1], 2))
            dif = e[0][0] - p  # e[0][0] is surface point evaluated at cuv(u,v)
            c1v = np.linalg.norm(dif)
            c1 = c1v < eps1  # |S(u,v) - p| < e1 (point coincidence)

            #  |Su(u,v)*(S(u,v) - P)|
            #  ----------------------  < e2 (cosine minima)
            #  |Su(u,v)| |S(u,v) - P|

            c2an = np.dot(e[1][0], dif)
            c2ad = np.linalg.norm(e[1][0]) * c1v
            c2av = c2an / c2ad
            c2a = c2av < eps2

            #  |Sv(u,v)*(S(u,v) - P)|
            #  ----------------------  < e2 (cosine minima)
            #  |Sv(u,v)| |S(u,v) - P|

            c2bn = np.dot(e[0][1], dif)
            c2bd = np.linalg.norm(e[0][1]) * c1v
            c2bv = c2bn / c2bd
            c2b = c2bv < eps2

            # exit if all tolerance are met,
            if (c1 and c2a and c2b):
                return cuv

            # otherwise, take a step
            ct = n(cuv, e, dif)

            #  correct for exceeding bounds
            if ct[0] < S.knotvector_u[0]:  # [ maxu - ( ct[0] - minu ), ct[1] ]
                if closedU():
                    ct = [S.knotvector_u[-1] - (S.knotvector_u[0] - ct[0]), ct[1]]
                # if closedU(): ct = [S.knotvector_u[0] + ct[0] % (S.knotvector_u[-1] - S.knotvector_u[0]), ct[1]]
                else:
                    ct = [S.knotvector_u[0] + eps_bspline, ct[1]]

            elif ct[0] > S.knotvector_u[-1]:
                if closedU():
                    ct = [S.knotvector_u[0] + (ct[0] - S.knotvector_u[-1]), ct[1]]
                # if closedU(): ct = [S.knotvector_u[-1] - ct[0] % (S.knotvector_u[-1] - S.knotvector_u[0]), ct[1]]
                else:
                    ct = [S.knotvector_u[-1] - eps_bspline, ct[1]]

            if ct[1] < S.knotvector_v[0]:
                if closedV():
                    ct = [ct[0], S.knotvector_v[-1] - (S.knotvector_v[0] - ct[1])]
                # if closedV(): ct = [ct[0], S.knotvector_v[0] + ct[1] % (S.knotvector_v[-1] - S.knotvector_v[0])]
                else:
                    ct = [ct[0], S.knotvector_v[0] + eps_bspline]

            elif ct[1] > S.knotvector_v[-1]:
                # if closedV(): ct = [ct[0], S.knotvector_v[0] + (ct[1] - S.knotvector_v[-1])]
                if closedV():
                    ct = [ct[0], S.knotvector_v[-1] - ct[1] % (S.knotvector_v[-1] - S.knotvector_v[0])]
                else:
                    ct = [ct[0], S.knotvector_v[-1] - eps_bspline]

            c3v0 = np.linalg.norm((ct[0] - cuv[0]) * e[1][0])
            c3v1 = np.linalg.norm((ct[1] - cuv[1]) * e[0][1])

            # if |(u* - u) C'(u)| < e1, halt
            if (c3v0 + c3v1 < eps1):
                return cuv
            cuv = ct
            i += 1

        if i == maxits:  # Newton-Raphson fails
            return (np.inf, np.inf)
        return cuv

    if curvatureTest:  # discrete grid summing curvatures along U and V axes
        localExtremaUV, localExtremaP, maxBound = localCurvatureExtrema(curvatureSearch())
    else:
        localExtremaUV, localExtremaP, maxBound = localNeighbourSearch(S)

    if len(localExtremaUV) == 0: return []

    extremaUV = []
    for luv, lp, mb in zip(localExtremaUV, localExtremaP, maxBound):
        mUV = NewtonRaphson(luv)

        if not np.inf in mUV:
            # find a maximum deviation to indicate when N-R fails
            if np.linalg.norm(lp - np.array(S.evaluate_single(mUV))) <= mb:
                extremaUV.append(mUV)
            else:
                mUV = (np.inf, np.inf)
        if np.inf in mUV:
            # Newton-Raphson fail, try hillclimb
            mUV = subSearchUV(luv[0], luv[1], eps_bspline)
            extremaUV.append(mUV)
        # print(S.evalpts[(S.sample_size_u * luv[0]) + luv[1]])

    # filter out identical values
    extremaUnique = []
    for m in extremaUV:
        # displacement between any 2 points < eps1
        # compare u, v individually to reduce comparison pool
        if not any([np.isclose(m, u, eps1).all() for u in extremaUnique]):
            extremaUnique.append(m)

    if (localExtrema and (
            len(extremaUnique) == 1)) or not localExtrema:  # if there is only one extrema, localextrema is moot
        if uv_xyz:
            return extremaUnique  # return single value in list for compatibility
        else:
            return [np.array(S.evaluate_single(extremaUnique[0])), ]
    else:
        if localExtrema:
            if uv_xyz:
                return extremaUnique
            else:
                return S.evaluate_list(extremaUnique)
        else:  # return single maxima
            extremaUniquePoint = S.evaluate_list(extremaUnique)
            dispsExtrema = [np.linalg.norm(np.array(e) - p) for e in extremaUniquePoint]
            if maxSearch:  # return single value in list for compatibility
                if uv_xyz:
                    return [extremaUnique[dispsExtrema.index(max(dispsExtrema))], ]
                else:
                    return [np.array(extremaUniquePoint[dispsExtrema.index(max(dispsExtrema))]), ]
            else:  # minsearch
                if uv_xyz:
                    return [extremaUnique[dispsExtrema.index(min(dispsExtrema))], ]
                else:
                    return [np.array(extremaUniquePoint[dispsExtrema.index(min(dispsExtrema))]), ]


def BsplineCurveExtremaDisp(C, p,
                            maxSearch=True,
                            localExtrema=False,
                            curvatureTest=False,
                            uv_xyz=True,
                            eps1=0.0001,
                            eps2=0.0005,
                            deltaFactor=4,
                            eps_bspline=1E-10):
    """
    Return parametric u,v parameters or cartesian point on a curve closest to a point based on Nurbs Book algorithm,
    Piegl & Tiller p.230. Same idea to find orthogonal tangent, but at maximum/minimum displacement from a centroid point
    Revised to test minima segment candidates using 3-point defined circle parameters, curvatureTest
    If Newton-Raphson fails, resorts to a subdivision bounding box search.

    :param maxSearch: return point furthest from point if true, nearest if false
    :param localExtrema: True for returning local maxima/minima, false for global maxima/minima point
    :param curvatureTest: True: seed N-R with a point taken from a minima curvature
                          False: use displacement comparison over sampled points
    :param uv_xyz: True: return u,v, parameter, False: return x, y, z point value
    :param eps1: Newton-Raphson halting condition
    :param eps2: Newton-Raphson halting condition
    :param deltaFactor: sampling factor for initial search
    :param eps_bspline: minimal point disp for sub-sample algorithm termination, see subSearchU()
    :return: list of Cartesian or parametric values
    """

    # def radiusCentre3points(p1, p2, p3):
    #     # radius + centre from 3 point circle
    #
    #     #collinear test
    #     coords = np.array([p1, p2, p3])
    #     coords -= coords[0]  # offset for collinear points to intersect the origin
    #     if np.linalg.matrix_rank(coords, tol=None) == 1:
    #         return np.inf, []
    #
    #     t = p2 - p1
    #     u = p3 - p1
    #     v = p3 - p2
    #
    #     w = np.cross(t, u)  # triangle normal
    #     wsl = np.linalg.norm(w)
    #     if (wsl < 10e-14):
    #         return False  # triangle area too small
    #
    #     wsl = np.dot(w, w)
    #     iwsl2 = 1. / (2. * wsl)
    #     tt = np.dot(t, t)
    #     uu = np.dot(u, u)
    #
    #     circCenter = p1 + (u * tt * (np.dot(u, v)) - t * uu * (np.dot(t, v))) * iwsl2
    #     circRadius = np.sqrt(tt * uu * (np.dot(v, v)) * iwsl2 * 0.5)
    #     # circAxis   = w / np.sqrt(wsl)
    #
    #     return circRadius, circCenter

    # eps1 = 0.0001   # Euclidean distance measure
    # eps2 = 0.0005   # zero cosine measure
    # eps_bspline = 1E-10

    C.delta = 1 / (len(C.knotvector) * deltaFactor)

    def closedC():
        return np.isclose(np.linalg.norm(np.array(C.ctrlpts[0]) - np.array(C.ctrlpts[-1])), eps_bspline)

    def localDispTest():
        # localDisp = [0] * len(pts)
        localExtremaU = []

        # get discrete curvature values and point displacement
        for i in range(1 - closedC(), len(pts) - 1 + closedC()):
            if (i == 0):
                p0 = np.array(pts[-2])
                p1 = np.array(pts[0])
                p2 = np.array(pts[1])

            elif (i == len(pts) - 1):
                p0 = np.array(pts[-2])
                p1 = np.array(pts[0])
                p2 = np.array(pts[1])

            else:
                p0 = np.array(pts[i - 1])
                p1 = np.array(pts[i])
                p2 = np.array(pts[i + 1])

            if i == 0 or i == 1 or i == len(pts) - 1:
                dp0 = np.linalg.norm(p0 - p)
                dp1 = np.linalg.norm(p1 - p)
                dp2 = np.linalg.norm(p2 - p)

            if (i > 1) and (i < len(pts) - 1):
                dp0 = dp1
                dp1 = dp2
                dp2 = np.linalg.norm(p2 - p)

            if maxSearch:
                if (dp1 >= dp2) and (dp1 >= dp0):
                    localExtremaU.append(kvs[i])
            else:  # minS
                if (dp1 <= dp2) and (dp1 <= dp0):
                    # where a point is orthogonal to a planar surface -> sphere radius test, or accept minima
                    localExtremaU.append(kvs[i])

            # localDisp[i] = dp1

        return localExtremaU

    def localCurvatureTest():
        discreteK = [np.inf] * len(pts)
        localExtremaU = []

        # get discrete curvature values and point displacement
        for i in range(1 - closedC(), len(pts) - 1 + closedC()):
            if (i == 0):
                p0 = np.array(pts[-2])
                p1 = np.array(pts[0])
                p2 = np.array(pts[1])

            elif (i == len(pts) - 1):
                p0 = np.array(pts[-2])
                p1 = np.array(pts[-1])
                p2 = np.array(pts[1])

            else:
                p0 = np.array(pts[i - 1])
                p1 = np.array(pts[i])
                p2 = np.array(pts[i + 1])

            pur, pucc = radiusCentre3points(p0, p1, p2)
            if not pur:
                break
            if (pur == np.inf):
                break
            else:
                # p is on the same side of the surface as the mean curvature
                if np.dot(p - p1, pucc - p1) / np.linalg.norm(pucc - p1) > 0:
                    discreteK[i] = pur
                else:
                    discreteK[i] = -pur

        # get local max curvature values
        for i in range(1 - closedC(), len(discreteK) - 1 + closedC()):
            if (i == 0):
                k0 = discreteK[-2]
                k1 = discreteK[0]
                k2 = discreteK[1]

            elif (i == len(discreteK) - 1):
                k0 = discreteK[-2]
                k1 = discreteK[-1]
                k2 = discreteK[0]

            else:
                k0 = discreteK[i - 1]
                k1 = discreteK[i]
                k2 = discreteK[i + 1]

            if (k0 != np.inf) and (k1 != np.inf) and (k2 != np.inf):
                if (np.abs(k1) <= np.abs(k2)) and (np.abs(k1) <= np.abs(k0)):
                    if (k1 > 0) and maxSearch:
                        localExtremaU.append(kvs[i])
                    elif (k1 < 0) and not maxSearch:
                        # where a point is orthogonal to a planar surface -> sphere radius test, or accept minima
                        localExtremaU.append(kvs[i])

        return localExtremaU

    def subSearchU(C, mpU, tol):
        # create simple subdivisions around a minimum point on curve as simple hillclimb search
        tolFlag = 0
        divU = C.delta
        dpulocal = np.linalg.norm(p - C.evaluate_single(mpU))

        while not tolFlag:
            divU /= 2
            subDivs = [mpU - divU, mpU, mpU + divU]
            pSubDivs = C.evaluate_list(subDivs)

            dpu_L = np.linalg.norm(p - pSubDivs[0])
            dpu_R = np.linalg.norm(p - pSubDivs[1])

            dispU = [dpu_L, dpulocal, dpu_R]

            if maxSearch:
                if np.abs(max(dispU) - dpulocal) >= tol:
                    maxUindex = dispU.index(max(dispU))
                    mpU = subDivs[maxUindex]
                    dpulocal = max(dispU)
                else:
                    tolFlag += 1
            else:
                if np.abs(min(dispU) - dpulocal) >= tol:
                    minUindex = dispU.index(min(dispU))
                    mpU = subDivs[minUindex]
                    dpulocal = min(dispU)
                else:
                    tolFlag += 1
        return mpU

    def NewtonRaphson(cu, maxits=5):

        def n(u2, e1, d):
            #   Newton's method: 	 u* = u - f / f'
            #   use product rule to form derivative, f':   f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
            #   d:  ( C(u) - p )
            return u2 - (np.dot(e1[1], d) / (np.dot(e1[2], d) + np.dot(e1[1], e1[1])))

        i = 0
        while (i < maxits):
            # Newton iteration to find orthogonal tangent and is max/min agnostic
            e = np.array(C.derivatives(cu, order=2))  # f(cu)
            dif = e[0] - p  # C(u) - p
            c1v = np.linalg.norm(dif)  # |C(u) - p|
            c2n = np.dot(e[1], dif)  # C'(u) * (C(u) - P)
            c2d = np.linalg.norm(e[1]) * c1v  # |C'(u)||C(u) - P|
            c2v = c2n / c2d  # |C'(u) * (C(u) - P)| / |C'(u)||C(u) - P|
            # c2v = np.dot(e[1], dif) / (np.linalg.norm(e[1]) * np.linalg.norm(dif))

            if (c1v < eps1) and (np.abs(c2v) < eps2):
                return cu
            # ct = n(cu, e, dif)                           #   u* = u - f / f'
            ct = cu - (np.dot(e[1], dif) / (np.dot(e[2], dif) + np.dot(e[1], e[1])))

            #  correct for exceeding bounds
            if ct < C.knotvector[0]:  # [ maxu - ( ct[0] - minu ), ct[1] ]
                if closedC():
                    ct = C.knotvector[0] + ct % (C.knotvector[-1] - C.knotvector[0])
                # if closedC(): ct = C.knotvector[-1] - (C.knotvector[0] - ct) # ct = C.knotvector[-1] - (ct - C.knotvector[0]) # NURBS book
                else:
                    ct = C.knotvector[0] + eps_bspline

            elif ct > C.knotvector[-1]:
                # if closedC(): ct = C.knotvector[0] + (ct - C.knotvector[-1])
                if closedC():
                    ct = C.knotvector[-1] - ct % (C.knotvector[-1] - C.knotvector[0])
                else:
                    ct = C.knotvector[-1] - eps_bspline

            c3v = np.linalg.norm(np.multiply(ct - cu, e[1]))
            if c3v < eps1:
                return cu

            cu = ct
            i += 1

        if i == maxits:  # Newton-Raphson fails
            return np.inf
        return cu

    pts, kvs = splineToPolyline(C, curveSampleFactor=deltaFactor)
    # original version fails on points of equal minimal curvature

    # 2-part selection for global maxima using both displacement and curve inflection
    # curve inflection determined from local curvature minima - prob doesn't work for self intersecting curves
    # get set of local minima curvature - use for curve characterisation
    # from set, find minima/maxima disp from point-of-interest to curve point tangent

    # curvature approach fails without C2 continuity at minima/maxima?
    # not much of an issue with STEP spline limitations

    if curvatureTest:  # discrete grid summing curvatures along U and V axes
        localExtremaU = localCurvatureTest()
    else:  # localDispTest()
        localExtremaU = localDispTest()

    if len(localExtremaU) == 0: return []

    extremaU = []
    for lu in localExtremaU:
        mU = NewtonRaphson(lu)
        if mU != np.inf:
            extremaU.append(mU)
        else:
            # Newton-Raphson fail, try hillclimb
            mU = subSearchU(C, lu, eps_bspline)
            extremaU.append(mU)

    # filter out identical values
    extremaUnique = []
    for m in extremaU:
        # displacement between any 2 points < eps1
        # compare u, v individually to reduce comparison pool
        if not any(np.isclose(m, u, eps1).all() for u in extremaUnique):
            extremaUnique.append(m)

    # if there is only one extrema, localextrema is moot
    if (localExtrema and (len(extremaUnique) == 1)) or not localExtrema:
        if uv_xyz:
            return extremaUnique  # return single value in list for compatibility
        else:
            return [np.array(C.evaluate_single(extremaUnique[0])), ]
    else:
        if localExtrema:
            if uv_xyz:
                return extremaUnique
            else:
                return [np.array(cp) for cp in C.evaluate_list(extremaUnique)]
        else:  # return single maxima
            extremaUniquePoint = C.evaluate_list(extremaUnique)
            dispsExtrema = [np.linalg.norm(np.array(e) - p) for e in extremaUniquePoint]
            if maxSearch:  # return single value in list for compatibility
                if uv_xyz:
                    return [extremaUnique[dispsExtrema.index(max(dispsExtrema))], ]
                else:
                    return [np.array(extremaUniquePoint[dispsExtrema.index(max(dispsExtrema))]), ]
            else:  # minsearch
                if uv_xyz:
                    return [extremaUnique[dispsExtrema.index(min(dispsExtrema))], ]
                else:
                    return [np.array(extremaUniquePoint[dispsExtrema.index(min(dispsExtrema))]), ]


def pointInArc(p, v1, v2, _refDir, _auxDir, _normDir, aCentreP, rotSym=True):
    """
    Return true if point p is within arc defined by vertex v1, v2 with centre point aCentreP
    :param p: cartesian point
    :param v1: STEP end arc vertex
    :param v2: STEP end arc vertex
    :param _refDir: reference direction vector
    :param _auxDir: auxilliary direction vector
    :param _normDir: normal direction vector
    :param aCentreP: arc centre point
    :param rotSym: bool, test whether radius is consistent for point
    :return bool
    """
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
        # print("pointInArc() transform error: " + str(max(zTest) - min(zTest)))
        # raise RuntimeError("error")

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


def pointOrderInArc(p, v1, v2, _refDir, _auxDir, _normDir, aCentreP, rotSym=True):
    """
    Return order of a list of cartesian points around a supplied axis
    within arc defined by vertex v1, v2 with centre point aCentreP
    :param p: cartesian point
    :param v1: STEP end arc vertex
    :param v2: STEP end arc vertex
    :param _refDir: reference direction vector
    :param _auxDir: auxilliary direction vector
    :param _normDir: normal direction vector
    :param aCentreP: arc centre point
    :param rotSym: bool, test whether radius is consistent for point
    :return sorted index list
    """

    # (note: if VOR, reverse v1, v2 entries)
    # determine the order of a list of cartesian points around a supplied axis
    # convert aCentreP, p, v1 & v2 to local coordinates with arcAxisDir as z
    # no requirement for equal radius (e.g. ellipse)

    #print("pointOrderInArc() undertested")

    if not isinstance(p, list):
        p = [p, ]
        print("pointOrderInArc() input point(s) in list")

    rotM = np.array([_refDir, _auxDir, _normDir])

    v1UV = np.matmul(rotM, v1 - aCentreP)
    v1UV_angle = np.arctan2(v1UV[1], v1UV[0])

    v2UV = np.matmul(rotM, v2 - aCentreP)
    v2UV_angle = np.arctan2(v2UV[1], v2UV[0])

    if np.abs(v1UV_angle) < eps_STEP_AP21:
        v1UV_angle = 0.0
    if (v1UV_angle - v2UV_angle) < eps_STEP_AP21:
        v2UV_angle = v1UV_angle + (np.pi * 2)

    # pUV_angle = []
    # if isinstance(p, list):
    #     for pp in p:
    #         ppUV = np.matmul(rotM, pp - aCentreP)
    #         ppUV_angle = np.arctan2(ppUV[1], ppUV[0])
    #         pUV_angle.append(ppUV_angle)
    # elif isinstance(p, np.ndarray):
    #     pUV = np.matmul(rotM, p - aCentreP)
    #     pUV_angle.append(np.arctan2(pUV[1], pUV[0]))

    pUV_angle = []
    for pp in p:
        # if pp is None:
        #     _1=1
        ppUV = np.matmul(rotM, pp - aCentreP)
        ppUV_angle = np.arctan2(ppUV[1], ppUV[0])
        ppUV_angle -= v1UV_angle
        if ppUV_angle < 0: ppUV_angle += (2 * np.pi)
        pUV_angle.append(ppUV_angle % (2 * np.pi))

    # # sanity check that Z-values are identical
    # if isinstance(p, list):
    #     zTest = [pp[2] for pp in pUV]
    #     zTest = zTest + [v1UV[2], v2UV[2]]
    # else:
    #     zTest = [pUV[2], v1UV[2], v2UV[2]]
    #
    # # if not np.isclose(max(zTest), min(zTest)):
    # planeTestFail = 0
    # if max(zTest) - min(zTest) > eps_STEP_AP21:
    #     planeTestFail += 1
    #     #print("pointInArc() transform error: " + str(max(zTest) - min(zTest)))
    #     #raise RuntimeError("error")

    # check if vertex order reversed before calling function

    # if (v1UV_angle < v2UV_angle): ???
    withinVertex1andVertex2 = [((pUV_a >= v1UV_angle) and (pUV_a <= v2UV_angle)) for pUV_a in pUV_angle]

    return withinVertex1andVertex2, np.argsort(pUV_angle), [pUV_a / (2 * np.pi) for pUV_a in pUV_angle]


def pointProjectAxis(p, axP, axN):  # todo is this also projection to surface via normal?
    '''
    Project point p to axis defined by point axP and normal axN
    :param p: point to project
    :param axP: axis point
    :param axN: axis normal vector
    :return cartesian point
    '''
    vpp = p - axP
    vpn = axP - axN
    return vpn * np.dot(vpp, vpn) / np.dot(vpn, vpn)


def pointCircleMinMaxDisp(p, centrePoint, normAxis, radius, interior=False):
    '''
    Get minima/maxima point on circle relative to point p
    :param p: point
    :param centrePoint: circle centre point
    :param normAxis: circle normal axis
    :param radius: circle radius
    :param interior: bool, if true and point is within interior, interior poinnt is returned,
                    otherwise nearest point on periphery
    :return minima, maxima point
    '''
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

    # centreNormal = False
    if np.isclose(ppCentreDisp, 0):
        # p lies on normal through circle centrePoint,
        return centrePoint, None

    pUVdir = (pUV - centrePoint) / ppCentreDisp

    nUVmin = centrePoint + radius * pUVdir
    nUVmax = centrePoint - radius * pUVdir

    # nUVmin = centrePoint + radius / (ppCentreDisp * (pUV - centrePoint))
    # nUVmax = centrePoint - radius / (ppCentreDisp * (pUV - centrePoint))

    if ppCentreDisp >= radius:  # and not interior:
        return nUVmin, nUVmax
    else:
        # if ppCentreDisp < eps_STEP_AP21:
        #     return centrePoint, centrePoint
        # else:
        centreOffset = pUV - centrePoint
        centreOffset = centreOffset / np.linalg.norm(centreOffset)
        nUVmax = centrePoint - radius * centreOffset
        #
        if interior:
            nUVmin = pUV
        else:
            nUVmin = centrePoint + radius * centreOffset

    return nUVmin, nUVmax,


def pointEllipseMinMaxDisp(p,
                           eCentre,
                           eLocalXaxis,
                           eLocalYaxis,
                           eNormalAxis,
                           eMajorRadius,
                           eMinorRadius,
                           interior=False):
    '''
    Get minima/maxima point on ellipse relative to point p
    Robert Nurnberg method, see http://wwwf.imperial.ac.uk/~rn/distance2ellipse.pdf
    Move origin to ellipse centre, and rotate point coordinate to local coordinates
    based on ellipse major & minor axis.
    Assuming direction_ratios of AXIS2_PLACEMENT_3D are normalised

    :param p: point
    :param eCentre: ellipse centre point
    :param eLocalXaxis: ellipse major axis
    :param eLocalYaxis: ellipse minor axis
    :param eNormalAxis: ellipse normal axis
    :param eMajorRadius: ellipse major radius
    :param eMinorRadius: ellipse minor radius
    :param interior: bool, if true and point is within interior, interior poinnt is returned,
                    otherwise nearest point on periphery
    :return minima, maxima point
    '''

    # eFeatureAxis = arrayTypeCast(eFeatureAxis)
    # eLocalXYaxis = arrayTypeCast(eLocalXYaxis)
    # eCentre = arrayTypeCast(eCentre)

    # cLocalYaxis = np.cross(eLocalXaxis - eCentre, eNormalAxis - eCentre) # STEP AP21 gives axes wrt local origin
    # eLocalYaxis = np.cross(eLocalXaxis, eNormalAxis)
    # eLocalYaxis = eLocalYaxis / np.linalg.norm(eLocalYaxis)  # unit normal vector

    rotM = np.array([eLocalXaxis, eLocalYaxis, eNormalAxis])
    pUV = np.matmul(rotM, p - eCentre)

    if np.linalg.norm(p - eCentre) < eps:
        print("should be 2x maxima/minima, write some code")  # todo: ellipse dual min/max
        # maxima will be +/- eMajorRadius points along major axis
        # minima +/- eMinorRadius points on minor axis

    # project to ellipse plane
    pUV = np.array([pUV[0], pUV[1], 0])

    theta = np.arctan2(eMajorRadius * pUV[1], eMinorRadius * pUV[0])
    i = 0
    n0 = pUV
    radiusConst = eMajorRadius ** 2 - eMinorRadius ** 2
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
            if (np.linalg.norm(n) > pUVdisp) and not interior:  # test if p < n, inside ellipse
                nUVmin = np.matmul(rotM.T, pUV) + eCentre
            else:
                nUVmin = np.matmul(rotM.T, n) + eCentre

            # furthest point on ellipse must be edge point on line that contains ellipse centre and closest point
            nUVmax = np.matmul(rotM.T, np.array([-n[0], -n[1], 0])) + eCentre

        n0 = n

    # if interior:
    #     nUVmin = pUV
    # else:
    #     nUVmin = centrePoint + radius * centreOffset

    return nUVmin, nUVmax


def intersectSegmentPlane(vertex0,
                          vertex1,
                          planeNorm,
                          planePoint,
                          precision=eps_STEP_AP21):
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
    return None  # segment is parallel to plane


def intersectPlaneCircle(pPoint, pNorm, cPoint, cNorm, cRadius, tol=eps_STEP_AP21):
    '''
    Get intersection of a plane and a circle
    4 cases of plane-circle intersection:

    1. intersect in 2 points (secant),
    2. intersect in 1 point (tangent),
    3. do not intersect, or
    4. coincide (circle.plane == plane).

    :param pPoint: plane point
    :param pNorm: plane norm
    :param cPoint: circle point
    :param cNorm: circle norm
    :param cRadius: circle radius
    :param tol:
    :return segment vertexes
    '''

    if np.abs(np.dot(pNorm, cNorm)) >= 1 - tol:
        return None

    d = np.cross(pNorm, cNorm)  # direction of intersection line
    # vector in plane 1 perpendicular to the direction of the intersection line
    v1 = np.cross(d, pNorm)
    p1 = pPoint + v1  # point on plane 1
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
            + cPoint[1] ** 2
            + cPoint[2] ** 2
            + l1[0] ** 2
            + l1[1] ** 2
            + l1[2] ** 2
            - 2.0 * (cPoint[0] * l1[0] + cPoint[1] * l1[1] + cPoint[2] * l1[2])
            - cRadius ** 2
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


def intersectArcPlane(planeNormal,
                      planePoint,
                      arcCentrePoint,
                      arcRefDir,
                      arcAuxDir,
                      arcNormal,
                      arcRadius,
                      arcVertex1,
                      arcVertex2):
    '''
    Get intersection segnemt between arc and plane,

    :param planeNormal:
    :param planePoint:
    :param arcCentrePoint:
    :param arcRefDir:
    :param arcAuxDir:
    :param arcNormal:
    :param arcRadius:
    :param arcVertex1:
    :param arcVertex2:
    :return intersection segment vertexes
    '''

    # planePoint = centroidProjAxisPoint for most calls

    # points within the edge
    validatedPoints = []

    # determine if circle and plane intersect, if so return intersection line/segment (cosine projection?)
    circleIntersectPoints = intersectPlaneCircle(planePoint, planeNormal, arcCentrePoint, arcNormal, arcRadius)
    if circleIntersectPoints is None:
        return []
    if len(circleIntersectPoints) == 0:
        return []  # no intersection
    if len(circleIntersectPoints) == 1:
        print("tangent exception")
        pass
    # check projected segment is between arc vertices,
    for cip in circleIntersectPoints:

        if pointInArc(cip, arcVertex1, arcVertex2, arcRefDir, arcAuxDir, arcNormal, arcCentrePoint, rotSym=True):
            validatedPoints.append(cip)

    return validatedPoints


def pointSegmentMinDisp(p, s0, s1):
    '''
    Get point on segment at minimum from point p.
    https://web.archive.org/web/20220121145748/http://geomalgorithms.com/index.html
    Final two boolean values indicate whether point is at or beyond endpoints of segments

    :param p: point
    :param s0: segment vertex 1
    :param s1: segment vertex 2
    :return point, bool, bool
    '''

    v = s1 - s0
    w = p - s0

    c1 = np.dot(w, v)
    if (c1 <= 0):
        return s0, np.linalg.norm(w), True, False

    c2 = np.dot(v, v)
    if (c2 <= c1):
        return s1, np.linalg.norm(p - s1), False, True

    b = c1 / c2
    pb = s0 + b * v
    return pb, np.linalg.norm(p - pb), False, False


def curveEnclosingRectangle(v1, v2, centreP, normV, radius):
    '''
    Get dimensions of box encompassing specified arc
    :param v1: vertex 1
    :param v2: vertex 2
    :param centreP: arc centre point
    :param normV: arc normal vector
    :param radius: rac radius
    :return list of box corners
    '''
    # in order to establish union of partial ellipse/circle, create box encompassing circle/ellipse
    # in the case of ellipse, use eMajorRadius
    # vertex midpoint
    midP = (v2 + v1) / 2
    midPv = np.cross(v2 - v1, normV)
    # midpoint - centroid vector, falls apart if coincident with centre point
    # midPv = midP - centreP
    midPv = midPv / np.linalg.norm(midPv)
    farP = centreP + midPv * radius
    orthoPv = np.cross(normV, midPv)
    orthoPv = orthoPv / np.linalg.norm(orthoPv)
    p1 = midP + radius * orthoPv
    p2 = midP - radius * orthoPv
    p3 = farP - radius * orthoPv
    p4 = farP + radius * orthoPv
    return [p1, p2, p3, p4]


# ================================================STEP stepcode.Part21 parsing related functions
def STEPlabelOffset(i, refStrList):
    '''account for missing text label string in early STEP versions'''
    if isinstance(refStrList[0], str):
        if refStrList[0]:
            if refStrList[0][0] == '#':
                return i - 1
            else:
                return i
        else:
            return i
    else:
        return i - 1


def cleanSubRefs(refStr):
    '''
    Flatten parameter sets and extract #d string references from stepcode parser
    :param refSt:
    :return list of reference strings
    '''

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


def ref2index(s): return int(s[1:]) - 1


def CP2point(pointRefString):
    '''
    Return np.array of STEP AP21 CARTESIAN_POINT, DIRECTION, VECTOR

    :param pointRefString: STEP reference string
    :return np.array of STEP AP21 [CARTESIAN_POINT, DIRECTION, VECTOR]
    '''

    # assert STEPobject[ref2index(pointRefString)].type_name == "CARTESIAN_POINT"
    # R3 = STEPobject[ref2index(pointRefString)].params[1]
    R3 = STEPobj[ref2index(pointRefString)].params[-1]
    return np.array([R3[0], R3[1], R3[2]])


def axis1Placement(strRef):
    '''
    Process STEP AXIS1_PLACEMENT entity
    axis: the direction of the axis.

    :param strRef: STEP reference string
    :return axisPoint, normDir, refDir
    '''
    # STEP ISO 10303-42 notes, EXPRESS descriptor
    # ENTITY Axis1_placement
    #   SUBTYPE OF (Axis_placement);
    #   axis : OPTIONAL Direction;

    # Attribute definitions:
    # axis: the direction of the axis.
    # The dimensionality of the Axis_placement_2d shall be 3.
    # The number of dimension ratios of the axis shall be 3.

    if STEPobj[ref2index(strRef)].type_name != 'AXIS1_PLACEMENT':
        print("axis2Placement3D assignment failure: " + STEPobj[ref2index(strRef)].type_name)
        return None

    axisPoint = None
    normDir = None
    subRefList = cleanSubRefs(STEPobj[ref2index(strRef)].params)

    if (STEPobj[ref2index(subRefList[0])].type_name == 'CARTESIAN_POINT'):
        axisPoint = CP2point(subRefList[0])

        if len(subRefList) > 1:
            if (STEPobj[ref2index(subRefList[1])].type_name == 'DIRECTION'):
                normDir = CP2point(subRefList[1])

    return axisPoint, normDir


def axis2Placement3D(strRef):
    '''
    Process STEP AXIS2_PLACEMENT_3D entity
    axis: the Direction that defines the second axis of the Axis_placement. (Y or V, not normal)
    ref_direction: the direction used to determine the direction of the local X axis. (or U)

    :param strRef: STEP reference string
    :return
    '''
    # STEP ISO 10303-42 notes, EXPRESS descriptor
    # process STEP AXIS2_PLACEMENT_3D entity
    # axis: the Direction that defines the second axis of the Axis_placement. (Y or V, not normal)
    # The value of this attribute need not be specified.
    # ref_direction: the direction used to determine the direction of the local X axis. (or U)
    # The value of this attribute need not be specified.
    # If axis or ref_direction is omitted, these directions are taken from the geometric coordinate system
    # If both axis and ref_direction are provided then the vector product of axis and ref_direction shall not be a null vector.

    if STEPobj[ref2index(strRef)].type_name != 'AXIS2_PLACEMENT_3D':
        print("axis2Placement3D assignment failure: " + STEPobj[ref2index(strRef)].type_name)
        return None

    axisPoint = None
    normDir = None
    refDir = None
    subRefList = cleanSubRefs(STEPobj[ref2index(strRef)].params)

    if (STEPobj[ref2index(subRefList[0])].type_name == 'CARTESIAN_POINT'):
        axisPoint = CP2point(subRefList[0])

        if len(subRefList) > 1:
            if (STEPobj[ref2index(subRefList[1])].type_name == 'DIRECTION'):
                normDir = CP2point(subRefList[1])

                # relying on order specified
                # https://www.steptools.com/stds/smrl/data/modules/elemental_geometric_shape/sys/4_info_reqs.htm#elemental_geometric_shape_arm.axis_placement
                if len(subRefList) > 2:
                    if (STEPobj[ref2index(subRefList[2])].type_name == 'DIRECTION'):
                        refDir = CP2point(subRefList[2])

    return axisPoint, normDir, refDir


# def axis2Placement3D(strRef):
#     '''
#     Process STEP AXIS2_PLACEMENT_3D entity
#     axis: the Direction that defines the second axis of the Axis_placement. (Y or V, not normal)
#     ref_direction: the direction used to determine the direction of the local X axis. (or U)
#
#     :param strRef: STEP reference string
#     :return axisPoint, auxDir, refDir
#     '''
#     # process STEP AXIS2_PLACEMENT_3D entity
#     assert STEPobj[ref2index(strRef)].type_name == 'AXIS2_PLACEMENT_3D'
#     # D = {}
#
#     if '#' in STEPobj[ref2index(strRef)].params[0]:  # check for name string
#         # not a robust approach as certain values are optional
#         offset = 0
#     else:
#         offset = 1
#
#     axisPoint = STEPobj[ref2index(strRef)].params[offset]
#     assert STEPobj[ref2index(axisPoint)].type_name == 'CARTESIAN_POINT'
#     axisPoint = CP2point(axisPoint)
#     # D['axisPoint'] = axisPoint
#
#     auxDir = STEPobj[ref2index(strRef)].params[offset + 1]
#     assert STEPobj[ref2index(auxDir)].type_name == 'DIRECTION'
#     auxDir = CP2point(auxDir)
#     # D['auxDir'] = auxDir
#
#     refDir = STEPobj[ref2index(strRef)].params[offset + 2]
#     assert STEPobj[ref2index(refDir)].type_name == 'DIRECTION'
#     refDir = CP2point(refDir)
#     # D['refDir'] = refDir
#
#     return axisPoint, auxDir, refDir


# search for cartesian points & check if points correspond with referenced curve/edge vertexes

# TRIMMED_CURVE
# ENTITY trimmed_curve
# SUBTYPE OF (bounded_curve);
# basis_curve : curve;
# trim_1 : SET[1:2] OF trimming_select;
# trim_2 : SET[1:2] OF trimming_select;
# sense_agreement : BOOLEAN;
# master_representation : trimming_preference;

# provides either parametric, point or angle end points to arc
# can be used to update vertex1/2


# local STEP (AdvancedFaceSurfaces) data assistant functions


def insideOutsideSurfaceTest(c, mPoint, AFSobj):
    '''
    Project ray from centroid (c) toward mPoint in order to determine whether this ray intersects the surface defined in AFSobj
    The angle of edge vertices of an intercepted face are calculated with respect to this point and summed
    Return true if point within surface edges ("Angle Sum Algorithm")
    No accommodation for holes in surface
    :param AFSobj: parsed surfaces data object
    :param c: estimated surface median point based on extrema points identified
    :param mPoint:
    :return bool True if within surface
    '''


    def vertexAngle(mPoint, v1, v2):
        # project vertex points to a plane defined by mPoint and sNorm
        m_v1 = v1 - mPoint
        p_v1 = v1 - np.dot(m_v1, sNorm) * sNorm

        m_v2 = v2 - mPoint
        p_v2 = v2 - np.dot(m_v2, sNorm) * sNorm

        mppv1 = p_v1 - mPoint
        mppv2 = p_v2 - mPoint

        sign = np.sign(np.cross(mppv1, mppv2).dot(sNorm))

        v1v2cos = np.dot(mppv1, mppv2) / (np.linalg.norm(mppv1) * np.linalg.norm(mppv2))  # cosine

        if np.abs(v1v2cos) < eps_STEP_AP21:
            print("insideOutsideSurfaceTest(): edge warning")

        # zero value implies mPoint lies on edge
        return (sign * np.arccos(np.clip(v1v2cos, -1, 1)))

    sNormDen = np.linalg.norm(mPoint - c)
    if sNormDen > eps_STEP_AP21:
        sNorm = (mPoint - c) / sNormDen
    else:
        # point is at centroid, find normal from provided surface parameters
        print("insideOutsideSurfaceTest(): median centroid at surface, surface normal substituted")
        # add offset to everything and try again?
        # mirror mPoint around surface normal? just use surface normal instead? -> requires earlier calc of axis2Placement3D
        # _1=1
        sNorm = AFSobj['ParsedSurface']['normDir']

    edgeTotalAngle = 0
    vertexLoopCount = 0
    for edge in AFSobj['outerBoundEdgeLoop']:
        v1 = edge['vertex1']
        v2 = edge['vertex2']

        # #if not np.isclose(v1, v2).all():
        # if edge['vertex1ref'] == edge['vertex2ref']: # loop, discard
        #     if edge['typeName'] == 'CIRCLE' or edge['typeName'] == 'B_SPLINE_CURVE_WITH_KNOTS' or edge['typeName'] == 'ELLIPSE':
        #         vertexLoopCount +=1 # should be discarded if point is not within loop?

        # catch spline and convert to polyline vertices
        if edge['typeName'] == 'B_SPLINE_CURVE_WITH_KNOTS':
            if edge['vertex1ref'] == edge['vertex2ref']:  # loop, discard
                vertexLoopCount += 1  # should be discarded if point is not within loop?
            else:
                splinePolyline, splinePolylineU = splineToPolyline(edge['BsplineKnotCurve'], curveSampleFactor=1)
                # find closest point on spline and ensure that this is not mPoint

                # nearestSplinePointU = splineCurveMinMaxPointDisp(edge['BsplineKnotCurve'], mPoint, maxSearch=False)

                # if mPoint[0]==0: # and mPoint[1]== -0.30397087 and mPoint[2]==-16.61759168:
                #     _1=1
                nearestSplinePointU = BsplineCurveExtremaDisp(edge['BsplineKnotCurve'],
                                                              mPoint,
                                                              maxSearch=False,
                                                              uv_xyz=True)

                if len(nearestSplinePointU) == 0: # usually means mPoint is closest to endpoints
                    #
                    if np.linalg.norm(mPoint - v1) <= np.linalg.norm(mPoint - v2):
                        nearestSplinePoint = v1
                        nearestSplinePointU = 0.
                    else:
                        nearestSplinePoint = v2
                        nearestSplinePointU = 1.
                else:
                    nearestSplinePoint = np.array(edge['BsplineKnotCurve'].evaluate_list(nearestSplinePointU))

                if np.linalg.norm(mPoint - nearestSplinePoint) < eps_STEP_AP21:
                    print('insideOutsideSurfaceTest(): spline edge warning - ignore polyline warnings')
                if not any([np.isclose(nearestSplinePointU, splU) for splU in splinePolylineU]):
                    # has to be ordered correctly
                    splinePolylineU = np.insert(splinePolylineU, splinePolylineU.searchsorted(nearestSplinePointU),
                                                nearestSplinePointU)
                    iu = np.argmin(np.abs(splinePolylineU - nearestSplinePointU))
                    splinePolyline = [np.array(spl) for spl in splinePolyline]
                    splinePolyline = np.insert(splinePolyline, iu, nearestSplinePoint, axis=0)
                for i in range(len(splinePolyline) - 1):
                    edgeTotalAngle += vertexAngle(mPoint, splinePolyline[i], splinePolyline[i + 1])
                # _1=1

        # account for circle or ellipse arcs, find nearest point on arc to mPoint,
        # add a point there, effectively converting arc into two segments
        elif (edge['typeName'] == 'CIRCLE'):
            if edge['vertex1ref'] == edge['vertex2ref']:  # loop, discard
                vertexLoopCount += 1  # should be discarded if point is not within loop?
            else:
                nearArcPoint, _ = pointCircleMinMaxDisp(mPoint,
                                                        edge['axisPoint'],
                                                        edge['normDir'],
                                                        edge['radius'])
                if pointInArc(nearArcPoint,
                              v1, v2,
                              edge['refDir'],
                              edge['auxDir'],
                              edge['normDir'],
                              edge['axisPoint'],
                              rotSym=False):
                    edgeTotalAngle += vertexAngle(mPoint, v1, nearArcPoint)
                    edgeTotalAngle += vertexAngle(mPoint, nearArcPoint, v2)
                    print("untested circle arc division")

        elif edge['typeName'] == 'ELLIPSE':
            if edge['vertex1ref'] == edge['vertex2ref']:  # loop, discard
                vertexLoopCount += 1  # should be discarded if point is not within loop?
            else:
                nearArcPoint, _ = pointEllipseMinMaxDisp(mPoint,
                                                         edge['axisPoint'],
                                                         edge['refDir'],
                                                         edge['auxDir'],
                                                         edge['normDir'],
                                                         edge['majorRadius'],
                                                         edge['minorRadius'],
                                                         interior=False)
                if pointInArc(nearArcPoint,
                              v1, v2,
                              edge['refDir'],
                              edge['auxDir'],
                              edge['normDir'],
                              edge['axisPoint'],
                              rotSym=False):
                    edgeTotalAngle += vertexAngle(mPoint, v1, nearArcPoint)
                    edgeTotalAngle += vertexAngle(mPoint, nearArcPoint, v2)
                    print("untested ellipse arc division")

        else:
            # project vertex points to a plane defined by mPoint and sNorm
            edgeTotalAngle += vertexAngle(mPoint, v1, v2)

    # surfaces that curve > 90deg can give false negative, assume inside/outside test is to detect nearest surface

    # exception for cylinder, conic, sphere
    if (AFSobj['SurfaceTypeName'] in ['CYLINDRICAL_SURFACE', ]) or (AFSobj['SurfaceTypeName'] in ['CONICAL_SURFACE', ]):
        # edges angle sum to zero from 2x loops and seam edge ---------------------needs testing
        if (np.abs(edgeTotalAngle % (2 * np.pi)) < eps_STEP_AP21) and (
                vertexLoopCount >= 1):  # (len(AFSobj['outerBoundEdgeLoop']) < 4):
            # only need to determine whether point lies within section of cone axis within cone as circle or ellipse loops at ends signify
            #   fully enclosed surface
            # vertex angles that sum to zero will still indicate a crossing, for either CW, or ACW direction configuration

            # if method is uniquely used for closest surface point on a plane encompassing projection and cylinder axis
            # can then ignore far-side intersection

            # ?? find nearest points on circles/ellipses/bsplines at each curved edge, then test if both points are to one side??
            curvedEdgeNearestPoint = []
            normDir = np.array(AFSobj['ParsedSurface']['normDir'])
            for edge in AFSobj['outerBoundEdgeLoop']:
                if edge['typeName'] == 'CIRCLE':
                    nearArcPoint, _ = pointCircleMinMaxDisp(mPoint,
                                                            edge['axisPoint'],
                                                            edge['normDir'],
                                                            edge['radius'])
                    curvedEdgeNearestPoint.append(np.array(nearArcPoint))

                elif edge['typeName'] == 'ELLIPSE':
                    nearArcPoint, _ = pointEllipseMinMaxDisp(mPoint,
                                                             edge['axisPoint'],
                                                             edge['refDir'],
                                                             edge['auxDir'],
                                                             edge['normDir'],
                                                             edge['majorRadius'],
                                                             edge['minorRadius'])
                    curvedEdgeNearestPoint.append(np.array(nearArcPoint))

                elif edge['typeName'] == 'B_SPLINE_CURVE_WITH_KNOTS':

                    # should go to a plane-curve intersection function

                    # define plane through axisPoint, mPoint, normDir
                    axisPoint = np.array(AFSobj['ParsedSurface']['axisPoint'])
                    # normDir = np.array(AFSobj['ParsedSurface']['normDir'])
                    planeNorm = np.cross((mPoint - axisPoint), normDir)
                    orthoPlane = np.cross(planeNorm, normDir)
                    mSign = np.dot(axisPoint - mPoint, orthoPlane) > 0
                    sPolyline, sPolylineU = splineToPolyline(edge['BsplineKnotCurve'], curveSampleFactor=1)
                    nearArcPointSet = []
                    U = []

                    for i in range(0, len(sPolyline) - 1):
                        nearArcPoint = intersectSegmentPlane(np.array(sPolyline[i]),
                                                             np.array(sPolyline[i + 1]),
                                                             planeNorm,
                                                             mPoint)
                        if nearArcPoint is not None:
                            nearArcPointSet.append(nearArcPoint)
                            U.append(i)

                    # determine intersection on mPoint side of plane
                    # test whether vector from mPoint to intersection is parallel with normDir
                    if len(nearArcPointSet) == 0:
                        break
                    elif len(nearArcPointSet) == 1:
                        nearArcPoint = nearArcPointSet[0]
                        iu = U[0]
                    else:
                        for i_naps, naps in enumerate(nearArcPointSet):
                            napsSign = np.dot(axisPoint - naps, orthoPlane) > 0
                            if napsSign == mSign:
                                # if np.isclose(np.cross(mPoint - naps, normDir), np.array([0., 0., 0.]), eps_STEP_AP21).all():
                                nearArcPoint = naps
                                iu = U[i_naps]

                    # binary search over u-interval
                    if iu + 1 > len(sPolylineU):
                        _1 = 1
                    hiU = sPolylineU[iu + 1]
                    hiP = sPolyline[iu + 1]
                    loU = sPolylineU[iu]
                    loP = sPolyline[iu]

                    while (hiU - loU) > eps_STEP_AP21:  # 1E-10: #==================================================
                        midU = loU + (hiU - loU) / 2
                        midP = np.array(edge['BsplineKnotCurve'].evaluate_single(midU))
                        nearArcPoint = intersectSegmentPlane(loP, midP, planeNorm, mPoint)
                        if nearArcPoint is not None:
                            hiP = midP
                            hiU = midU
                        else:
                            nearArcPoint = intersectSegmentPlane(midP, hiP, planeNorm, mPoint)
                            if nearArcPoint is not None:
                                loP = midP
                                loU = midU
                            else:
                                print("fail")

                    curvedEdgeNearestPoint.append(np.array(nearArcPoint))

            if len(curvedEdgeNearestPoint) > 2:
                print("insideOutsideSurfaceTest cylinder/conical fail")
                print(curvedEdgeNearestPoint)
                print(mPoint)
            if len(curvedEdgeNearestPoint) < 2:
                # print("insideOutsideSurfaceTest cylinder/conical fail")
                return False
            # project intersections to plane where axis is normal (normDir), through mPoint
            # if all points have same sign, all are on same side
            if not all([np.dot(mPoint - cenp, normDir) > 0 for cenp in curvedEdgeNearestPoint]):
                return True
            else:
                return False

    # if AFSobj['SurfaceTypeName'] in ['CONICAL_SURFACE',]:
    #     # edges angle sum to zero from 2x loops and seam edge ---------------------needs testing
    #     if (np.abs(edgeTotalAngle) < eps_STEP_AP21) and (vertexLoopCount == 2): # (len(AFSobj['outerBoundEdgeLoop']) < 4):
    #         # only need to determine whether point lies within section of cone axis within cone
    #         return True
    #     else:
    #         return False

    if AFSobj['SurfaceTypeName'] in ['PLANE', ]:  # full circle
        if (np.abs(edgeTotalAngle) < eps_STEP_AP21) and (len(AFSobj['outerBoundEdgeLoop']) < 2):
            return True
        elif ((np.abs(edgeTotalAngle) % (2 * np.pi)) < eps_STEP_AP21) and (len(AFSobj['outerBoundEdgeLoop']) > 1):
            return True
        else:
            return False

    # may have to take VOR into account depending on orientation of surface,
    # e.g. concave or convex spherical surface
    # assume that minima point is always closest to centroid?
    # if not edge['VOR']:
    if not AFSobj['SurfaceNormalOutwards']:
        if (edgeTotalAngle % (2 * np.pi)) < eps_STEP_AP21 and (edgeTotalAngle > 0):
            return True
        else:
            return False
    else:
        if (edgeTotalAngle % (2 * np.pi)) < eps_STEP_AP21 and (edgeTotalAngle < 0):
            return True
        else:
            return False



testDir = os.path.normpath(
    r"/media/foobert/Dell2HDD/STEP_test_files"
)
primitivesDir = os.path.normpath(
    testDir + r"primitives"
)

# filepath = primitivesDir + "Cube/unit_cube_inc8.0_blend.06.stp"

# filepath = primitivesDir + "Cylinder/unit_cyl.stp"
# filepath = primitivesDir + "Primitive_Cone-PartCone.step"

#filepath = testDir + "RR_STEP_test_1A.step" # pass.. I think, no representation available
#filepath = testDir + "RR_STEP_test_N.step" # fail, no trimmed_curve circle code================
# filepath = testDir + "00000001_1ffb81a71e5b402e966b9341_step_000.step"
# filepath = testDir + "00000010_b4b99d35e04b4277931f9a9c_step_000.step"
#filepath = testDir + os.sep + "LEA-M8F(AP203).STEP" # too large to debug?

#filepath = testDir + "9341_step_000.step" #FAIL - reverse list from spline=======================
#filepath = testDir + "revolve_spline.step" #PASS

# filepath = testDir + "/simpleSpindle1.step" #PASS
# filepath = testDir + "Maltese_cruciform.step" #PASS
# filepath = testDir + "TiltedConeAP203.step" # pass
filepath = "OffsetCone_AP214_noPlacement_noParametric.step"
#filepath = "OffsetCone-PartCone.step" #PASS2
#filepath = "TiltedCylinder4_AP214_PC.step" #PASS2
#filepath = "TiltedCylinder5-AP214_draft.step" #PASS2
#filepath = "TiltedCylinder4c-AP203-ESOP.step" #PASS2
#filepath = "Cylinder5_AP214_PC.step" #PASS2
#filepath = "DrexelBlendedCylinder_topOvalPlane.step" #PASS2
#filepath = "DrexelBlendedCylinder_curvedBlend.step" #PASS2
#filepath = "DrexelBlendedCylinder_midriffTube.step" #PASS2 - note disparity between ellipse minima and vertexes

# can a minima exist on an NURB surface edge? test all local NURB boundary minima against centroid disp??

#filepath = "TiltedCylinder2.step" #PASS2
#filepath = "TiltedCylinder3.step" #PASS2
#filepath = "primitives/Cube/unit_cube.stp" # pass2
#filepath = "Drexel_blended_ellipse_plain.step" # fail, no surfaces ellipse
#filepath = "Drexel_blended_cylinder.step" # pass2
#filepath = "DrexelBlendedCylinder_curvedBlend.step" # pass2
#filepath = "TiltedCylinder.step" # pass2
#filepath = "Synth_ellipse_plain.step" # pass (no Advanced Face Surface)

#filepath = primitivesDir + os.sep + "Cube/unit_cube.stp" #PASS2
# OFF format
#filepath = testDir + "testOFF.off"

filepath = testDir + os.sep + filepath
print(filepath)

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

def parseSTEPstr(s):
    # following operations unique to stepcode.Part21 Parser
    s = s.replace(";\n", ";;")
    s = s.replace("\n", "")
    s = s.replace(";;", ";\n")

    P = Parser()
    # model = P.parse(/home/foobert/.config/JetBrains/PyCharmCE2023.1/scratches/simple_schema.exp)
    model = P.parse(s)
    SO = model.sections[0].entities

    # fixing step files with non-contiguous line and #ref indices
    if ref2index(SO[-1].ref) is not len(SO):
        SO.sort(key=lambda x: int(x.ref[1:]))
        padded_SO = [0] * int(SO[-1].ref[1:])
        for ss in SO:
            padded_SO[ref2index(ss.ref)] = ss
        #SO = padded_SO

    return padded_SO


STEPobj = parseSTEPstr(STEPdata)


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


def segmentParse(c, ParsedEdge, calcExtrema=False):
    '''
    Parse STEP segment representation

    :param ParsedEdge: : parsed edge data object
    :param c: estimated surface median point based on extrema points identified
    :return
    '''

    # STEP ISO 10303-42 notes, EXPRESS descriptor
    # ENTITY Line
    #   SUBTYPE OF (Curve);
    #   point : Cartesian_point;
    #   line_direction : Direction;

    # if ParsedEdge.get('subedgeParams') is not None:
    #     edgeParams_ = ParsedEdge['subEdgeParams']
    # else:

    if ParsedEdge.get('lineVectorDir') is None:

        edgeParams_ = ParsedEdge['edgeParams']
        cleanEdgeParams = cleanSubRefs(edgeParams_)
        for cet in cleanEdgeParams:
            if (STEPobj[ref2index(cet)].type_name == 'CARTESIAN_POINT'):
                linePoint = CP2point(cet)
                ParsedEdge['linePoint'] = linePoint

            if (STEPobj[ref2index(cet)].type_name == 'VECTOR'):
                lineVectorDisp = STEPobj[ref2index(cet)].params[-1]
                ParsedEdge['lineVectorDisp'] = lineVectorDisp
                lineVectorDir = STEPobj[ref2index(cet)].params[-2]
                if (STEPobj[ref2index(lineVectorDir)].type_name == 'DIRECTION'):
                    lineVectorDir = CP2point(lineVectorDir)
                    ParsedEdge['lineVectorDir'] = lineVectorDir

    if calcExtrema:

        # v1 = ParsedEdge['vertex1']
        # v2 = ParsedEdge['vertex2']  # v1 always equals v2 for full circle

        minPoint, minPointCentroidDisp, v1ext, v2ext = pointSegmentMinDisp(c,
                                                                           ParsedEdge['vertex1'],
                                                                           ParsedEdge['vertex2'])

        # recalculated every centroid update
        ParsedEdge['pointFeature']['xyz'] = [ParsedEdge['vertex1'],
                                             ParsedEdge['vertex2']]
        ParsedEdge['pointFeature']['centroidDisp'] = [np.linalg.norm(c - ParsedEdge['vertex1']).item(),
                                                      np.linalg.norm(c - ParsedEdge['vertex2']).item()]
        ParsedEdge['pointFeature']['u'] = [0, 1]  # normalise to [0, 1]

        if not v1ext and not v2ext:
            ParsedEdge['pointFeature']['xyz'].insert(1, minPoint)
            ParsedEdge['pointFeature']['centroidDisp'].insert(1, minPointCentroidDisp)
            u = np.linalg.norm(minPoint - ParsedEdge['vertex1']) / np.linalg.norm(ParsedEdge['vertex2'] - ParsedEdge['vertex1'])
            ParsedEdge['pointFeature']['u'].insert(1, u)


# def segmentParse(c, ParsedEdge):
#     '''
#     Parse STEP segment representation
#
#     :param ParsedEdge: : parsed edge data object
#     :param c: estimated surface median point based on extrema points identified
#     :return
#     '''
#
#     # STEP ISO 10303-42 notes, EXPRESS descriptor
#     # ENTITY Line
#     #   SUBTYPE OF (Curve);
#     #   point : Cartesian_point;
#     #   line_direction : Direction;
#
#     # if ParsedEdge.get('subedgeParams') is not None:
#     #     edgeParams_ = ParsedEdge['subEdgeParams']
#     # else:
#     edgeParams_ = ParsedEdge['edgeParams']
#
#     cleanEdgeParams = cleanSubRefs(edgeParams_)
#     for cet in cleanEdgeParams:
#         if (STEPobj[ref2index(cet)].type_name == 'CARTESIAN_POINT'):
#             linePoint = CP2point(cet)
#             ParsedEdge['linePoint'] = linePoint
#
#         if (STEPobj[ref2index(cet)].type_name == 'VECTOR'):
#             lineVectorDisp = STEPobj[ref2index(cet)].params[-1]
#             ParsedEdge['lineVectorDisp'] = lineVectorDisp
#             lineVectorDir = STEPobj[ref2index(cet)].params[-2]
#             if (STEPobj[ref2index(lineVectorDir)].type_name == 'DIRECTION'):
#                 lineVectorDir = CP2point(lineVectorDir)
#                 ParsedEdge['lineVectorDir'] = lineVectorDir
#
#     minPoint, minPointCentroidDisp, v1ext, v2ext = pointSegmentMinDisp(c, ParsedEdge['vertex1'],
#                                                                        ParsedEdge['vertex2'])
#     # todo, if there are multiple minima, return list? store minpoint as list?
#
#     v1 = ParsedEdge['vertex1']
#     v2 = ParsedEdge['vertex2']  # v1 always equals v2 for full circle
#
#     # recalculated every centroid update
#     ParsedEdge['pointFeature']['xyz'] = [v1, v2]
#     ParsedEdge['pointFeature']['centroidDisp'] = [np.linalg.norm(centroid - v1).item(),
#                                                   np.linalg.norm(centroid - v2).item()]
#     ParsedEdge['pointFeature']['u'] = [0, 1]  # normalise to [0, 1]
#
#     if not v1ext and not v2ext:
#         ParsedEdge['pointFeature']['xyz'].insert(1, minPoint)
#         ParsedEdge['pointFeature']['centroidDisp'].insert(1, minPointCentroidDisp)
#         u = np.linalg.norm(minPoint - ParsedEdge['vertex1']) / np.linalg.norm(
#             ParsedEdge['vertex2'] - ParsedEdge['vertex1'])
#         ParsedEdge['pointFeature']['u'].insert(1, u)
#


def circleParse(c, ParsedEdge, calcExtrema=False):
    '''
    Parse STEP circle representation

    :param ParsedEdge: : parsed edge data object
    :param c: estimated surface median point based on extrema points identified
    :return
    '''

    # STEP ISO 10303-42 notes
    # axis: the Direction that defines the second axis of the Axis_placement. The value of this attribute need not be specified.
    # ref_direction: the direction used to determine the direction of the local X axis.
    # The value of this attribute need not be specified. If axis or ref_direction is omitted, these directions are taken from the geometric coordinate system

    # if ParsedEdge.get('subEdgeParams') is not None:
    #     edgeParams_ = ParsedEdge['subEdgeParams']
    # else:

    if ParsedEdge.get('normDir') is None:

        edgeParams_ = ParsedEdge['edgeParams']

        radius = edgeParams_[-1]
        ParsedEdge['radius'] = radius

        axisPoint, normDir, refDir = axis2Placement3D(edgeParams_[-2])
        ParsedEdge['axisPoint'] = axisPoint
        ParsedEdge['normDir'] = normDir
        ParsedEdge['refDir'] = refDir
        # axisDir, refDir guaranteed unit normals?
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        ParsedEdge['auxDir'] = auxDir

    if calcExtrema:

        # test if cylinder is rotationally symmetrical, if c is close to axisDir through axisPoint
        # rotSymDisp = np.linalg.norm((c - axisPoint) - np.dot((c - axisPoint), normDir) * normDir)

        # why doesn't the below work?-------------------------------------------------------------------------
        # rotSymDisp = np.linalg.norm(c - pointProjectAxis(c, axisPoint, normDir))

        # # project centroid to disc plane
        # # distance d of point p to the plane is dot product of plane unit normal, normAxis
        # # with the difference vector from the point in the plane P and S: d = (S − P) · u
        # # distance is signed: positive when S is on the side of the plane where u faces and
        # # negative on the other side. zero, S is in the plane

        # min/max point on arc wrt centroid
        minPoint, maxPoint = pointCircleMinMaxDisp(
            c,
            ParsedEdge['axisPoint'],
            ParsedEdge['normDir'],
            ParsedEdge['radius'],
            interior=False
        )

        # 'vertex1' always equals 'vertex2' for full circle

        # recalculated every centroid update
        ParsedEdge['pointFeature']['xyz'] = [ParsedEdge['vertex1'], ParsedEdge['vertex2']]
        ParsedEdge['pointFeature']['centroidDisp'] = [np.linalg.norm(c - ParsedEdge['vertex1']).item(),
                                                      np.linalg.norm(c - ParsedEdge['vertex2']).item()]
        ParsedEdge['pointFeature']['u'] = [0, 1]  # normalise to [0, 1]

        extrema = []
        if maxPoint is not None: extrema.append(np.array(maxPoint))
        if minPoint is not None: extrema.append(np.array(minPoint))

        if ParsedEdge['vertex1ref'] == ParsedEdge['vertex2ref']:
            # provide an extra vertex on opposite side of rotSym axis
            # to balance median point centroid search
            #oppVertex = [ParsedEdge['axisPoint'][0] - ParsedEdge['vertex1'][0], ParsedEdge['axisPoint'][1] - ParsedEdge['vertex1'][1], ParsedEdge['axisPoint'][2] - ParsedEdge['vertex1'][2]]
            oppVertex = ParsedEdge['axisPoint'] - (ParsedEdge['vertex1'] - ParsedEdge['axisPoint'])
            extrema.append(np.array(oppVertex))

        pointsInArc, pointOrderIndex, uArc = pointOrderInArc(extrema,
                                                             ParsedEdge['vertex1'],
                                                             ParsedEdge['vertex2'],
                                                             ParsedEdge['refDir'],
                                                             ParsedEdge['auxDir'],
                                                             ParsedEdge['normDir'],
                                                             ParsedEdge['axisPoint'])

        pointsInArc = [pointsInArc[poi] for poi in pointOrderIndex]
        extrema = [extrema[poi] for poi in pointOrderIndex]
        uArc = [uArc[poi] for poi in pointOrderIndex]
        # remove points that are not in arc
        extrema = [i for (i, v) in zip(extrema, pointsInArc) if v]

        # test whether extrema values already coincide with vertex values,
        for ex_i, ex in enumerate(extrema):
            if not np.isclose(ParsedEdge['vertex1'], ex).all() and not np.isclose(ParsedEdge['vertex2'], ex).all():
                if not any([np.isclose(uArc[ex_i], uu) for uu in ParsedEdge['pointFeature']['u']]):
                    ParsedEdge['pointFeature']['xyz'].insert(-1, ex)
                    ParsedEdge['pointFeature']['centroidDisp'].insert(-1, np.linalg.norm(ex - c).item())
                    ParsedEdge['pointFeature']['u'].insert(-1, uArc[ex_i].item())  # check u proximity?

        Cdisp = np.dot((c - ParsedEdge['axisPoint']), ParsedEdge['normDir'])
        Cproj = c - np.dot(Cdisp, ParsedEdge['normDir'])  # point projected to plane along orthogonal
        rotSymDisp = np.linalg.norm(Cproj - ParsedEdge['axisPoint'])
        # print("rotSymDisp: " + str(rotSymDisp))

        # centroid orthogonal to centrePoint
        if rotSymDisp < eps_STEP_AP21:
            # case of minPoint equidistant from circle centre - edge case
            # assign rotSymCentre feature point at centre

            ParsedEdge['rotSymFeature'] = dict()
            ParsedEdge['rotSymFeature']['rotSymCentre'] = [ParsedEdge['axisPoint'],]
            ParsedEdge['rotSymFeature']['rotSymRadius'] = [ParsedEdge['radius'],]
            ParsedEdge['rotSymFeature']['radiusCentroidDisp'] = \
                [np.sqrt( ParsedEdge['radius'] ** 2 + np.linalg.norm(ParsedEdge['axisPoint'] - c) ** 2).item(),]

            # # remove previously found point extrema
            # ParsedEdge['pointFeature']['xyz'] = []
            # ParsedEdge['pointFeature']['centroidDisp'] = []
            # ParsedEdge['pointFeature']['u'] = []

            # alternative is to provide an extra vertex on opposide side of rotSym axis
            # to balance median point centroid search


def ellipseParse(c, ParsedEdge, calcExtrema=False):
    '''
    Parse STEP ellipse representation

    :param ParsedEdge: : parsed edge data object
    :param c: estimated surface median point based on extrema points identified
    :return
    '''

    # STEP ISO 10303-42 notes
    # axis: the Direction that defines the second axis of the Axis_placement. The value of this attribute need not be specified.
    # ref_direction: the direction used to determine the direction of the local X axis.
    # The value of this attribute need not be specified. If axis or ref_direction is omitted, these directions are taken from the geometric coordinate system

    # first_semi_axis: half the length of the first diameter of the Ellipse.
    # second_semi_axis: half the length of the second diameter of the Ellipse.

    # if ParsedEdge.get('subEdgeParams') is not None:
    #     edgeParams_ = ParsedEdge['subEdgeParams']
    # else:

    if ParsedEdge.get('normDir') is None:

        edgeParams = ParsedEdge['edgeParams']

        majorRadius = edgeParams[-2]
        ParsedEdge['majorRadius'] = majorRadius

        minorRadius = edgeParams[-1]
        ParsedEdge['minorRadius'] = minorRadius

        axisPoint, normDir, refDir = axis2Placement3D(edgeParams[-3])
        ParsedEdge['axisPoint'] = axisPoint
        ParsedEdge['normDir'] = normDir
        ParsedEdge['refDir'] = refDir
        # axisDir, refDir guaranteed unit normals?
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        ParsedEdge['auxDir'] = auxDir

    if calcExtrema:

        minPoint, maxPoint = pointEllipseMinMaxDisp(
            c,
            ParsedEdge['axisPoint'],
            ParsedEdge['refDir'],
            ParsedEdge['auxDir'],
            ParsedEdge['normDir'],
            ParsedEdge['majorRadius'],
            ParsedEdge['minorRadius'],
            interior=False
        )

        # test if minPoint, maxPoint is on segment between vertex1 & vertex2

        # recalculated every centroid update
        ParsedEdge['pointFeature']['xyz'] = [ParsedEdge['vertex1'],
                                             ParsedEdge['vertex2']]
        ParsedEdge['pointFeature']['centroidDisp'] = [np.linalg.norm(c - ParsedEdge['vertex1']).item(),
                                                      np.linalg.norm(c - ParsedEdge['vertex2']).item()]
        ParsedEdge['pointFeature']['u'] = [0, 1]  # normalise to [0, 1]

        extrema = []
        if maxPoint is not None: extrema.append(np.array(maxPoint))
        if minPoint is not None: extrema.append(np.array(minPoint))

        if ParsedEdge['vertex1ref'] == ParsedEdge['vertex2ref']:
            # provide an extra vertex on opposite side of rotSym axis
            # to balance median point centroid search
            oppVertex = [c[0] - ParsedEdge['vertex1'][0],
                         c[1] - ParsedEdge['vertex1'][1],
                         c[2] - ParsedEdge['vertex1'][2]]
            extrema.append(np.array(oppVertex))

        pointsInArc, pointOrderIndex, uArc = pointOrderInArc(extrema,
                                                             ParsedEdge['vertex1'],
                                                             ParsedEdge['vertex2'],
                                                             ParsedEdge['refDir'],
                                                             ParsedEdge['auxDir'],
                                                             ParsedEdge['normDir'],
                                                             ParsedEdge['axisPoint'])
        pointsInArc = [pointsInArc[poi] for poi in pointOrderIndex]
        extrema = [extrema[poi] for poi in pointOrderIndex]
        uArc = [uArc[poi] for poi in pointOrderIndex]
        # drop points not in arc
        extrema = [i for (i, v) in zip(extrema, pointsInArc) if v]

        # test whether extrema values already coincide with vertex values,
        for ex_i, ex in enumerate(extrema):
            if not np.isclose(ParsedEdge['vertex1'], ex).all() and not np.isclose(ParsedEdge['vertex2'], ex).all():
                ParsedEdge['pointFeature']['xyz'].insert(-1, ex)
                ParsedEdge['pointFeature']['centroidDisp'].insert(-1, np.linalg.norm(ex - c).item())
                ParsedEdge['pointFeature']['u'].insert(-1, uArc[ex_i].item())  # check u proximity?



def boundedCurve(c, complexEntityList, calcExtrema=False):
    pass


# def complexEntityParse(complexEntityRef):
#     typeNameList = [cel.type_name for cel in STEPobj[complexEntityRef].params]
#     if typeNameList[0] == 'BOUNDED_SURFACE':
#         boundedSurface(np.array([0.,0.,0.]), ref2index(se2))
#
#         SurfaceClass['SurfaceTypeName'] = 'BOUNDED_SURFACE'
#         # NOTE    A surface, such as a cylinder or sphere, closed in one, or two, parametric directions is not a bounded surface.
#         #             SurfaceClass = {
#         #                 'SurfaceNormalOutwards': SurfaceNormalOutwards,
#         #                 'SurfaceTypeName': None,
#         #                 'SurfaceRef': None,
#         #                 'SurfaceParams': None,
#         #                 'EdgeLoopList': [],
#         #             }
#         #                         SurfaceClass['SurfaceTypeName'] = afTypeName
#         #                         SurfaceClass['SurfaceRef'] = se2
#         #                         SurfaceClass['SurfaceParams'] = STEPobj[
#         #                             ref2index(se2)
#         #                         ].params


def boundedSurface(AFSobj, c, calcExtrema=False):
    # todo untested
    # ENTITY bounded_surface
    # SUPERTYPE OF (ONEOF(b_spline_surface, rectangular_trimmed_surface, curve_bounded_surface, rectangular_composite_surface))
    # SUBTYPE OF (surface);

    #  assume provision of a list containing elements of rational weighted B spline, a complex entity
    #assert(STEPobj[complexEntityList].params[0].type_name == 'BOUNDED_SURFACE')

    # test whether AFSobj instance is simply a list where Bounded_Surface is the first reference,
    # or whether the list is AFSobj['SurfaceRef']

    if isinstance(AFSobj, list):
        complexEntity = AFSobj
    elif isinstance(AFSobj, dict):
        complexEntity = STEPobj[ref2index(AFSobj['SurfaceRef'])].params
    else:
        print('unknown complex entity formulation')

    typeNameList = [cel.type_name for cel in complexEntity]

    # ENTITY

    # SUPERTYPE OF (ONEOF(b_spline_surface, rectangular_trimmed_surface, curve_bounded_surface, rectangular_composite_surface))
    # SUBTYPE OF (surface)
    if typeNameList[0] != 'BOUNDED_SURFACE':
        print(" BoundedSurface type error")
        return

    if ('B_SPLINE_SURFACE' not in typeNameList) or \
            ('B_SPLINE_SURFACE_WITH_KNOTS' not in typeNameList) or \
            ('RATIONAL_B_SPLINE_SURFACE' not in typeNameList):
        print(" unhandled complex entity")

    ParsedSurface = {}

    if 'B_SPLINE_SURFACE' in typeNameList:
        fieldRef = typeNameList.index('B_SPLINE_SURFACE')

        # ENTITY b_spline_surface
        # SUPERTYPE OF (ONEOF(b_spline_surface_with_knots, uniform_surface, quasi_uniform_surface, bezier_surface)
        # ANDOR rational_b_spline_surface)
        # SUBTYPE OF (bounded_surface);
        # u_degree : INTEGER;
        # v_degree : INTEGER;
        # control_points_list : LIST [2:?] OF LIST [2:?] OF cartesian_point;
        # surface_form : b_spline_surface_form;
        # u_closed : LOGICAL;
        # v_closed : LOGICAL;
        # self_intersect : LOGICAL;
        # DERIVE
        # u_upper : INTEGER := SIZEOF(control_points_list) - 1;
        # v_upper : INTEGER := SIZEOF(control_points_list[1]) - 1;
        # control_points : ARRAY [0:u_upper] OF ARRAY [0:v_upper] OF
        # cartesian_point := make_array_of_array(control_points_list,0,u_upper,0,v_upper);
        # WHERE
        # WR1: (’GEOMETRY_SCHEMA.UNIFORM_SURFACE’ IN TYPEOF(SELF)) OR
        # (’GEOMETRY_SCHEMA.QUASI_UNIFORM_SURFACE’ IN TYPEOF(SELF)) OR
        # (’GEOMETRY_SCHEMA.BEZIER_SURFACE’ IN TYPEOF(SELF)) OR
        # (’GEOMETRY_SCHEMA.B_SPLINE_SURFACE_WITH_KNOTS’ IN TYPEOF(SELF));
        # get cartesian point references for eventual spline calculation
        BsplineSurfaceRefs = complexEntity[fieldRef].params

        ParsedSurface['surfaceUdegree'] = BsplineSurfaceRefs[0]
        ParsedSurface['surfaceVdegree'] = BsplineSurfaceRefs[1]
        ParsedSurface['controlPointsRefs'] = BsplineSurfaceRefs[2]

        # ParsedSurface['u_upper'] = len(ParsedSurface['controlPointsRefs']) - 1
        # ParsedSurface['v_upper'] = len(ParsedSurface['controlPointsRefs'][0]) - 1

        # controPointsListRefs structured as U control points length list of V control points length list

        controlPointsLenU = len(ParsedSurface['controlPointsRefs'])
        controlPointsLenV = len(ParsedSurface['controlPointsRefs'][0])

        # control_points_list : LIST [2:?] OF LIST [2:?] OF cartesian_point;
        controlPoints = []
        for cplU in ParsedSurface['controlPointsRefs']:
            for cplV in cplU:
                if (STEPobj[ref2index(cplV)].type_name == 'CARTESIAN_POINT'):
                    controlPoints.append(CP2point(cleanSubRefs(cplV)[0]))

        ParsedSurface['controlPoints'] = controlPoints
        ParsedSurface['surfaceForm'] = BsplineSurfaceRefs[3]
        ParsedSurface['closedU'] = 'T' in BsplineSurfaceRefs[4]
        ParsedSurface['closedV'] = 'T' in BsplineSurfaceRefs[5]
        ParsedSurface['selfIntersect'] = BsplineSurfaceRefs[6] # '.U.' ???

    if 'B_SPLINE_SURFACE_WITH_KNOTS' in typeNameList:
        fieldRef = typeNameList.index('B_SPLINE_SURFACE_WITH_KNOTS')

        BsplineSurfaceWithKnotsRefs = complexEntity[fieldRef].params
        ParsedSurface['knotUmultiplicities'] = BsplineSurfaceWithKnotsRefs[0]
        ParsedSurface['knotVmultiplicities'] = BsplineSurfaceWithKnotsRefs[1]
        ParsedSurface['knotsU'] = BsplineSurfaceWithKnotsRefs[2]
        ParsedSurface['knotsV'] = BsplineSurfaceWithKnotsRefs[3]
        ParsedSurface['knotSpec'] = BsplineSurfaceWithKnotsRefs[4]

    if 'RATIONAL_B_SPLINE_SURFACE' in typeNameList:
        fieldRef = typeNameList.index('RATIONAL_B_SPLINE_SURFACE')

        ParsedSurface['weightSpec'] = complexEntity[fieldRef].params[0]
        surfaceWeights = [item for sublist in ParsedSurface['weightSpec'] for item in sublist]

    # construct full U & V knotvectors
    STEPknotUvector = []
    for kml in range(len(ParsedSurface['knotUmultiplicities'])):
        for i in range(ParsedSurface['knotUmultiplicities'][kml]):
            STEPknotUvector.append(ParsedSurface['knotsU'][kml])
    ParsedSurface['STEPknotUvector'] = STEPknotUvector

    if (len(STEPknotUvector) - controlPointsLenU - ParsedSurface['surfaceUdegree']) != 1:
        print("maldefined B-spline surface (U) !")

    # p.106 https://quaoar.su/files/standards/Standard%20-%202003%20-%20ISO%2010303-42.pdf

    STEPknotVvector = []
    for kml in range(len(ParsedSurface['knotVmultiplicities'])):
        for i in range(ParsedSurface['knotVmultiplicities'][kml]):
            STEPknotVvector.append(ParsedSurface['knotsV'][kml])
    ParsedSurface['STEPknotVvector'] = STEPknotVvector

    # note _knotXvector values not used as Piegel+Tiller algorithms don't always work with explicit knots

    if (len(STEPknotVvector) - controlPointsLenV - ParsedSurface['surfaceVdegree']) != 1:
        print("maldefined B-spline surface (V) !")

    # BsplineKnotSurface = BSpline.Surface(normalize_kv=True)
    BsplineKnotSurface = NURBS.Surface()
    BsplineKnotSurface.degree_u = ParsedSurface['surfaceUdegree']
    BsplineKnotSurface.degree_v = ParsedSurface['surfaceVdegree']

    BsplineKnotSurface.ctrlpts_size_u = controlPointsLenU
    BsplineKnotSurface.ctrlpts_size_v = controlPointsLenV

    ctrlptsw = compatibility.combine_ctrlpts_weights(controlPoints, weights=surfaceWeights)
    BsplineKnotSurface.set_ctrlpts(ctrlptsw, controlPointsLenU, controlPointsLenV)

    # BsplineKnotSurface.ctrlpts2d = controlPoints
    # BsplineKnotSurface.knotvector_u = STEPknotUvector
    # BsplineKnotSurface.knotvector_v = STEPknotVvector

    BsplineKnotSurface.knotvector_u = utilities.generate_knot_vector(ParsedSurface['surfaceUdegree'],
                                                                     controlPointsLenU)
    BsplineKnotSurface.knotvector_v = utilities.generate_knot_vector(ParsedSurface['surfaceVdegree'],
                                                                     controlPointsLenV)

    ParsedSurface['BsplineKnotCurve'] = BsplineKnotSurface

    # No compelling reason to set evaluation delta & evaluate surface points as evaluation seems incorrect
    # BsplineKnotSurface.delta = 0.025  # this seems to be the minima delta under adaptive tesselation
    # BsplineKnotSurface.evaluate()

    if not calcExtrema:
        AFSobj['ParsedSurface'] = ParsedSurface
        return

    # local minima/maxima
    maxPointsUV = rationalSurfaceExtremaParam_5(BsplineKnotSurface,
                                                c,
                                                maxSearch=True,
                                                localExtrema=True,
                                                curvatureTest=False,
                                                uv_xyz=True)

    maxPoints = BsplineKnotSurface.evaluate_list(maxPointsUV)

    minPointsUV = rationalSurfaceExtremaParam_5(BsplineKnotSurface,
                                                c,
                                                maxSearch=False,
                                                localExtrema=True,
                                                curvatureTest=False,
                                                uv_xyz=True)

    minPoints = BsplineKnotSurface.evaluate_list(minPointsUV)

    # normalise u,v to [0, 1] - not required with geomdl generated knot vectors
    # knotUrange = STEPknotUvector[-1] - STEPknotUvector[0]
    # knotVrange = STEPknotVvector[-1] - STEPknotVvector[0]
    # maxPointsUV = [((mpu[0] - STEPknotUvector[0])/knotUrange, (mpu[1] - STEPknotVvector[0])/knotVrange ) for mpu in maxPointsUV]
    # minPointsUV = [((mpu[0] - STEPknotUvector[0])/knotUrange, (mpu[1] - STEPknotVvector[0])/knotVrange ) for mpu in minPointsUV]

    maxima = [True, ] * len(maxPoints) + [False, ] * len(minPoints)
    minima = [False, ] * len(maxPoints) + [True, ] * len(minPoints)

    maximaUV = np.array(maxPointsUV + minPointsUV)
    maximaPoints = np.array(maxPoints + minPoints)
    if len(maximaPoints) > 1:
        # create explicit dtype fields to permit sorting via u, then v
        maximaUV_field = maximaUV.ravel().view(dtype=[('u', maximaUV.dtype), ('v', maximaUV.dtype)])
        extremaUVindex = np.argsort(maximaUV_field, order=('u', 'v'))
        maximaUV = [maximaUV[euvi] for euvi in extremaUVindex]
        maximaPoints = [maximaPoints[euvi] for euvi in extremaUVindex]
        maxima = [maxima[euvi] for euvi in extremaUVindex]
        minima = [minima[euvi] for euvi in extremaUVindex]

        # because Newton-Raphson will converge on a tangent, separate max/min values may be identical
        maxima_truth = [np.allclose(maximaUV[mu + 1], maximaUV[mu], eps_STEP_AP21) for mu in
                        range(0, len(maximaUV) - 1)]
        maxima_truth = maxima_truth + [False, ]
        maximaUV = [mu for (mu, mt) in zip(maximaUV, maxima_truth) if not mt]
        maximaPoints = [mp for (mp, mt) in zip(maximaPoints, maxima_truth) if not mt]

    ParsedSurface['pointFeature']['xyz'] = [mp.tolist() for mp in maximaPoints]
    ParsedSurface['pointFeature']['uv'] = [muv.tolist() for muv in maximaUV]
    ParsedSurface['pointFeature']['centroidDisp'] = [np.linalg.norm(mp - c) for mp in maximaPoints]
    ParsedSurface['pointFeature']['maxima'] = maxima
    ParsedSurface['pointFeature']['minima'] = minima

    # detect rotationally symmetric spline surfaces, spline surface cylinder
    # get the median point of rings/arcs of control points
    spinePoints = np.array([np.array(s).mean(axis=0) for s in BsplineKnotSurface.ctrlpts2d])

    # test for deviation from mean point
    rowDelta = []
    for s in range(0, len(spinePoints) - 1):
        ctrlptsRowDisp = [np.linalg.norm(spinePoints[s] - sc) for sc in BsplineKnotSurface.ctrlpts2d[s]]
        rowDelta.append(max(ctrlptsRowDisp) - min(ctrlptsRowDisp))

    # low rowDelta values suggest rotational symmetry,
    # similar rowDelta values suggest consistent cross-section (with possible helicity):
    # consistentCrossSection = np.array(rowDelta).var()
    # should show up with edge maxima/minima

    # spinePointDisps = [np.linalg.norm(spinePoints[s] - BsplineKnotSurface.ctrlpts2d[s]) for s in range(0, len(BsplineKnotSurface.ctrlpts2d))]
    spinePointsMean = spinePoints.mean(axis=0)
    uu, dd, vv = np.linalg.svd(spinePoints - spinePointsMean)

    # vv[0] contains the first principal component, i.e. the direction
    # vector of the 'best fit' line in the least squares sense.
    # project mean control point rows to derived PCA spine and assess spine axis straightness
    spineDir = vv[0] - spinePointsMean
    spinePointDir = [s - spinePointsMean for s in spinePoints]
    spineIP = np.dot(spineDir, spineDir)
    projSpinePointDir = [np.dot(s, spineDir) / spineIP for s in spinePointDir]

    projSpinePoint = [spinePointsMean + pspd * spineDir for pspd in projSpinePointDir]
    # test for centroid disp from spine axis, giving minima rotation

    spineDeviations = [np.linalg.norm(projSpinePoint[n] - spinePoints[n]) for n in range(0, len(spinePoints - 1))]

    if (max(rowDelta) < eps_STEP_AP21) and (
            max(spineDeviations) < eps_STEP_AP21):  # rotSymLimit constant --------------------------------------------
        # ParsedSurface['rotSymFlag'] = 1
        ParsedSurface['rotSymFeature'] = dict()
        # ParsedSurface['rotSymFeature']['rotSymCentre'] = ?? centroid minima within
        ParsedSurface['rotSymFeature']['rotSymRadius'] = spinePointsMean

    AFSobj['ParsedSurface'] = ParsedSurface


def trimmedCurveParse(ParsedEdge, calcExtrema=False):
    pass
    # )
    # ENTITY trimmed_curve
    # SUBTYPE OF (bounded_curve);
    # basis_curve : curve;
    # trim_1 : SET[1:2] OF trimming_select;
    # trim_2 : SET[1:2] OF trimming_select;
    # sense_agreement : BOOLEAN;
    # master_representation : trimming_preference;
    # WHERE
    # WR1: (HIINDEX(trim_1) = 1) OR (TYPEOF(trim_1[1]) <> TYPEOF(trim_1[2]));
    # WR2: (HIINDEX(trim_2) = 1) OR (TYPEOF(trim_2[1]) <> TYPEOF(trim_2[2]));

    # ['#1951', [180.0], [360.0], '.T.', '.UNSPECIFIED.']
    # first value: reference to curve
    # TRIMMED_CURVE('',#17,(#22,PARAMETER_VALUE(0.)),(#23, PARAMETER_VALUE(6.28318530718)),.T.,.PARAMETER.)
    # TRIMMED_CURVE('',#17,(#96,PARAMETER_VALUE(1.570796326795)),(#97, PARAMETER_VALUE(7.853981633974)),.T.,.PARAMETER.)
    # TRIMMED_CURVE('',#89,(#93,PARAMETER_VALUE(0.)),(#94, PARAMETER_VALUE(10.)), .T.,.PARAMETER.)
    # TRIMMED_CURVE(#2197,(0.0),(1.0),.T.,.UNSPECIFIED.);

    # :
    # basis_curve: The curve to be trimmed. For curves with multiple representations any parameter values
    # given as trim_1 or trim_2 refer to the master representation of the basis_curve only.
    # trim_1: The first trimming point which may be specified as a cartesian point (point_1), as a real parameter value (parameter_1 = t1), or both.
    # c ISO 2003 — All rights reserved 69
    # ISO 10303-42:2003
    # trim_2: The second trimming point which may be specified as a cartesian point (point_2), as a real
    # parameter value (parameter_2 = t2 ), or both.
    # sense_agreement: Flag to indicate whether the direction of the trimmed curve agrees with or is opposed
    # to the direction of basis_curve.
    # — sense agreement = TRUE if the curve is being traversed in the direction of increasing parametric
    # value;
    # — sense agreement = FALSE otherwise. For an open curve, sense agreement = FALSE if t1 > t2. If
    # t2 > t1, sense agreement = TRUE. The sense information is redundant in this case but is essential
    # for a closed curve.
    # master_representation: Where both parameter and point are present at either end of the curve this
    # indicates the preferred form. Multiple representations provide the ability to communicate data in more
    # than one form, even though the data are expected to be geometrically identical. (See 4.3.9.)
    # NOTE 3 The master_representation attribute acknowledges the impracticality of ensuring that multiple forms
    # are indeed identical and allows the indication of a preferred form. This would probably be determined by the
    # creator of the data. All characteristics, such as parametrisation, domain, and results of evaluation, for an entity
    # having multiple representations, are derived from the master representation. Any use of the other representations
    # is a compromise for practical considerations.


def BsplineCurveWithKnotsParse(c, ParsedEdge, calcExtrema=False):
    '''
    Parse STEP BsplineCurveWithKnotsParse representation

    :param ParsedEdge: : parsed edge data object
    :param c: estimated surface median point based on extrema points identified
    :return
    '''

    # STEP ISO 10303-42 notes, EXPRESS descriptor
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

    # if ParsedEdge.get('subEdgeParams') is not None:
    #     edgeParams_ = ParsedEdge['subEdgeParams']
    # else:

    if ParsedEdge.get('knots') is None: # parse from STEP data

        edgeParams = ParsedEdge['edgeParams']

        if type(edgeParams[0]) == str:  # check for name string
            offset = 1
        else:
            offset = 0
        curveDegree = edgeParams[offset]
        controlPointsRefs = edgeParams[offset + 1]
        curveForm = edgeParams[offset + 2]
        closedCurve = 'T' in edgeParams[offset + 3]
        selfIntersect = 'T' in edgeParams[offset + 4]
        knotMultiplicities = edgeParams[offset + 5]
        knots = edgeParams[offset + 6]
        knotSpec = edgeParams[offset + 7]

        ParsedEdge['curveDegree'] = curveDegree
        ParsedEdge['controlPointsRefs'] = controlPointsRefs
        ParsedEdge['curveForm'] = curveForm
        ParsedEdge['closedCurve'] = closedCurve
        ParsedEdge['selfIntersect'] = selfIntersect
        ParsedEdge['knotMultiplicities'] = knotMultiplicities
        ParsedEdge['knots'] = knots
        ParsedEdge['knotSpec'] = knotSpec

        # extract control points
        controlPoints = []
        for cpl in controlPointsRefs:
            if (STEPobj[ref2index(cpl)].type_name == 'CARTESIAN_POINT'):
                controlPoint = CP2point(cleanSubRefs(cpl)[0])
                controlPoints.append(controlPoint)
        ParsedEdge['controlPoints'] = controlPoints

        knotvector = []
        for kml in range(len(knotMultiplicities)):
            for i in range(knotMultiplicities[kml]):
                knotvector.append(knots[kml])
        ParsedEdge['STEPknotvector'] = knotvector

        if (len(knotvector) - len(controlPoints) - curveDegree) != 1:
            print("maldefined B-spline!")

        #  rule for closed periodic curves is :
        #  number of input control points + curve degree = the total number of control points required.
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

        # BsplineKnotCurve = BSpline.Curve(normalize_kv=True)
        # BsplineKnotCurve.degree = curveDegree
        # BsplineKnotCurve.ctrlpts = [list(c) for c in controlPoints]
        # BsplineKnotCurve.knotvector = _knotvector

        BsplineKnotCurve = NURBS.Curve()
        BsplineKnotCurve.degree = curveDegree
        BsplineKnotCurve.set_ctrlpts(compatibility.combine_ctrlpts_weights(controlPoints, weights=None))

        # Auto-generate knot vector
        BsplineKnotCurve.knotvector = utilities.generate_knot_vector(BsplineKnotCurve.degree, len(BsplineKnotCurve.ctrlpts))
        # BsplineKnotCurve.evaluate()

        ParsedEdge['BsplineKnotCurve'] = BsplineKnotCurve

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

    if calcExtrema:

        maxPointsU = BsplineCurveExtremaDisp(ParsedEdge['BsplineKnotCurve'],
                                             c,
                                             localExtrema=True,
                                             # curvatureTest=True,
                                             maxSearch=True,
                                             uv_xyz=True)

        #maxPoints = [np.array(cp) for cp in BsplineKnotCurve.evaluate_list(maxPointsU)]
        maxPoints = ParsedEdge['BsplineKnotCurve'].evaluate_list(maxPointsU)

        minPointsU = BsplineCurveExtremaDisp(ParsedEdge['BsplineKnotCurve'],
                                             c,
                                             localExtrema=True,
                                             # curvatureTest=True,
                                             maxSearch=False,
                                             uv_xyz=True)

        #minPoints = [np.array(cp) for cp in BsplineKnotCurve.evaluate_list(minPointsU)]
        minPoints = ParsedEdge['BsplineKnotCurve'].evaluate_list(minPointsU)

        # normalise u to [0, 1] - unrequired where normalize_kv=True
        # knotRange = _knotvector[-1] - _knotvector[0]
        # maxPointsU = [(mpu - _knotvector[0])/knotRange for mpu in maxPointsU]
        # minPointsU = [(mpu - _knotvector[0])/knotRange for mpu in minPointsU]

        # exception to provided vertices in revolved surface
        if ParsedEdge.get('vertex1') is None:
            ParsedEdge['vertex1'] = np.array(ParsedEdge['BsplineKnotCurve'].evaluate_single(0.))
        if ParsedEdge.get('vertex2') is None:
            ParsedEdge['vertex2'] = np.array(ParsedEdge['BsplineKnotCurve'].evaluate_single(1.))

        extremaU = maxPointsU + minPointsU
        maximaPoints = maxPoints + minPoints
        if len(maximaPoints) > 0:
            extremaUindex = np.argsort(extremaU)
            extremaU = [extremaU[eui] for eui in extremaUindex]
            maximaPoints = [maximaPoints[eui] for eui in extremaUindex]

            # because Newton-Raphson will converge on a tangent, separate max/min values may be identical
            maxima_truth = [(extremaU[mu + 1] - extremaU[mu]) > eps_STEP_AP21 for mu in range(0, len(extremaU) - 1)]
            maxima_truth = maxima_truth + [True, ]
            extremaU = [mu for (mu, mt) in zip(extremaU, maxima_truth) if mt]
            maximaPoints = [np.array(mp) for (mp, mt) in zip(maximaPoints, maxima_truth) if mt]

            # if vertex1 centroidDisp is close to u[0] or vertex2 centroidDisp is close to u[-1], remove u[0], u[-1]
            # this metric is less relevant to a curvature based extrema: which is not used here
            centroidDisp = [np.linalg.norm(mp - c) for mp in maximaPoints]
            vertex1centroidDisp = np.linalg.norm(ParsedEdge['vertex1'] - c)
            vertex2centroidDisp = np.linalg.norm(ParsedEdge['vertex2'] - c)

            if np.abs(vertex1centroidDisp - centroidDisp[0]) < eps_STEP_AP21:
                extremaU.pop(0)
                maximaPoints.pop(0)
            else:
                centroidDisp.insert(0, vertex1centroidDisp)
            if np.abs(vertex2centroidDisp - centroidDisp[-1]) < eps_STEP_AP21:
                extremaU.pop(-1)
                maximaPoints.pop(-1)
            else:
                centroidDisp.append(vertex2centroidDisp)
        else:
            centroidDisp = []

        # maximaPoints.insert(0, ParsedEdge['vertex1']) # ParsedEdge['pointFeature']['xyz'][0])
        # maximaPoints.insert(-1, ParsedEdge['vertex2']) # ParsedEdge['pointFeature']['xyz'][-1])
        # ParsedEdge['pointFeature']['xyz'] = np.array([ParsedEdge['vertex1'], ] + maximaPoints + [ParsedEdge['vertex2'], ])
        ParsedEdge['pointFeature']['xyz'] = [ParsedEdge['vertex1'], ] + maximaPoints + [ParsedEdge['vertex2'], ]
        ParsedEdge['pointFeature']['u'] = [0., ] + [mu.item() for mu in extremaU] + [1.,]
        ParsedEdge['pointFeature']['centroidDisp'] = [cd.item() for cd in centroidDisp] # .. is this really necessary?

        # no point adding maxima & minima fields here

        # if len(maximaPoints) > 1:
        #     ParsedEdge['maximaPoints'] = maximaPoints
        #     maximaPointsCentroidDisp = [np.linalg.norm(np.array(e) - c) for e in maximaPoints]
        #     ParsedEdge['maximaPointsCentroidDisp'] = maximaPointsCentroidDisp
        #     maxPoint = maximaPoints[maximaPointsCentroidDisp.index(max(maximaPointsCentroidDisp))]
        #     ParsedEdge['maxPoint'] = maxPoint
        #     ParsedEdge['maxPointCentroidDisp'] = np.linalg.norm(maxPoint - c)
        # if len(maximaPoints) == 1:
        #     ParsedEdge['maxPoint'] = maximaPoints[0]
        #     ParsedEdge['maxPointCentroidDisp'] = np.linalg.norm(maxPoint - c)
        #     ParsedEdge['maximaPointsCentroidDisp'] = None
        #     ParsedEdge['maximaPoints'] = None
        # else:
        #     ParsedEdge['maximaPointsCentroidDisp'] = None
        #     ParsedEdge['maximaPoints'] = None
        #     ParsedEdge['maxPointCentroidDisp'] = None
        #     ParsedEdge['maxPoint'] = None

        # if len(minimaPoints) > 1:
        #     ParsedEdge['minimaPoints'] = minimaPoints
        #     minimaPointsCentroidDisp = [np.linalg.norm(np.array(e) - c) for e in minimaPoints]
        #     ParsedEdge['minimaPointsCentroidDisp'] = minimaPointsCentroidDisp
        #     minPoint = minimaPoints[minimaPointsCentroidDisp.index(min(minimaPointsCentroidDisp))]
        #     ParsedEdge['maxPoint'] = maxPoint
        #     ParsedEdge['maxPointCentroidDisp'] = np.linalg.norm(minPoint - c)
        # if len(minimaPoints) == 1:
        #     ParsedEdge['minPoint'] = minimaPoints[0]
        #     ParsedEdge['minPointCentroidDisp'] = np.linalg.norm(minPoint - c)
        #     ParsedEdge['minimaPointsCentroidDisp'] = None
        #     ParsedEdge['minimaPoints'] = None
        # else:
        #     ParsedEdge['minimaPointsCentroidDisp'] = None
        #     ParsedEdge['maximaPoints'] = None
        #     ParsedEdge['minPointCentroidDisp'] = None
        #     ParsedEdge['minPoint'] = None

        # minPointCentroidDisp = np.linalg.norm(minPoint - c)
        # ParsedEdge['minPointCentroidDisp'] = minPointCentroidDisp
        # ParsedEdge['minPoint'] = minPoint
        # maxPointCentroidDisp = np.linalg.norm(maxPoint - c)
        # ParsedEdge['maxPointCentroidDisp'] = maxPointCentroidDisp
        # ParsedEdge['maxPoint'] = maxPoint

        # v1 = ParsedEdge['vertex1']
        # v2 = ParsedEdge['vertex2']

        # # independent check for vertices,
        # ParsedEdge['vertex1extremaMin'] = False
        # ParsedEdge['vertex2extremaMin'] = False
        # ParsedEdge['vertex1extremaMax'] = False
        # ParsedEdge['vertex2extremaMax'] = False
        #
        # if ParsedEdge['minPointCentroidDisp'] is not None:
        #     if ParsedEdge['vertex1centroidDisp'] < ParsedEdge['minPointCentroidDisp']:
        #         ParsedEdge['minPoint'] = ParsedEdge['vertex1']
        #         ParsedEdge['minPointCentroidDisp'] = ParsedEdge['vertex1centroidDisp']
        #         ParsedEdge['vertex1extremaMin'] = True
        #
        # if ParsedEdge['minPointCentroidDisp'] is not None:
        #     if ParsedEdge['vertex2centroidDisp'] < ParsedEdge['minPointCentroidDisp']:
        #         ParsedEdge['minPoint'] = ParsedEdge['vertex2']
        #         ParsedEdge['minPointCentroidDisp'] = ParsedEdge['vertex2centroidDisp']
        #         ParsedEdge['vertex2extremaMin'] = True
        #
        # if ParsedEdge['maxPointCentroidDisp'] is not None:
        #     if ParsedEdge['vertex1centroidDisp'] > ParsedEdge['maxPointCentroidDisp']:
        #         ParsedEdge['maxPoint'] = ParsedEdge['vertex1']
        #         ParsedEdge['maxPointCentroidDisp'] = ParsedEdge['vertex1centroidDisp']
        #         ParsedEdge['vertex1extremaMax'] = True
        #
        # if ParsedEdge['maxPointCentroidDisp'] is not None:
        #     if ParsedEdge['vertex2centroidDisp'] > ParsedEdge['maxPointCentroidDisp']:
        #         ParsedEdge['maxPoint'] = ParsedEdge['vertex2']
        #         ParsedEdge['maxPointCentroidDisp'] = ParsedEdge['vertex2centroidDisp']
        #         ParsedEdge['vertex2extremaMax'] = True

        # #BsplineKnotCurve = BSpline.Curve(normalize_kv=False)
        # BsplineKnotCurve = NURBS.Curve()
        # BsplineKnotCurve.degree = curveDegree
        # #BsplineKnotCurve.ctrlpts = [list(c) for c in controlPoints]
        # #BsplineKnotCurve.knotvector = STEPknotvector
        # BsplineKnotCurve.set_ctrlpts(compatibility.combine_ctrlpts_weights(controlPoints, weights=None))

        # detect rotationally symmetric spline surfaces, spline surface cylinder
        # get the mean point of rings/arcs of control points
        splineMeanPoint = np.array(ParsedEdge['BsplineKnotCurve'].ctrlpts).mean(axis=0)

        # test for deviation from mean point
        ctrlptsDelta = [np.linalg.norm(splineMeanPoint - sc) for sc in ParsedEdge['BsplineKnotCurve'].ctrlpts]

        # low ctrlptsDelta values suggest rotational symmetry,
        if (max(np.abs(ctrlptsDelta)) < eps_STEP_AP21):  # rotSymLimit constant --------------------------------------------
            ParsedEdge['rotSymFlag'] = 1

        #ParsedEdge['BsplineKnotCurve'] = BsplineKnotCurve


def edgeSTEPparse(edgeInstance_):
    '''
    Parse common STEP representation

    :param edgeInstance_:  parsed common STEP data
    :return
    '''

    # STEP ISO 10303-42 notes
    # axis: the Direction that defines the second axis of the Axis_placement. The value of this attribute need not be specified.
    # ref_direction: the direction used to determine the direction of the local X axis.
    # The value of this attribute need not be specified. If axis or ref_direction is omitted, these directions are taken from the geometric coordinate system
    # orientation: a BOOLEAN flag. If TRUE, the topological orientation as used coincides with the orientation,
    # from start vertex to end vertex, of the edge_definition, if FALSE the vertices are reversed in order.
    # edge_start: the start vertex of the Oriented_edge. This is derived from the vertices of the edge_definition after taking account of the orientation.
    # edge_end: the end vertex of the Oriented_edge. This is derived from the vertices of the edge_definition after taking account of the orientation.

    # extract data common to all STEP edge instances
    ParsedEdge = {}
    if (STEPobj[ref2index(edgeInstance_)].type_name == 'ORIENTED_EDGE'):
        # if not hasattr(STEPobj[ref2index(edgeLoopInstance)].type_name, 'ORIENTED_EDGE'): ----------------------
        VOR = 'F' in STEPobj[ref2index(edgeInstance_)].params[-1]  # vertexOrderReversed
        ParsedEdge['VOR'] = VOR

    if (STEPobj[ref2index(edgeInstance_)].type_name == 'TRIMMED_CURVE'):
        trimmedCurveParamPrecedence = STEPobj[ref2index(edgeInstance_)].params[-1]
        VOR = 'F' in STEPobj[ref2index(edgeInstance_)].params[-2]  # vertexOrderReversed
        trimmedCurveParams = cleanSubRefs(STEPobj[ref2index(edgeInstance_)].params)
        if len(trimmedCurveParams) == 3: # curve ref & both vertex refs available
            edgeCurveRef = trimmedCurveParams[0]
            if (STEPobj[ref2index(trimmedCurveParams[1])].type_name == 'VERTEX_POINT'):
                vertex1 = STEPobj[ref2index(trimmedCurveParams[1])].params
                vertex1 = CP2point(cleanSubRefs(vertex1)[0])
        # trimmedCurveRef = STEPobj[ref2index(edgeInstance_)].params[-2]
        # assert STEPobj[ref2index(edgeCurveRef)].type_name == 'TRIMMED_CURVE'
        # # edge_geometry: the Curve defining the geometric shape of the edge.
        # # same_sense: a BOOLEAN variable giving the relationship between the topological sense
        # #             of the edge and the parametric sense of the curve.
        # trimmedCurveParams = STEPobj[ref2index(trimmedCurveRef)].params
        # # edgeSameSense = 'T' in edgeCurveParams[-1]; ParsedEdge['sameSense'] = edgeSameSense
        # trimmedCurveParams = cleanSubRefs(trimmedCurveParams)

    edgeCurveRef = STEPobj[ref2index(edgeInstance_)].params[-2]
    assert STEPobj[ref2index(edgeCurveRef)].type_name == 'EDGE_CURVE'
    # edge_geometry: the Curve defining the geometric shape of the edge.
    # same_sense: a BOOLEAN variable giving the relationship between the topological sense
    #             of the edge and the parametric sense of the curve.
    edgeCurveParams = STEPobj[ref2index(edgeCurveRef)].params
    # edgeSameSense = 'T' in edgeCurveParams[-1]; ParsedEdge['sameSense'] = edgeSameSense
    edgeCurveParams = cleanSubRefs(edgeCurveParams)

    if (STEPobj[ref2index(edgeCurveParams[0])].type_name == 'VERTEX_POINT'):
        vertex1 = STEPobj[ref2index(edgeCurveParams[0])].params
        vertex1 = CP2point(cleanSubRefs(vertex1)[0])
        if not VOR:
            ParsedEdge['vertex1'] = np.array(vertex1)
            ParsedEdge['vertex1ref'] = STEPobj[ref2index(edgeCurveParams[0])].ref
            # vertex1centroidDisp = np.linalg.norm(centroid - vertex1)
            # ParsedEdge['vertex1centroidDisp'] = vertex1centroidDisp
        else:
            ParsedEdge['vertex2'] = np.array(vertex1)
            ParsedEdge['vertex2ref'] = STEPobj[ref2index(edgeCurveParams[0])].ref
            # vertex2centroidDisp = np.linalg.norm(centroid - vertex1)
            # ParsedEdge['vertex2centroidDisp'] = vertex2centroidDisp

    if (STEPobj[ref2index(edgeCurveParams[1])].type_name == 'VERTEX_POINT'):
        vertex2 = STEPobj[ref2index(edgeCurveParams[1])].params
        vertex2 = CP2point(cleanSubRefs(vertex2)[0])
        if not VOR:
            ParsedEdge['vertex2'] = np.array(vertex2)
            ParsedEdge['vertex2ref'] = STEPobj[ref2index(edgeCurveParams[1])].ref
            # vertex2centroidDisp = np.linalg.norm(centroid - vertex2)
            # ParsedEdge['vertex2centroidDisp'] = vertex2centroidDisp
        else:
            ParsedEdge['vertex1'] = np.array(vertex2)
            ParsedEdge['vertex1ref'] = STEPobj[ref2index(edgeCurveParams[1])].ref
            # vertex2centroidDisp = np.linalg.norm(centroid - vertex2)
            # ParsedEdge['vertex1centroidDisp'] = vertex2centroidDisp

    # must be recalculated on every centroid update
    ParsedEdge['pointFeature'] = {}
    # ParsedEdge['pointFeature']['xyz'] = [ParsedEdge['vertex1'], ParsedEdge['vertex2']]
    # ParsedEdge['pointFeature']['centroidDisp'] = [np.linalg.norm(centroid - vertex1), np.linalg.norm(centroid - vertex2)]
    # ParsedEdge['pointFeature']['u'] = [0, 1]  # normalise to [0, 1]

    ParsedEdge['typeName'] = STEPobj[ref2index(edgeCurveParams[2])].type_name
    ParsedEdge['edgeRef'] = STEPobj[ref2index(edgeCurveParams[2])].ref
    ParsedEdge['edgeParams'] = STEPobj[ref2index(edgeCurveParams[2])].params

    if (ParsedEdge['typeName'] == 'SURFACE_CURVE') or (ParsedEdge['typeName'] == 'SEAM_CURVE'):

        #   SEAM_CURVE (not sure if this is worthwhile retaining)
        #   associated geometry in following list does not seem to be defined elsewhere

        #   Attribute	                Type	                                        Defined By
        #   name	                    label (STRING)	                                representation_item
        #   curve_3d	                curve (ENTITY)	                                surface_curve
        #   associated_geometry	        LIST OF pcurve_or_surface (SELECT)	            surface_curve
        #   master_representation	    preferred_surface_curve_representation (ENUM)	surface_curve

        # todo: just create a superTypeName flag?

        cleanEdgeParams = cleanSubRefs(ParsedEdge['edgeParams'])
        #edgeParamsName = STEPobj[ref2index(cleanEdgeParams[0])].type_name
        ParsedEdge['superTypeName'] = ParsedEdge['typeName']

        edgeRef = STEPobj[ref2index(cleanEdgeParams[0])].ref

        ParsedEdge['superEdgeRef'] = ParsedEdge['edgeRef']
        ParsedEdge['edgeRef'] = edgeRef
        ParsedEdge['superEdgeParams'] = ParsedEdge['edgeParams']
        ParsedEdge['edgeParams'] = STEPobj[ref2index(cleanEdgeParams[0])].params  # todo: change edgeParams
        ParsedEdge['superTypeName'] = ParsedEdge['typeName']
        ParsedEdge['typeName'] = STEPobj[ref2index(edgeRef)].type_name
        # if (ParsedEdge['typeName'] != 'SURFACE_CURVE') and (ParsedEdge['typeName'] != 'SEAM_CURVE'):
        #     edgeParse(centroid, edgeRef)

    return ParsedEdge

    # a proper object initiation would create extra terms to avoid messy ".get('x') is None" clauses
    # while scope is unsettled this is parked here for edge object, presumably surface requires the same

    # 'subTypeName'
    # 'subedgeParams'
    # 'normDir'
    # 'radius'
    # 'axisPoint'
    #
    # 'rotSymRadius'
    # 'rotSymCentre'
    # 'rotSymMin'
    # 'rotSymMax'


def edgeParse(c, ParsedEdge, calcExtrema=False):
    '''
    Process edge parsing according to STEP edge type

    :param ParsedEdge:  parsed STEP edge data
    :return
    '''
    # print("simpleEdgeParse edge type: " + ParsedEdge['typeName'])

    if (ParsedEdge['typeName'] == 'LINE'):
        segmentParse(c, ParsedEdge, calcExtrema)

    elif (ParsedEdge['typeName'] == 'TRIMMED_CURVE'): # identified here?
        #
        _1=1
        print("TRIMMED_CURVE - write code!")

    elif (ParsedEdge['typeName'] == 'SURFACE_CURVE') or (ParsedEdge['typeName'] == 'SEAM_CURVE'):

        #   SEAM_CURVE (not sure if this is worthwhile retaining)
        #   associated geometry in following list does not seem to be defined elsewhere

        #   Attribute	                Type	                                        Defined By
        #   name	                    label (STRING)	                                representation_item
        #   curve_3d	                curve (ENTITY)	                                surface_curve
        #   associated_geometry	        LIST OF pcurve_or_surface (SELECT)	            surface_curve
        #   master_representation	    preferred_surface_curve_representation (ENUM)	surface_curve

        # todo: just create a superTypeName flag?

        cleanEdgeParams = cleanSubRefs(ParsedEdge['edgeParams'])
        #edgeParamsName = STEPobj[ref2index(cleanEdgeParams[0])].type_name
        ParsedEdge['superTypeName'] = ParsedEdge['typeName']

        edgeRef = STEPobj[ref2index(cleanEdgeParams[0])].ref

        ParsedEdge['superEdgeRef'] = ParsedEdge['edgeRef']
        ParsedEdge['edgeRef'] = edgeRef
        ParsedEdge['superEdgeParams'] = ParsedEdge['edgeParams']
        ParsedEdge['edgeParams'] = STEPobj[ref2index(cleanEdgeParams[0])].params  # todo: change edgeParams
        ParsedEdge['superTypeName'] = ParsedEdge['typeName']
        ParsedEdge['typeName'] = STEPobj[ref2index(edgeRef)].type_name
        # if (ParsedEdge['typeName'] != 'SURFACE_CURVE') and (ParsedEdge['typeName'] != 'SEAM_CURVE'):
        #     # try this as an actual curve/line, as non-vertex minima or maxima will be winnowed out
        #     edgeParse(centroid, edgeRef)

        if (ParsedEdge['typeName'] != 'SURFACE_CURVE') and (ParsedEdge['typeName'] != 'SEAM_CURVE'):
            _1=1

        # try this as an actual curve/line, as non-vertex minima or maxima will be winnowed out
        edgeParse(centroid, edgeRef)

        # if ParsedEdge['edgeRef'] not in [pel['edgeRef'] for pel in ParsedEdgeLoop]:
        #     ParsedEdgeLoop.append(ParsedEdge)

        for cet in cleanEdgeParams[1:]:  # get associated surfaces, LIST OF pcurve_or_surface
            if STEPobj[ref2index(cet)].type_name == 'PCURVE':
                # extract surface data from second parameter term
                surfaceRef = STEPobj[ref2index(cet)].params[1]

                # surfaceRefCollection = []
                # for AFS2 in AFSobj:
                #     surfaceRefCollection.append(AFS2['SurfaceRef'])
                # if surfaceRef not in surfaceRefCollection:
                if surfaceRef not in [afs2['SurfaceRef'] for afs2 in AdvancedFaceSurfaces]: # remove once PCURVE exception is ruled out
                    _1 = 1
                    print('undetected surface! PCURVE')
            else:
                # todo: code for other surface option missing
                print('undetermined surfacetype failure line2824')
                _1 = 1

    elif (ParsedEdge['typeName'] == 'CIRCLE'):
        circleParse(c, ParsedEdge, calcExtrema)

    elif (ParsedEdge['typeName'] == 'ELLIPSE'):
        ellipseParse(c, ParsedEdge, calcExtrema)

    elif (ParsedEdge['typeName'] == 'B_SPLINE_CURVE_WITH_KNOTS'):
        BsplineCurveWithKnotsParse(c, ParsedEdge, calcExtrema)

    # TRIMMED_CURVE ?

    else:
        print("simpleEdgeParse unknown edge type: " + ParsedEdge['typeName'])

    return ParsedEdge


def revolvedSurfaceParse(AFSobj, c, calcExtrema=False):
    '''
    Parse a revolved surface, optionally extract surface minima points relative to a local centroid point.

    :param AFSobj: parsed surface data object
    :param c: estimated surface median point based on extrema points identified
    :param calcExtrema: bool additional calculation to determine surface extrema
    :return ParsedSurface object
    '''

    # A surface_of_revolution is the surface obtained by rotating a curve one complete revolution about an axis.
    # Will assume that any curvature will be in the plane of the axis.
    # Reposition curve on both sides of axis in plane of axis and local centroid to find maxima/minima, see cone

    # SURFACE_OF_REVOLUTION
    # ENTITY surface_of_revolution
    # SUBTYPE OF (swept_surface);
    # axis_position : axis1_placement;
    # DERIVE
    # axis_line : line := representation_item(’’)||
    # geometric_representation_item()|| curve()||
    # line(axis_position.location, representation_item(’’)||
    # geometric_representation_item()||
    # vector(axis_position.z, 1.0));

    if AFSobj.get('ParsedSurface') is not None:  # recalculation of max/min
        ParsedSurface = AFSobj['ParsedSurface']
        axisPoint = ParsedSurface['axisPoint']
        normDir = ParsedSurface['normDir']  # Z
        revolvedCurveRef = ParsedSurface['revolvedCurveRef']

    else:

        axisPoint, normDir = axis1Placement(AFSobj['SurfaceParams'][-1])
        revolvedCurveRef = AFSobj['SurfaceParams'][-2]

        ParsedSurface = {}
        ParsedSurface['axisPoint'] = axisPoint
        ParsedSurface['normDir'] = normDir  # Z
        ParsedSurface['revolvedCurveRef'] = revolvedCurveRef

        ParsedSurface['pointFeature'] = {}
        ParsedSurface['pointFeature']['xyz'] = []
        ParsedSurface['pointFeature']['centroidDisp'] = []
        ParsedSurface['pointFeature']['u'] = [] #-------------------------------------edge?
        ParsedSurface['pointFeature']['maxima'] = []
        ParsedSurface['pointFeature']['minima'] = []

    AFSobj['ParsedSurface'] = ParsedSurface
    if not calcExtrema:
        # return ParsedSurface
        return

    SurfaceProfile = {}
    SurfaceProfile['edgeRef'] = STEPobj[ref2index(revolvedCurveRef)].ref
    SurfaceProfile['edgeParams'] = STEPobj[ref2index(revolvedCurveRef)].params
    SurfaceProfile['typeName'] = STEPobj[ref2index(revolvedCurveRef)].type_name

    SurfaceProfile['pointFeature'] = {}
    SurfaceProfile['pointFeature']['xyz'] = []
    SurfaceProfile['pointFeature']['centroidDisp'] = []
    SurfaceProfile['pointFeature']['maxima'] = []
    SurfaceProfile['pointFeature']['minima'] = []
    SurfaceProfile['pointFeature']['u'] = []

    edgeParse(centroid, SurfaceProfile)

    #centroid = np.array([10., 10., 0.]) #--------------------------------------------------------------------TEST

    # test if centroid is close to axisDir through axisPoint
    rotSymDisp = np.linalg.norm((centroid - axisPoint) - np.dot((centroid - axisPoint), normDir) * normDir)

    # Cdisp = np.dot((centroid - axisPoint), normDir)
    # Cproj = centroid - np.dot(Cdisp, normDir)  # point projected to plane along orthogonal
    # rotSymDisp2 = np.linalg.norm(Cproj - axisPoint)
    # print("rotSymDisp: " + str(rotSymDisp))
    # print("rotSymDisp2: " + str(rotSymDisp2))
    # print("offset centroid revolvedSurfaceParse untested, line 3869")

    if rotSymDisp < eps_STEP_AP21:
        # no requirement to rotate centroid, surface maxima & minima are returned point maxima & minima from edgeParse
        # but reassigned as rotSym feature
        ParsedSurface['rotSymFeature'] = dict()
        ParsedSurface['rotSymFeature']['rotSymCentre'] = []
        ParsedSurface['rotSymFeature']['rotSymRadius'] = []

        if len(SurfaceProfile['pointFeature']['xyz']) > 0:
             for rsf in range(0, len(SurfaceProfile['pointFeature']['xyz'])):
                 ParsedSurface['rotSymFeature']['rotSymCentre'].append(pointProjectAxis(SurfaceProfile['pointFeature']['xyz'][rsf], axisPoint, normDir))
                 ParsedSurface['rotSymFeature']['rotSymRadius'].append(SurfaceProfile['pointFeature']['centroidDisp'][rsf])

        # reset SurfaceProfile['pointFeature'] to avoid being reclassed as general point features
        # SurfaceProfile['pointFeature']['xyz'] = []
        # SurfaceProfile['pointFeature']['centroidDisp'] = []
        # SurfaceProfile['pointFeature']['u'] = []
        # provide 180 degree points see circleParse


    else:
        # rotate centroid point to plane of axis & curve rather than rotating curve,
        # find minima/maxima point on curve wrt centroid,
        # again reliant on assumption that curvature is in the plane of axis

        # project centroid and a vertex to plane defined by normDir and axisPoint
        projectCentroid = centroid - np.dot(centroid - axisPoint, normDir) * normDir
        projectVertex1 = SurfaceProfile['vertex1'] - np.dot(SurfaceProfile['vertex1'] - axisPoint, normDir) * normDir
        projectVertex2 = SurfaceProfile['vertex2'] - np.dot(SurfaceProfile['vertex2'] - axisPoint, normDir) * normDir
        if np.linalg.norm(projectVertex1 - axisPoint) >= np.linalg.norm(projectVertex2 - axisPoint):
            projectCurve = projectVertex1
        else:
            projectCurve = projectVertex2

        centroidRot = np.dot(projectCentroid, projectCurve)/(np.linalg.norm(projectCurve - axisPoint) * np.linalg.norm(projectCentroid - axisPoint))
        centroidRotM = rotationMatrix(normDir, np.arccos(centroidRot))

        # θ = arccos[(a · b) / (| a | | b |)].


def sphereSurfaceParse(AFSobj, c, calcExtrema=False):
    '''
    Parse a spherical surface, optionally extract surface minima points relative to a local centroid point.

    :param AFSobj: parsed surface data object
    :param c: estimated surface median point based on extrema points identified
    :param calcExtrema: bool additional calculation to determine surface extrema
    :return ParsedSurface object
    '''

    # STEP ISO 10303-42 notes, EXPRESS descriptor
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

    if AFSobj.get('ParsedSurface') is not None:  # recalculation of max/min
        ParsedSurface = AFSobj['ParsedSurface']
        centrePoint = ParsedSurface['axisPoint']
        normDir = ParsedSurface['normDir']  # Z
        refDir = ParsedSurface['refDir']  # X
        auxDir = ParsedSurface['auxDir']  # Y
    else:
        centrePoint, normDir, refDir = axis2Placement3D(AFSobj['SurfaceParams'][1])
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        radius = AFSobj['SurfaceParams'][-1]

        ParsedSurface = {}
        ParsedSurface['radius'] = radius
        ParsedSurface['axisPoint'] = centrePoint
        ParsedSurface['normDir'] = normDir  # Z
        ParsedSurface['refDir'] = refDir  # X
        ParsedSurface['auxDir'] = auxDir

        ParsedSurface['pointFeature'] = {}
        ParsedSurface['pointFeature']['xyz'] = []
        ParsedSurface['pointFeature']['centroidDisp'] = []
        ParsedSurface['pointFeature']['maxima'] = []
        ParsedSurface['pointFeature']['minima'] = []

    AFSobj['ParsedSurface'] = ParsedSurface
    if not calcExtrema:  # not much point
        # return ParsedSurface
        return

    # sphere surface maxima & minima exist on a vector between centroid and sphere centre point
    # minima/maxima = centre point -/+ radius

    # "The local z-axis corresponds to the normal of planar and spherical surfaces and to the axis of cylindrical,
    # conical and toroidal surfaces."
    # Jaider Oussama et al Int. Journal of Engineering Research and Applications www.ijera.com
    # ISSN : 2248-9622, Vol. 4, Issue 5( Version 6), May 2014, pp.14-25
    # refDir, auxDir doesn't mean much in stock implementaions.

    sphereCentroidOffset = np.linalg.norm(centrePoint - c)
    if sphereCentroidOffset < eps_STEP_AP21:  # determine if centroid is coincident with sphere centre
        ParsedSurface['sphereFeatureFlag'] = True

        ParsedSurface['rotSymFeature'] = dict()
        ParsedSurface['rotSymFeature']['rotSymCentre'] = [centrePoint,]
        ParsedSurface['rotSymFeature']['rotSymRadius'] = [radius,]

    else:
        sphereCentreCentroidDir = (centrePoint - c) / sphereCentroidOffset

        # determine if maxPoint/minPoint are inside or outside an edge_loop defined boundary
        # CCW directionality of edge vertices determines outward surface normal
        # surface minPoint/maxPoint always orthonormal to c

        # sphere sector has at most one maximum and one minimum surface point.

        maxPoint = centrePoint + sphereCentreCentroidDir * radius
        minPoint = centrePoint - sphereCentreCentroidDir * radius

        if insideOutsideSurfaceTest(c, maxPoint, AFSobj):  # maxima not among edges
            ParsedSurface['pointFeature']['xyz'].append(maxPoint)
            ParsedSurface['pointFeature']['centroidDisp'].append(np.linalg.norm(maxPoint - c))
            ParsedSurface['pointFeature']['maxima'].append(False)
            ParsedSurface['pointFeature']['minima'].append(True)

        if insideOutsideSurfaceTest(c, minPoint, AFSobj):  # maxima not among edges
            ParsedSurface['pointFeature']['xyz'].append(minPoint)
            ParsedSurface['pointFeature']['centroidDisp'].append(np.linalg.norm(minPoint - c))
            ParsedSurface['pointFeature']['maxima'].append(True)
            ParsedSurface['pointFeature']['minima'].append(False)

    AFSobj['ParsedSurface'] = ParsedSurface


def toroidSurfaceParse(AFSobj, c, calcExtrema=False):
    '''
    Parse a toroidal surface, optionally extract surface minima points relative to a local centroid point.

    :param AFSobj: parsed surface data object
    :param c: estimated surface median point based on extrema points identified
    :param calcExtrema: bool additional calculation to determine surface extrema
    :return ParsedSurface object
    '''

    # STEP ISO 10303-42 notes, EXPRESS descriptor
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

    if AFSobj.get('ParsedSurface') is not None:  # recalculation of max/min
        ParsedSurface = AFSobj['ParsedSurface']
        axisPoint = ParsedSurface['axisPoint']
        normDir = ParsedSurface['normDir']  # Z
        refDir = ParsedSurface['refDir']  # X
        auxDir = ParsedSurface['auxDir']  # Y

        minorRadius = ParsedSurface['minorRadius']
        majorRadius = ParsedSurface['majorRadius']
    else:
        ParsedSurface = {}
        axisPoint, normDir, refDir = axis2Placement3D(AFSobj['SurfaceParams'][1])
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        minorRadius = AFSobj['SurfaceParams'][-1]
        majorRadius = AFSobj['SurfaceParams'][-2]

        ParsedSurface['axisPoint'] = axisPoint
        ParsedSurface['normDir'] = normDir  # Z
        ParsedSurface['refDir'] = refDir  # X
        ParsedSurface['auxDir'] = auxDir

        ParsedSurface['minorRadius'] = minorRadius
        ParsedSurface['majorRadius'] = majorRadius

        ParsedSurface['pointFeature'] = {}
        ParsedSurface['pointFeature']['xyz'] = []
        ParsedSurface['pointFeature']['centroidDisp'] = []
        ParsedSurface['pointFeature']['maxima'] = []
        ParsedSurface['pointFeature']['minima'] = []

    AFSobj['ParsedSurface'] = ParsedSurface
    if not calcExtrema:
        # return ParsedSurface
        return

    # c = np.array([1, 1, 0])
    # axisPoint = np.array([0, 0, 10])
    # axisDir = np.array([0, 0, 1])
    # refDir = np.array([1, 0, 0])
    # tMajorRadius = 10
    # tMinorRadius = 1

    # test if toroid is rotationally symmetrical, if c is close to normDir through axisPoint

    rotSymDisp = np.linalg.norm(
        (centroid - axisPoint) - np.dot((centroid - axisPoint), normDir) * normDir
    )

    # need a precision factor that reflects increasing accuracy arising from iterative calculations of centroid --------------------------------------<<<<<<<<<<<<<<<<<
    # does one find a centroid through axis-intersection?
    # also scale factor affecting numerical precision
    # also STEP seems to be E-6 decimal places at best
    # iterativeScaledPrecision = 1e-4  # for testing

    if rotSymDisp < eps_STEP_AP21:
        ParsedEdge['rotSymFeature'] = dict()
        ParsedEdge['rotSymFeature']['rotSymCentre'] = [axisPoint,]
        # ParsedEdge['rotSymFeature']['rotSymRadius'] = radius
        # ParsedSurface['rotSymCentre'] = axisPoint
    else:
        # create axes for circles based at points of majorMinPoint, majorMaxPoint
        # minima case, create new axis at majorMinPoint
        # get min, max of major radius
        majorMinPoint, majorMaxPoint = pointCircleMinMaxDisp(
            c,
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

        # point on minor circle cross-section at minimum disp from centroid
        minorMinPoint, _ = pointCircleMinMaxDisp(
            c,
            majorMinPoint,
            majorMinTangentDir,
            minorRadius,
        )

        # ParsedSurface['minPoint'] = minorMinPoint

        # same for maxima
        majorMaxRadialDir = (majorMaxPoint - axisPoint) / np.linalg.norm(
            majorMaxPoint - axisPoint
        )
        majorMaxTangentDir = np.cross(majorMaxRadialDir, auxDir)
        majorMaxTangentDir = majorMaxTangentDir / np.linalg.norm(
            majorMaxTangentDir
        )

        minorMaxPoint, _ = pointCircleMinMaxDisp(
            c,
            majorMaxPoint,
            majorMaxTangentDir,
            minorRadius,
        )

        # np.linalg.norm(centroid - tMinorMinPoint) < np.linalg.norm(centroid - tMinorMaxPoint)
        # ParsedSurface['maxPoint'] = minorMaxPoint

        # determine if maxPoint/minPoint are inside or outside an edge_loop defined boundary
        # CCW directionality of edge vertices determines outward surface normal
        # surface minPoint/maxPoint always orthonormal to centroid

        if insideOutsideSurfaceTest(c, minorMaxPoint, AFSobj):  # maxima is among edges
            ParsedSurface['pointFeature']['xyz'].append(minorMaxPoint)
            ParsedSurface['pointFeature']['centroidDisp'].append(np.linalg.norm(minorMaxPoint - c))
            ParsedSurface['pointFeature']['maxima'].append(False)
            ParsedSurface['pointFeature']['minima'].append(True)

        if insideOutsideSurfaceTest(c, minorMinPoint, AFSobj):  # minima is among edges
            ParsedSurface['pointFeature']['xyz'].append(minorMinPoint)
            ParsedSurface['pointFeature']['centroidDisp'].append(np.linalg.norm(minorMinPoint - c))
            ParsedSurface['pointFeature']['maxima'].append(True)
            ParsedSurface['pointFeature']['minima'].append(False)

    AFSobj['ParsedSurface'] = ParsedSurface
    # return ParsedSurface


def cylinderSurfaceParse(AFSobj, c, calcExtrema=False):
    '''
    Parse a cylindrical surface, optionally extract surface minima points relative to a local centroid point.
    calcExtrema flag, exit function without expensive minima calculations.

    Used during centroid iteration
    :param AFSobj: parsed surface data object
    :param c: estimated surface median point based on extrema points identified
    :param calcExtrema: bool additional calculation to determine surface extrema
    :return ParsedSurface object
    '''

    # STEP ISO 10303-42 notes, EXPRESS descriptor
    # ENTITY Cylindrical_surface
    #   SUBTYPE OF (Surface);
    #   position : Axis_placement_3d;
    #   radius : positive_length_measure;

    # position: an Axis_placement_3d that defines the location and orientation of the surface.
    # The axis of the Cylindrical_surface passes through the location and is normal to the plane of ref_direction and axis.
    # radius: the radius of the circular curve of intersection between the surface of the cylinder and a plane
    # perpendicular to the axis of the Cylindrical_surface.

    if AFSobj.get('ParsedSurface') is not None:  # recalculation of max/min
        ParsedSurface = AFSobj['ParsedSurface']
        axisPoint = ParsedSurface['axisPoint']
        normDir = ParsedSurface['normDir']  # Z
        refDir = ParsedSurface['refDir']  # X
        auxDir = ParsedSurface['auxDir']  # Y
        radius = ParsedSurface['radius']
    else:
        axisPoint, normDir, refDir = axis2Placement3D(AFSobj['SurfaceParams'][1])
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        radius = AFSobj['SurfaceParams'][-1]

        ParsedSurface = {}
        ParsedSurface['radius'] = radius
        ParsedSurface['axisPoint'] = axisPoint
        ParsedSurface[
            'normDir'] = normDir  # Z # FreeCAD STEP generator reverses normal wrt to cylinder end planes.. ignore?
        ParsedSurface['refDir'] = refDir  # X
        ParsedSurface['auxDir'] = auxDir

        ParsedSurface['pointFeature'] = {}
        ParsedSurface['pointFeature']['xyz'] = []
        ParsedSurface['pointFeature']['centroidDisp'] = []
        ParsedSurface['pointFeature']['maxima'] = []
        ParsedSurface['pointFeature']['minima'] = []


    if not calcExtrema:
        AFSobj['ParsedSurface'] = ParsedSurface
        return

    # test if cylinder is rotationally symmetrical, if centroid is close to axisDir through axisPoint
    rotSymDisp = np.linalg.norm(
        (c - axisPoint) - np.dot((c - axisPoint), normDir) * normDir
    )

    if rotSymDisp < eps_STEP_AP21:
        minSurfaceDispPoint = c  # for creating a minima plane intersection
        # centroidProjAxisPoint = centroid
        # maxima is circular edge centre, determined from vertices
        # minima is orthogonal to centroid
        # shouldn't there also be a radius measure?? check: not in doctorate schema

        # todo, update concept where the region between a ridge and groove is a sphere
        # todo: truth determined via wandering hillclimb point path, or grad

        # test if any edge crosses cylinder midsection
        # edges should all have been calculated prior to surface

        # circle edge cannot be declared groove or ridge prior to discovery of adjacent surfaces
        # rotationally-symmetric edge must be assigned indeterminate status with centroid & radius

        # insideOutsideSurfaceTest() doesn't work for objects curved more than 180 degrees

        # test if any edges intersect a plane through centroid and normal to axis
        # edges must contain a minimum radius value
        # minimum value of all cylinder edge must equal radius at (points)/2 = 0 if surface contains minimum radius feature

        edgeMinimaCount = 0
        for edge in AFSobj['outerBoundEdgeLoop']:
            for cd in edge['pointFeature']['centroidDisp']:
                if np.isclose(cd, radius, eps):
                    edgeMinimaCount += 1

        if edgeMinimaCount % 2 == 0:
            ParsedSurface['rotSymFeature'] = dict()
            ParsedSurface['rotSymFeature']['rotSymCentre'] = c
            ParsedSurface['rotSymFeature']['rotSymRadius'] = radius

        # if edgeMinimaCount == 1: can't be rotational minima if only one point

        # featureExtrema may also encompass radial symmetric feature,
        # ['xyCentre'] = ['xy']
        # ['centreCentroidDisp']??
        # ['radiusCentroidDisp']

        #         # project centroid to symmetric axis to find rotSymCentre---HUH???
        #         centroidProjAxisPoint = pointProjectAxis(c, axisPoint, normDir)
        #         ParsedSurface['rotSymCentre'] = centroidProjAxisPoint
        #         # minPointCentroidDisp/maxPointCentroidDisp exist where centroid is on rotSym axis
        #         ParsedSurface['minPointCentroidDisp'] = np.sqrt(np.linalg.norm(centroidProjAxisPoint - c)**2 + radius**2)

    else:
        # cylinder/conic maxima is a search of greatest disp of vertices from centroid
        # if centroid is outside full cylindrical surface, edges contain maxima/minima
        # circle edges are maxima where circle normal contains local centroid,
        # cylinder minima -> test for intersection with surface on orthogonal projection from axis to centroid
        # if centroid within surface, |centroid - projAxisDir| < cRadius

        # test:
        # axisPoint = np.array([1, 0, 0])
        # c = np.array([2, 2, 0])
        # axisDir = np.array([1.5, 0, 0])
        # cRadius = 1

        # orthogonal projection of centroid to cylinder axis
        # vCA = centroid - axisPoint
        # vNA = axisPoint - normDir
        # centroidProjAxisPoint = vNA*np.dot(vCA, vNA)/np.dot(vNA, vNA)

        centroidProjAxisPoint = pointProjectAxis(c, axisPoint, normDir)

        # centroidProjAxisPoint = axisPoint + (
        #             np.dot((centroid - axisPoint), normDir) / np.dot(normDir, normDir)) * normDir
        centroidAxisDir = centroidProjAxisPoint - c
        centroidAxisDisp = np.linalg.norm(centroidAxisDir)
        centroidAxisDir = centroidAxisDir / centroidAxisDisp

        smp = ((centroidAxisDisp - radius) * centroidAxisDir)

        if np.linalg.norm(c - smp) > np.linalg.norm(c + smp):
            surfaceMinPoint = c - smp
        else:
            surfaceMinPoint = c + smp

        # can't determine in advance whether rotSymMax or rotSymCentre => requires subsequent comparison to adjacent nodes

        if insideOutsideSurfaceTest(c, surfaceMinPoint, AFSobj):  # minima is among edges
            ParsedSurface['pointFeature']['xyz'].append(surfaceMinPoint)
            ParsedSurface['pointFeature']['centroidDisp'].append(np.linalg.norm(surfaceMinPoint - c))
            ParsedSurface['pointFeature']['maxima'].append(False)
            ParsedSurface['pointFeature']['minima'].append(True)

    AFSobj['ParsedSurface'] = ParsedSurface


def coneSurfaceParse(AFSobj, c, calcExtrema=False):
    '''
    Parse a conical surface, optionally extract surface minima points relative to a local centroid point.
    calcExtrema flag, exit function without expensive minima calculations (used during centroid iteration)

    :param AFSobj: parsed surface data object
    :param c: estimated surface median point based on extrema points identified
    :param calcExtrema: bool additional calculation to determine surface extrema
    :return ParsedSurface object
    '''

    # STEP ISO 10303-42 notes, EXPRESS descriptor
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
    # note 2 that the radial surface may define the top of a frustum rather than the base, need to check enclosing edges
    # https://quaoar.su/files/standards/Standard%20-%202003%20-%20ISO%2010303-42.pdf p88

    if AFSobj.get('ParsedSurface') is not None:  # recalculation of max/min
        ParsedSurface = AFSobj['ParsedSurface']
        axisPoint = ParsedSurface['axisPoint']
        normDir = ParsedSurface['normDir']  # Z
        refDir = ParsedSurface['refDir']  # X
        auxDir = ParsedSurface['auxDir']  # Y
        radius = ParsedSurface['radius']
        semiAngle = ParsedSurface['semiAngle']

    else:
        axisPoint, normDir, refDir = axis2Placement3D(AFSobj['SurfaceParams'][1])
        auxDir = np.cross(normDir, refDir)
        auxDir = auxDir / np.linalg.norm(auxDir)
        semiAngle = AFSobj['SurfaceParams'][-1]
        radius = AFSobj['SurfaceParams'][-2]

        ParsedSurface = {}
        ParsedSurface['axisPoint'] = axisPoint
        ParsedSurface['normDir'] = normDir  # Z
        ParsedSurface['refDir'] = refDir  # X
        ParsedSurface['auxDir'] = auxDir
        ParsedSurface['semiAngle'] = semiAngle
        ParsedSurface['radius'] = radius

        ParsedSurface['pointFeature'] = {}
        ParsedSurface['pointFeature']['xyz'] = []
        ParsedSurface['pointFeature']['centroidDisp'] = []
        ParsedSurface['pointFeature']['maxima'] = []
        ParsedSurface['pointFeature']['minima'] = []

        # for the purpose of extracting axisPoint, normDir for inside/outside point detection
        AFSobj['ParsedSurface'] = ParsedSurface

    AFSobj['ParsedSurface'] = ParsedSurface
    if not calcExtrema:
        # return ParsedSurface
        return

    # test if cone is rotationally symmetrical, if centroid is close to axisDir through axisPoint
    rotSymDisp = np.linalg.norm((c - axisPoint) - np.dot((c - axisPoint), normDir) * normDir)

    # rotSymDisp = pointProjectAxis(c, axisPoint, normDir)

    axisPointApexDisp = radius / np.tan(semiAngle)
    # find point at cone/triangle apex
    apexPoint = axisPoint + axisPointApexDisp * (-normDir)

    # circle edge cannot be declared groove or ridge prior to discovery of adjacent surfaces
    # rotationally-symmetric edge must be assigned indeterminate status with centroid & radius

    # test whether centroid is on the far side of the apex
    if np.linalg.norm(centroid - apexPoint) < np.linalg.norm(centroid - axisPoint):
        if insideOutsideSurfaceTest(c, apexPoint, AFSobj):
            ParsedEdge['rotSymFeature'] = dict()
            ParsedEdge['rotSymFeature']['rotSymCentre'] = apexPoint
            ParsedEdge['rotSymFeature']['rotSymRadius'] = 0

            # ParsedSurface['rotSymCentre'] = apexPoint
            # ParsedSurface['rotSymRadius'] = 0

    else:
        if rotSymDisp > eps_STEP_AP21:
            # get normal to plane containing centroid & cone axis, cross product of (apex - centroid) X (apex - coneLocation)
            minPlaneNorm = np.cross((apexPoint - c), (apexPoint - axisPoint))
            minPlaneNorm = minPlaneNorm / np.linalg.norm(minPlaneNorm)

            coneEdgeCW = np.dot(rotationMatrix(minPlaneNorm, semiAngle), (apexPoint - axisPoint))
            coneEdgeCW = coneEdgeCW / np.linalg.norm(coneEdgeCW)
            coneMinCW = pointProjectAxis(c, apexPoint, coneEdgeCW)

            coneEdgeCCW = np.dot(rotationMatrix(minPlaneNorm, -semiAngle), (apexPoint - axisPoint))
            coneEdgeCCW = coneEdgeCCW / np.linalg.norm(coneEdgeCCW)
            coneMinCCW = pointProjectAxis(c, apexPoint, coneEdgeCCW)

            if np.linalg.norm(c - coneMinCW) >= np.linalg.norm(c - coneMinCCW):
                coneMin = coneMinCCW
            else:
                coneMin = coneMinCW

            if insideOutsideSurfaceTest(c, coneMin, AFSobj):
                ParsedSurface['pointFeature']['xyz'].append(coneMin)
                ParsedSurface['pointFeature']['centroidDisp'].append(np.linalg.norm(coneMin - c))
                ParsedSurface['pointFeature']['maxima'].append(False)
                ParsedSurface['pointFeature']['minima'].append(True)

        else:
            # there is a band around a cone corresponding to a minima
            coneMinimaCentreDisp = np.linalg.norm(centroid - apexPoint) * np.cos(semiAngle)
            #coneMinimaCentrePoint = apexPoint + (coneMinimaCentreDisp * normDir)
            coneMinimaCentrePoint = apexPoint + (coneMinimaCentreDisp * np.cos(semiAngle) * normDir)

            # check is band within edges,
            rotSymCentroidDisp = np.linalg.norm(centroid - apexPoint) * np.sin(semiAngle)

            edgeMinimaCount = 0
            for edge in AFSobj['outerBoundEdgeLoop']:
                for cd in edge['pointFeature']['centroidDisp']:
                    if np.isclose(cd, rotSymCentroidDisp, eps_STEP_AP21):
                        edgeMinimaCount += 1

            if edgeMinimaCount % 2 == 0:
                ParsedSurface['rotSymFeature'] = dict()
                ParsedSurface['rotSymFeature']['rotSymCentre'] = [coneMinimaCentrePoint,]
                ParsedSurface['rotSymFeature']['rotSymRadius'] = [np.linalg.norm(
                    coneMinimaCentrePoint - apexPoint) * np.tan(semiAngle),]
                # ParsedSurface['rotSymRadius'] = np.linalg.norm(coneMinimaCentrePoint - apexPoint) * np.tan(semiAngle)
                # ParsedSurface['rotSymCentre'] = coneMinimaCentrePoint

    # else:
    #
    #     # project centroid to nearest side, i.e.
    #     # construct 2 rays through apex at half angles from cone axis within same plane as centroid
    #     # shortest orthogonal projection to either line is surface minima,
    #
    #     # project centroid to symmetric axis to find rotSymCentre
    #
    #     # cylinder/conic maxima is a search of greatest disp of vertices from centroid
    #     # if centroid is outside full cylindrical/conical surface radius,
    #     # test for intersection with surface on orthogonal projection from axis to centroid
    #     # if centroid within surface, |centroid - projAxisDir| < cRadius
    #
    #     # test:
    #     # axisPoint = np.array([1, 0, 0])
    #     # c = np.array([2, 2, 0])
    #     # axisDir = np.array([1.5, 0, 0])
    #     # cRadius = 1
    #
    #     # orthogonal projection of centroid to cylinder axis
    #     # vCA = centroid - axisPoint
    #     # vNA = axisPoint - normDir
    #     # centroidProjAxisPoint = vNA*np.dot(vCA, vNA)/np.dot(vNA, vNA)
    #
    #     #centroidProjAxisPoint = pointProjectAxis(c, axisPoint, normDir)-------------------------------------broken??
    #
    #     centroidProjAxisPoint = axisPoint + ( np.dot((centroid - axisPoint), normDir) / np.dot(normDir, normDir)) * normDir
    #     centroidAxisDir = centroidProjAxisPoint - c
    #     centroidAxisDisp = np.linalg.norm(centroidAxisDir)
    #     centroidAxisDir = centroidAxisDir / centroidAxisDisp
    #
    #     coneMinimaCentreDisp = axisPointApexDisp / np.cos(semiAngle)
    #     coneMinimaCentrePoint = apexPoint + (coneMinimaCentreDisp * normDir)
    #
    #     # calculate cone radius at min
    #     coneMinRadius = np.linalg.norm(apexPoint - coneMinimaCentrePoint) * np.tan(semiAngle)
    #
    #     smp1 = coneMinimaCentrePoint + (coneMinRadius * -centroidAxisDir)
    #     smp2 = coneMinimaCentrePoint - (coneMinRadius * -centroidAxisDir)
    #
    #     if np.linalg.norm(c - smp1) > np.linalg.norm(c - smp2):
    #         surfaceMinPoint = smp1
    #     else:
    #         surfaceMinPoint = smp2
    #
    #     if insideOutsideSurfaceTest(c, surfaceMinPoint, AFSobj):  # minima is among edges
    #         ParsedSurface['pointFeature']['xyz'].append(surfaceMinPoint)
    #         ParsedSurface['pointFeature']['centroidDisp'].append(np.linalg.norm(surfaceMinPoint - c))

    AFSobj['ParsedSurface'] = ParsedSurface


def BSplineSurfaceWithKnotsParse(AFSobj, c, calcExtrema=False):
    '''
    Parse a BSplineSurfaceWithKnotsParse surface, optionally extract surface minima points relative to a
    local centroid point.
    calcExtrema flag, exit function without expensive minima calculations (used during centroid iteration)

    :param AFSobj: parsed surface data object
    :param c: estimated surface median point based on extrema points identified
    :param STEPobj: object containing STEP data
    :param calcExtrema: bool additional calculation to determine surface extrema
    :return ParsedSurface object
    '''

    # STEP ISO 10303-42 notes, EXPRESS descriptor
    # ENTITY Surface_with_explicit_knots
    #   SUBTYPE OF (B_spline_surface);
    #   u_knot_multiplicities : LIST[2:?] OF INTEGER;
    #   u_knot_values : LIST[2:?] OF parameter_value;
    #   v_knot_multiplicities : LIST[2:?] OF INTEGER;
    #   v_knot_values : LIST[2:?] OF parameter_value;
    # WHERE
    #   WR1: SIZEOF(u_knot_multiplicities) = SIZEOF(u_knot_values);
    #   WR2: SIZEOF(v_knot_multiplicities) = SIZEOF(v_knot_values);

    # Attribute definitions:
    # u_knot_multiplicities: the list of multiplicities for the u_knot_values, these are in the first parametric direction for the surface.
    # u_knot_values: the list of knot values used to define the B-spline basis functions in the first parametric direction of the surface.
    # v_knot_multiplicities: the list of multiplicities for the v_knot_values, these are in the second parametric direction for the surface.
    # v_knot_values: the list of knot values used to define the B-spline basis functions in the second parametric direction of the surface.

    # u_multiplicities: The multiplicities of the knots in the u parameter direction.
    # v_multiplicities: The multiplicities of the knots in the v parameter direction.
    # u_knots: The list of the distinct knots in the u parameter direction.
    # v_knots: The list of the distinct knots in the v parameter direction.
    # knot_spec: The description of the knot type.
    # knot_u_upper: The number of distinct knots in the u parameter direction.
    # knot_v_upper: The number of distinct knots in the v parameter direction.
    # SELFnb_spline_surface.u_degree: Algebraic degree of basis functions in u.
    # SELFnb_spline_surface.v_degree: Algebraic degree of basis functions in v.
    # SELFnb_spline_surface.control_points_list: This is a list of lists of control points.
    # SELF\b_spline_surface.surface_form: Indicator of special surface types. (See 4.3.4.)
    # SELF\b_spline_surface.u_closed: Indication of whether the surface is closed in the u direction; this
    # is for information only.
    # SELFnb_spline_surface.v_closed: Indication of whether the surface is closed in the v direction; this is
    # for information only.
    # SELFnb_spline_surface.self_intersect: Flag to indicate whether, or not, surface is self-intersecting;
    # this is for information only.
    # SELFnb_spline_surface.u_upper: Upper index on control points in u direction.
    # SELFnb_spline_surface.v_upper: Upper index on control points in v direction.
    # SELFnb_spline_surface.control_points: Array (two-dimensional) of control points defining surface
    # geometry. This array is constructed from the control points list

    if AFSobj.get('ParsedSurface') is not None:
        ParsedSurface = AFSobj['ParsedSurface']
        surfaceUdegree = ParsedSurface['surfaceUdegree']
        surfaceVdegree = ParsedSurface['surfaceVdegree']
        controlPointsRefs = ParsedSurface['controlPointsRefs']
        surfaceForm = ParsedSurface['surfaceForm']
        closedU = ParsedSurface['closedU']
        closedV = ParsedSurface['closedV']
        selfIntersect = ParsedSurface['selfIntersect']
        knotUmultiplicities = ParsedSurface['knotUmultiplicities']
        knotVmultiplicities = ParsedSurface['knotVmultiplicities']
        knotsU = ParsedSurface['knotsU']
        knotsV = ParsedSurface['knotsV']
        knotSpec = ParsedSurface['knotSpec']

        controlPoints = ParsedSurface['controlPoints']
        STEPknotUvector = ParsedSurface['STEPknotUvector']
        STEPknotVvector = ParsedSurface['STEPknotVvector']
        BsplineKnotSurface = ParsedSurface['BsplineKnotCurve']

    else:
        if type(AFSobj['SurfaceParams'][0]) == str:  # check for name string
            offset = 1
        else:
            offset = 0

        surfaceUdegree = AFSobj['SurfaceParams'][offset]
        surfaceVdegree = AFSobj['SurfaceParams'][offset + 1]
        controlPointsRefs = AFSobj['SurfaceParams'][offset + 2]  # need to flatten
        surfaceForm = AFSobj['SurfaceParams'][offset + 3]
        closedU = 'T' in AFSobj['SurfaceParams'][offset + 4]
        closedV = 'T' in AFSobj['SurfaceParams'][offset + 5]
        selfIntersect = 'T' in AFSobj['SurfaceParams'][offset + 6]
        knotUmultiplicities = AFSobj['SurfaceParams'][offset + 7]
        knotVmultiplicities = AFSobj['SurfaceParams'][offset + 8]
        knotsU = AFSobj['SurfaceParams'][offset + 9]
        knotsV = AFSobj['SurfaceParams'][offset + 10]
        knotSpec = AFSobj['SurfaceParams'][offset + 11]

        ParsedSurface = {}
        ParsedSurface['pointFeature'] = {}
        ParsedSurface['pointFeature']['xyz'] = []
        ParsedSurface['pointFeature']['centroidDisp'] = []
        ParsedSurface['pointFeature']['maxima'] = []
        ParsedSurface['pointFeature']['minima'] = []

        ParsedSurface['surfaceUdegree'] = surfaceUdegree
        ParsedSurface['surfaceVdegree'] = surfaceVdegree
        ParsedSurface['controlPointsRefs'] = controlPointsRefs
        ParsedSurface['surfaceForm'] = surfaceForm
        ParsedSurface['closedU'] = closedU
        ParsedSurface['closedV'] = closedV
        ParsedSurface['selfIntersect'] = selfIntersect
        ParsedSurface['knotUmultiplicities'] = knotUmultiplicities
        ParsedSurface['knotVmultiplicities'] = knotVmultiplicities
        ParsedSurface['knotsU'] = knotsU
        ParsedSurface['knotsV'] = knotsV
        ParsedSurface['knotSpec'] = knotSpec

        # controPointsListRefs structured as U control points length list of V control points length list

        controlPointsLenU = len(controlPointsRefs)
        controlPointsLenV = len(controlPointsRefs[0])
        # controlPointsRefs = cleanSubRefs(controlPointsRefs)

        # extract control points
        # controlPoints = []
        # for cplU in controlPointsRefs:
        #     cpU = []
        #     for cplV in cplU:
        #         if (STEPobj[ref2index(cplV)].type_name == 'CARTESIAN_POINT'):
        #             cpV = CP2point(cleanSubRefs(cplV)[0])
        #             cpU.append(cpV)
        #     controlPoints.append(cpU)

        controlPoints = []
        for cplU in controlPointsRefs:
            for cplV in cplU:
                if (STEPobj[ref2index(cplV)].type_name == 'CARTESIAN_POINT'):
                    controlPoints.append(CP2point(cleanSubRefs(cplV)[0]))

        ParsedSurface['controlPoints'] = controlPoints

        # construct full U & V knotvectors
        STEPknotUvector = []
        for kml in range(len(knotUmultiplicities)):
            for i in range(knotUmultiplicities[kml]):
                STEPknotUvector.append(knotsU[kml])
        ParsedSurface['STEPknotUvector'] = STEPknotUvector

        if (len(STEPknotUvector) - controlPointsLenU - surfaceUdegree) != 1:
            print("maldefined B-spline surface (U) !")

        # p.106 https://quaoar.su/files/standards/Standard%20-%202003%20-%20ISO%2010303-42.pdf

        STEPknotVvector = []
        for kml in range(len(knotVmultiplicities)):
            for i in range(knotVmultiplicities[kml]):
                STEPknotVvector.append(knotsV[kml])
        ParsedSurface['STEPknotVvector'] = STEPknotVvector

        # note _knotXvector values not used as Piegel+Tiller algorithms don't always work with explicit knots

        if (len(STEPknotVvector) - controlPointsLenV - surfaceVdegree) != 1:
            print("maldefined B-spline surface (V) !")

        # BsplineKnotSurface = BSpline.Surface(normalize_kv=True)
        BsplineKnotSurface = NURBS.Surface()
        BsplineKnotSurface.degree_u = surfaceUdegree
        BsplineKnotSurface.degree_v = surfaceVdegree

        BsplineKnotSurface.ctrlpts_size_u = controlPointsLenU
        BsplineKnotSurface.ctrlpts_size_v = controlPointsLenV

        ctrlptsw = compatibility.combine_ctrlpts_weights(controlPoints, weights=None)
        BsplineKnotSurface.set_ctrlpts(ctrlptsw, controlPointsLenU, controlPointsLenV)

        # BsplineKnotSurface.ctrlpts2d = controlPoints
        # BsplineKnotSurface.knotvector_u = STEPknotUvector
        # BsplineKnotSurface.knotvector_v = STEPknotVvector

        BsplineKnotSurface.knotvector_u = utilities.generate_knot_vector(surfaceUdegree,
                                                                         controlPointsLenU)
        BsplineKnotSurface.knotvector_v = utilities.generate_knot_vector(surfaceVdegree,
                                                                         controlPointsLenV)

        ParsedSurface['BsplineKnotCurve'] = BsplineKnotSurface

        # No compelling reason to set evaluation delta & evaluate surface points as evaluation seems incorrect
        # BsplineKnotSurface.delta = 0.025  # this seems to be the minima delta under adaptive tesselation
        # BsplineKnotSurface.evaluate()

    if not calcExtrema:
        AFSobj['ParsedSurface'] = ParsedSurface
        # return ParsedSurface
        return

    # local minima/maxima
    maxPointsUV = rationalSurfaceExtremaParam_5(BsplineKnotSurface,
                                                c,
                                                maxSearch=True,
                                                localExtrema=True,
                                                curvatureTest=False,
                                                uv_xyz=True)

    maxPoints = BsplineKnotSurface.evaluate_list(maxPointsUV)

    minPointsUV = rationalSurfaceExtremaParam_5(BsplineKnotSurface,
                                                c,
                                                maxSearch=False,
                                                localExtrema=True,
                                                curvatureTest=False,
                                                uv_xyz=True)

    minPoints = BsplineKnotSurface.evaluate_list(minPointsUV)

    # normalise u,v to [0, 1] - not required with geomdl generated knot vectors
    # knotUrange = STEPknotUvector[-1] - STEPknotUvector[0]
    # knotVrange = STEPknotVvector[-1] - STEPknotVvector[0]
    # maxPointsUV = [((mpu[0] - STEPknotUvector[0])/knotUrange, (mpu[1] - STEPknotVvector[0])/knotVrange ) for mpu in maxPointsUV]
    # minPointsUV = [((mpu[0] - STEPknotUvector[0])/knotUrange, (mpu[1] - STEPknotVvector[0])/knotVrange ) for mpu in minPointsUV]

    maxima = [True, ] * len(maxPoints) + [False, ] * len(minPoints)
    minima = [False, ] * len(maxPoints) + [True, ] * len(minPoints)

    extremaUV = maxPointsUV + minPointsUV
    extremaUV = [euv.tolist() for euv in extremaUV]
    maximaPoints = maxPoints + minPoints
    if len(maximaPoints) > 1:
        # create explicit dtype fields to permit sorting via u, then v
        #extremaUV_field = extremaUV.ravel().view(dtype=[('u', extremaUV.dtype), ('v', extremaUV.dtype)]) # ravel() fails under Pypi 3.9
        extremaUV_field = np.array(extremaUV)
        extremaUV_field = extremaUV_field.ravel().view(dtype=[('u', float), ('v', float)])
        extremaUVindex = np.argsort(extremaUV_field, order=('u', 'v'))
        extremaUV = [extremaUV[euvi] for euvi in extremaUVindex]
        maximaPoints = [maximaPoints[euvi] for euvi in extremaUVindex]
        maxima = [maxima[euvi] for euvi in extremaUVindex]
        minima = [minima[euvi] for euvi in extremaUVindex]

        # because Newton-Raphson will converge on a tangent, separate max/min values may be identical
        maxima_truth = [np.allclose(extremaUV[mu + 1], extremaUV[mu], eps_STEP_AP21) for mu in
                        range(0, len(extremaUV) - 1)]
        maxima_truth = maxima_truth + [False, ]
        extremaUV = [mu for (mu, mt) in zip(extremaUV, maxima_truth) if not mt]
        maximaPoints = [mp for (mp, mt) in zip(maximaPoints, maxima_truth) if not mt]

    ParsedSurface['pointFeature']['xyz'] = maximaPoints
    ParsedSurface['pointFeature']['uv'] = extremaUV #= [euv.tolist() for euv in extremaUV]
    ParsedSurface['pointFeature']['centroidDisp'] = [np.linalg.norm(mp - c) for mp in maximaPoints]
    ParsedSurface['pointFeature']['maxima'] = maxima
    ParsedSurface['pointFeature']['minima'] = minima

    # detect rotationally symmetric spline surfaces, spline surface cylinder
    # get the median point of rings/arcs of control points
    spinePoints = np.array([np.array(s).mean(axis=0) for s in BsplineKnotSurface.ctrlpts2d])

    # test for deviation from mean point
    rowDelta = []
    for s in range(0, len(spinePoints) - 1):
        ctrlptsRowDisp = [np.linalg.norm(spinePoints[s] - sc) for sc in BsplineKnotSurface.ctrlpts2d[s]]
        rowDelta.append(max(ctrlptsRowDisp) - min(ctrlptsRowDisp))

    # low rowDelta values suggest rotational symmetry,
    # similar rowDelta values suggest consistent cross-section (with possible helicity):
    # consistentCrossSection = np.array(rowDelta).var()
    # should show up with edge maxima/minima

    # spinePointDisps = [np.linalg.norm(spinePoints[s] - BsplineKnotSurface.ctrlpts2d[s]) for s in range(0, len(BsplineKnotSurface.ctrlpts2d))]
    spinePointsMean = spinePoints.mean(axis=0)
    uu, dd, vv = np.linalg.svd(spinePoints - spinePointsMean)

    # vv[0] contains the first principal component, i.e. the direction
    # vector of the 'best fit' line in the least squares sense.
    # project mean control point rows to derived PCA spine and assess spine axis straightness
    spineDir = vv[0] - spinePointsMean
    spinePointDir = [s - spinePointsMean for s in spinePoints]
    spineIP = np.dot(spineDir, spineDir)
    projSpinePointDir = [np.dot(s, spineDir) / spineIP for s in spinePointDir]

    projSpinePoint = [spinePointsMean + pspd * spineDir for pspd in projSpinePointDir]
    # test for centroid disp from spine axis, giving minima rotation

    spineDeviations = [np.linalg.norm(projSpinePoint[n] - spinePoints[n]) for n in range(0, len(spinePoints - 1))]

    if (max(rowDelta) < eps_STEP_AP21) and (
            max(spineDeviations) < eps_STEP_AP21):  # rotSymLimit constant --------------------------------------------
        # ParsedSurface['rotSymFlag'] = 1
        ParsedSurface['rotSymFeature'] = dict()
        # ParsedSurface['rotSymFeature']['rotSymCentre'] = ?? centroid minima within
        ParsedSurface['rotSymFeature']['rotSymRadius'] = spinePointsMean

    AFSobj['ParsedSurface'] = ParsedSurface


def planeSurfaceParse_2(AFSobj, c, calcExtrema=False):
    '''
    Parse a planar surface, optionally extract surface minima points relative to a
    local centroid point.
    calcExtrema flag, exit function without expensive minima calculations (used during centroid iteration)

    :param AFSobj: parsed surface data object
    :param c: estimated surface median point based on extrema points identified
    :param calcExtrema: bool additional calculation to determine surface extrema
    :return ParsedSurface object
    '''

    # STEP ISO 10303-42 notes, EXPRESS descriptor
    # ENTITY Plane
    #   SUBTYPE OF (Surface);
    #   position : Axis_placement_3d;

    # with a planar surface, closest point is that point closest to orthogonal intersection of ray through
    # centroid with plane of surface

    if AFSobj.get('ParsedSurface') is not None:  # recalculation of max/min
        ParsedSurface = AFSobj['ParsedSurface']
        axisPoint = ParsedSurface['axisPoint']
        normDir = ParsedSurface['normDir']  # Z
        # refDir = ParsedSurface['refDir'] # X
        # auxDir = ParsedSurface['auxDir']  # Y

    else:
        ParsedSurface = {}
        axisPoint, normDir, refDir = axis2Placement3D(AFSobj['SurfaceParams'][STEPlabelOffset(1, AFSobj['SurfaceParams'])])
        ParsedSurface['axisPoint'] = axisPoint
        ParsedSurface['normDir'] = normDir  # Z
        ParsedSurface['refDir'] = refDir  # X

        ParsedSurface['pointFeature'] = {}
        ParsedSurface['pointFeature']['xyz'] = []
        ParsedSurface['pointFeature']['centroidDisp'] = []
        ParsedSurface['pointFeature']['maxima'] = []
        ParsedSurface['pointFeature']['minima'] = []

    if not calcExtrema:
        AFSobj['ParsedSurface'] = ParsedSurface
        # return ParsedSurface
        return

    # project centroid to plane defined by normDir and axisPoint
    surfaceMinPoint = c - np.dot(c - axisPoint, normDir) * normDir

    if insideOutsideSurfaceTest(c, surfaceMinPoint, AFSobj):  # minima is among edges
        # replace not append on every calculation
        ParsedSurface['pointFeature']['xyz'] = [surfaceMinPoint]
        ParsedSurface['pointFeature']['centroidDisp'] = [np.linalg.norm(surfaceMinPoint - c)]
        ParsedSurface['pointFeature']['maxima'] = [False, ]
        ParsedSurface['pointFeature']['minima'] = [True, ]

    AFSobj['ParsedSurface'] = ParsedSurface


def surfaceParse(AFSobj, c, calcExtrema=True):
    '''
    Process surface parsing according to STEP shape type
    calcExtrema flag, exit function without expensive minima calculations (used during centroid iteration)

    :param AFSobj: parsed surface data object
    :param c: estimated surface median point based on extrema points identified
    :param calcExtrema: bool additional calculation to determine surface extrema
    :return ParsedSurface object
    '''

    if AFSobj['SurfaceTypeName'] not in [
        'SPHERICAL_SURFACE',
        'TOROIDAL_SURFACE',
        'CYLINDRICAL_SURFACE',
        'CONICAL_SURFACE',
        'PLANE',
        'B_SPLINE_SURFACE_WITH_KNOTS',
        'SURFACE_OF_REVOLUTION',
        'BOUNDED_SURFACE'
    ]:
        print("Untreated surface type: " + AFSobj['SurfaceTypeName'])
        # SURFACE_OF_REVOLUTION
    else:
        if AFSobj['SurfaceTypeName'] == 'BOUNDED_SURFACE':
            boundedSurface(AFSobj, c, calcExtrema)

        if AFSobj['SurfaceTypeName'] == 'SURFACE_OF_REVOLUTION':
            revolvedSurfaceParse(AFSobj, c, calcExtrema)

        if AFSobj['SurfaceTypeName'] == 'SPHERICAL_SURFACE':
            sphereSurfaceParse(AFSobj, c, calcExtrema)

        if AFSobj['SurfaceTypeName'] == 'TOROIDAL_SURFACE':
            toroidSurfaceParse(AFSobj, c, calcExtrema)

        if (AFSobj['SurfaceTypeName'] == 'CYLINDRICAL_SURFACE'):
            cylinderSurfaceParse(AFSobj, c, calcExtrema)

        if (AFSobj['SurfaceTypeName'] == 'CONICAL_SURFACE'):  # still must test with simple STEP file
            coneSurfaceParse(AFSobj, c, calcExtrema)

        if AFSobj['SurfaceTypeName'] == 'PLANE':
            planeSurfaceParse_2(AFSobj, c, calcExtrema)

        if AFSobj['SurfaceTypeName'] == 'B_SPLINE_SURFACE_WITH_KNOTS':
            BSplineSurfaceWithKnotsParse(AFSobj, c, calcExtrema)


def parseStepcodePart21object():
    # get STEP surface instances and edge loop references to associated edges, parse within AdvancedFaceSurfaces list of features (not applicable to OFF format)
    AFSobj = []
    for se in STEPobj:
        if hasattr(se, 'type_name'):  # stepcode.Part21.ComplexEntity type missing attribute
            if (se.type_name == 'ADVANCED_FACE') or (se.type_name == 'FACE_SURFACE'):
                SurfaceNormalOutwards = 'T' in se.params[-1]
                afRefs = cleanSubRefs(se.params)
                SurfaceClass = {
                    'SurfaceNormalOutwards': SurfaceNormalOutwards,
                    'SurfaceTypeName': None,
                    'SurfaceRef': None,
                    'SurfaceParams': None,
                    'EdgeLoopList': [],
                }
                for se2 in afRefs: # should prolly assume [0] is FACE_*-----------------<<<<<<
                    if hasattr(STEPobj[ref2index(se2)], 'type_name'):
                        afTypeName = STEPobj[ref2index(se2)].type_name
                        # no risk of complex entities from here on in
                        if (afTypeName == 'FACE_OUTER_BOUND') or ( afTypeName == 'FACE_BOUND'): # get edge loop references
                            # todo assume 'FACE_INNER_BOUND' or 'HOLE' etc can run the same edge loop from here ----------------------<<<<<<
                            se3 = STEPobj[ref2index(se2)].params
                            se3 = cleanSubRefs(se3)[0]
                            if STEPobj[ref2index(se3)].type_name == 'EDGE_LOOP':
                                SurfaceClass['EdgeLoopList'] = STEPobj[ref2index(se3)].params[-1]
                        else:
                            # print(afTypeName)
                            SurfaceClass['SurfaceTypeName'] = afTypeName
                            SurfaceClass['SurfaceRef'] = se2
                            SurfaceClass['SurfaceParams'] = STEPobj[
                                ref2index(se2)
                            ].params

                    elif not hasattr(STEPobj[ref2index(se2)], 'type_name'): # complex entity

                        # several different flavours of complex entity - need a disambiguation function
                        # complex entity as a subset of advanced face surfaces
                        if STEPobj[ref2index(se2)].params[0].type_name == 'BOUNDED_SURFACE':
                            #complexEntityParse(ref2index(se2))
                            SurfaceClass['SurfaceTypeName'] = STEPobj[ref2index(se2)].params[0].type_name
                            SurfaceClass['SurfaceRef'] = se2
                        else:
                            print('unhandled complex entity @todo')
                AFSobj.append(SurfaceClass)

    return AFSobj


def getConvexHullPointSetSTEP():
    # get convex hull points of model in order to determine a reliable median centroid (not applicable to OFF format)
    # exclude Advanced_brep_shape_representation which contains a [0,0,0] origin point, e.g.
    # 10 = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#11,#15),#113);
    # 11 = AXIS2_PLACEMENT_3D('',#12,#13,#14);
    # 12 = CARTESIAN_POINT('',(0.,0.,0.));
    # 13 = DIRECTION('',(0.,0.,1.));
    # 14 = DIRECTION('',(1.,0.,-0.));

    outermostPoints = []
    ABSRpointRef = None
    for ABSR in STEPobj:
        if hasattr(ABSR, 'type_name'):  # stepcode.Part21.ComplexEntity type missing attribute
            if ABSR.type_name == 'ADVANCED_BREP_SHAPE_REPRESENTATION':
                cleanABSR = cleanSubRefs(ABSR.params)
                if STEPobj[ref2index(cleanABSR[0])].type_name == 'AXIS2_PLACEMENT_3D':
                    ABSRpointRef = cleanSubRefs(STEPobj[ref2index(cleanABSR[0])].params)[0]

        for vp in STEPobj:
            if hasattr(vp, 'type_name'):  # stepcode.Part21.ComplexEntity type missing attribute
                if vp.type_name == 'CARTESIAN_POINT':
                    if (len(vp.params[-1]) == 3):
                        if (vp.ref == ABSRpointRef):
                            #print(" unresolved ComplexEntity type attribute")
                            pass
                        else:
                            outermostPoints.append(np.array([vp.params[-1][0], vp.params[-1][1], vp.params[-1][2]]))

    return outermostPoints


#=========FORMAT stepcode.Part21 objects INTO USEFUL LOCAL DATA STRUCTURES WHERE INTERRELATIONSHIPS ARE RECORDED


AdvancedFaceSurfaces = parseStepcodePart21object()


# get convex hull points of model in order to determine a reliable median centroid (not applicable to OFF format)
outermostPoints = getConvexHullPointSetSTEP()

# OFF convex hull points is simply points set
# outermostPoints = [np.array(c) for c in OFF_entities.points]

centroid = medianPoint(outermostPoints)
# continue to ignore local maxima at both radial extents and NURB surface maxima?

def calculateEdgeMaxMin(AFSobj):
    # for all elements, surfaces and bounding edges, iterate maxima discovery relative to centroid
    # unique to stepcode.Part21 input
    for afs in AFSobj:
        EdgeLoops = afs['EdgeLoopList']
        ParsedEdgeLoop = []
        for edgeRefInstance in EdgeLoops:
            ParsedEdge = edgeSTEPparse(edgeRefInstance)
            ParsedEdgeLoop.append(edgeParse(centroid, ParsedEdge, calcExtrema=True))

        afs['outerBoundEdgeLoop'] = ParsedEdgeLoop
        # don't parse surface minima until centroid is stable
        surfaceParse(afs, centroid, calcExtrema=False)


# if not OFF
calculateEdgeMaxMin(AdvancedFaceSurfaces)


def centroidIterate(AFSobj, centroid, verbose=False):
    # iteratively calculate max/min points and median centroid point
    # extract max values for all surfaces to add to median vertex calculation for centroid
    # (should work with OFF)
    #lastCentroid2 = centroid # twin attractor issue------------------------------------------------------------------
    lastCentroid = centroid
    centroid = medianOfFeaturePoints(AFSobj)
    if len(centroid) == 3:
        # iterate centroid-features relationship
        while (np.linalg.norm(lastCentroid - centroid) > eps): # and (np.linalg.norm(lastCentroid2 - centroid) > eps):
            # recalculate surfaces & edges max/min points according to new centroid
            for afs in AFSobj:
                # update surface after edges
                for edge in afs['outerBoundEdgeLoop']:
                    edgeParse(centroid, edge, calcExtrema=True)
                # doesn't seem much point in parsing surface minima until centroid is stable, possible exception of bspline maxima, minimaFlag ?
                surfaceParse(afs, centroid, calcExtrema=False)

            #lastCentroid2 = lastCentroid
            lastCentroid = centroid
            centroid = medianOfFeaturePoints(AFSobj)
            # for afs in AFSobj:
            #     if afs['ParsedSurface']['maxPoint'] is not None:
            #         outermostPoints.append(afs['ParsedSurface']['maxPoint'])
            #         for edge in afs['outerBoundEdgeLoop']:
            #             outermostPoints.append(edge['maxPoint'])
            #
            # centroid = medianPoint(outermostPoints)
            if verbose:
                print('centroid update: ' + str(np.linalg.norm(lastCentroid - centroid)))
                print('centroid: ' + str(centroid))
    else:
        # case of no surfaces
        centroid = lastCentroid

    return centroid


# iteratively calculate max/min points and median centroid point
# extract max values for all surfaces to add to median vertex calculation for centroid
centroid = centroidIterate(AdvancedFaceSurfaces, centroid, verbose=True)


# expensive surface minima calculations after centroid is settled
for afs in AdvancedFaceSurfaces:
    surfaceParse(afs, centroid, calcExtrema=True)


# todo: search for all edges & surfaces related to point features (maybe symmetrical arcs too) get nearest points
# to originating points
# note that one can expect C2 continuity in STEP format, meaning that two point maxima will not be adjacent except
# where the surface or edge is of a fixed radius from the local centroid.


def parseAdjacentFeatures(AFSobj):  # getNearestFeatures(sourceRef, parsedFaceSet):
    """ given a STEP reference associated with a feature, find all associated objects (edges & surfaces) and return  """
    # seems there should be sets of associated references with every adjacent feature??
    # only vertex points have associated Ref

    # create structures that only represent adjacent features-
    # max/min @ vertex: return edges w/ second vertex for all edge set, surfaces that touch vertex
    # max/min @ edge: return vertex1 & 2 plus touching surfaces
    # max/min @ surface: return surface & its edges & vertices

    for afs_index, afs in enumerate(AFSobj):
        for se_index, se in enumerate(afs['outerBoundEdgeLoop']):
            for sRef in ['vertex1ref', 'vertex2ref']:
                if se.get(sRef):
                    ref1 = se[sRef]
                    # indexDictName = sRef + '_' + ref1 + '_featuresIndex'
                    indexDictName = sRef + '_surroundIndex'
                    se[indexDictName] = {}
                    # se[indexDictName][ref1] = (afs_index, se_index)================================================================
                    for afs2_index, afs2 in enumerate(AFSobj):
                        for te_index, te in enumerate(afs2['outerBoundEdgeLoop']):
                            te_strings = [et for et in te.values() if type(et) == str]
                            if ref1 in te_strings:  #
                                # add adjacent surfaces unnecessary as these are the first field of feature mapping
                                # for tRef in ['vertex1ref', 'vertex2ref']:#, 'edgeRef']:
                                for tRef in ['edgeRef']:
                                    if te.get(tRef):
                                        ref2 = te[tRef]
                                        if (ref2 != ref1):  # and ref2 not in edgeSource[indexDictName].keys():
                                            se[indexDictName][ref2] = (afs2_index, te_index)


    # instead of all purpose function, need to find sets of edges at a common vertex
    # (alongside set of vertices on all edges common to any one vertex)
    # sets of surfaces at a common vertex (do not exclude common ref (s, e) in search)
    # intersection of two surface sets to find common surfaces adjacent to common edge


def getRotSymFeatures(AFSobj):
    '''
    order rotational features to return unique rotSym features
    '''

    # rank rotationally symmetric edges as grooves/spherical/ridges based on adjoining surfaces
    # namely for the 2 surfaces adjoining an arc/circle, measure the displacement of the surface relative to this axis of rotation
    # project each point on adjoining surfaces to rotational axis to determine average displacement
    # Limit to adjoining surfaces to determine the extent of relevant points and edges on each side
    # (not check edges completely enclose rotational axis)
    # average the displacement to this axis of vertices & edges

    # test models: Maltese Cruciform
    # import FreeCAD as App
    # import Part
    # v = App.Vector
    # doc = App.newDocument()
    #
    # p0 = v(0,0,-20)
    # p1 = v(10,0,-20)
    # p2 = v(8,0,-10)
    # p3 = v(5,0,0)
    # p4 = v(8,0,10)
    # p5 = v(10,0,20)
    # p6 = v(0,0,20)
    #
    # L0=Part.makePolygon([p0, p1, p2, p3, p4, p5, p6])
    # S0=L0.revolve(App.Vector(0,0,0),App.Vector(0,0,1),360)
    # # Part.show(S0, "Surface")
    #
    # p7 = v(-20,0,0)
    # p8 = v(-20,10,0)
    # p9 = v(-10,8,0)
    # p10 = v(0,5,0)
    # p11 = v(10,8,0)
    # p12 = v(20,10,0)
    # p13 = v(20,0,0)
    #
    # L1=Part.makePolygon([p7, p8, p9, p10, p11, p12, p13])
    # S1=L1.revolve(App.Vector(0,0,0),App.Vector(1,0,0),360)
    # # Part.show(S1, "Surface")
    #
    # fusedX = S0.fuse(S1)
    # Part.show(fusedX, "Surface")

    # rather than searching for all features around a rotSym feature, restrict to rotSym features
    # recall that distinction between geometric surfaces is determinined by the most simple shapes permissible.
    # note that features must be tested to determine they are on the same axis.
    # multiple data sets created for multiple axes through centroid (e.g. sea-urchin shape, or throwing jacks)

    # import pickle
    # p_file = open('AFS_cruciform.pkl', 'wb')
    # pickle.dump(AFSobj, p_file)
    # p_file.close()

    # f = open('AFS_cruciform.pkl','rb')
    # AFSobj = pickle.load(f)

    allRotSymFeatures = dict()
    allRotSymFeatures['rotSymCentre'] = []
    allRotSymFeatures['rotSymRadius'] = []
    allRotSymFeatures['normDir'] = []

    for afs_index, afs in enumerate(AFSobj): # todo ------------------------------------------------afs_index not used?
        # get surface rotSymFeatures
        for edge in afs['outerBoundEdgeLoop']:
            if (edge.get('rotSymFeature') is not None):
                if (edge['rotSymFeature'].get('rotSymCentre') is not None) and \
                        (edge['rotSymFeature'].get('rotSymRadius') is not None):
                    if type(edge['rotSymFeature']['rotSymCentre']) == list:
                        for rsf in range(0, len(edge['rotSymFeature']['rotSymCentre'])):
                            allRotSymFeatures['rotSymCentre'].append(edge['rotSymFeature']['rotSymCentre'][rsf])
                            allRotSymFeatures['rotSymRadius'].append(edge['rotSymFeature']['rotSymRadius'][rsf])
                            allRotSymFeatures['normDir'].append(edge['normDir'])
                    else:
                        allRotSymFeatures['rotSymCentre'].append(edge['rotSymFeature']['rotSymCentre'])
                        allRotSymFeatures['rotSymRadius'].append(edge['rotSymFeature']['rotSymRadius'])
                        allRotSymFeatures['normDir'].append(edge['normDir'])

        # get surface rotSymFeatures
        if (afs['ParsedSurface'].get('rotSymFeature') is not None):
            if (afs['ParsedSurface']['rotSymFeature'].get('rotSymCentre') is not None) and \
                    (afs['ParsedSurface']['rotSymFeature'].get('rotSymRadius') is not None):
                if type(afs['ParsedSurface']['rotSymFeature']['rotSymCentre']) == list:
                    for rsf in range(0, len(afs['ParsedSurface']['rotSymFeature']['rotSymCentre'])):
                        allRotSymFeatures['rotSymCentre'].append(
                            afs['ParsedSurface']['rotSymFeature']['rotSymCentre'][rsf])
                        allRotSymFeatures['rotSymRadius'].append(
                            afs['ParsedSurface']['rotSymFeature']['rotSymRadius'][rsf])
                        allRotSymFeatures['normDir'].append(afs['ParsedSurface']['normDir'])
                else:
                    allRotSymFeatures['rotSymCentre'].append(afs['ParsedSurface']['rotSymFeature']['rotSymCentre'])
                    allRotSymFeatures['rotSymRadius'].append(afs['ParsedSurface']['rotSymFeature']['rotSymRadius'])
                    allRotSymFeatures['normDir'].append(afs['ParsedSurface']['normDir'])

    # get unique axis set (ignore reverse directions)
    if len(allRotSymFeatures['normDir']) < 1:
        return []

    normDirSet = [allRotSymFeatures['normDir'][0]]
    for arfsc in range(0, len(allRotSymFeatures['rotSymCentre'])):
        if not any([np.allclose(np.abs(allRotSymFeatures['normDir'][arfsc]), np.abs(nds), eps) for nds in normDirSet]):
            normDirSet.append(allRotSymFeatures['normDir'][arfsc])

    superAxisUniqueSet = []
    for ndc in normDirSet:
        axisUniqueSet = []
        for arfsc in range(0, len(allRotSymFeatures['rotSymCentre'])):
            if np.allclose(np.abs(allRotSymFeatures['normDir'][arfsc]), np.abs(ndc), eps):
                axisUniqueSet.append(arfsc)
        superAxisUniqueSet.append(axisUniqueSet)

    superAxesSet = []
    # need to find adjacent rotSym to establish maxima/minima
    # for sausIndex, saus in enumerate(superAxisUniqueSet):
    for saus in superAxisUniqueSet:
        axisSet = dict()
        # need to determine which rotSymCentre lie on the same side of centroid
        # 1. find point maxCentroidDisp with max disp from centroid
        # 2. sort all other points relative to maxDispCentroid

        # define centroid plane, normDir and centroid point, then find dot product
        # centreDisps = [pointProjectAxis(allRotSymFeatures['rotSymCentre'][s], centroid, normDirSet[sausIndex]) for s in saus]

        centreDisps = [np.linalg.norm(allRotSymFeatures['rotSymCentre'][s] - centroid) for s in saus]
        maxCentreDisps = saus[centreDisps.index(max(centreDisps))]
        axisDisps = [
            np.linalg.norm(allRotSymFeatures['rotSymCentre'][s] - allRotSymFeatures['rotSymCentre'][maxCentreDisps]) for
            s in saus]
        axisDispsSorted = np.argsort(axisDisps)

        saus = [saus[ads] for ads in axisDispsSorted]

        axisSet['centreCentroidDisp'] = [centreDisps[ads].item() for ads in axisDispsSorted]
        axisSet['rotSymCentre'] = [allRotSymFeatures['rotSymCentre'][s] for s in saus]
        axisSet['rotSymRadius'] = [allRotSymFeatures['rotSymRadius'][s] for s in saus]
        axisSet['rotSymNorm'] = [allRotSymFeatures['normDir'][s] for s in saus]
        # check neighbours for max/min determination (equality := min)
        # furthest rotSym centre is always a maxima ? ----------------------------------------------------------
        superAxesSet.append(axisSet)

    superAxesSet2 = []
    # remove duplicates
    for sas in superAxesSet:
        axisSet2 = dict()
        axisSet2['rotSymCentre'] = [sas['rotSymCentre'][0], ]
        axisSet2['rotSymRadius'] = [sas['rotSymRadius'][0], ]
        axisSet2['rotSymNorm'] = [sas['rotSymNorm'][0], ]
        axisSet2['centreCentroidDisp'] = [sas['centreCentroidDisp'][0]]
        for s in range(1, len(sas['rotSymCentre'])):
            if not (np.allclose(sas['rotSymCentre'][s], sas['rotSymCentre'][s - 1], eps) and \
                    np.isclose(sas['rotSymRadius'][s], sas['rotSymRadius'][s - 1], eps)):
                axisSet2['rotSymCentre'].append(sas['rotSymCentre'][s])
                axisSet2['rotSymRadius'].append(sas['rotSymRadius'][s])
                axisSet2['rotSymNorm'].append(sas['rotSymNorm'][s])
                axisSet2['centreCentroidDisp'].append(sas['centreCentroidDisp'][s])
        superAxesSet2.append(axisSet2)

    superAxesSet = superAxesSet2

    # get individual rotSym max/min values from relative radii
    # can assume outlying rotSym features are maxima, as face minima are points - adopt as convention?
    for sas in superAxesSet:
        sas['maxima'] = [False, ] * len(sas['rotSymCentre'])
        sas['minima'] = [False, ] * len(sas['rotSymCentre'])
        sas['maxima'][0] = True
        # sas['minima'][0] = False
        for s in range(1, len(sas['rotSymCentre']) - 1):
            if (sas['rotSymRadius'][s] > sas['rotSymRadius'][s - 1]) and (
                    sas['rotSymRadius'][s] > sas['rotSymRadius'][s + 1]):
                sas['maxima'][s] = True
                sas['minima'][s] = False
            elif (sas['rotSymRadius'][s] <= sas['rotSymRadius'][s - 1]) and (
                    sas['rotSymRadius'][s] <= sas['rotSymRadius'][s + 1]):
                sas['maxima'][s] = False
                sas['minima'][s] = True
        sas['maxima'][-1] = True
        # sas['minima'][-1] = False

    return superAxesSet


# for every local minima or maxima feature point, find the surrounding local minima points by way of geometrical features
# (i.e. opposing minima on a thin plate don't count; adjacency is not calculated by cartesian distance but surface topography)

# create similar dict structure listing local minima/maxima tuple (index1, index2) key and index list of surrounding local maxima/minima
# ignore overall minima within a surface and edge boundaries

# surface maxima/minima may be associated with related edge maxima/minima (note B-spline surface may have several max/min)
# find the minima/maxima associated with every edge


def getPointFeatures(AFSobj):
    '''
    order point features to return unique points and characteristics
    '''

    def vertexExtremaSearch(AFSobj, extremaType='maxima'):

        pointFeatures = dict()
        pointFeatures['xyz'] = []
        pointFeatures['centroidDisp'] = []
        pointFeatures['maxima'] = []
        pointFeatures['minima'] = []

        vertexX = ['vertex1', 'vertex2']
        vertexXref = ['vertex1ref', 'vertex2ref']
        vertexXref_surroundIndex = ['vertex1ref_surroundIndex', 'vertex2ref_surroundIndex']

        extrema = dict()
        for v in range(0, 2):  # vertex1, vertex2
            for afs_index, afs in enumerate(AFSobj): #----------------------------------afs_index unused?
                for edge in afs['outerBoundEdgeLoop']:
                    if (edge.get('superTypeNameFlag') is not None):
                        if edge['superTypeName'] == 'SEAM_CURVE':  # exclude SEAM_CURVE data
                            print("is this code reachable??")
                            break
                            # unnecessary?

                    elif not edge.get('pointFeature'):
                        # commence with searching local minima around edges associated with AFS surface
                        # for each edge, get 'pointFeatures' fields
                        print("missing pointFeature field")  # this should be checked elsewhere
                        break

                    elif ((v == 0) and (edge[vertexXref[v]] not in extrema.keys())) or (
                            (v == 1) and (edge[vertexXref[v]] not in extrema.keys())):
                        # get localMinima for all vertex and attached edges,
                        # retain extrema value in separate set, [surface]['outerBoundEdgeLoop'][edge]['pointFeature'][index]

                        # for each pointFeature element, get max/min
                        # for vertex1/2 find if relative max/min from adjacent pointFeature element - ignore surfaces for now

                        # u == 0, vertex1, get attached pointFeature
                        # surroundIndex should give address of adjacent edges
                        pfSurround = []
                        adjacentEdgeRefs = []
                        reversedOrder = []
                        for vseRef in edge[vertexXref_surroundIndex[v]].keys():
                            adjacentEdgeRefs.append(vseRef)
                            vse = edge[vertexXref_surroundIndex[v]][vseRef]
                            adjPF = AFSobj[vse[0]]['outerBoundEdgeLoop'][vse[1]]['pointFeature']
                            if AFSobj[vse[0]]['outerBoundEdgeLoop'][vse[1]][vertexXref[v]] == edge[vertexXref[v]]:
                                pfSurround.append(adjPF)
                                reversedOrder.append(False)
                            else:
                                reversePF = dict()
                                reversePF['u'] = adjPF['u'].copy()
                                reversePF['u'].reverse()
                                reversePF['xyz'] = adjPF['xyz'].copy()
                                # if not isinstance(reversePF['xyz'], list):
                                #     _1=1
                                reversePF['xyz'].reverse()
                                reversePF['centroidDisp'] = adjPF['centroidDisp'].copy()
                                reversePF['centroidDisp'].reverse()
                                pfSurround.append(reversePF)
                                reversedOrder.append(True)

                        # test vertexX with immediate neighbours
                        if extremaType == 'maxima':
                            if all([pfs['centroidDisp'][0] > pfs['centroidDisp'][1] + eps_STEP_AP21 for pfs in
                                    pfSurround]):
                                if not array3x1match(edge[vertexX[v]], pointFeatures['xyz']):
                                    pointFeatures['xyz'].append(edge[vertexX[v]].tolist())
                                    pointFeatures['centroidDisp'].append(pfSurround[0]['centroidDisp'][0])
                                    pointFeatures['maxima'].append(True)
                                    pointFeatures['minima'].append(False)

                        elif extremaType == 'minima':
                            if all([pfs['centroidDisp'][0] < pfs['centroidDisp'][1] - eps_STEP_AP21 for pfs in
                                    pfSurround]):  # local min
                                if not array3x1match(edge[vertexX[v]], pointFeatures['xyz']):
                                    pointFeatures['xyz'].append(edge[vertexX[v]].tolist())
                                    pointFeatures['centroidDisp'].append(pfSurround[0]['centroidDisp'][0])
                                    pointFeatures['maxima'].append(False)
                                    pointFeatures['minima'].append(True)

                        maxPFlen = max([len(pfs['u']) for pfs in pfSurround]) - 1
                        if maxPFlen > 1:
                            crawlindex = 1
                            while crawlindex < maxPFlen:
                                # should capture all local max/min before second vertex
                                for pfsIndex, pfs in enumerate(pfSurround):
                                    # 'centroidDisp' find local maxima/minima
                                    if crawlindex < len( pfs['centroidDisp']) - 1:
                                        pfcd0 = pfs['centroidDisp'][crawlindex - 1]
                                        pfcd1 = pfs['centroidDisp'][crawlindex]
                                        pfcd2 = pfs['centroidDisp'][crawlindex + 1]

                                        if reversedOrder[pfsIndex]:
                                            u = len(pfs['u']) - crawlindex - 1  # CHECK
                                        else:
                                            u = crawlindex

                                        if extremaType == 'maxima':  # local max
                                            if (pfcd1 > pfcd0 + eps_STEP_AP21) and (pfcd1 > pfcd2 + eps_STEP_AP21):
                                                if not array3x1match(pfs['xyz'][u], pointFeatures['xyz']):
                                                    pointFeatures['xyz'].append(pfs['xyz'][u].tolist())
                                                    pointFeatures['centroidDisp'].append(pfcd1)
                                                    pointFeatures['maxima'].append(True)
                                                    pointFeatures['minima'].append(False)

                                        elif extremaType == 'minima':  # local min
                                            if (pfcd1 < pfcd0 - eps_STEP_AP21) and (pfcd1 < pfcd2 - eps_STEP_AP21):
                                                if not array3x1match(pfs['xyz'][u], pointFeatures['xyz']):
                                                    pointFeatures['xyz'].append(pfs['xyz'][u].tolist())
                                                    pointFeatures['centroidDisp'].append(pfcd1)
                                                    pointFeatures['maxima'].append(False)
                                                    pointFeatures['minima'].append(True)

                                crawlindex += 1
        return pointFeatures

    maxPointFeatures = vertexExtremaSearch(AFSobj, extremaType='maxima')
    minPointFeatures = vertexExtremaSearch(AFSobj, extremaType='minima')

    # surfacePointFeature = maxPointFeatures | minPointFeatures # > Python 3.9, requires explicit merge function prior to 3.4
    surfacePointFeature = {**maxPointFeatures, **minPointFeatures}

    # should add surface point instances to min/maxPointFeatures here?

    # collate maxima minima extrema within all surfaces, note that local max/min implies that no global hierarchy is required
    surfacePointFeature = dict()  # todo: this should become an initiated class
    surfacePointFeature['xyz'] = []
    surfacePointFeature['centroidDisp'] = []
    surfacePointFeature['maxima'] = []
    surfacePointFeature['minima'] = []

    for mpf in [maxPointFeatures, minPointFeatures]:
        surfacePointFeature['xyz'] = surfacePointFeature['xyz'] + mpf['xyz']
        surfacePointFeature['centroidDisp'] = surfacePointFeature['centroidDisp'] + mpf['centroidDisp']
        surfacePointFeature['maxima'] = surfacePointFeature['maxima'] + mpf['maxima']
        surfacePointFeature['minima'] = surfacePointFeature['minima'] + mpf['minima']

    surfaceSphericalFeature = dict()
    surfaceSphericalFeature['rotSymCentre'] = []
    surfaceSphericalFeature['rotSymRadius'] = []

    for afs in AFSobj:
        if (afs.get('ParsedSurface') is not None):
            # ['ParsedSurface']['pointFeature'] only covers internal features
            for index_xyz, xyz in enumerate(afs['ParsedSurface']['pointFeature']['xyz']):
                if not array3x1match(xyz, surfacePointFeature['xyz']):
                    surfacePointFeature['xyz'].append(xyz)
                    surfacePointFeature['centroidDisp'].append(
                        afs['ParsedSurface']['pointFeature']['centroidDisp'][index_xyz])
                    surfacePointFeature['maxima'].append(afs['ParsedSurface']['pointFeature']['maxima'][index_xyz])
                    surfacePointFeature['minima'].append(afs['ParsedSurface']['pointFeature']['minima'][index_xyz])

                # for OBE in afs['outerBoundEdgeLoop']:
                #     for index_OBExyz, OBExyz in enumerate(OBE['pointFeature']['xyz']):
                #         if not array3x1match(OBExyz, surfacePointFeature['xyz']):
                #             surfacePointFeature['xyz'].append(OBExyz)
                #             surfacePointFeature['centroidDisp'].append(OBE['pointFeature']['centroidDisp'][index_OBExyz])
                #             surfacePointFeature['maxima'].append(OBE['pointFeature']['maxima'][index_OBExyz])
                #             surfacePointFeature['minima'].append(OBE['pointFeature']['minima'][index_OBExyz])

            if afs['SurfaceTypeName'] == 'SPHERICAL_SURFACE':
                if afs['ParsedSurface']['sphereFeatureFlag']:
                    surfaceSphericalFeature['rotSymCentre'] = afs['ParsedSurface']['pointFeature']['rotSymFeature'][
                        'rotSymCentre']
                    surfaceSphericalFeature['rotSymRadius'] = afs['ParsedSurface']['pointFeature']['rotSymFeature'][
                        'rotSymRadius']

    return surfacePointFeature, surfaceSphericalFeature


def translateShapeMatchFormat():
    # translate to shapeMatch friendly categories
    # based on data structure within STEPaffine, viz
    #         self.featureMaxPoints = []
    #         self.featureMinPoints = []  # arrange as a list of lists representing points
    #         self.rotSymRidgePoints = []
    #         self.rotSymGroovePoints = []
    #         self.featureMaxCurveDisps = []
    #         self.featureMaxCentres = []
    #         self.featureMinCurveDisps = []
    #         self.featureMinCentres = []
    #         self.featureSphereDisps = []
    #         self.spherePoints = []
    #         self.centrePoints = []

    shapeMatchFormat = dict()  # todo: this should become an initiated class
    shapeMatchFormat['featureMaxPoints'] = []
    shapeMatchFormat['featureMinPoints'] = []

    for spf_i, spf in enumerate(surfacePointFeature['xyz']):
        if surfacePointFeature['maxima'][spf_i]:
            if not array3x1match(spf, shapeMatchFormat['featureMaxPoints'], eps_STEP_AP21):
                shapeMatchFormat['featureMaxPoints'].append(spf)
        if surfacePointFeature['minima'][spf_i]:
            if not array3x1match(spf, shapeMatchFormat['featureMinPoints'], eps_STEP_AP21):
                shapeMatchFormat['featureMinPoints'].append(spf)

    shapeMatchFormat['featureSphereDisps'] = []
    shapeMatchFormat['spherePoints'] = []

    for ssf_i, ssf in enumerate(surfaceSphericalFeature['rotSymCentre']):
        if not array3x1match(ssf, shapeMatchFormat['spherePoints'], eps_STEP_AP21):
            shapeMatchFormat['spherePoints'].append(ssf)
            shapeMatchFormat['featureSphereDisps'].append(surfaceSphericalFeature['rotSymRadius'][ssf_i])

    shapeMatchFormat['rotSymRidgePoints'] = []
    shapeMatchFormat['rotSymGroovePoints'] = []
    shapeMatchFormat['featureMaxCurveDisps'] = []
    shapeMatchFormat['featureMaxCentres'] = []
    shapeMatchFormat['featureMaxCurveNorms'] = []
    shapeMatchFormat['featureMinCurveDisps'] = []
    shapeMatchFormat['featureMinCentres'] = []
    shapeMatchFormat['featureMinCurveNorms'] = []
    shapeMatchFormat['featureSeamCurveDisps'] = [] # neither maxima nor minima, e.g. cone minima radius
    shapeMatchFormat['featureSeamCentres'] = []
    shapeMatchFormat['featureSeamCurveNorms'] = []

    for rotSymAxis in rotSymFeatures: # per axis
        for rsf_i, rsf in enumerate(rotSymAxis['rotSymCentre']): # todo style rewrite with zip()?
            if rotSymAxis['maxima'][rsf_i]:
                if not array3x1match(rsf, shapeMatchFormat['featureMaxCentres'], eps_STEP_AP21):
                    shapeMatchFormat['featureMaxCentres'].append(rsf)
                    shapeMatchFormat['featureMaxCurveDisps'].append(rotSymAxis['rotSymRadius'][rsf_i])
                    shapeMatchFormat['featureMaxCurveNorms'].append(rotSymAxis['rotSymNorm'][rsf_i])
            if rotSymAxis['minima'][rsf_i]:
                if not array3x1match(rsf, shapeMatchFormat['featureMinCentres'], eps_STEP_AP21):
                    shapeMatchFormat['featureMinCentres'].append(rsf)
                    shapeMatchFormat['featureMinCurveDisps'].append(rotSymAxis['rotSymRadius'][rsf_i])
                    shapeMatchFormat['featureMinCurveNorms'].append(rotSymAxis['rotSymNorm'][rsf_i])
            if not rotSymAxis['maxima'][rsf_i] and not rotSymAxis['minima'][rsf_i]:
                if not array3x1match(rsf, shapeMatchFormat['featureMinCentres'], eps_STEP_AP21):
                    shapeMatchFormat['featureSeamCentres'].append(rsf)
                    shapeMatchFormat['featureSeamCurveDisps'].append(rotSymAxis['rotSymRadius'][rsf_i])
                    shapeMatchFormat['featureSeamCurveNorms'].append(rotSymAxis['rotSymNorm'][rsf_i])

    # problem of feature points on edges that form part of surfaces of revolution,
    # e.g. spline minima/maxima not at revolved surface boundary
    # note not consistently marked as "seam_curve", e.g. FreeCAD 7.6 STEP engine

    featureMaxPointValid = np.ones(len(shapeMatchFormat['featureMaxPoints']))
    featureMinPointValid = np.ones(len(shapeMatchFormat['featureMinPoints']))

    for fmaxc, fmaxr, fmaxn in zip(shapeMatchFormat['featureMaxCentres'],
                                   shapeMatchFormat['featureMaxCurveDisps'],
                                   shapeMatchFormat['featureMaxCurveNorms']):

        for ifmaxp, fmaxp in enumerate(shapeMatchFormat['featureMaxPoints']):
            if (np.dot(fmaxp - fmaxc, fmaxn) < eps_STEP_AP21) and \
                    (np.abs(np.linalg.norm(fmaxp - fmaxc) - fmaxr) < eps_STEP_AP21):
                featureMaxPointValid[ifmaxp] = 0

        for ifminp, fminp in enumerate(shapeMatchFormat['featureMinPoints']):
            if (np.dot(fminp - fmaxc, fmaxn) < eps_STEP_AP21) and \
                    (np.abs(np.linalg.norm(fminp - fmaxc) - fmaxr) < eps_STEP_AP21):
                featureMinPointValid[ifminp] = 0

    for fminc, fminr, fminn in zip(shapeMatchFormat['featureMinCentres'],
                                   shapeMatchFormat['featureMinCurveDisps'],
                                   shapeMatchFormat['featureMinCurveNorms']):

        for ifmaxp, fmaxp in enumerate(shapeMatchFormat['featureMaxPoints']):
            if (np.dot(fmaxp - fminc, fminn) < eps_STEP_AP21) and \
                    (np.abs(np.linalg.norm(fmaxp - fminc) - fminr) < eps_STEP_AP21):
                featureMaxPointValid[ifmaxp] = 0

        for ifminp, fminp in enumerate(shapeMatchFormat['featureMinPoints']):
            # dot product orthogonal vectors == 0
            if (np.abs(np.dot(fminp - fminc, fminn)) < eps_STEP_AP21) and \
                    (np.abs(np.linalg.norm(fminp - fminc) - fminr) < eps_STEP_AP21):
                featureMinPointValid[ifminp] = 0

    for fseamc, fseamr, fseamn in zip(shapeMatchFormat['featureSeamCentres'],
                                   shapeMatchFormat['featureSeamCurveDisps'],
                                   shapeMatchFormat['featureSeamCurveNorms']):

        for ifmaxp, fmaxp in enumerate(shapeMatchFormat['featureMaxPoints']):
            if (np.dot(fmaxp - fseamc, fseamn) < eps_STEP_AP21) and \
                    (np.abs(np.linalg.norm(fmaxp - fseamc) - fseamr) < eps_STEP_AP21):
                featureMaxPointValid[ifmaxp] = 0

        for ifminp, fminp in enumerate(shapeMatchFormat['featureMinPoints']):
            if (np.dot(fminp - fseamc, fseamn) < eps_STEP_AP21) and \
                (np.abs(np.linalg.norm(fminp - fseamc) - fseamr) < eps_STEP_AP21):
                featureMinPointValid[ifminp] = 0

    shapeMatchFormat['featureMaxPoints'] = [fmaxp for fmaxp, fmaxpv in zip(shapeMatchFormat['featureMaxPoints'], featureMaxPointValid) if fmaxpv]
    shapeMatchFormat['featureMinPoints'] = [fminp for fminp, fminpv in zip(shapeMatchFormat['featureMinPoints'], featureMinPointValid) if fminpv]

    return shapeMatchFormat


# todo: looking for 'superTypeName', e.g TiltedCylinder.step, eventually delete
for afs in AdvancedFaceSurfaces:
    for edge in afs['outerBoundEdgeLoop']:
        if (edge.get('superTypeNameFlag') is not None):
            if edge['superTypeName'] == 'SEAM_CURVE':  # exclude SEAM_CURVE data
                print("is this code reachable??")


parseAdjacentFeatures(AdvancedFaceSurfaces)

# Note that rotSymFeatures returns fields consistent with rotationally symmetrical objects along consistent axes.
# Relative maxima/minima would be lost if rotational features divided between individual surface
rotSymFeatures = getRotSymFeatures(AdvancedFaceSurfaces)

surfacePointFeature, surfaceSphericalFeature = getPointFeatures(AdvancedFaceSurfaces)

shapeMatchFeatures = translateShapeMatchFormat()
shapeMatchFeatures['centroid'] = centroid

#     def __init__(self):
#         self.objectHandle = (
#             ""  # pointer or string used to reference model in CAD model-space
#         )
#         self.filepath = ""
#         self.name = ""
#         self.generationTime = 0
#         self.surfaceStatus = ""
#         self.insertionPoint = Point(0.0, 0.0, 0.0)
#         self.rotation = 0.0
#         self.rotationAxis = None
#         self.rotationMatrix = None
#         self.scale = 1.0

#         self.featureMaxPoints = []
#         self.featureMinPoints = []  # arrange as a list of lists representing points

#         self.surfacePoints = []

#         self.centroid = None
#         self.rotSymRidgePoints = []
#         self.rotSymGroovePoints = []

#         self.featureMaxCurveDisps = []
#         self.featureMaxCentres = []

#         self.featureMinCurveDisps = []
#         self.featureMinCentres = []

#         self.featureSphereDisps = []
#         self.spherePoints = []

#         self.centrePoints = []

# export to shapeMatch format



pprint(shapeMatchFeatures)

#pass  # stop here for now
_1=1

# The standard CSG primitives are the cone, eccentric_cone, cylinder, sphere, torus,
# block, right_angular_wedge, ellipsoid, tetrahedron and pyramid.

#        b_spline_surface,
#        b_spline_surface_with_knots,
#        bezier_surface,
#        conical_surface,

###        curve_bounded_surface,

#        cylindrical_surface,
###        degenerate_toroidal_surface,
###        offset_surface,
###        quasi_uniform_surface,
#        rational_b_spline_surface,
###        rectangular_composite_surface,
#        rectangular_trimmed_surface,
#        spherical_surface,
#        surface,
###        surface_of_linear_extrusion,
#        surface_of_revolution,
###        surface_replica,
###        swept_surface,
#        toroidal_surface,
###        uniform_surface,

#
# class STEPmodel(object):
#     """class and methods to encapsulate objects and methods unique to FreeCAD environment"""
#
#     def __init__(self):
#         self.objectHandle = (
#             ""  # pointer or string used to reference model in CAD model-space
#         )
#         self.filepath = ""
#         self.name = ""
#         self.generationTime = 0
#         self.surfaceStatus = ""
#         self.insertionPoint = Point(0.0, 0.0, 0.0)
#         self.rotation = 0.0
#         self.rotationAxis = None
#         self.rotationMatrix = None
#         self.scale = 1.0
#         self.featureMaxPoints = []
#         self.featureMinPoints = []  # arrange as a list of lists representing points
#         self.surfacePoints = []
#         self.centroid = None
#         self.rotSymRidgePoints = []
#         self.rotSymGroovePoints = []
#         self.featureMaxCurveDisps = []
#         self.featureMaxCentres = []
#         self.featureMinCurveDisps = []
#         self.featureMinCentres = []
#         self.featureSphereDisps = []
#         self.spherePoints = []
#         self.centrePoints = []
#
#     # noinspection PyTypeChecker
#     def importSTEPmodel(
#         self,
#         filepath="",
#         insertionPoint=Point(0.0, 0.0, 0.0),
#         scale=1,
#         rotation=0.0,
#         rotationAxis=[],
#     ):
#         """
#         FreeCAD dependent
#
#         :param filepath: string
#         :param insertionPoint: {x,y,z} Cartesian tuple, not equivalent to centroid
#         :param scale: floating point scalar
#         :param rotation: floating point scalar for Rhino degree value (0.0 <= rotation <= 360.0)
#         :param rotationAxis:
#         """
#
#         # get default values from existing object or replace default initialisation values with given parameters
#         if (filepath == "") and not (self.filepath == ""):
#             filepath = self.filepath
#         elif self.filepath == "":
#             self.filepath = filepath
#
#         if (insertionPoint == Point(0.0, 0.0, 0.0)) and not (
#             self.insertionPoint == Point(0.0, 0.0, 0.0)
#         ):
#             insertionPoint = self.insertionPoint
#         elif self.insertionPoint == Point(0.0, 0.0, 0.0):
#             self.insertionPoint = insertionPoint
#
#         if (scale == 1) and not (self.scale == 1):
#             scale = self.scale
#         elif self.scale == 1:
#             self.scale = scale
#
#         if (rotation == 0.0) and not (self.rotation == 0.0):
#             rotation = self.rotation
#         elif self.rotation == 0.0:
#             self.rotation = rotation
#
#         if (rotationAxis == []) and not (self.rotationAxis is None):
#             rotationAxis = self.rotationAxis
#         elif self.rotationAxis is None:
#             self.rotationAxis = rotationAxis
#
#         if self.name == "":
#             self.name = os.path.basename(filepath)
#
#         # def angleAxisTranslate2AffineMat(rotAxis, theta, disp):
#         #     """
#         #     Affine transformation matrix from 3x3 rotation matrix and displacement vector.
#         #     3x3 rotation matrix derived from CCW angle theta around axis
#         #     From: http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle & transform3d
#         #
#         #     :param rotAxis: axis of rotation vector axis emanating from origin
#         #     :param theta: scalar angle of rotation
#         #     :param disp: displacement from origin point, presumed (0,0,0)
#         #     :return: 4x4 matrix representing affine transformation
#         #     """
#         #
#         #     s = sin(theta); c = cos(theta); C = (1 - c)
#         #
#         #     xs = rotAxis[0] * s; ys = rotAxis[1] * s; zs = rotAxis[2] * s
#         #
#         #     xC = rotAxis[0] * C; yC = rotAxis[1] * C; zC = rotAxis[2] * C
#         #
#         #     xyC = rotAxis[0] * yC; yzC = rotAxis[1] * zC; zxC = rotAxis[2] * xC
#         #
#         #     xxC = rotAxis[0] * xC; yyC = rotAxis[1] * yC; zzC = rotAxis[2] * zC
#         #
#         #     return F_app.Matrix(xxC + c, xyC - zs, zxC + ys, disp[0],
#         #                         xyC + zs, yyC + c, yzC - xs, disp[1],
#         #                         zxC - ys, yzC + xs, zzC + c, disp[2],
#         #                         0,        0,        0,       1      )
#
#         try:
#             # filepath sanity check
#             checkPath = glob.glob(filepath)
#             if len(checkPath) != 1:
#                 raise Exception("pathname failure")
#             else:
#                 self.objectHandle = Part_app.read(filepath)
#
#             if not self.objectHandle.isValid():
#                 raise Exception("imported STEP file fails geometry checks")
#
#             if self.objectHandle.check():
#                 print("imported STEP file passes detailed geometry shape.check()")
#
#             if not self.objectHandle.isClosed():
#                 raise Exception("imported STEP file fails closed geometry checks")
#
#             # uniform scaling
#             if scale - 1.0 >= eps:
#                 scalingMatrix = F_app.Matrix()
#                 scalingMatrix.scale(scale, scale, scale)
#                 self.objectHandle = self.objectHandle.transformGeometry(scalingMatrix)
#
#             # shape rotation
#             # if math.fabs(rotation) > eps:
#             if not rotationAxis:
#                 # rotationAxis = F_app.Vector(0.0, 0.0, 1.0)  # Z-axis
#                 rotationAxis = [0.0, 0.0, 0.0]
#
#             self.rotationMatrix = angleAxis2RotMat(rotationAxis, math.radians(rotation))
#             self.objectHandle.Placement.Matrix = F_app.Matrix(
#                 self.rotationMatrix[0][0],
#                 self.rotationMatrix[0][1],
#                 self.rotationMatrix[0][2],
#                 self.insertionPoint.x,
#                 self.rotationMatrix[1][0],
#                 self.rotationMatrix[1][1],
#                 self.rotationMatrix[1][2],
#                 self.insertionPoint.y,
#                 self.rotationMatrix[2][0],
#                 self.rotationMatrix[2][1],
#                 self.rotationMatrix[2][2],
#                 self.insertionPoint.z,
#                 0,
#                 0,
#                 0,
#                 1,
#             )
#
#         except Exception as e:
#             raise RuntimeError("Model insertion failure")
#
#     def rayModelIntersection(self, V1, V2):
#         """
#         Returns the vertex that is at the intersection of the face
#         and the line that passes through vertex1 and vertex2
#
#         :param V1: origin of ray
#         :param V2: direction of ray
#         :return:
#         """
#
#         # Define a faraway extent that allows FreeCAD to approximate an infinite ray for intersection points.
#         FreeCADEXTENTS = 100  # calculation speed sensitive to this scaling value
#         # not required for rayModelIntersection3()
#
#         # input variable check, vertices must be of Point type
#         for v in [V1, V2]:
#             if v.__class__.__name__ == "Shape":
#                 v = Point(v.Vertexes.Point.x, v.Vertexes.Point.y, v.Vertexes.Point.z)
#
#         # offset V1 to cartesian origin, (0, 0, 0), move V2 accordingly
#         _V2 = Point(V2.x - V1.x, V2.y - V1.y, V2.z - V1.z)
#         lenV = sqrt(_V2.x**2 + _V2.y**2 + _V2.z**2)
#         unitV = Point(_V2.x / lenV, _V2.y / lenV, _V2.z / lenV)
#         farV2 = Point(
#             FreeCADEXTENTS * unitV.x, FreeCADEXTENTS * unitV.y, FreeCADEXTENTS * unitV.z
#         )
#         farV1 = Point(-farV2.x, -farV2.y, -farV2.z)
#         farV2 = Point(farV2.x + V1.x, farV2.y + V1.y, farV2.z + V1.z)
#         farV1 = Point(farV1.x + V1.x, farV1.y + V1.y, farV1.z + V1.z)
#
#         # """
#         # Redefine second vertex as the second vertex projected to the surface
#         # of a large sphere centered at the first vertex.
#         #
#         # :param V1: origin of ray
#         # :param V2: direction of ray, to be moved to far model extents
#         # :param sphereR: radius of far model extents
#         # :return: redefined V2
#         # """
#
#         ray = Part_app.makeLine(
#             (farV1.x, farV1.y, farV1.z), (farV2.x, farV2.y, farV2.z)
#         )
#         intersections = ray.common(self.objectHandle)
#         # if intersections.Vertexes:
#         #     [print(i.Point.x, i.Point.y, i.Point.z) for i in intersections.Vertexes]
#         intersections = [
#             Point(i.Point.x, i.Point.y, i.Point.z) for i in intersections.Vertexes
#         ]
#         # if vertex1 in intersections:
#         #     intersections.remove(vertex1)
#         if intersections:
#             [
#                 self.surfacePoints.append(i)
#                 for i in intersections
#                 if i not in self.surfacePoints
#             ]
#
#         return intersections
#
#     def rayModelIntersection2(self, V1, V2):
#         """
#         Returns the vertex that is at the intersection of the face
#         and the line that passes through vertex1 and vertex2
#
#         :param V1: origin of ray
#         :param V2: direction of ray
#         :return:
#         """
#
#         # Define a faraway extent that allows FreeCAD to approximate an infinite ray for intersection points.
#         FreeCADEXTENTS = 300  # calculation speed sensitive to this scaling value
#         # not required for rayModelIntersection3()
#
#         # input variable check, vertices must be of Point type
#         # for v in [V1, V2]:
#         #     if v.__class__.__name__ == 'Shape':
#         #         v = Point(v.Vertexes.Point.x, v.Vertexes.Point.y, v.Vertexes.Point.z)
#
#         # offset V1 to cartesian origin, (0, 0, 0), move V2 accordingly
#         _V2 = Point(V2.x - V1.x, V2.y - V1.y, V2.z - V1.z)
#         lenV = sqrt(_V2.x**2 + _V2.y**2 + _V2.z**2)
#         unitV = Point(_V2.x / lenV, _V2.y / lenV, _V2.z / lenV)
#         farV2 = Point(
#             FreeCADEXTENTS * unitV.x, FreeCADEXTENTS * unitV.y, FreeCADEXTENTS * unitV.z
#         )
#         # farV1 = Point(-farV2.x, -farV2.y, -farV2.z)
#         farV2 = Point(farV2.x + V1.x, farV2.y + V1.y, farV2.z + V1.z)
#         # farV1 = Point(farV1.x + V1.x, farV1.y + V1.y, farV1.z + V1.z)
#
#         # """
#         # Redefine second vertex as the second vertex projected to the surface
#         # of a large sphere centered at the first vertex.
#         #
#         # :param V1: origin of ray
#         # :param V2: direction of ray, to be moved to far model extents
#         # :param sphereR: radius of far model extents
#         # :return: redefined V2
#         # """
#
#         intersectSet = []
#         ray = Part_app.makeLine((V1.x, V1.y, V1.z), (farV2.x, farV2.y, farV2.z))
#         intersections = ray.common(self.objectHandle)
#         if intersections.Vertexes:
#             for iv in intersections.Vertexes:
#                 if not (
#                     (abs(iv.Point.x - V1.x) < eps)
#                     and (abs(iv.Point.y - V1.y) < eps)
#                     and (abs(iv.Point.z - V1.z) < eps)
#                 ):
#                     if (
#                         (abs(iv.Point.x - farV2.x) < eps)
#                         and (abs(iv.Point.y - farV2.y) < eps)
#                         and (abs(iv.Point.z - farV2.z) < eps)
#                     ):
#                         print(
#                             "WARNING: far model intersection constant likely to have been exceeded by model extents"
#                         )
#                     else:
#                         intersectSet.append(Point(iv.Point.x, iv.Point.y, iv.Point.z))
#
#             [
#                 self.surfacePoints.append(i)
#                 for i in intersectSet
#                 if i not in self.surfacePoints
#             ]
#
#         return intersectSet
#
#         # face.Surface.intersect(edge.Curve)[0][0] should look like
#         # face.Shape.Surface.intersect(edge.Curve)[0][0] if working with objects
#
#     def rayModelIntersection3(self, V1, V2):
#         """
#         Returns the vertex that is at the intersection of the face
#         and the line that passes through vertex1 and vertex2
#
#         :param V1: origin of ray
#         :param V2: direction of ray
#         :return:
#         """
#
#         ray = Part_app.makeLine(V1, V2)
#         # extract the curve from the ray, then test intersections with shape faces
#         allFaceSets = [ray.Curve.intersect(f.Surface) for f in self.objectHandle.Faces]
#         intersectPoints = [x for face in allFaceSets for x in face if x]
#         intersectPoints = [x for face in intersectPoints for x in face]
#         intersectPoints = [Point(p.X, p.Y, p.Z) for p in intersectPoints]
#         [
#             self.surfacePoints.append(i)
#             for i in intersectPoints
#             if i not in self.surfacePoints
#         ]
#         return intersectPoints
#
#     def printModelAlignment(self):
#         """
#         Print insertion point and rotation of model to stdout.
#
#         :return: stdout display
#         """
#         print(
#             "Model insertion point: x: {:f}, y: {:f}, z: {:f}".format(
#                 self.insertionPoint.x, self.insertionPoint.y, self.insertionPoint.z
#             )
#         )
#         print("Model insertion scale: {:f}".format(self.scale))
#         print("Model insertion rotation (single axis?): {:f}".format(self.rotation))
#         if self.rotationAxis:
#             print(
#                 "Model insertion rotation axis: [{:f}, {:f}, {:f}]".format(
#                     self.rotationAxis[0], self.rotationAxis[1], self.rotationAxis[2]
#                 )
#             )
#         else:
#             print("Model insertion rotation axis: []")
#
#     def featureClean(self):
#         """
#         remove spurious NaN values from feature values,
#         (malformed tests)
#         """
#
#         def cleanPointNaN(P):
#             Pclean = [
#                 p for p in P if not any([np.isnan(p.x), np.isnan(p.y), np.isnan(p.z)])
#             ]
#             if len(Pclean) != len(Pclean):
#                 print("NaN cleaned")
#             return Pclean
#
#         def cleanNaN(P):
#             Pclean = [p for p in P if not np.isnan(p)]
#             if len(Pclean) != len(Pclean):
#                 print("NaN cleaned")
#             return Pclean
#
#         # cleanPointNaN = lambda P: [p for p in P if not any([np.isnan(p.x), np.isnan(p.y), np.isnan(p.z)])]
#         # cleanNaN = lambda P: [p for p in P if not np.isnan(p)]
#
#         self.featureMaxPoints = cleanPointNaN(self.featureMaxPoints)
#         self.featureMinPoints = cleanPointNaN(self.featureMinPoints)
#         self.surfacePoints = cleanPointNaN(self.surfacePoints)
#         self.rotSymRidgePoints = cleanPointNaN(self.rotSymRidgePoints)
#         self.rotSymGroovePoints = cleanPointNaN(self.rotSymGroovePoints)
#         self.featureMaxCurveDisps = cleanNaN(self.featureMaxCurveDisps)
#         self.featureMaxCentres = cleanPointNaN(self.featureMaxCentres)
#         self.featureMinCurveDisps = cleanNaN(self.featureMinCurveDisps)
#         self.featureMinCentres = cleanPointNaN(self.featureMinCentres)
#         self.featureSphereDisps = cleanNaN(self.featureSphereDisps)



# *****SURFACE*****

# Swept_surface - swept_curve - Curve

# Swept_surface - Extruded_surface - extrusion_axis - Direction

# Swept_surface - Surface_of_revolution - axis_direction - Direction
#                                       - axis_point - Cartesian_point

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

# Toroidal_surface - position - Axis_placement_3d
#                  - radius - Positive_length_measure
#                  - minor_radius - Positive_length_measure
# ( if toroid is complete at outer/inner radius, increment rut-symm maxima/minima )

# Spherical_surface - position - Axis_placement_3d
#                   - radius - Positive_length_measure
# ( point on surface cotangent to centre and origin point)

# Bounded_surface - B_spline_surface - U_degree - integer
#                                    - V_degree - integer
#                                    - control_points - Cartesian_point
#                                    - U_closed - boolean
#                                    - V_closed - boolean

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


# STEP Part21 decomposition:
# https://www.mbx-if.org/documents/AP203e2_html/AP203e2.htm

# root node, CLOSED_SHELL (OPEN_SHELL), Connected_Face_Set SUPERTYPE OF (ONEOF (Closed_Shell, Open_Shell))

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


# for finding nearest vertices, can use KD-tree, e.g.
# from scipy import spatial # numpy < 1.24
# A = np.random.random((10,3))*100
# pt = [6, 30]  # <-- the point to find
# A[spatial.KDTree(A).query(pt)[1]] # <-- the nearest point
# distance,index = spatial.KDTree(A).query(pt)

# strategy for planar polygonal surfaces is to divide into non-overlapping triangles where corners are polygon vertices
# and test for nearest/furthest point on triangle

# import and build pymesh locally, see https://github.com/PyMesh/PyMesh/blob/main/README.md


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


# def planeMinMaxPoint(p, V, minmax='min'):
#     import pymesh
#     def pointTriangleMinDistance(TRI, P):
#         # return distance between a point and triangle
#         # SYNTAX
#         #   dist = pointTriangleDistance(TRI,P)
#         #   [dist,PP0] = pointTriangleDistance(TRI,P)
#         #
#         # DESCRIPTION
#         #   Calculate the distance of a given point P from a triangle TRI.
#         #   Point P is a row vector of the form 1x3. The triangle is a matrix
#         #   formed by three rows of points TRI = [P1;P2;P3] each of size 1x3.
#         #   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
#         #   to the triangle TRI.
#         #   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
#         #   closest point PP0 to P on the triangle TRI.
#         #
#         # Author: Gwolyn Fischer
#         # Release: 1.0
#         # Release date: 09/02/02
#
#         # The algorithm is based on
#         # "David Eberly, 'Distance Between Point and Triangle in 3D',
#         # Geometric Tools, LLC, (1999)"
#         # http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
#         #
#         #        ^t
#         #  \     |
#         #   \reg2|
#         #    \   |
#         #     \  |
#         #      \ |
#         #       \|
#         #        *P2
#         #        |\
#         #        | \
#         #  reg3  |  \ reg1
#         #        |   \
#         #        |reg0\
#         #        |     \
#         #        |      \ P1
#         # -------*-------*------->s
#         #        |P0      \
#         #  reg4  | reg5    \ reg6
#         # rewrite triangle in normal form
#         B = TRI[0, :]
#         E0 = TRI[1, :] - B
#         # E0 = E0/np.sqrt(sum(E0.^2)); %normalize vector
#         E1 = TRI[2, :] - B
#         # E1 = E1/np.sqrt(sum(E1.^2)); %normalize vector
#         D = B - P
#         a = np.dot(E0, E0)
#         b = np.dot(E0, E1)
#         c = np.dot(E1, E1)
#         d = np.dot(E0, D)
#         e = np.dot(E1, D)
#         f = np.dot(D, D)
#
#         # print "{0} {1} {2} ".format(B,E1,E0)
#         det = a * c - b * b
#         s = b * e - c * d
#         t = b * d - a * e
#
#         # Teribble tree of conditionals to determine in which region of the diagram
#         # shown above the projection of the point into the triangle-plane lies.
#         if (s + t) <= det:
#             if s < 0.0:
#                 if t < 0.0:
#                     # region4
#                     if d < 0:
#                         t = 0.0
#                         if -d >= a:
#                             s = 1.0
#                             sqrDistance = a + 2.0 * d + f
#                         else:
#                             s = -d / a
#                             sqrDistance = d * s + f
#                     else:
#                         s = 0.0
#                         if e >= 0.0:
#                             t = 0.0
#                             sqrDistance = f
#                         else:
#                             if -e >= c:
#                                 t = 1.0
#                                 sqrDistance = c + 2.0 * e + f
#                             else:
#                                 t = -e / c
#                                 sqrDistance = e * t + f
#
#                                 # of region 4
#                 else:
#                     # region 3
#                     s = 0
#                     if e >= 0:
#                         t = 0
#                         sqrDistance = f
#                     else:
#                         if -e >= c:
#                             t = 1
#                             sqrDistance = c + 2.0 * e + f
#                         else:
#                             t = -e / c
#                             sqrDistance = e * t + f
#                             # of region 3
#             else:
#                 if t < 0:
#                     # region 5
#                     t = 0
#                     if d >= 0:
#                         s = 0
#                         sqrDistance = f
#                     else:
#                         if -d >= a:
#                             s = 1
#                             sqrDistance = a + 2.0 * d + f
#                             # GF 20101013 fixed typo d*s ->2*d
#                         else:
#                             s = -d / a
#                             sqrDistance = d * s + f
#                 else:
#                     # region 0
#                     invDet = 1.0 / det
#                     s = s * invDet
#                     t = t * invDet
#                     sqrDistance = (
#                             s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f
#                     )
#         else:
#             if s < 0.0:
#                 # region 2
#                 tmp0 = b + d
#                 tmp1 = c + e
#                 if tmp1 > tmp0:  # minimum on edge s+t=1
#                     numer = tmp1 - tmp0
#                     denom = a - 2.0 * b + c
#                     if numer >= denom:
#                         s = 1.0
#                         t = 0.0
#                         sqrDistance = a + 2.0 * d + f
#                         # GF 20101014 fixed typo 2*b -> 2*d
#                     else:
#                         s = numer / denom
#                         t = 1 - s
#                         sqrDistance = (
#                                 s * (a * s + b * t + 2 * d) + t * (b * s + c * t + 2 * e) + f
#                         )
#
#                 else:  # minimum on edge s=0
#                     s = 0.0
#                     if tmp1 <= 0.0:
#                         t = 1
#                         sqrDistance = c + 2.0 * e + f
#                     else:
#                         if e >= 0.0:
#                             t = 0.0
#                             sqrDistance = f
#                         else:
#                             t = -e / c
#                             sqrDistance = e * t + f
#                             # of region 2
#             else:
#                 if t < 0.0:
#                     # region6
#                     tmp0 = b + e
#                     tmp1 = a + d
#                     if tmp1 > tmp0:
#                         numer = tmp1 - tmp0
#                         denom = a - 2.0 * b + c
#                         if numer >= denom:
#                             t = 1.0
#                             s = 0
#                             sqrDistance = c + 2.0 * e + f
#                         else:
#                             t = numer / denom
#                             s = 1 - t
#                             sqrDistance = (
#                                     s * (a * s + b * t + 2.0 * d)
#                                     + t * (b * s + c * t + 2.0 * e)
#                                     + f
#                             )
#
#                     else:
#                         t = 0.0
#                         if tmp1 <= 0.0:
#                             s = 1
#                             sqrDistance = a + 2.0 * d + f
#                         else:
#                             if d >= 0.0:
#                                 s = 0.0
#                                 sqrDistance = f
#                             else:
#                                 s = -d / a
#                                 sqrDistance = d * s + f
#                 else:
#                     # region 1
#                     numer = c + e - b - d
#                     if numer <= 0:
#                         s = 0.0
#                         t = 1.0
#                         sqrDistance = c + 2.0 * e + f
#                     else:
#                         denom = a - 2.0 * b + c
#                         if numer >= denom:
#                             s = 1.0
#                             t = 0.0
#                             sqrDistance = a + 2.0 * d + f
#                         else:
#                             s = numer / denom
#                             t = 1 - s
#                             sqrDistance = (
#                                     s * (a * s + b * t + 2.0 * d)
#                                     + t * (b * s + c * t + 2.0 * e)
#                                     + f
#                             )
#
#         # account for numerical round-off error
#         if sqrDistance < 0:
#             sqrDistance = 0
#
#         dist = np.sqrt(sqrDistance)
#
#         PP0 = B + s * E0 + t * E1
#         return dist, PP0
#
#     disp = []
#     vertices = []  # np.array([])
#     #p = arrayTypeCast(p)
#     for v in V:
#         #v = arrayTypeCast(v)
#         vertices.append(v)
#         if minmax == 'max':
#             d = np.linalg.norm(v - p)
#             disp.append(d)
#
#     if minmax == 'max':
#         maxDisp = max(disp)
#         maxPoint = vertices[disp.index(maxDisp)]
#         return maxPoint, maxDisp
#
#     tri = pymesh.triangle()
#     tri.points = vertices  # np.array(vertices)
#     tri.split_boundary = False
#     tri.verbosity = 0
#     tri.run()  # Execute triangle
#     minPointIntersect = []
#     for t in tri.faces:
#         tt = np.array([vertices[t[0]], vertices[t[1]], vertices[t[2]]])
#         d, pd = pointTriangleMinDistance(tt, p)
#         disp.append(d)
#         minPointIntersect.append(pd)
#
#     minDisp = min(disp)
#     minPoint = minPointIntersect[disp.index(minDisp)]
#     return minPoint, minDisp


# def checkFilePath(filePath):
#     if os.path.isfile(filePath) and os.access(filePath, os.R_OK):
#         return filePath
#     else:
#         raise Exception("{0} is not a readable filepath".format(filePath))

