import numpy as np
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from geomdl import BSpline
from geomdl import evaluators
from geomdl import NURBS
from geomdl import compatibility
from geomdl import operations
from geomdl import knotvector
from geomdl import utilities
from geomdl import helpers
from geomdl import exchange
# eps machine precision constant
eps = np.finfo(float).eps
# from geomdl.visualization import VisMPL
from geomdl.visualization import VisPlotly

eps_STEP_AP21 = 1e-6  # STEP precision seems to end here

def set_aspect_equal_3d(ax):
    """Fix equal aspect bug for 3D plots."""

    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    from numpy import mean
    xmean = mean(xlim)
    ymean = mean(ylim)
    zmean = mean(zlim)

    plot_radius = max([abs(lim - mean_)
                       for lims, mean_ in ((xlim, xmean),
                                           (ylim, ymean),
                                           (zlim, zmean))
                       for lim in lims])

    ax.set_xlim3d([xmean - plot_radius, xmean + plot_radius])
    ax.set_ylim3d([ymean - plot_radius, ymean + plot_radius])
    ax.set_zlim3d([zmean - plot_radius, zmean + plot_radius])

def radiusCentre3points_3(p1, p2, p3):
    # radius + centre from 3 point circle
    D21x = p2[0] - p1[0] #P2x-P1x
    D21y = p2[1] - p1[1] #P2y-P1y
    D21z = p2[2] - p1[2] #P2z-P1z
    D31x = p3[0] - p1[0] #P3x-P1x
    D31y = p3[1] - p1[1] #P3y-P1y
    D31z = p3[2] - p1[2] #P3z-P1z

    F2 = 1/2*(D21x**2+D21y**2+D21z**2)
    F3 = 1/2*(D31x**2+D31y**2+D31z**2)

    M23xy = D21x*D31y-D21y*D31x
    M23yz = D21y*D31z-D21z*D31y
    M23xz = D21z*D31x-D21x*D31z

    F23x = F2*D31x-F3*D21x
    F23y = F2*D31y-F3*D21y
    F23z = F2*D31z-F3*D21z

    Cx = (M23xy*F23y-M23xz*F23z)/(M23xy**2+M23yz**2+M23xz**2) + p1[0] #P1x
    Cy = (M23yz*F23z-M23xy*F23x)/(M23xy**2+M23yz**2+M23xz**2) + p1[1] #P1y
    Cz = (M23xz*F23x-M23yz*F23y)/(M23xy**2+M23yz**2+M23xz**2) + p1[2] #P1z

    return np.array([Cx, Cy, Cz])

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

def splineToPolyline(C, curveSampleFactor=2):
    numSamples = C.ctrlpts_size * C.degree * curveSampleFactor
    span = (C.knotvector[-1] - C.knotvector[0]) / (numSamples - 1)
    kvs = np.array([C.knotvector[0] + (span * i) for i in range(0, numSamples)])
    pts = C.evaluate_list(kvs)
    return (pts, kvs)


def BsplineCurveExtremaDisp_2(C, p,
                            maxSearch=True,
                            localExtrema=False,
                            curvatureTest=False,
                            uv_xyz=True,
                            eps1=0.0001,
                            eps2=0.0005,
                            deltaFactor=4,
                            eps_bspline = 1E-10):

    # based on Nurbs Book, Piegl & Tiller p.230
    # same idea to find orthogonal tangent, but at maximum/minimum displacement from a centroid point
    # revised to test minima segment candidates using 3-point defined circle parameters

    # def radiusCentre3points2(p1, p2, p3):
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
    # minU = C.knotvector[0]
    # maxU = C.knotvector[-1]

    C.delta = 1 / len(C.knotvector) * deltaFactor

    def closedC():
        return np.isclose(np.linalg.norm(np.array(C.ctrlpts[0]) - np.array(C.ctrlpts[-1])), eps_bspline)

    def localDispTest():
        #localDisp = [0] * len(pts)
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

            #localDisp[i] = dp1

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
            dif = e[0] - p                      # C(u) - p
            c1v = np.linalg.norm(dif)           # |C(u) - p|
            c2n = np.dot(e[1], dif)             # C'(u) * (C(u) - P)
            c2d = np.linalg.norm(e[1]) * c1v    # |C'(u)||C(u) - P|
            c2v = c2n / c2d                     # |C'(u) * (C(u) - P)| / |C'(u)||C(u) - P|
            # c2v = np.dot(e[1], dif) / (np.linalg.norm(e[1]) * np.linalg.norm(dif))

            if (c1v < eps1) and (np.abs(c2v) < eps2):
                return cu
            # ct = n(cu, e, dif)                           #   u* = u - f / f'
            ct = cu - (np.dot(e[1], dif) / (np.dot(e[2], dif) + np.dot(e[1], e[1])))

            #  correct for exceeding bounds
            if ct < C.knotvector[0]:  # [ maxu - ( ct[0] - minu ), ct[1] ]
                if closedC(): ct = C.knotvector[0] + ct % (C.knotvector[-1] - C.knotvector[0])
                #if closedC(): ct = C.knotvector[-1] - (C.knotvector[0] - ct) # ct = C.knotvector[-1] - (ct - C.knotvector[0]) # NURBS book
                else: ct = C.knotvector[0] + eps_bspline

            elif ct > C.knotvector[-1]:
                #if closedC(): ct = C.knotvector[0] + (ct - C.knotvector[-1])
                if closedC(): ct = C.knotvector[-1] - ct % (C.knotvector[-1] - C.knotvector[0])
                else: ct = C.knotvector[-1] - eps_bspline

            c3v = np.linalg.norm(np.multiply(ct - cu, e[1]))
            if c3v < eps1:
                return cu

            cu = ct
            i += 1

        if i == maxits:  # Newton-Raphson fails
            return np.inf
        return cu

    # numSamples = C.ctrlpts_size * C.degree * CSampleFactor
    # span = (C.knotvector[-1] - C.knotvector[0]) / (numSamples - 1)
    # kvs = np.array([C.knotvector[0] + (span * i) for i in range(0, numSamples)])
    # pts = C.evaluate_list(kvs)

    #pts, kvs = splineToPolyline(C, curveSampleFactor=delta) todo 1/C.delta?
    pts, kvs = splineToPolyline(C, curveSampleFactor=deltaFactor)

    # original version fails on points of equal minimal curvature

    # 2-part selection for global maxima using both displacement and curve inflection
    # curve inflection determined from local curvature minima - prob doesn't work for self intersecting curves
    # get set of local minima curvature - use for curve characterisation
    # from set, find minima/maxima disp from point-of-interest to curve point tangent

    # curvature approach fails without C2 continuity at minima/maxima?
    # not much of an issue with STEP spline limitations
    #
    # localDisp = [0] * len(pts)
    # localExtremaU = []
    #
    # # get discrete curvature values and point displacement
    # for i in range(1-closedC(), len(pts) - 1 + closedC()):
    #     if (i == 0):
    #         p0 = np.array(pts[-1])
    #         p1 = np.array(pts[0])
    #         p2 = np.array(pts[1])
    #
    #     elif (i == len(pts)-1):
    #         p0 = np.array(pts[-2])
    #         p1 = np.array(pts[-1])
    #         p2 = np.array(pts[0])
    #
    #     else:
    #         p0 = np.array(pts[i - 1])
    #         p1 = np.array(pts[i])
    #         p2 = np.array(pts[i + 1])
    #
    #     if i==0 or i==1 or i==len(pts)-1:
    #         dp0 = np.linalg.norm(p0 - p)
    #         dp1 = np.linalg.norm(p1 - p)
    #         dp2 = np.linalg.norm(p2 - p)
    #
    #     if (i>1) and (i<len(pts)-1):
    #         dp0 = dp1
    #         dp1 = dp2
    #         dp2 = np.linalg.norm(p2 - p)
    #
    #     if maxSearch:
    #         if (dp1 >= dp2) and (dp1 >= dp0):
    #             localExtremaU.append(kvs[i])
    #     else:  # minS
    #         if (dp1 <= dp2) and (dp1 <= dp0):
    #             # where a point is orthogonal to a planar surface -> sphere radius test, or accept minima
    #             localExtremaU.append(kvs[i])
    #
    #     localDisp[i] = dp1
    #
    #     if curvatureTest:
    #         pur, pucc = radiusCentre3points_2(p0, p1, p2)
    #         if (pur == np.inf):
    #             break
    #         else:
    #             # p is on the same side of the surface as the mean curvature
    #             if np.dot(p - p1, pucc - p1)/np.linalg.norm(pucc - p1) > 0:
    #                 discreteK[i] = pur
    #             else:
    #                 discreteK[i] = -pur

    if curvatureTest:  # discrete grid summing curvatures along U and V axes
        localExtremaU = localCurvatureTest()
    else: # localDispTest()
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

    if (localExtrema and (len(extremaUnique) == 1)) or not localExtrema:  # if there is only one extrema, localextrema is moot
        if uv_xyz: return extremaUnique # return single value in list for compatibility
        else: return [np.array(C.evaluate_single(extremaUnique[0])),]
    else:
        if localExtrema:
            if uv_xyz: return extremaUnique
            else: return C.evaluate_list(extremaUnique)
        else: # return single maxima
            extremaUniquePoint = C.evaluate_list(extremaUnique)
            dispsExtrema = [np.linalg.norm(np.array(e) - p) for e in extremaUniquePoint]
            if maxSearch: # return single value in list for compatibility
                if uv_xyz: return [extremaUnique[dispsExtrema.index(max(dispsExtrema))],]
                else: return [np.array(extremaUniquePoint[dispsExtrema.index(max(dispsExtrema))]),]
            else: # minsearch
                if uv_xyz: return [extremaUnique[dispsExtrema.index(min(dispsExtrema))],]
                else: return [np.array(extremaUniquePoint[dispsExtrema.index(min(dispsExtrema))]),]

    # if len(extremaU) == 0: #return []
    #     # minima/maxima at curve ends
    #     if all([localDisp[1] < ll for ll in localDisp[2:-2]]): # u = 0
    #         if uv_xyz: return [0.,]
    #         else: return C.evaluate_single(0.)
    #     elif all([localDisp[-2] > ll for ll in localDisp[1:-1:-3]]): # u = 1
    #         if uv_xyz: return [1.,]
    #         else: return C.evaluate_single(1.)
    #     else:
    #         print("Bspline assumption fail")

    # if (localExtrema and (len(extremaUnique) == 1)) or not localExtrema:  # if there is only one extrema, localextrema is moot
    #     if uv_xyz: return [extremaUnique,] # return single value in list for compatibility
    #     else: return [np.array(C.evaluate_single(extremaUnique)),]
    # else:
    #     if localExtrema:
    #         if uv_xyz: return extremaUnique
    #         else: return C.evaluate_list(extremaUnique)
    #     else: # return single maxima
    #         extremaUniquePoint = C.evaluate_list(extremaUnique)
    #         dispsExtrema = [np.linalg.norm(np.array(e) - p) for e in extremaUniquePoint]
    #         if maxSearch: # return single value in list for compatibility
    #             if uv_xyz: return [extremaUnique[dispsExtrema.index(max(dispsExtrema))],]
    #             else: return [np.array(extremaUniquePoint[dispsExtrema.index(max(dispsExtrema))]),]
    #         else: # minsearch
    #             if uv_xyz: return [extremaUnique[dispsExtrema.index(min(dispsExtrema))],]
    #             else: return [np.array(extremaUniquePoint[dispsExtrema.index(min(dispsExtrema))]),]


#======================================================================================
# filepath = testDir + "/DrexelBlendedCylinder_topOvalPlane.step"

knotvector = [0.785398163398, 1.570796326795, 1.570796326795, 1.570796326795, 1.733091285395, 1.733091285395,
              2.284183032662, 2.284183032662, 2.628615374704, 2.628615374704, 2.863120266837, 2.863120266837,
              3.034908145413, 3.034908145413, 3.120802084701, 3.120802084701, 3.163749054345, 3.163749054345,
              3.206696023989, 3.206696023989, 3.378483902565, 3.378483902565, 3.550271781141, 3.550271781141,
              4.031867678959, 4.031867678959, 4.513463576777, 4.513463576777, 5.201166769448, 5.201166769448,
              5.630981264867, 5.630981264867, 5.899615324504, 5.899615324504, 6.067511611777, 6.067511611777,
              6.192645594127, 6.192645594127, 6.255212585302, 6.255212585302, 6.317779576477, 6.317779576477,
              6.442913558827, 6.442913558827, 6.568047541177, 6.568047541177, 7.068583470577, 7.068583470577,
              7.853981633974, 7.853981633974, 7.853981633974, 8.016276592574]

controlPoints = [[0.44, 0. , 2.  ], [ 0.44 , -0.07041708,  2. ], [ 0.4390651 , -0.14097727,  2. ], [ 0.43074754, -0.4513205 ,  2. ],
                     [ 0.41375243, -0.69145119,  2. ], [ 0.36260374, -1.07735828,  2. ], [ 0.33738697, -1.22490475,  2. ], [ 0.27883076, -1.46917541,  2. ],
                     [ 0.25043654, -1.56872952,  2. ], [ 0.18594868, -1.73343788,  2. ], [ 0.15500922, -1.80159469,  2. ], [ 0.0873894 , -1.89162136,  2. ],
                     [ 0.06126296, -1.92130129,  2. ], [ 0.00908672, -1.941482  ,  2. ], [-0.01092797, -1.94116291,  2. ], [-0.04561823, -1.9269163 ,  2. ],
                     [-0.0607161 , -1.91519705,  2. ], [-0.12831229, -1.85096302,  2. ], [-0.16256897, -1.78200904,  2. ], [-0.22253043, -1.64537717,  2. ],
                     [-0.24665442, -1.57412589,  2. ], [-0.32637208, -1.30107955,  2. ], [-0.36166841, -1.09502398,  2. ], [-0.41366057, -0.67914515,  2. ],
                     [-0.42881919, -0.4691928 ,  2. ], [-0.44559926,  0.03914896,  2. ], [-0.43902399,  0.33862537,  2. ], [-0.39798116,  0.82312728,  2. ],
                     [-0.37496332,  1.0086884 ,  2. ], [-0.31888405,  1.30657205,  2. ], [-0.29276722,  1.42122432,  2. ], [-0.23672487,  1.6036873 ,  2. ],
                     [-0.21214317,  1.67259416,  2. ], [-0.15835089,  1.78758639,  2. ], [-0.13212057,  1.83594557,  2. ], [-0.07969348,  1.89862881,  2. ],
                     [-0.0598763 ,  1.91842336,  2. ], [-0.01127545,  1.94296363,  2. ], [0.01874379, 1.94084753, 2. ], [0.08921129, 1.89947525, 2. ],
                     [0.12038053, 1.85115693, 2. ], [0.1744747 , 1.75703212, 2. ], [0.19651641, 1.70718636, 2. ], [0.29456594, 1.45276485, 2. ],
                     [0.336495  , 1.23856467, 2. ], [0.41977928, 0.68578778, 2. ], [0.44      , 0.34077118, 2. ], [0.44, 0.  , 2.  ]]

#BsplineKnotCurve = BSpline.Curve(normalize_kv=False)
BsplineKnotCurve = NURBS.Curve()
BsplineKnotCurve.degree = 3

BsplineKnotCurve.set_ctrlpts(compatibility.combine_ctrlpts_weights(controlPoints, weights=None))

#BsplineKnotCurve.knotvector = knotvector

# Auto-generate knot vector
BsplineKnotCurve.knotvector  = utilities.generate_knot_vector(BsplineKnotCurve.degree, len(BsplineKnotCurve.ctrlpts))
BsplineKnotCurve.evaluate()

localCentroid = np.array([-0.03,       -0.02851818,  1.        ])
localCentroid = np.array([-0.03002627, -1.01507624,  2.        ])
localCentroid = np.array([-0.03, -1.,    2.  ])

# maxPointU = splineCurveMinMaxPointDisp(BsplineKnotCurve, localCentroid, maxSearch=True)
# maxPoint_1 = np.array(BsplineKnotCurve.evaluate_single(maxPointU))

maxPoints_c = BsplineCurveExtremaDisp_2(BsplineKnotCurve,
                                     localCentroid,
                                     maxSearch=True,
                                     localExtrema=False,
                                     curvatureTest=True,
                                     uv_xyz=False)

maxPoints_s = BsplineCurveExtremaDisp_2(BsplineKnotCurve,
                                     localCentroid,
                                     maxSearch=True,
                                     localExtrema=False,
                                     curvatureTest = False,
                                     uv_xyz=False)

maxPoints_l = BsplineCurveExtremaDisp_2(BsplineKnotCurve,
                                     localCentroid,
                                     maxSearch=False,
                                     localExtrema=False,
                                     curvatureTest=True,
                                     uv_xyz=False) # should return null for ellipse spline

# multiple local minima (should return nearest orthogonal tangents)
minPoints_m = BsplineCurveExtremaDisp_2(BsplineKnotCurve,
                                     localCentroid,
                                     maxSearch=False,
                                     localExtrema=True,
                                     curvatureTest=False,
                                     uv_xyz=False)

# multiple local maxima (should return major radius extrema)
maxPoints_m = BsplineCurveExtremaDisp_2(BsplineKnotCurve,
                                     localCentroid,
                                     maxSearch=True,
                                     localExtrema=True,
                                     curvatureTest=False,
                                     uv_xyz=False)

# # Plot the control point polygon and the evaluated curve
# vis_comp = VisPlotly.VisCurve3D()
# BsplineKnotCurve.vis = vis_comp
# BsplineKnotCurve.render()

# Visualize data and evaluated points together
# import numpy as np
# import matplotlib.pyplot as plt
evalpts = np.array(BsplineKnotCurve.evalpts)
fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot(evalpts[:, 0], evalpts[:, 1], evalpts[:, 2])
ax.scatter([mm[0] for mm in maxPoints_m], [mm[1] for mm in maxPoints_m], [mm[2] for mm in maxPoints_m], color="red")
ax.scatter([mm[0] for mm in minPoints_m], [mm[1] for mm in minPoints_m], [mm[2] for mm in minPoints_m], color="blue")
ax.scatter(localCentroid[0], localCentroid[1], localCentroid[2], color="green")

set_aspect_equal_3d(ax)
plt.show()
_1=1
pass