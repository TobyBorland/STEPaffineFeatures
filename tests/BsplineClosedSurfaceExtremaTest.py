from geomdl import BSpline
from geomdl import evaluators
import numpy as np
from geomdl import NURBS
from geomdl import compatibility
#from geomdl import knotvector
from geomdl import utilities
from geomdl import helpers
import copy

# eps machine precision constant
eps = np.finfo(float).eps

eps_STEP_AP21 = 1e-6  # STEP precision seems to end here

def FreeCADpointSyntax(A, color=(0., 0., 0.)):
    '''print the instructions to create a FreeCAD model of a list of points'''
    print('import FreeCAD as App')
    print('import Draft')
    for i, xyz in enumerate(A):
        print('P{:} = Draft.make_point( {:1.8f}, {:1.8f}, {:1.8f}, color=({:1.8f}, {:1.8f}, {:1.8f}))'.format(i, xyz[0], xyz[1], xyz[2], color[0], color[1], color[2]))
    print('App.ActiveDocument.recompute()')
    #print('setview()')


def radiusCentre3points_2(p1, p2, p3):
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


def rationalSurfaceExtremaParam_4(S, p,
                                  maxSearch=True,
                                  localExtrema=False,
                                  curvatureTest=False,
                                  uv_xyz=True,
                                  eps1=0.0001,
                                  eps2=0.0005,
                                  deltaFactor=4,
                                  eps_bspline = 1E-10):
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

    #S.delta = delta  # set evaluation delta
    # tradeoff between sample density and feature identification

    #deltaFactor = 5
    #deltaFactor = 4

    # S.delta_u = 1 / (S.ctrlpts_size_u * deltaFactor)
    # S.delta_v = 1 / (S.ctrlpts_size_v * deltaFactor)

    S.delta_u = 1 / len(S.knotvector_u) * deltaFactor
    S.delta_v = 1 / len(S.knotvector_v) * deltaFactor

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

        #test1 = S.evaluate_single(uvArray[0]) #(0.0, 0.0))
        #test2 = S.evaluate_single(uvArray[207]) #(0.0, 1.0))
        #minSample=[]

        for CC, uv in neighbouringUV(len(uSamples), len(vSamples), closedU(), closedV()):

            pNW = pArray[CC[0]]
            pN  = pArray[CC[1]]
            pNE = pArray[CC[2]]
            pE  = pArray[CC[3]]
            pC =  pArray[CC[4]]
            pW  = pArray[CC[5]]
            pSW = pArray[CC[6]]
            pS  = pArray[CC[7]]
            pSE = pArray[CC[8]]

            dpuv_NW = np.linalg.norm(p - pNW)
            dpuv_N  = np.linalg.norm(p - pN)
            dpuv_NE = np.linalg.norm(p - pNE)
            dpuv_E  = np.linalg.norm(p - pE)
            dpuv    = np.linalg.norm(p - pC)
            dpuv_W  = np.linalg.norm(p - pW)
            dpuv_SW = np.linalg.norm(p - pSW)
            dpuv_S  = np.linalg.norm(p - pS)
            dpuv_SE = np.linalg.norm(p - pSE)

            #testSample = [pNW, pN, pNE, pW, pC, pE, pSW, pS, pSE]
            # if ((dpuv_S <= dpuv_NW) and (dpuv_S <= dpuv_N) and (dpuv_S <= dpuv_NE) and
            #         (dpuv_S <= dpuv_W) and (dpuv_S <= dpuv_E) and
            #         (dpuv_S <= dpuv_SW) and (dpuv_S <= dpuv) and (dpuv_S <= dpuv_SE)):
            #     minSample.append(pS)

            if maxSearch:
                if ((dpuv >= dpuv_NW) and (dpuv >= dpuv_N) and (dpuv >= dpuv_NE) and
                        (dpuv >= dpuv_W) and (dpuv >= dpuv_E) and
                        (dpuv >= dpuv_SW) and (dpuv >= dpuv_S) and (dpuv >= dpuv_SE)):

                    #b_NW = np.linalg.norm(pC - pNW)
                    b_N = np.linalg.norm(pC - pN)
                    #b_NE = np.linalg.norm(pC - pNE)
                    b_E = np.linalg.norm(pC - pE)

                    b_W = np.linalg.norm(pC - pW)
                    #b_SW = np.linalg.norm(pC - pSW)
                    b_S = np.linalg.norm(pC - pS)
                    #b_SE = np.linalg.norm(pC - pSE)

                    maxBound.append(max([b_N, b_E, b_W, b_S]))
                    extremaUV.append(uvArray[CC[4]]) # uv[0]*S.sample_size_v + uv[1]
                    extremaP.append(pArray[CC[4]])

            else:  # minS
                if ((dpuv <= dpuv_NW) and (dpuv <= dpuv_N) and (dpuv <= dpuv_NE) and
                        (dpuv <= dpuv_W) and (dpuv <= dpuv_E) and
                        (dpuv <= dpuv_SW) and (dpuv <= dpuv_S) and (dpuv <= dpuv_SE)):
                    # where a point is orthogonal to a planar surface -> sphere radius test, or accept minima


                    #b_NW = np.linalg.norm(pC - pNW)
                    b_N = np.linalg.norm(pC - pN)
                    #b_NE = np.linalg.norm(pC - pNE)
                    b_E = np.linalg.norm(pC - pE)

                    b_W = np.linalg.norm(pC - pW)
                    #b_SW = np.linalg.norm(pC - pSW)
                    b_S = np.linalg.norm(pC - pS)
                    #b_SE = np.linalg.norm(pC - pSE)

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
            pur, pucc = radiusCentre3points_2(pArray[CC[3]], puv, pArray[CC[5]])
            pvr, pvcc = radiusCentre3points_2(pArray[CC[7]], puv, pArray[CC[1]])

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
                    if (all((np.abs(discreteK[ssn]) >= discreteK[(S.sample_size_u * uv[0]) + uv[1]]) for ssn in surroundNodes) and
                            (discreteK[(S.sample_size_u * uv[0]) + uv[1]] > 0)):
                        # this version excludes negative curves from surrounding curvature values
                        # if all(np.abs(discreteK[ssn]) > discreteK[S.sample_size_u * u + v] or
                        # (discreteK[ssn]<0) for ssn in surroundNodes) and
                        # (discreteK[S.sample_size_u * u + v] > 0):

                        pC = np.array(S.evaluate_single(uv))
                        localExtremaP(pC)
                        maxBound.append(max([np.linalg.norm(pC - np.array(S.evaluate_single(uvk))) for uvk in surroundNodes]))
                        leu.append(uv)
                else:
                    if all((np.abs(discreteK[ssn]) >= discreteK[(S.sample_size_u * uv[0]) + uv[1]]) for ssn in
                           surroundNodes) and (discreteK[(S.sample_size_u * uv[0]) + uv[1]] < 0):

                        pC = np.array(S.evaluate_single(uv))
                        localExtremaP(pC)
                        maxBound.append(max([np.linalg.norm(pC - np.array(S.evaluate_single(uvk))) for uvk in surroundNodes]))
                        leu.append(uv)

        return leu

    def subSearchUV(mpU, mpV, tol):
        # create simple subdivisions around a minimum point as simple hillclimb search
        tolFlag = 0
        divU = S.delta_u
        divV = S.delta_v
        dpuvlocal = np.linalg.norm(p - np.array(S.evaluate_single((mpU, mpV,))))
        #dpuvlocal = np.inner(p - np.array(S.evaluate_single((mpU, mpV,))))

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
                      #dpuvlocal,
                      dpuv_E,
                      dpuv_SW,
                      dpuv_S,
                      dpuv_SE]

            if maxSearch:
                if np.abs(max(dispUV) - dpuvlocal) >= tol:
                    #np.abs(max(dispUV) - dpuvlocal) >= tol:
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
                return (np.inf, np.inf) # derivatives sometimes > |2|

            #test = S.derivatives(1.0, 0.0, 0)

            e = np.array(S.derivatives(cuv[0], cuv[1], 2))
            dif = e[0][0] - p # e[0][0] is surface point evaluated at cuv(u,v)
            c1v = np.linalg.norm(dif)
            c1 = c1v < eps1 #  |S(u,v) - p| < e1 (point coincidence)

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
            if ct[0] < S.knotvector_u[0]: # [ maxu - ( ct[0] - minu ), ct[1] ]
                if closedU(): ct = [S.knotvector_u[-1] - (S.knotvector_u[0] - ct[0]), ct[1]]
                #if closedU(): ct = [S.knotvector_u[0] + ct[0] % (S.knotvector_u[-1] - S.knotvector_u[0]), ct[1]]
                else: ct = [S.knotvector_u[0] + eps_bspline, ct[1]]

            elif ct[0] > S.knotvector_u[-1]:
                if closedU(): ct = [S.knotvector_u[0] + (ct[0] - S.knotvector_u[-1]), ct[1]]
                #if closedU(): ct = [S.knotvector_u[-1] - ct[0] % (S.knotvector_u[-1] - S.knotvector_u[0]), ct[1]]
                else: ct = [S.knotvector_u[-1] - eps_bspline, ct[1]]

            if ct[1] < S.knotvector_v[0]:
                if closedV(): ct = [ct[0], S.knotvector_v[-1] - (S.knotvector_v[0] - ct[1])]
                #if closedV(): ct = [ct[0], S.knotvector_v[0] + ct[1] % (S.knotvector_v[-1] - S.knotvector_v[0])]
                else: ct = [ct[0], S.knotvector_v[0] + eps_bspline]

            elif ct[1] > S.knotvector_v[-1]:
                #if closedV(): ct = [ct[0], S.knotvector_v[0] + (ct[1] - S.knotvector_v[-1])]
                if closedV(): ct = [ct[0], S.knotvector_v[-1] - ct[1] % (S.knotvector_v[-1] - S.knotvector_v[0])]
                else: ct = [ct[0], S.knotvector_v[-1] - eps_bspline]

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

    if (localExtrema and (len(extremaUnique) == 1)) or not localExtrema:  # if there is only one extrema, localextrema is moot
        if uv_xyz: return extremaUnique # return single value in list for compatibility
        else: return [np.array(S.evaluate_single(extremaUnique[0])),]
    else:
        if localExtrema:
            if uv_xyz: return extremaUnique
            else: return S.evaluate_list(extremaUnique)
        else: # return single maxima
            extremaUniquePoint = S.evaluate_list(extremaUnique)
            dispsExtrema = [np.linalg.norm(np.array(e) - p) for e in extremaUniquePoint]
            if maxSearch: # return single value in list for compatibility
                if uv_xyz: return [extremaUnique[dispsExtrema.index(max(dispsExtrema))],]
                else: return [np.array(extremaUniquePoint[dispsExtrema.index(max(dispsExtrema))]),]
            else: # minsearch
                if uv_xyz: return [extremaUnique[dispsExtrema.index(min(dispsExtrema))],]
                else: return [np.array(extremaUniquePoint[dispsExtrema.index(min(dispsExtrema))]),]


localCentroid = np.array([-9.77160901e-03, -5.00044450e-13,  0.00000000e+00])
localCentroid = np.array([0., 0.,  1.97])

# surfaceUdegree = 3
# surfaceVdegree = 3
# closedU = False
# closedV = True
# selfIntersect = False
# knotUmultiplicities = [4, 4]
# knotVmultiplicities = [1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1]
# knotsU = [-1.0, 1.0]
# knotsV = [-0.392699081699, -2.449293598295e-16, 0.392699081699, 0.785398163397, 1.178097245096, 1.570796326795, 1.963495408494, 2.356194490192, 2.94524311274, 3.534291735289, 3.926990816987, 4.319689898686, 4.712388980385, 5.105088062083, 5.497787143782, 5.890486225481, 6.28318530718, 6.675884388879]

_knotUvector = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
_knotVvector = [1.408501368196, 1.570796326795, 1.570796326795, 1.570796326795, 2.356194490192, 2.356194490192, 
                2.856730419592, 2.856730419592, 2.981864401942, 2.981864401942, 3.106998384292, 3.106998384292, 
                3.169565375467, 3.169565375467, 3.232132366642, 3.232132366642, 3.357266348992, 3.357266348992, 
                3.525162636266, 3.525162636266, 3.793796695902, 3.793796695902, 4.223611191322, 4.223611191322, 
                4.911314383992, 4.911314383992, 5.39291028181, 5.39291028181, 5.874506179628, 5.874506179628, 
                6.046294058204, 6.046294058204, 6.218081936781, 6.218081936781, 6.261028906425, 6.261028906425, 
                6.303975876069, 6.303975876069, 6.389869815357, 6.389869815357, 6.561657693933, 6.561657693933, 
                6.796162586065, 6.796162586065, 7.140594928107, 7.140594928107, 7.691686675375, 7.691686675375, 
                7.853981633974, 7.853981633974, 7.853981633974, 8.639379797371]

controlPointsList = [
                    [5.00000000e-01, 4.89842542e-16, 1.94000000e+00], [0.5 , 0.34334628, 1.94],
                    [0.47970362, 0.69078459, 1.94], [0.39544685, 1.25001608, 1.94], [0.35363397, 1.46634455, 1.94],
                    [0.25173233, 1.73076161, 1.94], [0.22875243, 1.78300213, 1.94], [0.17014475, 1.88498053, 1.94],
                    [0.14147652, 1.93836645, 1.94], [0.03817818, 1.9990141 , 1.94], [-0.02525655, 2.00311907, 1.94],
                    [-0.09998333, 1.96538689, 1.94], [-0.12277512, 1.94064082, 1.94], [-0.18402992, 1.86740293, 1.94],
                    [-0.21140698, 1.81577046, 1.94], [-0.26822354, 1.69431319, 1.94], [-0.29350409, 1.62318021, 1.94],
                    [-0.35104497, 1.43583831, 1.94], [-0.37754281, 1.31929509, 1.94], [-0.43441625, 1.01719295, 1.94],
                    [-0.45762795, 0.8298337 , 1.94], [-0.4990324 , 0.34106292, 1.94], [-0.50564302, 0.03948469, 1.94],
                    [-0.48873301, -0.47279363, 1.94], [-0.4734661 , -0.68443663, 1.94], [-0.42093597, -1.10461872, 1.94],
                    [-0.38541925, -1.31292279, 1.94], [-0.30373212, -1.59271487, 1.94], [-0.2789671 , -1.66608308, 1.94],
                    [-0.21601647, -1.80952625, 1.94], [-0.18303939, -1.88172746, 1.94], [-0.098698 , -1.96187381, 1.94],
                    [-0.07931412, -1.97794075, 1.94], [-0.02281914, -2.00114215, 1.94], [ 0.01938063, -2.00183212, 1.94],
                    [ 0.10560777, -1.96848121, 1.94], [ 0.13224001, -1.93181438, 1.94], [ 0.20923139, -1.82931073, 1.94],
                    [ 0.24071689, -1.75812759, 1.94], [ 0.30781136, -1.58676168, 1.94], [ 0.33673456, -1.48501092, 1.94],
                    [ 0.39638483, -1.23617637, 1.94], [ 0.42187291, -1.08683116, 1.94], [ 0.4735693 , -0.69679181, 1.94],
                    [ 0.49067666, -0.4547692 , 1.94], [ 0.49905809, -0.14204254, 1.94], [ 0.5 , -0.0709492, 1.94 ],
                    [5.00000000e-01, 4.89842542e-16, 1.94000000e+00],

                    [5.00000000e-01, 4.89842542e-16, 1.95570796e+00], [0.5 , 0.34318685, 1.95570796], [0.47966684, 0.69102865, 1.95570796],
                    [0.39547029, 1.24986054, 1.95570796], [0.35344318, 1.46683961, 1.95570796], [0.25178002, 1.73063785, 1.95570796],
                    [0.22861271, 1.78324524, 1.95570796], [0.17028446, 1.88473742, 1.95570796], [0.14012147, 1.93916201, 1.95570796],
                    [0.0388557 , 1.99861632, 1.95570796], [-0.02606528, 2.00271071, 1.95570796], [-0.0991746 , 1.96579525, 1.95570796],
                    [-0.1229573 , 1.940423 , 1.95570796], [-0.18366557, 1.86783855, 1.95570796], [-0.21148693, 1.81559954, 1.95570796],
                    [-0.26811626, 1.69454252, 1.95570796], [-0.29353978, 1.62306403, 1.95570796], [-0.35098788, 1.43602419, 1.95570796],
                    [-0.37756172, 1.31919464, 1.95570796], [-0.43438599, 1.01735367, 1.95570796], [-0.45763656, 0.82973205, 1.95570796],
                    [-0.49901862, 0.34122555, 1.95570796], [-0.50563829, 0.03934136, 1.95570796], [-0.48873633, -0.47269325, 1.95570796],
                    [-0.47344945, -0.68456984, 1.95570796], [-0.42095262, -1.10448551, 1.95570796], [-0.38532937, -1.31323062, 1.95570796],
                    [-0.30376417, -1.59260507, 1.95570796], [-0.27887457, -1.66629393, 1.95570796], [-0.216109 , -1.80931541, 1.95570796],
                    [-0.18221003, -1.88251557, 1.95570796], [-0.09890534, -1.96167678, 1.95570796], [-0.07863915, -1.97821795, 1.95570796],
                    [-0.02349411, -2.00086495, 1.95570796], [ 0.02008333, -2.00156033, 1.95570796], [ 0.10420237, -1.9690248 , 1.95570796],
                    [ 0.13243341, -1.9315569 , 1.95570796], [ 0.20884459, -1.8298257 , 1.95570796], [ 0.24078512, -1.75795331, 1.95570796],
                    [ 0.30771822, -1.58699958, 1.95570796], [ 0.336762 , -1.48489647, 1.95570796], [ 0.39634453, -1.23634447, 1.95570796],
                    [ 0.42188595, -1.08673276, 1.95570796], [ 0.47354843, -0.69694924, 1.95570796], [ 0.49067972, -0.45465521, 1.95570796],
                    [ 0.4990572 , -0.14207611, 1.95570796], [ 0.5 , -0.07091625, 1.95570796], [5.00000000e-01, 4.89842542e-16, 1.95570796e+00],

                    [4.86998832e-01, 4.89046476e-16, 1.98699883e+00], [0.48699883, 0.342868 , 1.98699883], [0.46673723, 0.68957982, 1.98699883],
                    [0.38266109, 1.24761249, 1.98699883], [0.34093014, 1.46315448, 1.98699883], [0.23974395, 1.72571507, 1.98699883],
                    [0.21706105, 1.77725323, 1.98699883], [0.15929167, 1.87777298, 1.98699883], [0.13082889, 1.92954148, 1.98699883],
                    [0.03362825, 1.9866091 , 1.98699883], [-0.02182267, 1.9902884 , 1.98699883], [-0.09169706, 1.95500636, 1.98699883],
                    [-0.11334884, 1.9316463 , 1.98699883], [-0.17296406, 1.86036874, 1.98699883], [-0.19987048, 1.80974882, 1.98699883],
                    [-0.25612534, 1.68949231, 1.98699883], [-0.28118298, 1.61901445, 1.98699883], [-0.33844554, 1.43257872, 1.98699883],
                    [-0.36482281, 1.31658841, 1.98699883], [-0.42154875, 1.01526977, 1.98699883], [-0.44469901, 0.82843134, 1.98699883],
                    [-0.4860363 , 0.34045341, 1.98699883], [-0.49263473, 0.03948361, 1.98699883], [-0.47574886, -0.47206357, 1.98699883],
                    [-0.4605154 , -0.68322344, 1.98699883], [-0.40808518, -1.10260628, 1.98699883], [-0.37266948, -1.31020261, 1.98699883],
                    [-0.29134814, -1.5887418 , 1.98699883], [-0.26678432, -1.66149098, 1.98699883], [-0.20438887, -1.80366908, 1.98699883],
                    [-0.17159547, -1.87466718, 1.98699883], [-0.09036418, -1.95185812, 1.98699883], [-0.07235017, -1.96674586, 1.98699883],
                    [-0.01990501, -1.98828408, 1.98699883], [ 0.01679873, -1.98889098, 1.98699883], [ 0.09670157, -1.95798619, 1.98699883],
                    [ 0.1224248 , -1.92323386, 1.98699883], [ 0.19767559, -1.82304756, 1.98699883], [ 0.22881526, -1.75286481, 1.98699883],
                    [ 0.2954256 , -1.58273541, 1.98699883], [ 0.3241739 , -1.48163681, 1.98699883], [ 0.38362096, -1.23364992, 1.98699883],
                    [ 0.40902358, -1.08482771, 1.98699883], [ 0.46061824, -0.69555585, 1.98699883], [ 0.47768933, -0.45407891, 1.98699883],
                    [ 0.48605889, -0.14179493, 1.98699883], [ 0.48699883, -0.07085037, 1.98699883], [4.86998832e-01, 4.89046476e-16, 1.98699883e+00],

                    [4.55707963e-01, 4.87130526e-16, 2.00000000e+00], [0.45570796, 0.34144534, 2. ], [0.43546743, 0.68709594, 2. ],
                    [0.35192856, 1.24156264, 2. ], [0.31002992, 1.45632001, 2. ], [0.2109719 , 1.71335835, 2. ],
                    [0.18868458, 1.76383105, 2. ], [0.13340877, 1.86001193, 2. ], [0.1028943 , 1.90965694, 2. ],
                    [0.0238317 , 1.95607551, 2. ], [-0.01493569, 1.95871229, 2. ], [-0.0703763 , 1.93071838, 2. ],
                    [-0.09097223, 1.90962753, 2. ], [-0.1457104 , 1.84418108, 2. ], [-0.17224094, 1.79496496, 2. ],
                    [-0.22682497, 1.67828019, 2. ], [-0.25158964, 1.60879053, 2. ], [-0.3080243 , 1.42505025, 2. ],
                    [-0.33424087, 1.30990294, 2. ], [-0.39052806, 1.01091489, 2. ], [-0.41359665, 0.82488302, 2. ],
                    [-0.45473415, 0.33926352, 2. ], [-0.46131868, 0.03923685, 2. ], [-0.4445046, -0.4701355, 2. ],
                    [-0.42931762, -0.68053046, 2. ], [-0.37718462, -1.09753587, 2. ], [-0.3418306, -1.3041801, 2. ],
                    [-0.26159733, -1.57899248, 2. ], [-0.23730552, -1.65079797, 2. ], [-0.17656149, -1.78921303, 2. ],
                    [-0.14263981, -1.85901713, 2. ], [-0.07065974, -1.927417 , 2. ], [-0.05443979, -1.94027447, 2. ],
                    [-0.01404107, -1.95686544, 2. ], [ 0.01178166, -1.95728162, 2. ], [ 0.0728724 , -1.93365296, 2. ],
                    [ 0.09913126, -1.90214387, 2. ], [ 0.16920455, -1.80885073, 2. ], [ 0.20028696, -1.73990163, 2. ],
                    [ 0.26545723, -1.57345033, 2. ], [ 0.29398994, -1.47332114, 2. ], [ 0.35283258, -1.22785565, 2. ],
                    [ 0.37812037, -1.07983828, 2. ], [ 0.42941245, -0.69284936, 2. ], [ 0.44643695, -0.45222337, 2. ],
                    [ 0.45477123, -0.14125616, 2. ], [ 0.45570796, -0.07055639, 2. ], [4.55707963e-01, 4.87130526e-16, 2.00000000e+00],

                    [4.40000000e-01, 4.86168722e-16, 2.00000000e+00], [0.44 , 0.34077118, 2. ], [0.41977928, 0.68578778, 2. ],
                    [0.336495 , 1.23856467, 2. ], [0.29456594, 1.45276485, 2. ], [0.19651641, 1.70718636, 2. ],
                    [0.1744747 , 1.75703212, 2. ], [0.12038053, 1.85115693, 2. ], [0.08921129, 1.89947525, 2. ],
                    [0.01874379, 1.94084753, 2. ], [-0.01127545, 1.94296363, 2. ], [-0.0598763 , 1.91842336, 2. ],
                    [-0.07969348, 1.89862881, 2. ], [-0.13212057, 1.83594557, 2. ], [-0.15835089, 1.78758639, 2. ],
                    [-0.21214317, 1.67259416, 2. ], [-0.23672487, 1.6036873 , 2. ], [-0.29276722, 1.42122432, 2. ],
                    [-0.31888405, 1.30657205, 2. ], [-0.37496332, 1.0086884 , 2. ], [-0.39798116, 0.82312728, 2. ],
                    [-0.43902399, 0.33862537, 2. ], [-0.44559926, 0.03914896, 2. ], [-0.42881919, -0.4691928 , 2. ],
                    [-0.41366057, -0.67914515, 2. ], [-0.36166841, -1.09502398, 2. ], [-0.32637208, -1.30107955, 2. ],
                    [-0.24665442, -1.57412589, 2. ], [-0.22253043, -1.64537717, 2. ], [-0.16256897, -1.78200904, 2. ],
                    [-0.12831229, -1.85096302, 2. ], [-0.0607161 , -1.91519705, 2. ], [-0.04561823, -1.9269163 , 2. ],
                    [-0.01092797, -1.94116291, 2. ], [ 0.00908672, -1.941482 , 2. ], [ 0.06126296, -1.92130129, 2. ],
                    [ 0.0873894 , -1.89162136, 2. ], [ 0.15500922, -1.80159469, 2. ], [ 0.18594868, -1.73343788, 2. ],
                    [ 0.25043654, -1.56872952, 2. ], [ 0.27883076, -1.46917541, 2. ], [ 0.33738697, -1.22490475, 2. ],
                    [ 0.36260374, -1.07735828, 2. ], [ 0.41375243, -0.69145119, 2. ], [ 0.43074754, -0.4513205 , 2. ],
                    [ 0.4390651 , -0.14097727, 2. ], [ 0.44 , -0.07041708, 2. ], [4.40000000e-01, 4.86168722e-16, 2.00000000e+00]
                    ]

#
# # Since we have no information on knot vectors, let's auto-generate them
# n_knotvector_u = utils.generate_knot_vector(p_degree_u, p_size_u)
# n_knotvector_v = utils.generate_knot_vector(p_degree_v, p_size_v)
#
#
# #
# # Import surface to NURBS-Python
# #
#
# # Create a NURBS surface instance
# surf = NURBS.Surface()
#
# # Fill the surface object
# surf.degree_u = p_degree_u
# surf.degree_v = p_degree_v
# surf_set_ctrlpts(n_ctrlptsw, p_size_u, p_size_v)
# surf.knotvector_u = n_knotvector_u
# surf.knotvector_v = n_knotvector_v
#
# # Set knot vectors
# surf.knotvector_u = utilities.generate_knot_vector(surf.degree_u, surf.ctrlpts_size_u)
# surf.knotvector_v = utilities.generate_knot_vector(surf.degree_v, surf.ctrlpts_size_v)


#_knotVvector = [kv - min(_knotVvector) for kv in _knotVvector]

BsplineKnotSurface = NURBS.Surface()
#BsplineKnotSurface = BSpline.Surface(normalize_kv=False)

ctrlptsw = compatibility.combine_ctrlpts_weights(controlPointsList, weights=None)

BsplineKnotSurface.degree_u = 3
BsplineKnotSurface.degree_v = 3
#BsplineKnotSurface.ctrlpts_size_u = 5
#BsplineKnotSurface.ctrlpts_size_v = 48
#BsplineKnotSurface.ctrlpts2d = controlPointsList
#BsplineKnotSurface.set_ctrlpts(controlPointsList, 5, 48)
BsplineKnotSurface.set_ctrlpts(ctrlptsw, 5, 48)
BsplineKnotSurface.set_ctrlpts(ctrlptsw, BsplineKnotSurface.ctrlpts_size_u, BsplineKnotSurface.ctrlpts_size_v)

BsplineKnotSurface.knotvector_u = utilities.generate_knot_vector(BsplineKnotSurface.degree_u, BsplineKnotSurface.ctrlpts_size_u)
BsplineKnotSurface.knotvector_v = utilities.generate_knot_vector(BsplineKnotSurface.degree_v, BsplineKnotSurface.ctrlpts_size_v)

# BsplineKnotSurface.knotvector_u = _knotUvector
# BsplineKnotSurface.knotvector_v = _knotVvector

#combine_ctrlpts_weights(ctrlpts, weights=None)
#
# if all([type(wcc)==np.ndarray for wc in controlPointsList for wcc in wc]):
#     wctrlpts = np.ones((BsplineKnotSurface.ctrlpts_size_v * BsplineKnotSurface.ctrlpts_size_u, 4))
#     wctrlpts.tolist()
#     for uindex in range(0, BsplineKnotSurface.ctrlpts_size_u):
#         for vindex in range(0, BsplineKnotSurface.ctrlpts_size_v):
#             wctrlpts[(uindex * BsplineKnotSurface.ctrlpts_size_v) + vindex][:3] = controlPointsList[uindex][vindex][:]
#
#     #wctrlpts.reshape(BsplineKnotSurface.ctrlpts_size_u, BsplineKnotSurface.ctrlpts_size_v, 4)

# # STEP only uses unitary weights

# Set control points
#BsplineKnotSurface.ctrlpts2d = ctrlpts # BSpline
# won't take numpy arrays
#BsplineKnotSurface.set_ctrlpts(wctrlpts.tolist(), BsplineKnotSurface.ctrlpts_size_u, BsplineKnotSurface.ctrlpts_size_v) #NURBS homogeneous weighted control points

#BsplineKnotSurface = BSpline.Surface()
#BsplineKnotSurface = BSpline.Surface(normalize_kv=False)

# Evaluates n-th order surface derivatives at the given (u, v) parameter pair.
# SKL[0][0] will be the surface point itself
# SKL[0][1] will be the 1st derivative w.r.t. v
# SKL[2][1] will be the 2nd derivative w.r.t. u and 1st derivative w.r.t. v

#BsplineKnotSurface.evaluator = evaluators.SurfaceEvaluator2()
BsplineKnotSurface.evaluator = evaluators.SurfaceEvaluatorRational()

# deltaFactor = 4
#
# # Set evaluation delta
# BsplineKnotSurface.delta_u = 1/(BsplineKnotSurface.ctrlpts_size_u * deltaFactor)
# BsplineKnotSurface.delta_v = 1/(BsplineKnotSurface.ctrlpts_size_v * deltaFactor)
# #BsplineKnotSurface.delta = 0.025  # this seems to be the minima delta under adaptive tesselation
#
# # Evaluate surface points
# BsplineKnotSurface.evaluate()

localCentroid = np.array([0.0, 0.0, 1.94])
localCentroid = np.array([0.0, 0.0, 0.])


minPointsUV = rationalSurfaceExtremaParam_4(BsplineKnotSurface,
                                            localCentroid,
                                            maxSearch=False,
                                            localExtrema=True,
                                            curvatureTest=True,
                                            uv_xyz=True)

minPoints = BsplineKnotSurface.evaluate_list(minPointsUV)


maxPointsUV = rationalSurfaceExtremaParam_4(BsplineKnotSurface,
                                            localCentroid,
                                            maxSearch=True,
                                            localExtrema=True,
                                            curvatureTest=False,
                                            uv_xyz=True)

maxPoints = BsplineKnotSurface.evaluate_list(maxPointsUV)


# Import and use Matplotlib's colormaps
from matplotlib import cm
import matplotlib.colors as mcolors

colors = [(1,0,0,c) for c in np.linspace(0,1,100)]
cmapred = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=5)
colors = [(0,0,1,c) for c in np.linspace(0,1,100)]
cmapblue = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=5)


from geomdl.visualization import VisMPL

vis_config = VisMPL.VisConfig(alpha=0.4, display_axes=False, display_ctrlpts=False)
BsplineKnotSurface.vis = VisMPL.VisSurface(vis_config)
BsplineKnotSurface.vis.mconf['others']='points'
BsplineKnotSurface.vis.mconf['alpha']=0.1

BsplineKnotSurface.render()#colormap=cm.cool)
pass
_1=1

NonMaxPointsL = []
MaxPointL=([np.array([-0.11418219,  1.92267258,  1.98242641]), localCentroid])

pNonMaxPoints = dict(points=NonMaxPointsL, name="NonMaxPoints", color="green", size=1)
#pSearchPoints = dict(points=SearchPointsL, name="pSearchPoints", color="black", size=1)
pMaxPoints = dict(points=MaxPointL, name="MaxPoint", color="red", size=1)

#BsplineKnotSurface.render(extras=[pMinPoints, pSearchPoints], colormap=cm.cool)
BsplineKnotSurface.render(extras=[pNonMaxPoints, pMaxPoints], colormap=cm.cool)
BsplineKnotSurface.render(extras=[pNonMaxPoints, pMaxPoints], colormap=cm.cool)


pass
_1=1