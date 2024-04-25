import verb
import numpy as np
from geomdl import BSpline
from geomdl import utilities
#from geomdl.visualization import VisMPL
from geomdl import exchange
from geomdl import operations

import numpy as np
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from matplotlib.patches import Ellipse, Circle
import mpl_toolkits.mplot3d.art3d as art3d

def curvature3points(p1, p2, p3):
    # Menger Curvature = 4*triangleArea/(sideLength1*sideLength2*sideLength3)
    # use numpy or verb methods?
    # 16 A^2 = (a+b+c) (b+c-a) (a+c-b) (a+b-c)

    # S1 = verb.verb_core_Vec.norm(verb.verb_core_Trig.sub(p1, p2))
    # S2 = verb.verb_core_Vec.norm(verb.verb_core_Trig.sub(p2, p3))
    # S3 = verb.verb_core_Vec.norm(verb.verb_core_Trig.sub(p1, p3))

    S1 = np.linalg.norm(p1, p2)
    S2 = np.linalg.norm(p2, p3)
    S3 = np.linalg.norm(p1, p3)

    SP = (S1 + S2 + S3) / 2  #  semi-perimeter

    # triangle area via Herons formula
    TA = np.sqrt(SP * (SP - S1) * (SP - S2) * (SP - S3))

    # curvature via Menger curvature formula
    return (4 * TA) / (S1 * S2 * S3)


# centre closer to p or centroid than 3 points =>


def radiusCentre3points(p1, p2, p3):
    # radius + centre from 3 point circle

    p1 = np.array([p1[0], p1[1], p1[2]])
    p2 = np.array([p2[0], p2[1], p2[2]])
    p3 = np.array([p3[0], p3[1], p3[2]])

    t = p2 - p1
    u = p3 - p1
    v = p3 - p2

    w = np.cross(t, u)    # triangle normal
    wsl = np.linalg.norm(w)
    if (wsl<10e-14): return False # triangle area too small (additionally check points for colinearity)

    wsl = np.dot(w, w)
    iwsl2 = 1. / (2. * wsl)
    tt = np.dot(t, t)
    uu = np.dot(u, u)

    circCenter = p1 + (u * tt * (np.dot(u, v)) - t * uu * (np.dot(t, v))) * iwsl2
    circRadius = np.sqrt(tt * uu * (np.dot(v, v)) * iwsl2 * 0.5)
    #circAxis   = w / np.sqrt(wsl)

    return circRadius, [circCenter[0], circCenter[1], circCenter[2]]


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

        print("dua: " + repr(dua) + "    dcc: " + repr(dcc))

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
        print("geomdl_cu: " + repr(cu))
        #print("geomdl_e: " + repr(curve.derivatives(cu, order=2)))
        e = np.array(curve.derivatives(cu, order=2) )#f(cu)
        dif = e[0] - p                               #   C(u) - p
        c1v = np.linalg.norm(dif)                    #   |C(u) - p|
        c2n = np.dot(e[1], dif)                      #   C'(u) * (C(u) - P)
        c2d = np.linalg.norm(e[1]) * c1v             #   |C'(u)||C(u) - P|
        c2v = c2n / c2d                              #   |C'(u) * (C(u) - P)| / |C'(u)||C(u) - P|
        #c2v = np.dot(e[1], dif) / (np.linalg.norm(e[1]) * np.linalg.norm(dif))
        c1 = c1v < eps1
        c2 = np.abs(c2v) < eps2
        if c1 and c2:
            return cu #, ccDump, dccDump
        #ct = n(cu, e, dif)                           #   u* = u - f / f'
        ct = cu - (np.dot(e[1], dif) / (np.dot(e[2], dif) + np.dot(e[1], e[1])))

        #print("geomdl_dif: " + repr(dif))
        #print("geomdl_c1v: " + repr(c1v))
        #print("geomdl_c2n: " + repr(c2n))
        #print("geomdl_c2v: " + repr(c2v))
        #print("geomdl_c1: " + repr(c1))
        #print("geomdl_c2: " + repr(c2))
        print("geomdl_ct: " + repr(ct))

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
            return cu #, ccDump, dccDump

        cu = ct
        i += 1

    return cu #, ccDump, dccDump


    # while (_g < _g1):
    #     i = _g
    #     _g = (_g + 1)
    #     u0 = python_internal_ArrayImpl._get((pts[i] if i >= 0 and i < len(pts) else None), 0)
    #     u1 = python_internal_ArrayImpl._get(python_internal_ArrayImpl._get(pts, (i + 1)), 0)
    #     p0 = (pts[i] if i >= 0 and i < len(pts) else None)[1:None]
    #     p1 = python_internal_ArrayImpl._get(pts, (i + 1))[1:None]
    #     proj = verb_core_Trig.segmentClosestPoint(p,p0,p1,u0,u1)
    #     d = verb_core_Vec.norm(verb_core_Vec.sub(p,proj.pt))
    #     if (d < _hx_min): #,_____________<<<<
    #         _hx_min = d
    #         u = proj.u

    # verb_core_Trig.segmentClosestPoint = function(pt,segpt0,segpt1,u0,u1) {
    # 	var dif = verb_core_Vec.sub(segpt1,segpt0);
    # 	var l = verb_core_Vec.norm(dif);
    # 	if(l < verb_core_Constants.EPSILON) return { u : u0, pt : segpt0};
    # 	var o = segpt0;
    # 	var r = verb_core_Vec.mul(1 / l,dif);
    # 	var o2pt = verb_core_Vec.sub(pt,o);
    # 	var do2ptr = verb_core_Vec.dot(o2pt,r);
    # 	if(do2ptr < 0) return { u : u0, pt : segpt0}; else if(do2ptr > l) return { u : u1, pt : segpt1};
    # 	return { u : u0 + (u1 - u0) * do2ptr / l, pt : verb_core_Vec.add(o,verb_core_Vec.mul(do2ptr,r))};
    # };

    # def segmentClosestPoint (p, s0, s1, u0, u1):
    # 	v = s1 - s0
    # 	dv = np.linalg.norm(v)
    # 	if dv < verb_EPSILON:
    # 	    return u0, s0
    # 	nv = v / dv
    # 	w = p - s0
    # 	os = np.dot(w, nv)
    # 	if os < 0:
    # 	    return u0, s0
    # 	elif os > dv:
    # 	    return u1, s1
    # 	return u0 + (u1 - u0) * os / dv, s0 + (os * nv)

    # s0 = np.array([-1, 1, 0])
    # s1 = np.array([1, 1, 0])
    # p = np.array([0, 0, 0])

    # def segmentClosestPoint2(p, s0, s1, u0, u1):
    #     #   https://web.archive.org/web/20220121145748/http://geomalgorithms.com/index.html
    #     v = s1 - s0
    #     w = p - s0
    #
    #     c1 = np.dot(w, v)
    #     if (c1 <= 0):
    #         return u0, s0#, np.linalg.norm(w)
    #
    #     c2 = np.dot(v, v)
    #     if (c2 <= c1):
    #         return u1, s1#, np.linalg.norm(p - s1)
    #
    #     b = c1 / c2
    #     pb = s0 + (b * v)
    #     return u0 + ((u1 - u0) * b), pb#, np.linalg.norm(p - pb)

    # for i in range(1, len(pts) - 2):
    #     #u0 = python_internal_ArrayImpl._get((pts[i] if i >= 0 and i < len(pts) else None), 0)
    #     u0 = kvs[i]
    #     #u1 = python_internal_ArrayImpl._get(python_internal_ArrayImpl._get(pts, (i + 1)), 0)
    #     u1 = kvs[i+1]
    #     #p0 = (pts[i] if i >= 0 and i < len(pts) else None)[1:None]
    #     p0 = np.array(pts[i])
    #     #p1 = python_internal_ArrayImpl._get(pts, (i + 1))[1:None]
    #     p1 = np.array(pts[i + 1])
    #     #proj = verb_core_Trig.segmentClosestPoint(p,p0,p1,u0,u1)
    #     proj = segmentClosestPoint2(p, p0, p1, u0, u1)
    #     #d = verb_core_Vec.norm(verb_core_Vec.sub(p, proj.pt))
    #     d = np.linalg.norm(p - proj[1])
    #     #     if (d < _hx_min): #,_____________<<<<
    #     if (d < d_min):
    #     #         _hx_min = d
    #         d_min = d
    #     #         u = proj.u
    #         u = proj[0]

    # for i in range(1, len(pts) - 2):
    #     u1 = kvs[i]
    #     p0 = np.array(pts[i - 1])
    #     p1 = np.array(pts[i])
    #     p2 = np.array(pts[i + 1])
    #     dcc, cc = radiusCentre3points2(p0, p1, p2)
    #     # du0_p = np.linalg.norm(p - p0)
    #     # du1_p = np.linalg.norm(p - p1)
    #     # du2_p = np.linalg.norm(p - p2)
    #     dcc_p = np.linalg.norm(p - cc)
    #     # dua = (du0_p + du1_p + du2_p) / 3
    #     dua = (np.linalg.norm(p - p0) + np.linalg.norm(p - p1) + np.linalg.norm(p - p2)) / 3 # ---------------------------
    #     # centre of 3 point circle is nearer or further than sample point barycentre from p
    #     orthoSign = dua > dcc_p  # closer means maximum curve
    #     #print("dua: " + repr(dua) + "    dcc: " + repr(dcc))
    #
    #     if (dcc < d_min) and orthoSign:
    #         d_min = dcc
    #         u = u1


# def _curveFurthestPoint(curve,p): #TODO: find multiple maxima
#     d_min = np.inf
#     u = 0.0
#     numSamples = curve.ctrlpts_size * curve.degree * 2
#     span = (curve.knotvector[-1] - curve.knotvector[0]) / (numSamples - 1)
#     kvs = np.array([curve.knotvector[0] + (span * i) for i in range(0, numSamples)])
#     pts = curve.evaluate_list(kvs)
#
#     #ccDump = []
#     #dccDump = []
#
#     for i in range(1, len(pts) - 2):
#         u1 = kvs[i]
#         p0 = np.array(pts[i - 1])
#         p1 = np.array(pts[i])
#         p2 = np.array(pts[i + 1])
#         dcc, cc = radiusCentre3points2(p0, p1, p2)
#         if (cc == np.inf):
#             # curvature method failure, use displacement minima , dcc_p = np.inf
#             du_p = [np.linalg.norm(p - i) for i in pts]
#             #u = pts[du_p.index(max(du_p))][0]
#             u = kvs[np.where(du_p == max(du_p))[0][0]]
#             break
#             # du_p = [du0_p, du1_p, du2_p]
#             # if (max(du_p) < d_min):
#             #     d_min = max(du_p)
#             #     offset = du_p.index(min(du_p)) - 1
#             #     u = pts[i+offset][0]
#         else:
#             dcc_p = np.linalg.norm(p - cc)
#             dua = (np.linalg.norm(p - p0) + np.linalg.norm(p - p1) + np.linalg.norm(p - p2)) / 3
#             # centre of 3 point circle is nearer or further than sample point barycentre from p
#             orthoSign = dua > dcc_p  # closer means maximum curve
#             # print("dua: "+repr(dua)+"    dcc: "+repr(dcc))
#             if (dcc < d_min) and orthoSign:  # dcc > d_max, initialised as 0.0 for minima
#                 d_min = dcc
#                 u = u1
#
#
#     for i in range(1, len(pts) - 2):
#         u1 = kvs[i]
#         p0 = np.array(pts[i - 1])
#         p1 = np.array(pts[i])
#         p2 = np.array(pts[i + 1])
#         dcc, cc = radiusCentre3points2(p0, p1, p2)
#         #ccDump.append(cc)
#         #dccDump.append(dcc)
#
#         # du0_p = np.linalg.norm(p - p0)
#         # du1_p = np.linalg.norm(p - p1)
#         # du2_p = np.linalg.norm(p - p2)
#         dcc_p = np.linalg.norm(p - cc)
#         # dua = (du0_p + du1_p + du2_p) / 3
#         dua = (np.linalg.norm(p - p0) + np.linalg.norm(p - p1) + np.linalg.norm(p - p2)) / 3 # ---------------------------
#         # centre of 3 point circle is nearer or further than sample point barycentre from p
#         orthoSign = dua > dcc_p  # closer means maximum curve
#         #print("dua: " + repr(dua) + "    dcc: " + repr(dcc))
#
#         if (dcc < d_min) and orthoSign:
#             d_min = dcc
#             u = u1
#
#     eps1 = 0.0001   # Euclidean distance measure
#     eps2 = 0.0005   # zero cosine measure
#     verb_EPSILON = 1E-10
#     minu = curve.knotvector[0]
#     maxu = curve.knotvector[-1]
#
#     curveClosed = np.linalg.norm(np.array(curve.ctrlpts[0]) - np.array(curve.ctrlpts[-1]))
#     closed = False
#     if (curveClosed < verb_EPSILON):
#         if (curveClosed*curveClosed < verb_EPSILON):
#             closed = True
#     cu = u
#
#     # def n(u2,e1,d):
#     #     #   Newton's method: 	 u* = u - f / f'
#     #     #   use product rule to form derivative, f':   f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
#     #     #   d:  ( C(u) - p )
#     #
#     #     #f1 = verb.verb_core_Vec.dot(e1[1], d)       # C'(u) * ( C(u) - p )
#     #     #s0 = verb.verb_core_Vec.dot(e1[2], d)       # C"(u) * ( C(u) - p )
#     #     #s1 = verb.verb_core_Vec.dot(e1[1], e1[1])   # C'(u) * C'(u)
#     #     f1 = np.dot(e1[1], d)       # C'(u) * ( C(u) - p )
#     #     s0 = np.dot(e1[2], d)       # C"(u) * ( C(u) - p )
#     #     s1 = np.dot(e1[1], e1[1])   # C'(u) * C'(u)
#     #     df = s0 + s1
#     #     return u2 - f1 / df
#     #     #return (u2 - np.dot(e1[1], d)) / (np.dot(e1[2], d) + np.dot(e1[1], e1[1]))
#
#     def n(u2, e1, d):
#         #   Newton's method: 	 u* = u - f / f'
#         #   use product rule to form derivative, f':   f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
#         #   d:  ( C(u) - p )
#         return u2 - (np.dot(e1[1], d) / (np.dot(e1[2], d) + np.dot(e1[1], e1[1])))
#
#     maxits = 5
#     i = 0
#     while(i < maxits):
#         print("geomdl_cu: " + repr(cu))
#         #print("geomdl_e: " + repr(curve.derivatives(cu, order=2)))
#         e = np.array(curve.derivatives(cu, order=2) )#f(cu)
#         dif = e[0] - p                               #   C(u) - p
#         c1v = np.linalg.norm(dif)                    #   |C(u) - p|
#         c2n = np.dot(e[1], dif)                      #   C'(u) * (C(u) - P)
#         c2d = np.linalg.norm(e[1]) * c1v             #   |C'(u)||C(u) - P|
#         c2v = c2n / c2d                              #   |C'(u) * (C(u) - P)| / |C'(u)||C(u) - P|
#         #c2v = np.dot(e[1], dif) / (np.linalg.norm(e[1]) * np.linalg.norm(dif))
#         c1 = c1v < eps1
#         c2 = np.abs(c2v) < eps2
#         if c1 and c2:
#             return cu #, ccDump, dccDump
#         ct = n(cu, e, dif)                           #   u* = u - f / f'
#
#         #print("geomdl_dif: " + repr(dif))
#         #print("geomdl_c1v: " + repr(c1v))
#         #print("geomdl_c2n: " + repr(c2n))
#         #print("geomdl_c2v: " + repr(c2v))
#         #print("geomdl_c1: " + repr(c1))
#         #print("geomdl_c2: " + repr(c2))
#         print("geomdl_ct: " + repr(ct))
#
#         if ct < minu:
#             if closed:
#                 #ct = maxu - (ct - minu) # NURBS book
#                 ct = maxu - (minu - ct)
#             else:
#                 ct = minu
#         elif ct > maxu:
#             if closed:
#                 ct = minu + (ct - maxu)
#             else:
#                 ct = maxu
#
#         c3v = np.linalg.norm(np.multiply(ct - cu, e[1]))
#         if c3v < eps1:
#             return cu #, ccDump, dccDump
#
#         cu = ct
#         i += 1
#
#     return cu #, ccDump, dccDump

# def rationalCurveClosestParam(curve,p):
#     _hx_min = Math.POSITIVE_INFINITY
#     u = 0.0
#     pts = verb_eval_Tess.rationalCurveRegularSample(curve,(len(curve.controlPoints) * curve.degree),True)
#     _g = 0
#     _g1 = (len(pts) - 1)
#     while (_g < _g1):
#         i = _g
#         _g = (_g + 1)
#         u0 = python_internal_ArrayImpl._get((pts[i] if i >= 0 and i < len(pts) else None), 0)
#         u1 = python_internal_ArrayImpl._get(python_internal_ArrayImpl._get(pts, (i + 1)), 0)
#         p0 = (pts[i] if i >= 0 and i < len(pts) else None)[1:None]
#         p1 = python_internal_ArrayImpl._get(pts, (i + 1))[1:None]
#         proj = verb_core_Trig.segmentClosestPoint(p,p0,p1,u0,u1)
#         d = verb_core_Vec.norm(verb_core_Vec.sub(p,proj.pt))
#         if (d < _hx_min): #,_____________<<<<
#             _hx_min = d
#             u = proj.u
#     maxits = 5
#     i = 0
#     e = None
#     eps1 = 0.0001
#     eps2 = 0.0005
#     dif = None
#     minu = (curve.knots[0] if 0 < len(curve.knots) else None)
#     maxu = verb_core_ArrayExtensions.last(curve.knots)
#     closed = (verb_core_Vec.normSquared(verb_core_Vec.sub((curve.controlPoints[0] if 0 < len(curve.controlPoints) else None),verb_core_ArrayExtensions.last(curve.controlPoints))) < verb_core_Constants.EPSILON)
#     cu = u
#     def _hx_local_0(u):
#         return verb_eval_Eval.rationalCurveDerivatives(curve,u,2)
#     f = _hx_local_0
#     def _hx_local_1(u,e,d):
#         f = verb_core_Vec.dot((e[1] if 1 < len(e) else None),d)
#         s0 = verb_core_Vec.dot((e[2] if 2 < len(e) else None),d)
#         s1 = verb_core_Vec.dot((e[1] if 1 < len(e) else None),(e[1] if 1 < len(e) else None))
#         df = (s0 + s1)
#         return (u - ((f / df)))
#     n = _hx_local_1
#     while (i < maxits):
#         e = f(cu)
#         dif = verb_core_Vec.sub((e[0] if 0 < len(e) else None),p)
#         c1v = verb_core_Vec.norm(dif)
#         c2n = verb_core_Vec.dot((e[1] if 1 < len(e) else None),dif)
#         c2d = (verb_core_Vec.norm((e[1] if 1 < len(e) else None)) * c1v)
#         c2v = (c2n / c2d)
#         c1 = (c1v < eps1)
#         c2 = (Reflect.field(Math,"fabs")(c2v) < eps2)
#         if (c1 and c2):
#             return cu
#         ct = n(cu,e,dif)
#         if (ct < minu):
#             ct = ((maxu - ((ct - minu))) if closed else minu)
#         elif (ct > maxu):
#             ct = ((minu + ((ct - maxu))) if closed else maxu)
#         c3v = verb_core_Vec.norm(verb_core_Vec.mul((ct - cu),(e[1] if 1 < len(e) else None)))
#         if (c3v < eps1):
#             return cu
#         cu = ct
#         i = (i + 1)
#     return cu



def rationalCurveFurthestParam(curve,p): #TODO: find multiple maxima
    # based on verb & Nurbs Book, Piegl & Tiller p.230
    # same idea to find orthogonal tangent, but at maximum displacement from a centroid point
    # revised to test minima segment candidates using 3-point defined circle parameters
    d_min = np.inf
    u = 0.0
    pts = verb.verb_eval_Tess.rationalCurveRegularSample(
        curve, (len(curve.controlPoints) * curve.degree * 2), True
    )

    ccDump = []
    dccDump = []
    for i in range(1, len(pts)-2):
        #u0 = pts[i-1][0]
        u1 = pts[i][0]
        #u2 = pts[i+1][0]
        p0 = pts[i-1][1:]
        p1 = pts[i][1:]
        p2 = pts[i+1][1:]
        dcc, cc = radiusCentre3points(p0, p1, p2)
        ccDump.append(cc)
        dccDump.append(dcc)

        du0_p = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, p0))
        du1_p = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, p1))
        du2_p = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, p2))
        dcc_p = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, cc))
        dua = (du0_p + du1_p + du2_p)/3
        # centre of 3 point circle is nearer or further than sample point barycentre from p
        orthoSign = dua > dcc_p # closer means maximum curve

        #print("dua: "+repr(dua)+"    dcc: "+repr(dcc))

        if (dcc < d_min) and orthoSign:
            d_min = dcc
            u = u1

    # for i in range(len(pts) - 1):
    #     u0 = pts[i][0]
    #     u1 = pts[i + 1][0]
    #     p0 = pts[i][1:]
    #     p1 = pts[i + 1][1:]
    #     proj = verb.verb_core_Trig.segmentClosestPoint(p, p0, p1, u0, u1)
    #     d1 = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, proj.pt))
    #     d2 = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, p0))
    #     d3 = verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, p1))
    #
    #     print("u: "+repr(u)+"    disp: "+repr(max([d_max, d1, d2, d3])))
    #
    #     if (d3 > d2) and (d3 > d1) and (d3 > d_max) and (u1 < 1.0):
    #         d_max = d3
    #         u = u1
    #     if (d2 > d3) and (d2 > d_max) and (u0 > 0.0):
    #         d_max = d2
    #         u = u0
    #         #print(u)
    #     if (d3 - d2 < verb.verb_core_Constants.EPSILON) and (d3 > d_max):
    #         d_max = d1
    #         u = proj.u
    #
    #     # if (d1 < min):
    #     #     min = d1
    #     #     u = proj.u
    #
    # print("u: "+repr(u))


    #   solve:  C'(u) * ( C(u) - P ) = 0 = f(u)    C(u) is the curve, p is the point, * is a dot product
    #
    #   use Newton-Raphson method:  u* = u - f / f'
    #
    #   use the product rule to form the derivative, f':   	f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
    #
    #   What is the conversion criteria? (Piegl & Tiller suggest)
    #
    #    |C(u) - p| < e1
    #
    #    |C'(u)*(C(u) - P)|
    #    ------------------  < e2
    #    |C'(u)| |C(u) - P|
    #
    #     1) first check 2 & 3
    #     2) if at least one of these is not, compute new value, otherwise halt
    #     3) ensure the parameter stays within range
    #    			* if not closed, don't allow outside of range a-b
    #    			* if closed (e.g. circle), allow to move back to beginning
    #     4)  if |(u* - u)C'(u)| < e1, halt

    # ?? discard u if u = 0.0 or u = 1.0 as cannot determine a tangent from end point.
    # compare against vertex point measurement

    eps1 = 0.0001   # Euclidean distance measure
    eps2 = 0.0005   # zero cosine measure
    minu = curve.knots[0]
    maxu = verb.verb_core_ArrayExtensions.last(curve.knots)
    closed = verb.verb_core_Vec.normSquared(verb.verb_core_Vec.sub(curve.controlPoints[0], verb.verb_core_ArrayExtensions.last(curve.controlPoints))) < verb.verb_core_Constants.EPSILON
    cu = u

    def f(u1):
        #print("verb_e: "+repr(verb.verb_eval_Eval.rationalCurveDerivatives(curve, u1,2)))
        return verb.verb_eval_Eval.rationalCurveDerivatives(curve, u1,2)

    def n(u2,e1,d):
        #   Newton's method: 	 u* = u - f / f'
        #   use product rule to form derivative, f':   f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
        #   d:  ( C(u) - p )

        f1 = verb.verb_core_Vec.dot(e1[1], d)       # C'(u) * ( C(u) - p )
        s0 = verb.verb_core_Vec.dot(e1[2], d)       # C"(u) * ( C(u) - p )
        s1 = verb.verb_core_Vec.dot(e1[1], e1[1])   # C'(u) * C'(u)
        df = s0 + s1
        return u2 - f1 / df

    maxits = 5
    i = 0
    while(i < maxits):
        print("verb_cu: " + repr(cu))
        e = f(cu)
        dif = verb.verb_core_Vec.sub(e[0], p)       #   C(u) - p
        c1v = verb.verb_core_Vec.norm(dif)          #   |C(u) - p|
        c2n = verb.verb_core_Vec.dot(e[1], dif)     #   C'(u) * (C(u) - P)
        c2d = verb.verb_core_Vec.norm(e[1]) * c1v   #   |C'(u)||C(u) - P|
        c2v = c2n / c2d                             #   |C'(u) * (C(u) - P)| / |C'(u)||C(u) - P|
        c1 = c1v < eps1
        c2 = np.abs(c2v) < eps2
        if c1 and c2:
            return cu #, ccDump, dccDump
        ct = n(cu, e, dif)                          #   u* = u - f / f'

        #print("verb_dif: " + repr(dif))
        #print("verb_c1v: " + repr(c1v))
        #print("verb_c2n: " + repr(c2n))
        #print("verb_c2v: " + repr(c2v))
        #print("verb_c1: " + repr(c1))
        #print("verb_c2: " + repr(c2))
        print("verb_ct: " + repr(ct))

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

        c3v = verb.verb_core_Vec.norm(verb.verb_core_Vec.mul(ct - cu, e[1]))
        if c3v < eps1:
            return cu #, ccDump, dccDump
        print(ct)
        cu = ct
        i += 1

    return cu #, ccDump, dccDump


curveDegree = 3

k = [0.785398163397, 1.570796326795, 1.570796326795, 1.570796326795, 1.733091285395, 1.733091285395, 2.284183032662, 2.284183032662, 2.628615374704, 2.628615374704, 2.863120266837, 2.863120266837, 3.034908145413, 3.034908145413, 3.120802084701, 3.120802084701, 3.163749054345, 3.163749054345, 3.206696023989, 3.206696023989, 3.378483902565, 3.378483902565, 3.550271781141, 3.550271781141, 4.031867678959, 4.031867678959, 4.513463576777, 4.513463576777, 5.201166769448, 5.201166769448, 5.630981264867, 5.630981264867, 5.899615324504, 5.899615324504, 6.067511611777, 6.067511611777, 6.192645594127, 6.192645594127, 6.255212585302, 6.255212585302, 6.317779576477, 6.317779576477, 6.442913558827, 6.442913558827, 6.568047541177, 6.568047541177, 7.068583470577, 7.068583470577, 7.853981633974, 7.853981633974, 7.853981633974, 8.016276592574]
#k = k[1:-1]
cp = [[0.44, 0.  , 2.  ],
      [ 0.44      , -0.07041708,  2.        ],
      [ 0.4390651 , -0.14097727,  2.        ],
      [ 0.43074754, -0.4513205 ,  2.        ],
      [ 0.41375243, -0.69145119,  2.        ],
      [ 0.36260374, -1.07735828,  2.        ],
      [ 0.33738697, -1.22490475,  2.        ],
      [ 0.27883076, -1.46917541,  2.        ],
      [ 0.25043654, -1.56872952,  2.        ],
      [ 0.18594868, -1.73343788,  2.        ],
      [ 0.15500922, -1.80159469,  2.        ],
      [ 0.0873894 , -1.89162136,  2.        ],
      [ 0.06126296, -1.92130129,  2.        ],
      [ 0.00908672, -1.941482  ,  2.        ],
      [-0.01092797, -1.94116291,  2.        ],
      [-0.04561823, -1.9269163 ,  2.        ],
      [-0.0607161 , -1.91519705,  2.        ],
      [-0.12831229, -1.85096302,  2.        ],
      [-0.16256897, -1.78200904,  2.        ],
      [-0.22253043, -1.64537717,  2.        ],
      [-0.24665442, -1.57412589,  2.        ],
      [-0.32637208, -1.30107955,  2.        ],
      [-0.36166841, -1.09502398,  2.        ],
      [-0.41366057, -0.67914515,  2.        ],
      [-0.42881919, -0.4691928 ,  2.        ],
      [-0.44559926,  0.03914896,  2.        ],
      [-0.43902399,  0.33862537,  2.        ],
      [-0.39798116,  0.82312728,  2.        ],
      [-0.37496332,  1.0086884 ,  2.        ],
      [-0.31888405,  1.30657205,  2.        ],
      [-0.29276722,  1.42122432,  2.        ],
      [-0.23672487,  1.6036873 ,  2.        ],
      [-0.21214317,  1.67259416,  2.        ],
      [-0.15835089,  1.78758639,  2.        ],
      [-0.13212057,  1.83594557,  2.        ],
      [-0.07969348,  1.89862881,  2.        ],
      [-0.0598763 ,  1.91842336,  2.        ],
      [-0.01127545,  1.94296363,  2.        ],
      [0.01874379, 1.94084753, 2.        ],
      [0.08921129, 1.89947525, 2.        ],
      [0.12038053, 1.85115693, 2.        ],
      [0.1744747 , 1.75703212, 2.        ],
      [0.19651641, 1.70718636, 2.        ],
      [0.29456594, 1.45276485, 2.        ],
      [0.336495  , 1.23856467, 2.        ],
      [0.41977928, 0.68578778, 2.        ],
      [0.44      , 0.34077118, 2.        ],
      [0.44, 0.  , 2.  ]]

km = [1, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 1]
k2 = [0.785398163397,
      1.570796326795,
      1.733091285395,
      2.284183032662,
      2.628615374704,
      2.863120266837,
      3.034908145413,
      3.120802084701,
      3.163749054345,
      3.206696023989,
      3.378483902565,
      3.550271781141,
      4.031867678959,
      4.513463576777,
      5.201166769448,
      5.630981264867,
      5.899615324504,
      6.067511611777,
      6.192645594127,
      6.255212585302,
      6.317779576477,
      6.442913558827,
      6.568047541177,
      7.068583470577,
      7.853981633974,
      8.016276592574]
cp2 = [[0.44, 0.  , 2.  ],
       [ 0.44      , -0.07041708,  2.        ],
       [ 0.4390651 , -0.14097727,  2.        ],
       [ 0.43074754, -0.4513205 ,  2.        ],
       [ 0.41375243, -0.69145119,  2.        ],
       [ 0.36260374, -1.07735828,  2.        ],
       [ 0.33738697, -1.22490475,  2.        ],
       [ 0.27883076, -1.46917541,  2.        ],
       [ 0.25043654, -1.56872952,  2.        ],
       [ 0.18594868, -1.73343788,  2.        ],
       [ 0.15500922, -1.80159469,  2.        ],
       [ 0.0873894 , -1.89162136,  2.        ],
       [ 0.06126296, -1.92130129,  2.        ],
       [ 0.00908672, -1.941482  ,  2.        ],
       [-0.01092797, -1.94116291,  2.        ],
       [-0.04561823, -1.9269163 ,  2.        ],
       [-0.0607161 , -1.91519705,  2.        ],
       [-0.12831229, -1.85096302,  2.        ],
       [-0.16256897, -1.78200904,  2.        ],
       [-0.22253043, -1.64537717,  2.        ],
       [-0.24665442, -1.57412589,  2.        ],
       [-0.32637208, -1.30107955,  2.        ],
       [-0.36166841, -1.09502398,  2.        ],
       [-0.41366057, -0.67914515,  2.        ],
       [-0.42881919, -0.4691928 ,  2.        ],
       [-0.44559926,  0.03914896,  2.        ],
       [-0.43902399,  0.33862537,  2.        ],
       [-0.39798116,  0.82312728,  2.        ],
       [-0.37496332,  1.0086884 ,  2.        ],
       [-0.31888405,  1.30657205,  2.        ],
       [-0.29276722,  1.42122432,  2.        ],
       [-0.23672487,  1.6036873 ,  2.        ],
       [-0.21214317,  1.67259416,  2.        ],
       [-0.15835089,  1.78758639,  2.        ],
       [-0.13212057,  1.83594557,  2.        ],
       [-0.07969348,  1.89862881,  2.        ],
       [-0.0598763 ,  1.91842336,  2.        ],
       [-0.01127545,  1.94296363,  2.        ],
       [0.01874379, 1.94084753, 2.        ],
       [0.08921129, 1.89947525, 2.        ],
       [0.12038053, 1.85115693, 2.        ],
       [0.1744747 , 1.75703212, 2.        ],
       [0.19651641, 1.70718636, 2.        ],
       [0.29456594, 1.45276485, 2.        ],
       [0.336495  , 1.23856467, 2.        ],
       [0.41977928, 0.68578778, 2.        ],
       [0.44      , 0.34077118, 2.        ],
       [0.44, 0.  , 2.  ]]

# kv = []
# for kml in range(len(km)):
#     for i in range(km[kml]):
#         kv.append(k2[kml])

def knotMultToVector(k, km):
    kv = []
    for kml in range(len(km)):
        for i in range(km[kml]):
            kv.append(k[kml])
    return kv

v2 = BSpline.Curve(normalize_kv=False)
v2.degree = 3
# Set control points
v2.ctrlpts = cp2
kv=knotMultToVector(k2, km)
# Set knot vector
v2.knotvector = kv

# from geomdl import multi
# # Create a curve container
# mcrv = multi.CurveContainer(v2)
#
# # Import Matplotlib visualization module
# from geomdl.visualization import VisMPL
#
# # Set the visualization component of the curve container
# mcrv.vis = VisMPL.VisCurve3D()
#
# # Plot the curves in the curve container
# mcrv.render()

# spline nearest/furthest point test using pboyer verb library
v1 = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(3, kv, cp2)
p = [ 4.,  -10.,  0.]
#maxInterceptU, ccD, dccD = rationalCurveFurthestParam(v1._data, p)
maxInterceptU = rationalCurveFurthestParam(v1._data, p)
maxIntercept1 = verb.verb_eval_Eval.dehomogenize(verb.verb_eval_Eval.curvePoint(v1._data, maxInterceptU))
#maxIntercept = v.closestParam(maxInterceptU)
print("verb maxIntercept: "+repr(maxIntercept1))

# repeat for geomdl NURBS-python library
maxInterceptU = splineCurveMinMaxPointDisp(v2, p)
maxIntercept2 = v2.evaluate_single(maxInterceptU)
print("geom maxIntercept: "+repr(maxIntercept2))

# def rationalCurveRegularSampleRange2(curve, start, end, numSamples, includeU):
#     if (numSamples < 1):
#         numSamples = 2
#     span = (end - start) / (numSamples - 1)
#     p = [start + (span * i) for i in range(0, numSamples)]
#
# # geomdl: evaluate_list(param_list)
# # Evaluates the curve for an input range of parameters.
#
# # Start evaluating from u=0.1 to u=0.5
# v2.evaluate(start=0.1, stop=0.5)
#
# # Get the evaluated points
# curve_points = v2.evalpts

#v = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(3, k, cp)

# verb_geom_NurbsCurve.byKnotsControlPointsWeights = function(degree,knots,controlPoints,weights) {
# 	return new verb_geom_NurbsCurve(new verb_core_NurbsCurveData(degree,knots.slice(),verb_eval_Eval.homogenize1d(controlPoints,weights)));
# };


#v = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(3, k[1:-1], cp[1:-1])

#vl = v.length()

# fecking plot in verb

#     //Tessellate a curve at a given tolerance
#     //
#     //**params**
#     //
#     //* The tolerance at which to sample the curve
#     //
#     //**returns**
#     //
#     //* A point represented as an array
#
#     public function tessellate(tolerance : Float = null) : Array<Point> {
#         return Tess.rationalCurveAdaptiveSample( _data, tolerance, false );
#     }



#print(v.length())
minIntercept1 = v1.closestPoint(p)
#minInterceptU = verb.verb_eval_Eval.curvePoint(v._data, minIntercept)
minInterceptU = v1.closestParam(minIntercept1)
print(minIntercept1)

# invert centroid to test minima
p2 = [ 4.,  10.,  0.]
invIntercept = v1.closestPoint(p2)
invInterceptU = v1.closestParam(invIntercept)
#print(invIntercept)

# use geoml representation for matplotlib display
curve = BSpline.Curve()
curve.degree = 3
curve.ctrlpts = cp
curve.knotvector = k
curve.delta = 0.01  # evaluation delta

# tangent vectors
curvetan = []

# Evaluate curve tangent at u = minInterceptU
ct1 = operations.tangent(curve, minInterceptU, normalize=True)
curvetan.append(ct1)

# Evaluate curve tangent at u = maxInterceptU
ct2 = operations.tangent(curve, maxInterceptU, normalize=True)
curvetan.append(ct2)

# Evaluate curve tangent at u = invInterceptU
ct3 = operations.tangent(curve, invInterceptU, normalize=True)
curvetan.append(ct3)

# # Evaluate curve tangent at u = 0.6
# ct4 = operations.tangent(curve, 0.6, normalize=True)
# curvetan.append(ct4)
#
# # Evaluate curve tangent at u = 0.8
# ct5 = operations.tangent(curve, 0.8, normalize=True)
# curvetan.append(ct5)
#
# # Evaluate curve tangent at u = 1.0
# ct6 = operations.tangent(curve, 1.0, normalize=True)
# curvetan.append(ct6)

#
# Control Points, Curve and Tangent Vector Plotting
#

# Arrange control points and evaluated curve points for plotting
ctrlpts = np.array(curve.ctrlpts)
curvepts = np.array(curve.evalpts)

# Convert tangent list into a NumPy array
ctarr = np.array(curvetan)

minIntercept_radius=verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, minIntercept))
minIntercept_circ = plt.Circle((p[0], p[1]), radius=minIntercept_radius, angle=0, ec='gray', fc=None, fill=False)

maxIntercept_radius=verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, maxIntercept))
maxIntercept_circ = plt.Circle((p[0], p[1]), radius=maxIntercept_radius, angle=0, ec='red', fc=None, fill=False)

invIntercept_radius=verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p2, invIntercept))
invIntercept_circ = plt.Circle((p2[0], p2[1]), radius=invIntercept_radius, angle=0, ec='cyan', fc=None, fill=False)

endIntercept_radius=verb.verb_core_Vec.norm(verb.verb_core_Vec.sub(p, verb.verb_eval_Eval.curvePoint(v._data, 1.0)))
endIntercept_circ = plt.Circle((p[0], p[1]), radius=endIntercept_radius, angle=0, ec='green', fc=None, fill=False)


fig, ax = plt.subplots(figsize=(10.67, 8), dpi=96)
ax.set_aspect(1)
#yaxis = plt.plot((-1, 25), (0, 0), "k-")  # y-axis line
cppolygon, = plt.plot(ctrlpts[:, 0], ctrlpts[:, 1], color='black', linestyle='-.', marker='o', markersize='3')  # control points polygon

# for c in range(len(ccD)-1):
#     C = plt.Circle((ccD[c][0], ccD[c][1]), radius=dccD[c], angle=0, ec='gray', fc=None, fill=False)
#     ax.add_artist(C)

ax.add_artist(minIntercept_circ)
ax.add_artist(maxIntercept_circ)
ax.add_artist(invIntercept_circ)
ax.add_artist(endIntercept_circ)

p_point, = plt.plot(p[0], p[1], color='red', linestyle='-.', marker='o', markersize='5')  #
p_to_minIntercept, = plt.plot((p[0], minIntercept[0]), (p[1], minIntercept[1]), color='gray', linestyle='-.', marker='o', markersize='5')  #
p_to_maxIntercept, = plt.plot((p[0], maxIntercept[0]), (p[1], maxIntercept[1]), color='red', linestyle='-.', marker='o', markersize='5')  #
p2_to_invIntercept, = plt.plot((p2[0], invIntercept[0]), (p2[1], invIntercept[1]), color='cyan', linestyle='-.', marker='o', markersize='5')  #

cppolygon, = plt.plot(ctrlpts[:, 0], ctrlpts[:, 1], color='gray', linestyle='-.', marker='o', markersize='3')  # control points polygon

curveplt, = plt.plot(curvepts[:, 0], curvepts[:, 1], color='green', linestyle='-')  # evaluated curve points
tanline = plt.quiver(ctarr[:, 0, 0], ctarr[:, 0, 1], ctarr[:, 1, 0], ctarr[:, 1, 1], color='blue', angles='xy', scale_units='xy', scale=1, width=0.003)  # tangents
#tanlinekey = plt.quiverkey(tanline, 23.75, -14.5, 35, "Tangent Vectors", coordinates='data', labelpos='W')
plt.legend([cppolygon, curveplt], ["Control Points", "Evaluated Curve"])
plt.axis([-5, 10, -10, 10])
#ax.add_patch(circ)

plt.show()

_1=1

import splipy
## using cubic splines
k = 4  # order
p = k - 1 # degree

# Create a clamped curve
cpts1 = np.array(
    [[0, 0, 0], [1, 0, 0], [1, 1, 1], [0, 1, 4], [-4, 3, 5], [4, 3, -2]])  # control points
k1 = np.array([0, 0, 0, 0, 1 / 3, 2 / 3, 1, 1, 1, 1])  # knot vector
basis1 = splipy.BSplineBasis(k, k1)
curve1 = splipy.Curve(basis1, cpts1)

# Unclamped version
n = len(cpts1) - 1
u = k1.copy()  # knot vector for unclamped curve
new_cpts = cpts1.copy() # controlpoints for unclamped curve

# unclamp left
for i in range(0, p - 1):
    u[p - i - 1] = u[p - i] - (u[n - i + 1] - u[n - i])
    k = p - 1
    for j in range(i, -1, -1):
        alfa = (u[p] - u[k]) / (u[p + j + 1] - u[k])
        new_cpts[j] = (new_cpts[j] - alfa * new_cpts[j + 1]) / (1 - alfa)
        k -= 1
u[0] = u[1] - (u[n - p + 2] - u[n - p + 1]) # set first knot

# unclamp right
for i in range(0, p - 1):
    u[n + i + 2] = u[n + i + 1] + (u[p + i + 1] - u[p + i])
    for j in range(i, -1, -1):
        alfa = (u[n + 1] - u[n - j]) / (u[n - j + i + 2] - u[n - j])
        new_cpts[n - j] = (new_cpts[n - j] - (1 - alfa) * new_cpts[n - j - 1]) / alfa
u[n + p + 1] = u[n + p] + (u[2 * p] - u[2 * p - 1]) # set last knot

# create unclamped bspline from above calculated knots, u, and controlpoints, new_cpts
basis2 = splipy.BSplineBasis(k, u)
curve2 = splipy.Curve(basis2, new_cpts)

# compare values
t1 = np.linspace(curve1.start(0), curve1.end(0), 1000)
t2 = np.linspace(curve2.start(0), curve2.end(0), 1000)
print(np.allclose(curve1.evaluate(t1), curve2.evaluate(t2))) # printing False, see visualization below

#
#
#     //Determine the closest point on a NURBS surface to a given point. *This is an experimental method and not hightly reliable.*
#     //
#     //**params**
#     //
#     //* The NURBS surface
#     //* The point to which we're trying to locate the closest point on the surface
#     //
#     //**returns**
#     //
#     //* The closest point on the surface, bounded by the parametric range of the surface
#
#     public static function rationalSurfaceClosestPoint( surface : NurbsSurfaceData, p : Point ) : Point {
#         var uv = Analyze.rationalSurfaceClosestParam( surface, p );
#         return Eval.rationalSurfacePoint( surface, uv[0], uv[1] );
#     }
#
#     //Determine the closest parameters on a NURBS surface to a given point. *This is an experimental method and not hightly reliable.*
#     //
#     //**params**
#     //
#     //* The NURBS surface
#     //* The point to which we're trying to locate the closest parameters on the surface
#     //
#     //**returns**
#     //
#     //* The closest parameters on the surface, bounded by the parametric domain of the surface
#
#     public static function rationalSurfaceClosestParam( surface : NurbsSurfaceData, p : Point ) : UV {
#
#         //for surfaces, we try to minimize the following:
#         //
#         //f = Su(u,v) * r = 0
#         //g = Sv(u,v) * r = 0
#         //
#         //  where r = S(u,v) - P
#         //
#         //Again, this requires newton iteration, but this time our objective function is vector valued
#         //
#         //    J d = k
#         //
#         //      d =   [ u* - u, v* - v ]
#         //		k = - [ f(u,v), g(u,v) ]
#         //		J =
#         //          |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
#         //		     Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
#         //
#         //
#         // 	we have similar halting conditions:
#         //
#         //  point coincidence
#         //
#         //		|S(u,v) - p| < e1
#         //
#         //  cosine
#         //
#         //   |Su(u,v)*(S(u,v) - P)|
#         //   ----------------------  < e2
#         //   |Su(u,v)| |S(u,v) - P|
#         //
#         //   |Sv(u,v)*(S(u,v) - P)|
#         //   ----------------------  < e2
#         //   |Sv(u,v)| |S(u,v) - P|
#         //
#         //  1) first check 2 & 3
#         // 	2) if at least one of these is not, compute new value, otherwise halt
#         // 	3) ensure the parameter stays within range
#         // 			* if not closed, don't allow outside of range a-b
#         // 			* if closed (e.g. circle), allow to move back to beginning
#         //  4)  if |(u* - u)C'(u)| < e1, halt
#         //
#
#         var maxits = 5;
#         var i = 0
#         , e
#         , eps1 = 0.0001
#         , eps2 = 0.0005
#         , dif
#         , minu = surface.knotsU[0]
#         , maxu = surface.knotsU.last()
#         , minv = surface.knotsV[0]
#         , maxv = surface.knotsV.last()
#         , closedu = isRationalSurfaceClosed(surface)
#         , closedv = isRationalSurfaceClosed(surface, false)
#         , cuv;
#
#         //todo: divide surface instead of a full on tessellation
#
#         //approximate closest point with tessellation
#         var tess = Tess.rationalSurfaceAdaptive( surface, new AdaptiveRefinementOptions() );
#
#         var dmin = Math.POSITIVE_INFINITY;
#
#         for ( i in 0...tess.points.length ){
#             var x = tess.points[i];
#             var d = Vec.normSquared( Vec.sub( p, x ) );
#
#             if ( d < dmin ){
#                 dmin = d;
#                 cuv = tess.uvs[i];
#             }
#         }
#
#         function f(uv : UV) : Array<Array<Point>> {
#             return Eval.rationalSurfaceDerivatives( surface, uv[0], uv[1], 2 );
#         }
#
#         function n(uv : UV, e : Array<Array<Point>>, r : Array<Float>) : UV {
#
#             //f = Su(u,v) * r = 0
#             //g = Sv(u,v) * r = 0
#
#             var Su = e[1][0];
#             var Sv = e[0][1];
#
#             var Suu = e[2][0];
#             var Svv = e[0][2];
#
#             var Suv = e[1][1];
#             var Svu = e[1][1];
#
#             var f = Vec.dot( Su, r );
#             var g = Vec.dot( Sv, r );
#
#             var k = [-f, -g];
#
#             var J00 = Vec.dot( Su, Su ) + Vec.dot( Suu, r );
#             var J01 = Vec.dot( Su, Sv ) + Vec.dot( Suv, r );
#             var J10 = Vec.dot( Su, Sv ) + Vec.dot( Svu, r );
#             var J11 = Vec.dot( Sv, Sv ) + Vec.dot( Svv, r );
#
#             var J = [ [ J00, J01 ], [ J10, J11 ] ];
#
#             //    d =   [ u* - u, v* - v ]
#             //		k = - [ f(u,v), g(u,v) ]
#             //		J =
#             //          |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
#             //		     Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
#             //
#
#             var d = Mat.solve( J, k );
#
#             return Vec.add( d, uv );
#
#         }
#
#         while( i < maxits ){
#
#             e = f(cuv);
#             dif = Vec.sub(e[0][0], p );
#
#             //  point coincidence
#             //
#             //		|S(u,v) - p| < e1
#             var c1v = Vec.norm( dif );
#
#             //
#             //  cosine
#             //
#             //   |Su(u,v)*(S(u,v) - P)|
#             //   ----------------------  < e2
#             //   |Su(u,v)| |S(u,v) - P|
#             //
#             //   |Sv(u,v)*(S(u,v) - P)|
#             //   ----------------------  < e2
#             //   |Sv(u,v)| |S(u,v) - P|
#             //
#             var c2an = Vec.dot( e[1][0], dif);
#             var c2ad = Vec.norm( e[1][0] ) * c1v;
#
#             var c2bn = Vec.dot( e[0][1], dif);
#             var c2bd = Vec.norm( e[0][1] ) * c1v;
#
#             var c2av = c2an / c2ad;
#             var c2bv = c2bn / c2bd;
#
#             var c1 = c1v < eps1;
#             var c2a = c2av < eps2;
#             var c2b = c2bv < eps2;
#
#             //if all of the tolerance are met, we're done
#             if (c1 && c2a && c2b){
#                 return cuv;
#             }
#
#             //otherwise, take a step
#             var ct = n(cuv, e, dif);
#
#             //correct for exceeding bounds
#             if ( ct[0] < minu ){
#                 ct = closedu ? [ maxu - ( ct[0] - minu ), ct[1] ] : [ minu + Constants.EPSILON, ct[1] ];
#             } else if (ct[0] > maxu){
#                 ct = closedu ? [ minu + ( ct[0] - maxu ), ct[1] ] : [ maxu - Constants.EPSILON, ct[1] ];
#             }
#
#             if ( ct[1] < minv ){
#                 ct = closedv ? [ ct[0], maxv - ( ct[1] - minv ) ] : [ ct[0], minv + Constants.EPSILON ];
#             } else if (ct[1] > maxv){
#                 ct = closedv ? [ ct[0], minv + ( ct[0] - maxv ) ] : [ ct[0], maxv - Constants.EPSILON ];
#             }
#
#             //if |(u* - u) C'(u)| < e1, halt
#             var c3v0 =  Vec.norm( Vec.mul(ct[0] - cuv[0], e[1][0] ) );
#             var c3v1 =  Vec.norm( Vec.mul(ct[1] - cuv[1], e[0][1] ) );
#
#             if (c3v0 + c3v1 < eps1) {
#                 return cuv;
#             }
#
#             cuv = ct;
#             i++;
#
#         }
#
#         return cuv;
#
#     }

_1=1

#     //Determine the closest point on a NURBS curve to a given point.
#     //
#     //**params**
#     //
#     //* The NURBS curve
#     //* The point to which we're trying to locate the closest point on the curve
#     //
#     //**returns**
#     //
#     //* The closest point on the surface, bounded by the parametric domain of the surface
#
#     public static function rationalCurveClosestPoint( curve : NurbsCurveData, p : Point ) : Point {
#         return Eval.rationalCurvePoint( curve, rationalCurveClosestParam(curve, p));
#     }
#
#     //Determine the closest parameters on a NURBS curve to a given point.
#     //
#     //**params**
#     //
#     //* The NURBS curve
#     //* The point to which we're trying to locate the closest parameter on the curve
#     //
#     //**returns**
#     //
#     //* The closest parameter on the curve, bounded by the parametric domain of the curve
#
#     public static function rationalCurveClosestParam( curve : NurbsCurveData, p : Point ) : Float {
#
#         //  We want to solve:
#         //
#         //   C'(u) * ( C(u) - P ) = 0 = f(u)
#         //
#         //  C(u) is the curve, p is the point, * is a dot product
#         //
#         // We'll use newton's method:
#         //
#         // 	 u* = u - f / f'
#         //
#         // We use the product rule in order to form the derivative, f':
#         //
#         //	f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
#         //
#         // What is the conversion criteria? (Piegl & Tiller suggest)
#         //
#         // |C(u) - p| < e1
#         //
#         // |C'(u)*(C(u) - P)|
#         // ------------------  < e2
#         // |C'(u)| |C(u) - P|
#         //
#         //  1) first check 2 & 3
#         // 	2) if at least one of these is not, compute new value, otherwise halt
#         // 	3) ensure the parameter stays within range
#         // 			* if not closed, don't allow outside of range a-b
#         // 			* if closed (e.g. circle), allow to move back to beginning
#         //  4)  if |(u* - u)C'(u)| < e1, halt
#         //
#
#         var min = Math.POSITIVE_INFINITY;
#         var u = 0.0;
#
#         var pts = Tess.rationalCurveRegularSample( curve, curve.controlPoints.length * curve.degree, true );
#
#         for ( i in 0...pts.length-1){
#
#             var u0 = pts[i][0];
#             var u1 = pts[i+1][0];
#
#             var p0 = pts[i].slice(1);
#             var p1 = pts[i+1].slice(1);
#
#             var proj = Trig.segmentClosestPoint( p, p0, p1, u0, u1 );
#             var d = Vec.norm( Vec.sub( p, proj.pt ) );
#
#             if ( d < min ){
#                 min = d;
#                 u = proj.u;
#             }
#         }
#
#         var maxits = 5
#         , i = 0
#         , e
#         , eps1 = 0.0001
#         , eps2 = 0.0005
#         , dif
#         , minu = curve.knots[0]
#         , maxu = curve.knots.last()
#         , closed = Vec.normSquared( Vec.sub( curve.controlPoints[0], curve.controlPoints.last() ) ) < Constants.EPSILON
#         , cu = u;
#
#         function f(u : Float) : Array<Point> {
#             return Eval.rationalCurveDerivatives( curve, u, 2 );
#         }
#
#         function n(u : Float, e: Array<Point>, d : Array<Float>) : Float {
#             //   C'(u) * ( C(u) - P ) = 0 = f(u)
#             var f = Vec.dot( e[1], d );
#
#             //	f' = C"(u) * ( C(u) - p ) + C'(u) * C'(u)
#             var s0 = Vec.dot( e[2], d )
#             , s1 = Vec.dot( e[1], e[1] )
#             , df = s0 + s1;
#
#             return u - f / df;
#         }
#
#         while( i < maxits ){
#
#             e = f( cu );
#             dif = Vec.sub( e[0], p );
#
#             // |C(u) - p| < e1
#             var c1v = Vec.norm( dif );
#
#             //C'(u) * (C(u) - P)
#             // ------------------ < e2
#             // |C'(u)| |C(u) - P|
#             var c2n = Vec.dot( e[1], dif);
#             var c2d = Vec.norm( e[1] ) * c1v;
#
#             var c2v = c2n / c2d;
#
#             var c1 = c1v < eps1;
#             var c2 = Math.abs(c2v) < eps2;
#
#             //if both tolerances are met
#             if (c1 && c2){
#                 return cu;
#             }
#
#             var ct = n(cu, e, dif);
#
#             //are we outside of the bounds of the curve?
#             if ( ct < minu ){
#                 ct = closed ? maxu - ( ct - minu ) : minu;
#             } else if (ct > maxu){
#                 ct = closed ? minu + ( ct - maxu ) : maxu;
#             }
#
#             //will our next step force us out of the curve?
#             var c3v = Vec.norm( Vec.mul(ct - cu, e[1] ) );
#
#             if (c3v < eps1) {
#                 return cu;
#             }
#
#             cu = ct;
#             i++;
#
#         }
#
#         return cu;
#
#     }