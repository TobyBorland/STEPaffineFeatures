import verb
import numpy as np
from geomdl import BSpline
from geomdl import utilities
#from geomdl.visualization import VisMPL
from geomdl import exchange
from geomdl import operations

# c1 = verb.verb_geom_Circle([0, 0, 0],  [1, 0, 0], [0, 1, 0],2)
#
# c2 = verb.verb_geom_Circle([0, 0, 0],  [1, 0, 0], [0, 1, 1],2)
#
# i = verb.verb_geom_Intersect.curves(c1, c2)

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


# knots = [0.0, 0.000233103228888006, 0.000466206457776013, 0.000699309686664019, 0.000932412915552025]
# knotWeights = [0.0, 0.000233103228888006, 0.000466206457776013, 0.000699309686664019, 0.000932412915552025]
# controlPoints = [[0.006858469252432, -0.0222699757367726, 0.106086187509177],
#                  [0.00686049915998826, -0.0222673534192968, 0.106164963690966],
#                  [0.00686931060539193, -0.0222554352354012, 0.106241262581358],
#                  [0.00690606323291456, -0.022214856700426, 0.106388430796954],
#                  [0.00693379834238463, -0.0221865453457744, 0.106456749382545],
#                  [0.00701623177466876, -0.0221231537578462, 0.106575616747196],
#                  [0.00706942236875169, -0.0220895511309412, 0.106623608851608],
#                  [0.00720010164787025, -0.022035851920572, 0.106692489108923],
#                  [0.00727614476396148, -0.0220159757367721, 0.106713624376135],
#                  [0.00735574211737335, -0.0220159757367718, 0.106713624376135]]
# knotMultiplicities = [4, 2, 2, 2, 4]
#
# knots = [0.0, 0.0, 0.0, 0.0,
#          0.000233103228888006, 0.000233103228888006,
#          0.000466206457776013, 0.000466206457776013,
#          0.000699309686664019, 0.000699309686664019,
#          0.000932412915552025, 0.000932412915552025, 0.000932412915552025, 0.000932412915552025
#          ]
#
# p = [-0.012790131663638498, -0.022726639030997087, -0.022726639030997087]
#v = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(curveDegree, knots, controlPoints)

# p = [ 0.,  0.,  -0.25]
#
# cp = [[-1., -1.,  0.],
#       [-1.,  1.,  0.],
#       [ 1.,  1.,  0.],
#       [ 1., -1.,  0.]]
#
# #p = [ 6.,  -1.,  0.] # converges
#
# cp = [[0.,  0.,  0.],
#       [0.,  2.,  0.],
#       [2.,  2.,  0.],
#       [2.,  0.,  0.]]
#
# cp = [[0.,  0.,  0.],
#       [1.,  2.,  0.],
#       [3.,  2.,  0.],
#       [5.,  0.,  0.],
#       [7.,  0.,  0.],
#       [8.,  2.,  0.]]
#
# cp = [[0.,  0.,  0.],
#       [0.,  4.,  0.],
#       [4.,  4.,  0.],
#       [4.,  0.,  0.],
#       [8.,  0.,  0.],
#       [8.,  4.,  0.]]

#k = [0, 0, 0, 0, 1, 1, 1, 1]
#k = [0, 0.111111, 0.222222, 0.333333, 0.444444, 0.555555, 0.666666, 0.777777, 0.888888, 1.]

# Orthogonal tangents work for both maxima and minima => Newton-Raphson will work if correct subdivision of curve selected.

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

        print("dua: "+repr(dua)+"    dcc: "+repr(dcc))

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
    #   use newton's method:  u* = u - f / f'
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
        e = f(cu)
        dif = verb.verb_core_Vec.sub(e[0], p)       #   C'(u) - p
        c1v = verb.verb_core_Vec.norm(dif)          #   |C'(u) - p|
        c2n = verb.verb_core_Vec.dot(e[1], dif)     #   C'(u) * (C(u) - P)
        c2d = verb.verb_core_Vec.norm(e[1]) * c1v   #   |C'(u)||C(u) - P|
        c2v = c2n / c2d                             #   |C'(u) * (C(u) - P)| / |C'(u)||C(u) - P|
        c1 = c1v < eps1
        c2 = np.abs(c2v) < eps2
        if c1 and c2:
            return cu, ccDump, dccDump
        ct = n(cu, e, dif)                          #   u* = u - f / f'

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
            return cu, ccDump, dccDump
        print(ct)
        cu = ct
        i += 1

    return cu, ccDump, dccDump


# 4. Calculating the closest distance between a point and a surface
# To calculate the closest distance of a point to a surface we use an iterative technique. We start with
# a point on the surface and compute the vector to the specified point. Using dot products with vectors v1
# and v2 we estimate increments in the local coordinates. We iterate until the distance from the specified
# point to the surface does not change. For surfaces which are very curved we need to limit the increments
# to a maximum value using a damping coefficient (about 0.5 to 0.75) otherwise the iteration may diverge.
# Algorithm 2 Algorithm for computing minimum distance to a surface
# Require: Coordinates of point xP , NURBS of surface,
# Max. number of iterations I , damping damp, Tolerance.
# 1: Set max. increments: maxξ = 1,maxη = 1
# 2: Select starting point coordinates ξ0,η0 and compute point on surface xS0
# 3: Compute R0 = kxP −xS0k
# 4: for i = 0 to I do
# 5:    Compute unit vectors v1 and v2 in ξ- and η-directions.
# 6:    Lξ = |v1|,Lη = |v2|
# 7:    4ξ = v1 ·(xP −xSi)/Lξ
# 8:    4η = v2 ·(xP −xSi)/Lη
# 9:    ξi+1 = ξi +4ξ
# 10:   if ξi+1 > maxξ then
# 11:       ξi+1 = maxξ
# 12:   else if ξi+1 < 0 then
# 13:       ξi+1 = 0
# 14:   end if
# 15:   ηi+1 = ηi +4η
# 16:   if ηi+1 > maxη then
# 17:        ηi+1 = maxη
# 18:   else if ηi+1 < 0 then
# 19:       ηi+1 = 0
# 20:   end if
# 21:   Compute new point xSi+1(ξi+1,ηi+1) on surface.
# 22:   Compute distance Ri+1 = kxP −xSi+1k
# 23:   if Ri+1 −Ri <Tolerance then
# 24:       Exit loop
# 25:   end if
# 26:   Set maxξ = maxξ ∗ damp
# 27:   Set maxη = maxη ∗ damp
# 28: end for
# 29: return Minimum distance R, local and global coordinates of point on surface.


#   for surfaces, we try to minimize the following:
#
#   f = Su(u,v) * r = 0
#   g = Sv(u,v) * r = 0
#
#   where r = S(u,v) - P
#
#   Again, this requires newton iteration, but this time our objective function is vector valued
#
#   J d = k
#
#   d = [ u* - u, v* - v ]
#   k = - [ f(u,v), g(u,v) ]
#   J =
#     |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
#   		     Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
#
#    	we have similar halting conditions:
#
#     point coincidence
#
#   		|S(u,v) - p| < e1
#
#     cosine
#
#      |Su(u,v)*(S(u,v) - P)|
#      ----------------------  < e2
#      |Su(u,v)| |S(u,v) - P|
#
#      |Sv(u,v)*(S(u,v) - P)|
#      ----------------------  < e2
#      |Sv(u,v)| |S(u,v) - P|
#
#    1) first check 2 & 3
#    2) if at least one of these is not, compute new value, otherwise halt
#    3) ensure the parameter stays within range
#    	* if not closed, don't allow outside of range a-b
#    	* if closed (e.g. circle), allow to move back to beginning
#    4)  if |(u* - u)C'(u)| < e1, halt

curveDegree = 3
cp = [[0.,  2.,  0.],
      [1.,  4.,  0.],
      [3.,  4.,  0.],
      [5.,  0.,  0.],
      [7.,  0.,  0.],
      [8.,  2.,  0.]]

cp = [[0.,  0.,  0.],
      [1.,  4.,  0.],
      [3.,  4.,  0.],
      [5.,  0.,  0.],
      [7.,  0.,  0.],
      [8.,  4.,  0.]]

k = [0, 0, 0, 0, 0.444444, 0.555555, 1, 1, 1, 1.]  # clamped at ends


cp = [[41.50152415, 19.62887329,  4.88969889],
      [41.50151013, 19.6288766,   2.96937888],
      [41.50151013, 19.6288766,  -0.11885608],
      [41.50151013, 19.6288766,  -3.3774912 ],
      [41.50151013, 19.6288766,  -4.71580029],
      [41.50152108, 19.62887415, -4.8861958 ]] # [41.50152415, 19.62887329, 4.88969889]

k = [1.24647562, 1.24647562, 1.24647562, 1.24647562, 7.00745369, 10.51118053, 11.02238096, 11.02238096, 11.02238096, 11.02238096]


k = [0.785398163397, 1.570796326795, 1.570796326795, 1.570796326795, 1.733091285395, 1.733091285395, 2.284183032662, 2.284183032662, 2.628615374704, 2.628615374704, 2.863120266837, 2.863120266837, 3.034908145413, 3.034908145413, 3.120802084701, 3.120802084701, 3.163749054345, 3.163749054345, 3.206696023989, 3.206696023989, 3.378483902565, 3.378483902565, 3.550271781141, 3.550271781141, 4.031867678959, 4.031867678959, 4.513463576777, 4.513463576777, 5.201166769448, 5.201166769448, 5.630981264867, 5.630981264867, 5.899615324504, 5.899615324504, 6.067511611777, 6.067511611777, 6.192645594127, 6.192645594127, 6.255212585302, 6.255212585302, 6.317779576477, 6.317779576477, 6.442913558827, 6.442913558827, 6.568047541177, 6.568047541177, 7.068583470577, 7.068583470577, 7.853981633974, 7.853981633974, 7.853981633974, 8.016276592574]

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

v = verb.verb_geom_NurbsCurve.byKnotsControlPointsWeights(3, k, cp)

p = [ 4.,  -10.,  0.]
maxInterceptU, ccD, dccD = rationalCurveFurthestParam(v._data, p)
maxIntercept = verb.verb_eval_Eval.dehomogenize(verb.verb_eval_Eval.curvePoint(v._data, maxInterceptU))
#maxIntercept = v.closestParam(maxInterceptU)
print(maxIntercept)

#print(v.length())
minIntercept = v.closestPoint(p)
#minInterceptU = verb.verb_eval_Eval.curvePoint(v._data, minIntercept)
minInterceptU = v.closestParam(minIntercept)
print(minIntercept)

# invert centroid to test minima
p2 = [ 4.,  10.,  0.]
invIntercept = v.closestPoint(p2)
invInterceptU = v.closestParam(invIntercept)
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



#
# # Plot the control point polygon and the evaluated curve
# curve.vis = VisMPL.VisCurve2D()
# curve.render()
# # Plot the control point polygon and the evaluated curve
# curve.vis = VisMPL.VisCurve2D()
# curve.render()

#
# import splipy as sp
# #import numpy as np
# import matplotlib.pyplot as plt
# import splipy.curve_factory as cf
#
# basis = sp.BSplineBasis(order=3, knots=k)
# t = np.linspace(0,1, 150)
# curve = sp.Curve(basis, cp)
#
# t = np.linspace(0,4,150) # 150 visualization points on our parametric domain [0,4]
# x = curve.evaluate(t)    # compute (x,y)-coordinates of the curve evaluation
# print(x.shape)           # 2 components at 150 evaluation points, this prints (150,2)
#
# plt.plot(x[:,0], x[:,1])
# plt.show()

# // Create the nodes needed for the B-Spline curve.
# SoSeparator *
# makeCurve()
# {
#    SoSeparator *curveSep = new SoSeparator();
#    curveSep->ref();
#
#    // Set the draw style of the curve.
#    SoDrawStyle *drawStyle  = new SoDrawStyle;
#    drawStyle->lineWidth = 4;
#    curveSep->addChild(drawStyle);
#
#    // Define the NURBS curve including the control points
#    // and a complexity.
#    SoComplexity  *complexity = new SoComplexity;
#    SoCoordinate3 *controlPts = new SoCoordinate3;
#    SoNurbsCurve  *curve      = new SoNurbsCurve;
#    complexity->value = 0.8;
#    controlPts->point.setValues(0, 7, pts);
#    curve->numControlPoints = 7;
#    curve->knotVector.setValues(0, 10, knots);
#    curveSep->addChild(complexity);
#    curveSep->addChild(controlPts);
#    curveSep->addChild(curve);
#
#    curveSep->unrefNoDelete();
#    return curveSep;
# }


#
#     //Determine the derivatives of a NURBS curve at a given parameter
#     //
#     //**params**
#     //
#     //* NurbsCurveData object representing the curve - the control points are in homogeneous coordinates
#     //* parameter on the curve at which the point is to be evaluated
#     //* number of derivatives to evaluate
#     //
#     //**returns**
#     //
#     //* a point represented by an array of length (dim)
#
#     public static function rationalCurveDerivatives( curve : NurbsCurveData, u : Float, numDerivs : Int = 1  ) : Array<Point> {
#
#         var ders = curveDerivatives( curve, u, numDerivs )
#         , Aders = rational1d(ders)
#         , wders = weight1d(ders)
#         , k = 0
#         , i  = 0
#         , CK = [];
#
#         for (k in 0...numDerivs+1) {
#             var v = Aders[k];
#
#             for (i in 1...k+1) {
#                 Vec.subMulMutate( v, Binomial.get(k, i) * wders[i], CK[k-i] );
#             }
#
#             Vec.mulMutate( 1/wders[0], v );
#             CK.push( v ); //demogenize
#         }
#
#         return CK;
#
#     }

#
#     //Determine the derivatives of a non-uniform, non-rational B-spline curve at a given parameter
#     // (corresponds to algorithm 3.1 from The NURBS book, Piegl & Tiller 2nd edition)
#     //
#     //**params**
#     //
#     //* integer number of basis functions - 1 = knots.length - degree - 2
#     //* NurbsCurveData object representing the curve
#     //* parameter on the curve at which the point is to be evaluated
#     //
#     //**returns**
#     //
#     //* a point represented by an array of length (dim)
#
#     public static function curveDerivativesGivenN( n : Int, curve : NurbsCurveData, u : Float, numDerivs : Int ) : Array<Point> {
#
#         var degree = curve.degree
#         , controlPoints = curve.controlPoints
#         , knots = curve.knots;
#
#         if ( !areValidRelations( degree, controlPoints.length, knots.length ) ) {
#             throw 'Invalid relations between control points, knot vector, and n';
#         }
#
#         var dim = controlPoints[0].length
#         , du = numDerivs < degree ? numDerivs : degree
#         , CK = Vec.zeros2d( numDerivs+1, dim )
#         , knotSpan_index = knotSpanGivenN( n, degree, u, knots )
#         , nders = derivativeBasisFunctionsGivenNI( knotSpan_index, u, degree, du, knots )
#         , k = 0
#         , j = 0;
#
#         for (k in 0...du+1) {
#             for (j in 0...degree+1){
#                 Vec.addMulMutate( CK[k], nders[k][j], controlPoints[ knotSpan_index - degree + j ] );
#             }
#         }
#         return CK;
#     }

a=1

#
#     //Compute the derivatives at a point on a NURBS surface
#     //
#     //**params**
#     //
#     //* NurbsSurfaceData object representing the surface
#     //* number of derivatives to evaluate
#     //* u parameter at which to evaluate the derivatives
#     //* v parameter at which to evaluate the derivatives
#     //
#     //**returns**
#     //
#     //* a point represented by an array of length (dim)
#
#     public static function rationalSurfaceDerivatives( 	surface : NurbsSurfaceData,
#                                                         u : Float,
#                                                         v : Float,
#                                                         numDerivs : Int = 1) : Array<Array<Array<Float>>> {
#
#         var ders = surfaceDerivatives( surface, u, v, numDerivs)
#         , Aders = rational2d(ders)
#         , wders = weight2d(ders)
#         , SKL = new Array<Array<Array<Float>>>()
#         , dim = Aders[0][0].length;
#
#         for (k in 0...numDerivs+1){
#             SKL.push( new Array<Array<Float>>() );
#
#             for (l in 0...numDerivs-k+1){
#                 var v = Aders[k][l];
#
#                 for (j in 1...l+1){
#                     Vec.subMulMutate( v, Binomial.get(l, j) * wders[0][j], SKL[k][l-j] );
#                 }
#
#                 for (i in 1...k+1){
#                     Vec.subMulMutate( v, Binomial.get(k, i) * wders[i][0], SKL[k-i][l] );
#
#                     var v2 = Vec.zeros1d(dim);
#
#                     for (j in 1...l+1){
#                         Vec.addMulMutate( v2, Binomial.get(l, j) * wders[i][j], SKL[k-i][l-j] );
#                     }
#
#                     Vec.subMulMutate( v, Binomial.get(k, i), v2 );
#                 }
#
#                 Vec.mulMutate(1 / wders[0][0], v );
#                 SKL[k].push( v ); //demogenize
#             }
#         }
#
#         return SKL;
#     }

# // The control points for this curve
# float pts[7][3] = {
#    { 4.0, -6.0,  6.0},
#    {-4.0,  1.0,  0.0},
#    {-1.5,  5.0, -6.0},
#    { 0.0,  2.0, -2.0},
#    { 1.5,  5.0, -6.0},
#    { 4.0,  1.0,  0.0},
#    {-4.0, -6.0,  6.0}};
#
# // The knot vector
# float knots[10] = {1, 2, 3, 4, 5, 5, 6, 7, 8, 9};
#


# 4. Calculating the closest distance between a point and a surface
# To calculate the closest distance of a point to a surface we use an iterative technique. We start with
# a point on the surface and compute the vector to the specified point. Using dot products with vectors v1
# and v2 we estimate increments in the local coordinates. We iterate until the distance from the specified
# point to the surface does not change. For surfaces which are very curved we need to limit the increments
# to a maximum value using a damping coefficient (about 0.5 to 0.75) otherwise the iteration may diverge.
# Algorithm 2 Algorithm for computing minimum distance to a surface
# Require: Coordinates of point xP , NURBS of surface,
# Max. number of iterations I , damping damp, Tolerance.
# 1: Set max. increments: maxξ = 1,maxη = 1
# 2: Select starting point coordinates ξ0,η0 and compute point on surface xS0
# 3: Compute R0 = kxP −xS0k
# 4: for i = 0 to I do
# 5:    Compute unit vectors v1 and v2 in ξ- and η-directions.
# 6:    Lξ = |v1|,Lη = |v2|
# 7:    4ξ = v1 ·(xP −xSi)/Lξ
# 8:    4η = v2 ·(xP −xSi)/Lη
# 9:    ξi+1 = ξi +4ξ
# 10:   if ξi+1 > maxξ then
# 11:       ξi+1 = maxξ
# 12:   else if ξi+1 < 0 then
# 13:       ξi+1 = 0
# 14:   end if
# 15:   ηi+1 = ηi +4η
# 16:   if ηi+1 > maxη then
# 17:        ηi+1 = maxη
# 18:   else if ηi+1 < 0 then
# 19:       ηi+1 = 0
# 20:   end if
# 21:   Compute new point xSi+1(ξi+1,ηi+1) on surface.
# 22:   Compute distance Ri+1 = kxP −xSi+1k
# 23:   if Ri+1 −Ri <Tolerance then
# 24:       Exit loop
# 25:   end if
# 26:   Set maxξ = maxξ ∗ damp
# 27:   Set maxη = maxη ∗ damp
# 28: end for
# 29: return Minimum distance R, local and global coordinates of point on surface.

# verb_eval_Eval.rationalCurveDerivatives = function(curve,u,numDerivs) {
# 	if(numDerivs == null) numDerivs = 1;
# 	var ders = verb_eval_Eval.curveDerivatives(curve,u,numDerivs);
# 	var Aders = verb_eval_Eval.rational1d(ders);
# 	var wders = verb_eval_Eval.weight1d(ders);
# 	var k = 0;
# 	var i = 0;
# 	var CK = [];
# 	var _g1 = 0;
# 	var _g = numDerivs + 1;
# 	while(_g1 < _g) {
# 		var k1 = _g1++;
# 		var v = Aders[k1];
# 		var _g3 = 1;
# 		var _g2 = k1 + 1;
# 		while(_g3 < _g2) {
# 			var i1 = _g3++;
# 			verb_core_Vec.subMulMutate(v,verb_core_Binomial.get(k1,i1) * wders[i1],CK[k1 - i1]);
# 		}
# 		verb_core_Vec.mulMutate(1 / wders[0],v);
# 		CK.push(v);
# 	}
# 	return CK;
# };

#   for surfaces, we try to minimize the following:
#
#   f = Su(u,v) * r = 0
#   g = Sv(u,v) * r = 0
#
#     where r = S(u,v) - P
#
#   Again, this requires newton iteration, but this time our objective function is vector valued
#
#       J d = k
#
#         d =   [ u* - u, v* - v ]
#   		k = - [ f(u,v), g(u,v) ]
#   		J =
#             |Su|^2   +  Suu * r       Su*Sv  +  Suv * r
#   		     Su*Sv   +  Svu * r      |Sv|^2  +  Svv * r
#
#
#    	we have similar halting conditions:
#
#     point coincidence
#
#   		|S(u,v) - p| < e1
#
#     cosine
#
#      |Su(u,v)*(S(u,v) - P)|
#      ----------------------  < e2
#      |Su(u,v)| |S(u,v) - P|
#
#      |Sv(u,v)*(S(u,v) - P)|
#      ----------------------  < e2
#      |Sv(u,v)| |S(u,v) - P|
#
#     1) first check 2 & 3
#                       	2) if at least one of these is not, compute new value, otherwise halt
#    	3) ensure the parameter stays within range
#    			* if not closed, don't allow outside of range a-b
#    			* if closed (e.g. circle), allow to move back to beginning
#     4)  if |(u* - u)C'(u)| < e1, halt
#
