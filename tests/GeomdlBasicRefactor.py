
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

# The rational shapes have some minor differences between the non-rational ones.
# This class is designed to operate with weighted control points (Pw) as described in The NURBS Book by Piegl and Tiller.

def rationalCurveRegularSampleRange(curve, start, end, numSamples, includeU):
    if (numSamples < 1):
        numSamples = 2
    p = []
    span = (((end - start)) / ((numSamples - 1)))
    u = 0
    _g = 0
    _g1 = numSamples
    while (_g < _g1):
        i = _g
        _g = (_g + 1)
        u = (start + ((span * i)))
        if includeU:
            x = ([u] + verb_eval_Eval.rationalCurvePoint(curve, u))
            p.append(x)
        else:
            x1 = verb_eval_Eval.rationalCurvePoint(curve, u)
            p.append(x1)
    return p

def rationalCurveRegularSampleRange2(curve, start, end, numSamples, includeU):
    if (numSamples < 1):
        numSamples = 2
    span = (end - start) / (numSamples - 1)
    p = [start + (span * i) for i in range(0, numSamples)]
    if includeU:
        p = [[u] + rationalCurvePoint(curve, u) for u in p]
    else:
        p = [rationalCurvePoint(curve, u) for u in p]
    return p

# rationalCurveRegularSampleRange(curve : NurbsCurveData, start : Float, end : Float, numSamples : Int, includeU : Bool) : Array<Point>
# Sample a range of a NURBS curve at equally spaced parametric intervals
#
# params:
# NurbsCurveData object
# start parameter for sampling
# end parameter for sampling
# integer number of samples
# whether to prefix the point with the parameter
# returns:
# an dictionary of parameter - point pairs


# geomdl: evaluate_list(param_list)
# Evaluates the curve for an input range of parameters.

def rationalCurvePoint(curve,u):
    return verb_eval_Eval.dehomogenize(verb_eval_Eval.curvePoint(curve,u))

# rationalCurvePoint(curve : NurbsCurveData, u : Float) : Point
# Compute a point on a NURBS curve (see rationalCurveRegularSampleRange)
#
# params:
# integer degree of curve
# array of nondecreasing knot values
# 2d array of homogeneous control points, where each control point is an array of length (dim+1) and form (wi*pi, wi)
# parameter on the curve at which the point is to be evaluated
# returns:
# a point represented by an array of length (dim)

# ignore for STEP BSPLINE_WITH_KNOTS
def dehomogenize(homoPoint):
    dim = len(homoPoint)
    point = []
    wt = python_internal_ArrayImpl._get(homoPoint, (dim - 1))
    l = (len(homoPoint) - 1)
    _g = 0
    _g1 = l
    while (_g < _g1):
        i = _g
        _g = (_g + 1)
        point.append(((homoPoint[i] if i >= 0 and i < len(homoPoint) else None) / wt))
    return point

# dehomogenize(homoPoint : Point) : Point
# Dehomogenize a point
# params:
# a point represented by an array (wi*pi, wi) with length (dim+1)
# returns:
# a point represented by an array pi with length (dim)

# divide weighted points by weight



def curvePoint(curve,u):
    n = ((len(curve.knots) - curve.degree) - 2)
    return verb_eval_Eval.curvePointGivenN(n,curve,u)

# curvePoint(curve : NurbsCurveData, u : Float)
# Compute a point on a non-uniform, non-rational b-spline curve
# params:
# NurbsCurveData object representing the curve
# parameter on the curve at which the point is to be evaluated
# returns:
# a point represented by an array of length (dim)

def curvePointGivenN(n,curve,u):
    degree = curve.degree
    controlPoints = curve.controlPoints
    knots = curve.knots
    if (not verb_eval_Eval.areValidRelations(degree,len(controlPoints),len(knots))):
        raise haxe_Exception.thrown("Invalid relations between control points, knot Array, and n")
    knotSpan_index = verb_eval_Eval.knotSpanGivenN(n,degree,u,knots)
    basis_values = verb_eval_Eval.basisFunctionsGivenKnotSpanIndex(knotSpan_index,u,degree,knots)
    position = verb_core_Vec.zeros1d(len((controlPoints[0] if 0 < len(controlPoints) else None)))
    _g = 0
    _g1 = (degree + 1)
    while (_g < _g1):
        j = _g
        _g = (_g + 1)
        position = verb_core_Vec.add(position,verb_core_Vec.mul((basis_values[j] if j >= 0 and j < len(basis_values) else None),python_internal_ArrayImpl._get(controlPoints, ((knotSpan_index - degree) + j))))
    return position

# curvePointGivenN(n : Int, curve : NurbsCurveData, u : Float) : Point
# Compute a point on a non-uniform, non-rational b-spline curve (corresponds to algorithm 3.1 from The NURBS book, Piegl & Tiller 2nd edition)
# params:
# integer number of basis functions - 1 = knots.length - degree - 2
# NurbsCurveData object representing the curve
# parameter on the curve at which the point is to be evaluated
# returns:
# a point represented by an array of length (dim)

# geomdl: evaluate_single method will return the point evaluated as the specified parameter.

def knotSpanGivenN(n,degree,u,knots):
    if knots[n + 1] is None:
        _1=1
    print("edit line 6889")

    if (u > ((python_internal_ArrayImpl._get(knots, (n + 1)) - verb_core_Constants.EPSILON))):
        return n
    if (u < (((knots[degree] if degree >= 0 and degree < len(knots) else None) + verb_core_Constants.EPSILON))):
        return degree
    low = degree
    high = (n + 1)
    mid = Math.floor((((low + high)) / 2))
    while ((u < (knots[mid] if mid >= 0 and mid < len(knots) else None)) or ((u >= python_internal_ArrayImpl._get(knots, (mid + 1))))):
        if (u < (knots[mid] if mid >= 0 and mid < len(knots) else None)):
            high = mid
        else:
            low = mid
        mid = Math.floor((((low + high)) / 2))
    return mid

# knotSpanGivenN(n : Int, degree : Int, u : Float, knots : Array<Float>) : Int
# Find the span on the knot Array knots of the given parameter (corresponds to algorithm 2.1 from The NURBS book, Piegl & Tiller 2nd edition)
# params:
# integer number of basis functions - 1 = knots.length - degree - 2
# integer degree of function
# parameter
# array of nondecreasing knot values
# returns:
# the index of the knot span

# def findspan(p,u,U):  # removed n from the input parameters.
# 	# determine the index of the knot span of u in U
# 	# returns the index
# 	m = len(U) - 1
# 	# check knot vector is pinned. this is a basic check, not an exhaustive check.
# 	if U[p] != U[0] or U[m-p] !=U[m]:
# 		return 'knot vector is not pinned for the given degree'
# 	n = m - p - 1 # makes n redundant in the function input
#
# 	# special case to catch the end point of the knot vector
# 	if u == U[n+1]:	# u==U[m-p] = last knot value
# 		return n 	#assign the last non-zero width knot span
# 	# special case to catch the start point of the knot vector
# 	# Piegl doesn't address u = uo = u1 ... = up, and the while loop below doesn't catch it either?
# 	if u == U[p]:	# u==U[m-p] = first knot value
# 		return p 	# assign the first non-zero width knot span
# 	# setup binary search
# 	low = p
# 	high = n + 1
# 	mid = ( low + high ) / 2
# 	while u < U[mid] or u >= U[mid+1]:
# 		if u <U[mid]:
# 			high = mid
# 		else:
# 			low = mid
# 		mid = ( low + high ) / 2
# 	index = mid
# 	return index

def basisFunctionsGivenKnotSpanIndex(knotSpan_index,u,degree,knots):
    basisFunctions = verb_core_Vec.zeros1d((degree + 1))
    left = verb_core_Vec.zeros1d((degree + 1))
    right = verb_core_Vec.zeros1d((degree + 1))
    saved = 0
    temp = 0
    python_internal_ArrayImpl._set(basisFunctions, 0, 1.0)
    _g = 1
    _g1 = (degree + 1)
    while (_g < _g1):
        j = _g
        _g = (_g + 1)
        python_internal_ArrayImpl._set(left, j, (u - python_internal_ArrayImpl._get(knots, ((knotSpan_index + 1) - j))))
        python_internal_ArrayImpl._set(right, j, (python_internal_ArrayImpl._get(knots, (knotSpan_index + j)) - u))
        saved = 0.0
        _g2 = 0
        _g3 = j
        while (_g2 < _g3):
            r = _g2
            _g2 = (_g2 + 1)
            temp = ((basisFunctions[r] if r >= 0 and r < len(basisFunctions) else None) / ((python_internal_ArrayImpl._get(right, (r + 1)) + python_internal_ArrayImpl._get(left, (j - r)))))
            python_internal_ArrayImpl._set(basisFunctions, r, (saved + ((python_internal_ArrayImpl._get(right, (r + 1)) * temp))))
            saved = (python_internal_ArrayImpl._get(left, (j - r)) * temp)
        python_internal_ArrayImpl._set(basisFunctions, j, saved)
    return basisFunctions

# basisFunctionsGivenKnotSpanIndex(knotSpan_index : Int, u : Float, degree : Int, knots : KnotArray) : Array<Float>
# Compute the non-vanishing basis functions (corresponds to algorithm 2.2 from The NURBS book, Piegl & Tiller 2nd edition)
# params:
# Number, integer knot span index
# Number, float parameter
# Number, integer degree of function
# array of nondecreasing knot values
# returns:
# list of non-vanishing basis functions

# def basisfuns(i,p,u,U):
# 	#computes the non vanishing basis functions of degree p at u on knot span i
# 	N = [0.0] * (p+1)		# for a given knot span, there are at most p+1 non zero basis functions of degree p
# 	N[0] = 1.0			# assuming i is the correct index for u
# 	left = [0] * (p+1)
# 	right = [0] * (p+1)
# 	for j in range(1, p+1):
# 		left[j] = u - U[i+1-j]
# 		right[j] = U[i+j] - u
# 		saved = 0.0
# 		for r in range(j):
# 			temp = N[r]/(right[r+1] + left[j-r])
# 			N[r] = saved + right[r+1] * temp
# 			saved = left[j-r] * temp
# 		N[j] = saved
# 	return N

def areValidRelations(degree,num_controlPoints,knots_length):
    return ((((num_controlPoints + degree) + 1) - knots_length) == 0)


def rationalCurveDerivatives(curve,u,numDerivs = None):
    if (numDerivs is None):
        numDerivs = 1
    ders = verb_eval_Eval.curveDerivatives(curve,u,numDerivs)
    Aders = verb_eval_Eval.rational1d(ders)
    wders = verb_eval_Eval.weight1d(ders)
    k = 0
    i = 0
    CK = []
    _g = 0
    _g1 = (numDerivs + 1)
    while (_g < _g1):
        k = _g
        _g = (_g + 1)
        v = (Aders[k] if k >= 0 and k < len(Aders) else None)
        _g2 = 1
        _g3 = (k + 1)
        while (_g2 < _g3):
            i = _g2
            _g2 = (_g2 + 1)
            v = verb_core_Vec.sub(v,verb_core_Vec.mul((verb_core_Binomial.get(k,i) * (wders[i] if i >= 0 and i < len(wders) else None)),python_internal_ArrayImpl._get(CK, (k - i))))
        x = verb_core_Vec.mul((1 / (wders[0] if 0 < len(wders) else None)),v)
        CK.append(x)
    return CK

class verb_core_Binomial:
    _hx_class_name = "verb.core.Binomial"
    __slots__ = ()
    _hx_statics = ["memo", "get", "get_no_memo", "memo_exists", "get_memo", "memoize"]

    @staticmethod
    def get(n, k):
        if (k == 0.0):
            return 1.0
        if ((n == 0) or ((k > n))):
            return 0.0
        if (k > ((n - k))):
            k = (n - k)
        if verb_core_Binomial.memo_exists(n, k):
            return verb_core_Binomial.get_memo(n, k)
        r = 1
        n_o = n
        _g = 1
        _g1 = (k + 1)
        while (_g < _g1):
            d = _g
            _g = (_g + 1)
            if verb_core_Binomial.memo_exists(n_o, d):
                n = (n - 1)
                r = verb_core_Binomial.get_memo(n_o, d)
                continue
            r1 = n
            n = (n - 1)
            r = (r * r1)
            r = (r / d)
            verb_core_Binomial.memoize(n_o, d, r)
        return r

    @staticmethod
    def get_no_memo(n, k):
        if (k == 0):
            return 1
        if ((n == 0) or ((k > n))):
            return 0
        if (k > ((n - k))):
            k = (n - k)
        r = 1
        n_o = n
        _g = 1
        _g1 = (k + 1)
        while (_g < _g1):
            d = _g
            _g = (_g + 1)
            r1 = n
            n = (n - 1)
            r = (r * r1)
            r = (r / d)
        return r

    @staticmethod
    def memo_exists(n, k):
        if (n in verb_core_Binomial.memo.h):
            return (k in verb_core_Binomial.memo.h.get(n, None).h)
        else:
            return False

    @staticmethod
    def get_memo(n, k):
        return verb_core_Binomial.memo.h.get(n, None).h.get(k, None)

    @staticmethod
    def memoize(n, k, val):
        if (not (n in verb_core_Binomial.memo.h)):
            verb_core_Binomial.memo.set(n, haxe_ds_IntMap())
        verb_core_Binomial.memo.h.get(n, None).set(k, val)

    verb_core_Binomial._hx_class = verb_core_Binomial
    _hx_classes["verb.core.Binomial"] = verb_core_Binomial

def curveDerivatives(crv,u,numDerivs):
    n = ((len(crv.knots) - crv.degree) - 2)
    return verb_eval_Eval.curveDerivativesGivenN(n,crv,u,numDerivs)

def curveDerivativesGivenN(n,curve,u,numDerivs):
    degree = curve.degree
    controlPoints = curve.controlPoints
    knots = curve.knots
    if (not verb_eval_Eval.areValidRelations(degree,len(controlPoints),len(knots))):
        raise haxe_Exception.thrown("Invalid relations between control points, knot vector, and n")
    dim = len((controlPoints[0] if 0 < len(controlPoints) else None))
    du = (numDerivs if ((numDerivs < degree)) else degree)
    CK = verb_core_Vec.zeros2d((numDerivs + 1),dim)
    knotSpan_index = verb_eval_Eval.knotSpanGivenN(n,degree,u,knots)
    nders = verb_eval_Eval.derivativeBasisFunctionsGivenNI(knotSpan_index,u,degree,du,knots)
    k = 0
    j = 0
    _g = 0
    _g1 = (du + 1)
    while (_g < _g1):
        k = _g
        _g = (_g + 1)
        _g2 = 0
        _g3 = (degree + 1)
        while (_g2 < _g3):
            j = _g2
            _g2 = (_g2 + 1)
            python_internal_ArrayImpl._set(CK, k, verb_core_Vec.add((CK[k] if k >= 0 and k < len(CK) else None),verb_core_Vec.mul(python_internal_ArrayImpl._get((nders[k] if k >= 0 and k < len(nders) else None), j),python_internal_ArrayImpl._get(controlPoints, ((knotSpan_index - degree) + j)))))
    return CK

#     def deriv(self, u, k, d=None, rtype='Vector', domain='local'):
#         """
#         Compute the *k* -th derivative at point *u*.
#
#         :param float u: Parametric point.
#         :param int k: Derivative to return (0 <= k <= d).
#         :param int d: Highest derivative to compute. If *d* = *None*, then only
#             the *k* th derivative will be computed.
#         :param str rtype: Option to return a NumPy array or a Vector instance
#             (rtype = 'Vector' or 'ndarray').
#         :param str domain: Option to use local (0 <= u <= 1) or global
#             (a <= u <= b) domain ('local', 'l', 'global', 'g').
#
#         :return: Curve *k* -th derivative.
#         :rtype: :class:`.Vector` or ndarray
#         """
#         if self._cp is None:
#             return None
#         if d is None:
#             d = k
#         if is_local_domain(domain):
#             u = self.local_to_global_param(u)
#         der = rat_curve_derivs(self._n, self._p, self._uk,
#                                self._cpw, u, d)
#         if is_array_type(rtype):
#             return der[k]
#         else:
#             p0 = Point(der[0])
#             return Vector(der[k], p0)

# curveDerivativesGivenN(n : Int, curve : NurbsCurveData, u : Float, numDerivs : Int) : Array<Point>
# Determine the derivatives of a non-uniform, non-rational B-spline curve at a given parameter (corresponds to algorithm 3.1 from The NURBS book, Piegl & Tiller 2nd edition)
#
# params:
# integer number of basis functions - 1 = knots.length - degree - 2
# NurbsCurveData object representing the curve
# parameter on the curve at which the point is to be evaluated
# returns:
# a point represented by an array of length (dim)

# geomdl: derivatives(u, order=0, **kwargs)
# Evaluates n-th order curve derivatives at the given parameter value.

def derivativeBasisFunctionsGivenNI(knotSpan_index,u,p,n,knots):
    ndu = verb_core_Vec.zeros2d((p + 1),(p + 1))
    left = verb_core_Vec.zeros1d((p + 1))
    right = verb_core_Vec.zeros1d((p + 1))
    saved = 0.0
    temp = 0.0
    python_internal_ArrayImpl._set((ndu[0] if 0 < len(ndu) else None), 0, 1.0)
    _g = 1
    _g1 = (p + 1)
    while (_g < _g1):
        j = _g
        _g = (_g + 1)
        python_internal_ArrayImpl._set(left, j, (u - python_internal_ArrayImpl._get(knots, ((knotSpan_index + 1) - j))))
        python_internal_ArrayImpl._set(right, j, (python_internal_ArrayImpl._get(knots, (knotSpan_index + j)) - u))
        saved = 0.0
        _g2 = 0
        _g3 = j
        while (_g2 < _g3):
            r = _g2
            _g2 = (_g2 + 1)
            python_internal_ArrayImpl._set((ndu[j] if j >= 0 and j < len(ndu) else None), r, (python_internal_ArrayImpl._get(right, (r + 1)) + python_internal_ArrayImpl._get(left, (j - r))))
            temp = (python_internal_ArrayImpl._get((ndu[r] if r >= 0 and r < len(ndu) else None), (j - 1)) / python_internal_ArrayImpl._get((ndu[j] if j >= 0 and j < len(ndu) else None), r))
            python_internal_ArrayImpl._set((ndu[r] if r >= 0 and r < len(ndu) else None), j, (saved + ((python_internal_ArrayImpl._get(right, (r + 1)) * temp))))
            saved = (python_internal_ArrayImpl._get(left, (j - r)) * temp)
        python_internal_ArrayImpl._set((ndu[j] if j >= 0 and j < len(ndu) else None), j, saved)
    ders = verb_core_Vec.zeros2d((n + 1),(p + 1))
    a = verb_core_Vec.zeros2d(2,(p + 1))
    s1 = 0
    s2 = 1
    d = 0.0
    rk = 0
    pk = 0
    j1 = 0
    j2 = 0
    _g = 0
    _g1 = (p + 1)
    while (_g < _g1):
        j = _g
        _g = (_g + 1)
        python_internal_ArrayImpl._set((ders[0] if 0 < len(ders) else None), j, python_internal_ArrayImpl._get((ndu[j] if j >= 0 and j < len(ndu) else None), p))
    _g = 0
    _g1 = (p + 1)
    while (_g < _g1):
        r = _g
        _g = (_g + 1)
        s1 = 0
        s2 = 1
        python_internal_ArrayImpl._set((a[0] if 0 < len(a) else None), 0, 1.0)
        _g2 = 1
        _g3 = (n + 1)
        while (_g2 < _g3):
            k = _g2
            _g2 = (_g2 + 1)
            d = 0.0
            rk = (r - k)
            pk = (p - k)
            if (r >= k):
                python_internal_ArrayImpl._set((a[s2] if s2 >= 0 and s2 < len(a) else None), 0, (python_internal_ArrayImpl._get((a[s1] if s1 >= 0 and s1 < len(a) else None), 0) / python_internal_ArrayImpl._get(python_internal_ArrayImpl._get(ndu, (pk + 1)), rk)))
                d = (python_internal_ArrayImpl._get((a[s2] if s2 >= 0 and s2 < len(a) else None), 0) * python_internal_ArrayImpl._get((ndu[rk] if rk >= 0 and rk < len(ndu) else None), pk))
            if (rk >= -1):
                j1 = 1
            else:
                j1 = -rk
            if ((r - 1) <= pk):
                j2 = (k - 1)
            else:
                j2 = (p - r)
            _g4 = j1
            _g5 = (j2 + 1)
            while (_g4 < _g5):
                j = _g4
                _g4 = (_g4 + 1)
                python_internal_ArrayImpl._set((a[s2] if s2 >= 0 and s2 < len(a) else None), j, (((python_internal_ArrayImpl._get((a[s1] if s1 >= 0 and s1 < len(a) else None), j) - python_internal_ArrayImpl._get((a[s1] if s1 >= 0 and s1 < len(a) else None), (j - 1)))) / python_internal_ArrayImpl._get(python_internal_ArrayImpl._get(ndu, (pk + 1)), (rk + j))))
                d = (d + ((python_internal_ArrayImpl._get((a[s2] if s2 >= 0 and s2 < len(a) else None), j) * python_internal_ArrayImpl._get(python_internal_ArrayImpl._get(ndu, (rk + j)), pk))))
            if (r <= pk):
                python_internal_ArrayImpl._set((a[s2] if s2 >= 0 and s2 < len(a) else None), k, (-python_internal_ArrayImpl._get((a[s1] if s1 >= 0 and s1 < len(a) else None), (k - 1)) / python_internal_ArrayImpl._get(python_internal_ArrayImpl._get(ndu, (pk + 1)), r)))
                d = (d + ((python_internal_ArrayImpl._get((a[s2] if s2 >= 0 and s2 < len(a) else None), k) * python_internal_ArrayImpl._get((ndu[r] if r >= 0 and r < len(ndu) else None), pk))))
            python_internal_ArrayImpl._set((ders[k] if k >= 0 and k < len(ders) else None), r, d)
            temp = s1
            s1 = s2
            s2 = temp
    acc = p
    _g = 1
    _g1 = (n + 1)
    while (_g < _g1):
        k = _g
        _g = (_g + 1)
        _g2 = 0
        _g3 = (p + 1)
        while (_g2 < _g3):
            j = _g2
            _g2 = (_g2 + 1)
            _hx_local_2 = (ders[k] if k >= 0 and k < len(ders) else None)
            _hx_local_3 = j
            _hx_local_4 = (_hx_local_2[_hx_local_3] if _hx_local_3 >= 0 and _hx_local_3 < len(_hx_local_2) else None)
            python_internal_ArrayImpl._set(_hx_local_2, _hx_local_3, (_hx_local_4 * acc))
            (_hx_local_2[_hx_local_3] if _hx_local_3 >= 0 and _hx_local_3 < len(_hx_local_2) else None)
        acc = (acc * ((p - k)))
    return ders

# derivativeBasisFunctionsGivenNI(knotIndex : Int, u : Float, p : Int, n : Int, knots : KnotArray) : Array<Array<Float>>
# Compute the non-vanishing basis functions and their derivatives (corresponds to algorithm 2.3 from The NURBS book, Piegl & Tiller 2nd edition)
#
# params:
# integer knot span index
# float parameter
# integer degree
# integer number of basis functions - 1 = knots.length - degree - 2
# array of nondecreasing knot values
# returns:
# 2d array of basis and derivative values of size (n+1, p+1) The nth row is the nth derivative and the first row is made up of the basis function values.

# STEP doesn't provide weights, so no NURBS, or non-rational or BSpline, or all weights assumed as 1.0

# verb_eval_Eval.homogenize1d = function(controlPoints,weights) {
# 	var rows = controlPoints.length;
# 	var dim = controlPoints[0].length;
# 	var homo_controlPoints = [];
# 	var wt = 0.0;
# 	var ref_pt = [];
# 	var weights1;
# 	if(weights != null) weights1 = weights; else weights1 = verb_core_Vec.rep(controlPoints.length,1.0);
# 	var _g = 0;
# 	while(_g < rows) {
# 		var i = _g++;
# 		var pt = [];
# 		ref_pt = controlPoints[i];
# 		wt = weights1[i];
# 		var _g1 = 0;
# 		while(_g1 < dim) {
# 			var k = _g1++;
# 			pt.push(ref_pt[k] * wt);
# 		}
# 		pt.push(wt);
# 		homo_controlPoints.push(pt);
# 	}
# 	return homo_controlPoints;