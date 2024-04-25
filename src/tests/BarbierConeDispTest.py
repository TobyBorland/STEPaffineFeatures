#https://liris.cnrs.fr/Documents/Liris-1297.pdf
# https://github.com/erich666/jgt-code/tree/master/Volume_09/Number_2/Barbier2004

from math import sqrt
import numpy as np


# class Axis:
#     def __init__(self, a, b):
#         self.a = a
#         self.b = b
#         self.axis = None
#         self.length = None
#         self.quadric = [0, 0, 0]
#         self.linear = [0, 0]

# class Cone(Axis):
#     def __init__(self, a, b, ra, rb):
#         super().__init__(a, b)
#         self.ra = ra
#         self.rb = rb
#         self.rrb = None
#         self.rra = None
#         self.conelength = None
#         self.side = None
#
#     def R(self, vector):
#         pass
#
#     def Set(self, vector1, vector2):
#         pass
#
#     def R(self, value):
#         pass

# class SphereCone(Axis):
#     def __init__(self, a, b, ra, rb):
#         super().__init__(a, b)
#         self.ra = ra
#         self.rb = rb
#         self.rrb = None
#         self.rra = None
#         self.ha = None
#         self.hb = None
#         self.hrb = None
#         self.hra = None
#         self.conelength = None
#         self.side = None
#
#     def R(self, vector):
#         pass
#
#     def Set(self, vector1, vector2):
#         pass
#
#     def R(self, value):
#         pass

# class Cylinder(Axis):
#     def __init__(self, a, b):
#         super().__init__(a, b)
#         self.r = [0, 0]


#   \class Vector vector.h
#   \brief This class implements a vector structure of three doubles.
#
#   Most binary operators have been overloaded as expected.
#   Destructive operators, such as addition and subtraction
#   have been implemented. Destructive operators += and -=
#   behave as one could expect.
#
#   Operators *= and /= behave in a specific way however,
#   scaling vector coordinates by the coordinates of the
#   argument vector.
#
#   The cross product of two vectors is defined by the operator /.

# def Coplanar(a, b, c, epsilon):
#     #   Check if three vectors are coplanar.
#     #   compute the cross product of a and b, and the dot product with c.
#     #   Compare the result with a given tolerance.
#     s = abs((a / b) * c) / ( np.linalg.norm(a) * np.linalg.norm(b) * np.linalg.norm(c))
#     return s < epsilon

# def Normalize(u):
#     #   Normalize a vector, computing the inverse of its norm and scaling
#     #   the components.
#     #   This function does not check if the vector is null,
#     #   which might resulting in floating point errors.
#     u *= 1.0 / np.linalg.norm(u)

# def Sine(u, v):
#     #  Returns the positive sine of two vectors, computes the
#     #  cross product of the vectors and normalizes the result.
#     return np.linalg.norm(u / v) / np.sqrt((u * u) * (v * v))

# def Cosine(u, v):
#     #  Returns the positive cosine of two vectors, computes the
#     #   dot product of the normalized vectors.
#     return (u*v)/np.sqrt((u*u)*(v*v))

# def Aligned(u, v):
#     #  Returns alignment boolean.
#     #  computes the cosine of the two vectors, and checks for unity.
#     c = Cosine(u, v)
#     c *= c
#     return (c > (1.0-0.0001))

# def Swap(a, b):
#     #  Swap two vectors.
#     t = a
#     a = b
#     b = t
#
# def Coplanar2(t, u, v, w):
#     #  Checks if four points are coplanar.
#     return Coplanar(u-t, v-t, w-t)

# def orthogonal(u):
#     #  Returns a new vector orthogonal to the argument vector.
#     a = [np.abs(x) for x in u]
#     i = 0
#     j = 1
#     if a[0] > a[1]:
#         if a[2] > a[1]:
#             j = 2
#     else:
#         i = 1
#         j = 2
#         if a[0] > a[2]:
#             j = 0
#     a = np.zeros(1, 3) #[0.0] * 3
#     a[i] = u[j]
#     a[j] = -u[i]
#     return a

# def __str__(u):
#     return '(' + str(u[0]) + ',' + str(u[1]) + ',' + str(u[2]) + ')'
#
# def min(a, b):
#     return a if a < b else b
#
# def max(a, b):
#     return a if a > b else b
#
# def min(a, b):
#     return a if a < b else b
#
# def max(a, b):
#     return a if a > b else b
#
# def min(a, b, c):
#     return a if a < b else (a if a < c else c)
#
# def max(a, b, c):
#     return a if a > b else (a if a > c else c)

# class Vector:
#     def __init__(self):
#         self.x = 0.0
#         self.y = 0.0
#         self.z = 0.0
#
#     def __init__(self, a):
#         self.x = a
#         self.y = a
#         self.z = a
#
#     def __init__(self, a, b, c):
#         self.x = a
#         self.y = b
#         self.z = c
#
#     def __init__(self, coords):
#         self.x = coords[0]
#         self.y = coords[1]
#         self.z = coords[2]
#
#     def __getitem__(self, i):
#         if i == 0:
#             return self.x
#         elif i == 1:
#             return self.y
#         elif i == 2:
#             return self.z
#
#     def __setitem__(self, i, value):
#         if i == 0:
#             self.x = value
#         elif i == 1:
#             self.y = value
#         elif i == 2:
#             self.z = value
#
#     def __pos__(self):
#         return Vector(self.x, self.y, self.z)
#
#     def __neg__(self):
#         return Vector(-self.x, -self.y, -self.z)
#
#     def __iadd__(self, other):
#         self.x += other.x
#         self.y += other.y
#         self.z += other.z
#         return self
#
#     def __isub__(self, other):
#         self.x -= other.x
#         self.y -= other.y
#         self.z -= other.z
#         return self
#
#     def __imul__(self, other):
#         self.x *= other.x
#         self.y *= other.y
#         self.z *= other.z
#         return self
#
#     def __idiv__(self, other):
#         self.x /= other.x
#         self.y /= other.y
#         self.z /= other.z
#         return self
#
#     def __imul__(self, scalar):
#         self.x *= scalar
#         self.y *= scalar
#         self.z *= scalar
#         return self
#
#     def __idiv__(self, scalar):
#         self.x /= scalar
#         self.y /= scalar
#         self.z /= scalar
#         return self
#
#     def __gt__(self, other):
#         return self.x > other.x and self.y > other.y and self.z > other.z
#
#     def __lt__(self, other):
#         return self.x < other.x and self.y < other.y and self.z < other.z
#
#     def __ge__(self, other):
#         return self.x >= other.x and self.y >= other.y and self.z >= other.z
#
#     def __le__(self, other):
#         return self.x <= other.x and self.y <= other.y and self.z <= other.z
#
#     def __add__(self, other):
#         return Vector(self.x + other.x, self.y + other.y, self.z + other.z)
#
#     def __sub__(self, other):
#         return Vector(self.x - other.x, self.y - other.y, self.z - other.z)
#
#     def __mul__(self, other):
#         return self.x * other.x + self.y * other.y + self.z * other.z
#
#     def __mul__(self, scalar):
#         return Vector(self.x * scalar, self.y * scalar, self.z * scalar)
#
#     class Vector:
#         def __init__(self, x, y, z):
#             self.x = x
#             self.y = y
#             self.z = z
#
#         def __mul__(self, other):
#             if isinstance(other, Vector):
#                 return self.x * other.x + self.y * other.y + self.z * other.z
#             elif isinstance(other, (int, float)):
#                 return Vector(self.x * other, self.y * other, self.z * other)
#
#         def __rmul__(self, other):
#             return self.__mul__(other)
#
#         def __truediv__(self, other):
#             if isinstance(other, Vector):
#                 return Vector(self.y * other.z - self.z * other.y, self.z * other.x - self.x * other.z,
#                               self.x * other.y - self.y * other.x)
#             elif isinstance(other, (int, float)):
#                 return Vector(self.x / other, self.y / other, self.z / other)
#
#         def __eq__(self, other):
#             return self.x == other.x and self.y == other.y and self.z == other.z
#
#         def __ne__(self, other):
#             return not self.__eq__(other)
#
#         def __getitem__(self, index):
#             if index == 0:
#                 return self.x
#             elif index == 1:
#                 return self.y
#             elif index == 2:
#                 return self.z
#
#         def __abs__(self):
#             return Vector(abs(self.x), abs(self.y), abs(self.z))
#
#         def Norm(self):
#             return np.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
#
#         def Normalized(self):
#             return self * (1.0 / self.Norm())
#
#         def NormInfinity(self):
#             return max(np.abs(self.x), np.abs(self.y), np.abs(self.z))
#
#         def Abs(u):
#             return [u[0] if u[0] > 0.0 else -u[0], u[1] if u[1] > 0.0 else -u[1], u[2] if u[2] > 0.0 else -u[2]]
#
#         def min(a, b):
#             return [a[0] if a[0] < b[0] else b[0], a[1] if a[1] < b[1] else b[1], a[2] if a[2] < b[2] else b[2]]
#
#         def max(a, b):
#             return [a[0] if a[0] > b[0] else b[0], a[1] if a[1] > b[1] else b[1], a[2] if a[2] > b[2] else b[2]]

class Axis:
    def __init__(self, a, b):
        self.a = a
        self.b = b
        self.axis = b - a
        self.length = np.sqrt(self.axis.x**2 + self.axis.y**2 + self.axis.z**2)
        self.axis /= self.length

class Cylinder:
    #   class implements a simple data-structure to define
    #   a simple cylinder characterized by its end vertices and radius

    def __init__(self, a, b, r):
        #   \brief Creates a generic cylinder given vertices and radius.
        #   \param a, b End vertices of the axis.
        #   \param r Radius of the cylinder.
        self.Axis = (a, b)
        self.r = [r, r * r]
        self.c = (a + b) * 0.5
        self.h = 0.5 * Axis.length

    def R(self, p):
        #   distance between a point in space and the cylinder.
        n = p - self.c
        y = self.Axis * n
        y = np.abs(y)
        yy = y * y
        xx = n * n - yy
        e = 0.0

        # cylinder
        if y < self.h:
            if xx > self.r[1]:
                x = np.sqrt(xx)
                x -= self.r[0]
                e = x * x

        # ends of cylinder
        else:
            y -= self.h
            yy = y * y

            # inside disc
            if xx < self.r[1]:
                e = yy
            else:
                x = np.sqrt(xx)
                x -= self.r[0]
                e = yy + x * x

        return e

    def Set(self, o, d):
        #   \brief Computes the pre-processing equations.
        #   \param o, d Ray origin and direction (which should be normalized).
        pa = self.c - o
        dx = d * self.Axis
        pax = pa * self.Axis
        dpa = d * pa
        quadric = [pa * pa - pax * pax,
                   2.0 * (dx * pax - dpa),
                   1.0 - dx*dx]
        linear = [-pax, dx]

    def R(self, t):
        #   \brief Compute the distance between a points on a line and a cylinder.
        #   The member function Cylinder::Set() should be called for pre-processing.
        #   \param t Parameter of the point on the line.
        y = linear[1] * t + linear[0]
        xx = (quadric[2] * t + quadric[1]) * t + quadric[0]
        y = np.abs(y)
        yy = y * y
        e = 0.0

        #  Cylinder
        if y < self.h:
            if xx > self.r[1]:
                x = np.sqrt(xx)
                x -= self.r[0]
                e = x * x
            # else
            #   e = 0.0

        else:
            # ends of cylinder
            y -= self.h
            yy = y * y

            if xx < self.r[1]:
                # inside disc
                e = yy
            else:
                x = sqrt(xx)
                x -= self.r[0]
                e = yy + x * x

        return e

class Cone:
    #   This class implements a truncated cone primitive.
    def __init__(self, a, b, ra, rb):
        '''
        Cone skeletal element, assume radius ra is greater than rb.
        :param a: end vertex of cone
        :param b: end vertex of cone
        :param ra: radius at a
        :param rb: radius at b
        '''
        self.axis = b - a
        self.ra = ra
        self.rb = rb
        self.rrb = rb * rb
        self.rra = ra * ra

        #  Compute the length of side of cone, i.e. its slant height
        self.conelength = np.sqrt((rb - ra) * (rb - ra) + self.axis.length * self.axis.length)

        #  Line segment
        self.side = np.array([rb - ra, self.axis.length, 0.0])
        self.side /= self.conelength

    def R(self, p):
        # Compute the distance between a point in space and a cone.
        #  Compute revolution coordinates
        n = p - self.a
        y = np.dot(n, self.axis)
        yy = y * y

        #  Squared radial distance to axis: postpone square root evaluation only when needed
        xx = np.dot(n, n) - yy
        e = 0.0

        #  Distance to large cap
        if y < 0.0:
            #  Disk : distance to plane cap
            if xx < self.rra:
                e = yy
            else:
                #  Distance to plane circle
                x = np.sqrt(xx) - self.ra
                e = x * x + yy

        #  Small cylinder test (optimization)
        elif xx < self.rrb:
            if y > self.length:
                e = y - self.length
                e *= e
                # inside cone
                # else
                #   e = 0.0
        else:
            #  Evaluate radial distance to axis
            x = np.sqrt(xx)
            #  Change frame, so that point is now on large cap
            x -= self.ra

            #  Distance to large cap
            if y < 0.0:
                #  Disk : distance to plane cap
                if x < 0.0:
                    e = yy
                else:
                    #  Distance to plane circle
                    e = x * x + yy

            else:
                #  Compute coordinates in the new rotated frame
                #  Postpone some computation that may not be needed in the following case
                ry = x * self.side[0] + y * self.side[1]
                if ry < 0.0:
                    e = x * x + y * y
                else:
                    rx = x * self.side[1] - y * self.side[0]
                    if ry > self.conelength:
                        ry -= self.conelength
                        e = rx * rx + ry * ry
                    else:
                        if rx > 0.0:
                            e = rx * rx
                    # else
                    #   e = 0.0
        return e

    def Set(self, o, d):
        #  Computes the pre-processing equations.
        # \param o, d Ray origin and direction (which should be normalized).
        pa = [a - o for a, o in zip(a, o)]
        dx = sum(d[i] * self.axis[i] for i in range(len(d)))
        pax = sum(pa[i] * self.axis[i] for i in range(len(pa)))
        dpa = sum(d[i] * pa[i] for i in range(len(d)))
        quadric = [sum(pa[i] * pa[i] for i in range(len(pa))) - pax * pax,
                   2.0 * (dx * pax - dpa),
                   1.0 - dx * dx]
        linear = [-pax, dx]

    def R(self, t):
        # Compute the distance between a points on a line and a cone.
        # The member function Cone::Set() should be called for pre-processing.
        # \param t  Parameter of the point on the line.

        y = self.linear[1] * t + self.linear[0]
        xx = (self.quadric[2] * t + self.quadric[1]) * t + self.quadric[0]
        yy = y * y
        e = 0.0

        #  Distance to large cap
        if y < 0.0:
            #  Disk : distance to plane cap
            if xx < self.rra:
                e = yy
            else:
                #  Distance to plane circle
                x = sqrt(xx) - self.ra
                e = x * x + yy
        #  Small cylinder test (optimization)
        elif xx < self.rrb:
            if y > self.length:
                e = y - self.length
                e *= e
                # inside cone
                # else
                #   e = 0.0
        else:
            #  Evaluate radial distance to axis
            x = np.sqrt(xx)
            # Change frame, so that point is now on large cap
            x -= self.ra
            # Distance to large cap
            if y < 0.0:
                # Disk : distance to plane cap
                if x < 0.0:
                    e = yy
                else:
                    # Distance to plane circle
                    e = x * x + yy
            else:
                #  Compute coordinates in the new rotated frame
                #  Postpone some computation that may not be needed in the following case
                ry = x * self.side[0] + y * self.side[1]
                if ry < 0.0:
                    e = x * x + yy
                else:
                    rx = x * self.side[1] - y * self.side[0]
                    if ry > self.conelength:
                        ry -= self.conelength
                        e = rx * rx + ry * ry
                    else:
                        if rx > 0.0:
                            e = rx * rx
        return e



