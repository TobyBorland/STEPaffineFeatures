from math import cos, sin
import numpy as np


def EllipseLine(A, B, C, rx, ry, theta):
    # adapted from https://www.xarg.org/book/computer-graphics/line-segment-ellipse-intersection/
    # C: x,y centre point
    # line segment between points A & B
    # rx, ry: ellipse axes along local ellipse x, y, axes
    # theta: ellipse rotation about centre point
    #
    #     A = Vector2(A).sub(C).rotate(-theta);
    #     B = Vector2(B).sub(C).rotate(-theta);
    A = A[:2]
    B = B[:2]
    C = C[:2]

    rot = np.array([[cos(-theta), -sin(-theta)], [sin(-theta), cos(-theta)]])
    unrot = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])

    A = np.dot(rot, A - C)
    B = np.dot(rot, B - C)

    #     rx *= rx;
    #     ry *= ry;

    rxrx = rx * rx
    ryry = ry * ry

    ret = []
    v = B - A
    #v = v/np.linalg.norm(v)

    #     const a = rx * v.y * v.y + ry * v.x * v.x;
    #     const b = 2 * (rx * A.y * v.y + ry * A.x * v.x);
    #     const c = rx * A.y * A.y + ry * A.x * A.x - rx * ry;

    a = rxrx * v[1] * v[1] + ryry * v[0] * v[0]
    b = 2 * (rxrx * A[1] * v[1] + ryry * A[0] * v[0])
    c = rxrx * A[1] * A[1] + ryry * A[0] * A[0] - rxrx * ryry

    D = b * b - 4 * a * c  # Discriminant

    if D >= 0:
        sqrtD = np.sqrt(D)
        t1 = (-b + sqrtD) / (2 * a)
        t2 = (-b - sqrtD) / (2 * a)

        if 0 <= t1 and t1 <= 1:
            ret.append(np.dot(unrot, t1 * v + A) + C)

        if 0 <= t2 and t2 <= 1 and np.abs(t1 - t2) > 1e-16:  # eps
            ret.append(np.dot(unrot, t2 * v + A) + C)

    if len(ret) == 0:
        return None
    else:
        return ret


linePoint = np.array([np.sqrt(2), 1, 0])
lineDir = np.array([0, 1, 0])

rx = 2
ry = 4
theta = 3.141592/5
C = np.array([0, 0, 0])

# create 2 points A & B from orthogonal line intersection, of disp max(rx, ry)
A = linePoint + 2 * rx * lineDir
B = linePoint - 2 * rx * lineDir

dunno = EllipseLine(A, B, C, rx, ry, theta)

A = A[:2]



