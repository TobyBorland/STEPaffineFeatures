# Point on ellipse (including interior points) closest to target point "p".
import numpy as np

from collections import namedtuple  # PyCharm bug will flag this as type error
Point = namedtuple("Point", "x y z")
# eps machine precision constant
eps = np.finfo(float).eps

p = np.array([8., 8., 20.])
eCentre = np.array([0.0, 0.0, 0.0])
eCentre = np.array([5.0, 0.0, 0.0])
eFeatureAxis = np.array([1.0, 0.0, 0.2])
eLocalXaxis = np.array([0.0, 1.0, 0.2])
eMajorRadius = 20.
eMinorRadius = 10.

# #Coord3d local_coord = new Coord3d(this.Center, this._v1, this._v2);
#     # /// Initializes coordinate system using origin point and two vectors.
#     # /// <param name="p">Origin of the coordinate system.</param>
#     # /// <param name="v1">Vector oriented along the X axis.</param>
#     # /// <param name="v2">Vector in the XY plane.</param>
#     public Coord3d(Point3d p, Vector3d v1, Vector3d v2, string name = "")
#     {
#         v1 = v1.ConvertToGlobal().Normalized;
#
#             v = _coord.Axes.Transpose() * v
#             v._coord = Coord3d.GlobalCS
#
#         Vector3d v3 = v1.Cross(v2).Normalized;
#         v2 = v3.Cross(v1).Normalized;
#
#         _origin = p.ConvertToGlobal();
#         _axes = new Matrix3d(v1, v2, v3);
#     }

# public Coord3d(Point3d p, double[] d1, double[] d2, string name = "")
#         {
#             Vector3d v1 = new Vector3d(d1);
#             Vector3d v2 = new Vector3d(d2);
#             if (v1.IsParallelTo(v2))
#             {
#                 throw new Exception("Vectors are parallel");
#             }
#
#             v1 = v1.Normalized;
#             Vector3d v3 = v1.Cross(v2).Normalized;
#             v2 = v3.Cross(v1).Normalized;
#
#             _origin = p.ConvertToGlobal();
#             _axes = new Matrix3d(v1, v2, v3);

#         }

#
#     # /// <summary>
#     # /// Convert vector to global coordinate system
#     # /// </summary>
#     public Vector3d ConvertToGlobal()
#     {
#         if (_coord == null || object.ReferenceEquals(_coord, Coord3d.GlobalCS))
#         {
#             return this.Copy();
#         }
#         else
#         {
#             Vector3d v = this.Copy();
#             v = _coord.Axes.Transpose() * v;
#             v._coord = Coord3d.GlobalCS;
#             return v;
#         }
#     }

#p = p.ConvertTo(local_coord);


# #def create_rotation_matrix_from_points(local_origin, local_x_dir, local_xy_plane):
#     # """
#     # Create a 4x4 homogeneous rotation matrix
#     # :param local_origin: np_array of x,y,z wrt global cs
#     # :param local_x_dir:  np_array of x,y,z point. vector from local_origin to this point defines X
#     # :param local_xy_plane: point in xy plane
#     # :return: rot_matrix: homogeneous rotation matrix
#     # """
#     local_x_vec = local_x_dir - local_origin
#     local_unit_x = unit_vector(local_x_vec) / np.linalg.norm(xyz_array)
#
#     local_xy_vec = local_xy_plane - local_origin
#     local_z_vec = np.cross(local_x_vec, local_xy_vec)
#     local_unit_z = unit_vector(local_z_vec)
#
#     local_y_vec = np.cross(local_z_vec, local_x_vec)
#     local_unit_y = unit_vector(local_y_vec)
#
#     # print('local_unit_x = {}'.format(local_unit_x))
#     # print('local_unit_y = {}'.format(local_unit_y))
#     # print('local_unit_z = {}'.format(local_unit_z))
#
#     rot_matrix = np.zeros((4, 4))
#     rot_matrix[3, 3] = 1.0
#     rot_matrix[0:3, 0] = local_unit_x
#     rot_matrix[0:3, 1] = local_unit_y
#     rot_matrix[0:3, 2] = local_unit_z
#     determinant = np.linalg.det(rot_matrix)
#     assert np.isclose(determinant, 1.0)
#
#     # print('rot_matrix = \n{}'.format(rot_matrix))
#     # print('determinant = {}\n'.format(determinant))
#
#     return rot_matrix

#def pointEllipseMinDisp(p, eCentre, eLocalYaxis, eLocalXaxis, eMajorRadius, eMinorRadius):
    # Algorithm by Dr. Robert Nurnberg
    # http://wwwf.imperial.ac.uk/~rn/distance2ellipse.pdf
    # Does not work for interior points

    # Move origin to ellipse centre, and rotate point coordinate to local coordinates based on ellipse major & minor axis.
    # Assuming direction_ratios of AXIS2_PLACEMENT_3D are normalised

    # eLocalXaxis = arrayTypeCast(eLocalXaxis)
    # eLocalYaxis = arrayTypeCast(eLocalYaxis)
    # eCentre = arrayTypeCast(eCentre)


# eNorm = np.cross(eLocalXaxis, eLocalYaxis)
# U = eLocalXaxis/np.linalg.norm(eLocalXaxis)
# V = np.cross(eNorm, eLocalXaxis)# - eCentre)
# V = V/np.linalg.norm(V)
# eNorm = eNorm / np.linalg.norm(eNorm)

# local axes assumed to start from origin

# U = eLocalXaxis # - eCentre
# U = U / np.linalg.norm(U)
#
# V = eLocalYaxis # - eCentre
# V = V / np.linalg.norm(V)
#
# eNorm = np.cross(U, V)
# eNorm = eNorm / np.linalg.norm(eNorm)

#         Vector3d v3 = v1.Cross(v2).Normalized;
#         v2 = v3.Cross(v1).Normalized;

#             Vector3d v1 = new Vector3d(d1);
#             Vector3d v2 = new Vector3d(d2);
#             if (v1.IsParallelTo(v2))
#             {
#                 throw new Exception("Vectors are parallel");
#             }
#
#             v1 = v1.Normalized;
#             Vector3d v3 = v1.Cross(v2).Normalized;
#             v2 = v3.Cross(v1).Normalized;

# https://nurnberg.maths.unitn.it/distance2ellipse.pdf

U = eFeatureAxis # - eCentre
U = U / np.linalg.norm(U)

eNorm = np.cross(U, eLocalXaxis)
eNorm = eNorm / np.linalg.norm(eNorm)

V = np.cross(eNorm, U)
V = V / np.linalg.norm(V)

rotM = np.array([U, V, eNorm])
pUV = np.matmul(rotM, p - eCentre)

# project to ellipse plane
pUV = np.array([pUV[0], pUV[1], 0])

if (np.isclose(p[0], 0) and np.isclose(p[1], 0)):   # Center point
   print(eCentre)

theta = np.arctan2(eMajorRadius * pUV[1], eMinorRadius * pUV[0])
i = 0
n0 = pUV
radiusConst = eMajorRadius**2 - eMinorRadius**2
pUVdisp = np.linalg.norm(pUV)
strike = 0

while (i < 100) and not strike:
    i += 1
    f = radiusConst * np.cos(theta) * np.sin(theta) - \
        pUV[0] * eMajorRadius * np.sin(theta) + \
        pUV[1] * eMinorRadius * np.cos(theta)

    f_ = radiusConst * (np.cos(theta)**2 - np.sin(theta)**2) - \
         pUV[0] * eMajorRadius * np.cos(theta) - \
         pUV[1] * eMinorRadius * np.sin(theta)

    theta = theta - f / f_

    n = np.array([eMajorRadius * np.cos(theta), eMinorRadius * np.sin(theta), 0])

    if np.allclose(n0, n): # tolerance = 1E-12
        strike += 1
        if np.linalg.norm(n) > pUVdisp: # test if p < n, inside ellipse
            #print(pUV - eCentre)
            nXY = np.matmul(rotM.T, pUV) + eCentre
        else:
            #print(np.matmul(rotM.T, n) + eCentre)
            nXY = np.matmul(rotM.T, n) + eCentre
    n0 = n

# furthest point on ellipse must be edge point on line that contains ellipse centre and closest point
nn = np.array([-n[0], -n[1], 0])
nnXY = np.matmul(rotM.T, nn) + eCentre

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import mpl_toolkits.mplot3d.art3d as art3d


def rotation_matrix(d):
    """
    Calculates a rotation matrix given a vector d. The direction of d
    corresponds to the rotation axis. The length of d corresponds to
    the sin of the angle of rotation.

    Variant of: http://mail.scipy.org/pipermail/numpy-discussion/2009-March/040806.html
    """
    sin_angle = np.linalg.norm(d)

    if sin_angle == 0:
        return np.identity(3)

    d /= sin_angle

    eye = np.eye(3)
    ddt = np.outer(d, d)
    skew = np.array([[    0,  d[2],  -d[1]],
                  [-d[2],     0,  d[0]],
                  [d[1], -d[0],    0]], dtype=np.float64)

    M = ddt + np.sqrt(1 - sin_angle**2) * (eye - ddt) + sin_angle * skew
    return M

def pathpatch_2d_to_3d(pathpatch, z = 0, normal = 'z'):
    """
    Transforms a 2D Patch to a 3D patch using the given normal vector.

    The patch is projected into they XY plane, rotated about the origin
    and finally translated by z.
    """
    if type(normal) is str: #Translate strings to normal vectors
        index = "xyz".index(normal)
        normal = np.roll((1.0,0,0), index)

    normal /= np.linalg.norm(normal)

    path = pathpatch.get_path()
    trans = pathpatch.get_patch_transform()

    path = trans.transform_path(path)

    pathpatch.__class__ = art3d.PathPatch3D
    pathpatch._code3d = path.codes #Copy the codes
    pathpatch._facecolor3d = pathpatch.get_facecolor #Get the face color

    verts = path.vertices #Get the vertices in 2D

    # U = np.cross(normal, np.array([0., 1., 0.]))
    # U = U / np.linalg.norm(U)
    #
    # V = np.cross(normal, U)
    # V = V / np.linalg.norm(V)
    #
    # rotM = np.array([U, V, normal])
    # A = np.array([np.matmul(rotM, np.array([x, y, 0])) + np.array([0, 0, z]) for x, y in verts])

    d = np.cross(normal, (0, 0, 1))
    M = rotation_matrix(d)

    pathpatch._segment3d = np.array([np.dot(M, (x, y, 0)) + (0, 0, z) for x, y in verts])
    #pathpatch._segment3d = A

    #return(A)

ax = plt.figure().add_subplot(projection='3d')
ax.set_xlim(-20, 20)
ax.set_ylim(-20, 20)
ax.set_zlim(-20, 20)

ax.plot3D((p[0],), (p[1],), (p[2],), label="p", marker="x", markersize=5)
ax.plot3D((pUV[0], n[0]), (pUV[1], n[1]), (pUV[2], n[2]), label="n", marker="o", markersize=5, color='red')
#ax.plot3D((p[0], pUV[0]), (p[1], pUV[1]), (p[2], pUV[2]), 'gray')
ax.plot3D((nXY[0], nnXY[0]), (nXY[1], nnXY[1]), (nXY[2], nnXY[2]), 'blue')
ax.plot3D((nXY[0], p[0]), (nXY[1], p[1]), (nXY[2], p[2]), 'gray')
ax.plot3D((eCentre[0], 20*eNorm[0]), (eCentre[1], 20*eNorm[1]), (eCentre[2], 20*eNorm[2]), 'green')

ellps = Ellipse((eCentre[0], eCentre[1]), width=2*eMajorRadius, height=2*eMinorRadius, angle=0, ec='black', fc=None, fill=False)
ax.add_patch(ellps)
#art3d.pathpatch_2d_to_3d(ellps, z=eCentre[2], zdir=(eNorm[0], eNorm[1], eNorm[2]))
pathpatch_2d_to_3d(ellps, z=eCentre[2], normal=(eNorm[0], eNorm[1], eNorm[2]))
#ellps2 = pathpatch_2d_to_3d(ellps, z=eCentre[2], normal=(eNorm[0], eNorm[1], eNorm[2]))
#ax.plot(ellps2[:,0],ellps2[:,1],ellps2[:,2])
ellps = Ellipse((eCentre[0], eCentre[1]), width=2*eMajorRadius, height=2*eMinorRadius, angle=0, ec='red', fc=None, fill=False)
ax.add_patch(ellps)
#art3d.pathpatch_2d_to_3d(ellps, z=eCentre[2], zdir=(eNorm[0], eNorm[1], eNorm[2]))
pathpatch_2d_to_3d(ellps, z=eCentre[2], normal=(eCentre[0], eCentre[1], eCentre[2]+1))
#ellps3 = pathpatch_2d_to_3d(ellps, z=eCentre[2], normal=(eCentre[0], eCentre[1], eCentre[2]+1))
#ax.plot(ellps3[:,0],ellps3[:,1],ellps3[:,2])
plt.show()
