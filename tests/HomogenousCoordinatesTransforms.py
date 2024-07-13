# Functions to convert points etc from one coordinate system to another
#
# These functions are based on a 4x4 homogeneous matrices or
# homogeneous coordinates introduced by August Ferdinand Möbius
# See: https://en.wikipedia.org/wiki/Homogeneous_coordinates
#
# Other work I found useful in developing this code includes:
# University of Illinois Robotics course http://motion.cs.illinois.edu/RoboticSystems/CoordinateTransformations.html

import numpy as np


def unit_vector(xyz_array):
    """
    Compute unit vector of vector from origin to specified xyz point
    :param xyz_array: numpy array of length 3. x = 1st element, y = 2nd element, z = 3rd element
    :return: xyz_unit : unit vector (ie ijk) in specified direction
    """
    num_xyz = len(xyz_array)
    if num_xyz != 3:
        raise ValueError('Input array to "unit_vector" must have length of 3.')

    xyz_unit = xyz_array / np.linalg.norm(xyz_array)

    return xyz_unit


def create_rotation_matrix_from_points(local_origin, local_x_dir, local_xy_plane):
    """
    Create a 4x4 homogeneous rotation matrix
    :param local_origin: np_array of x,y,z wrt global cs
    :param local_x_dir:  np_array of x,y,z point. vector from local_origin to this point defines X
    :param local_xy_plane: point in xy plane
    :return: rot_matrix: homogeneous rotation matrix
    """
    local_x_vec = local_x_dir - local_origin
    local_unit_x = unit_vector(local_x_vec)

    local_xy_vec = local_xy_plane - local_origin
    local_z_vec = np.cross(local_x_vec, local_xy_vec)
    local_unit_z = unit_vector(local_z_vec)

    local_y_vec = np.cross(local_z_vec, local_x_vec)
    local_unit_y = unit_vector(local_y_vec)

    # print('local_unit_x = {}'.format(local_unit_x))
    # print('local_unit_y = {}'.format(local_unit_y))
    # print('local_unit_z = {}'.format(local_unit_z))

    rot_matrix = np.zeros((4, 4))
    rot_matrix[3, 3] = 1.0
    rot_matrix[0:3, 0] = local_unit_x
    rot_matrix[0:3, 1] = local_unit_y
    rot_matrix[0:3, 2] = local_unit_z
    determinant = np.linalg.det(rot_matrix)
    assert np.isclose(determinant, 1.0)

    # print('rot_matrix = \n{}'.format(rot_matrix))
    # print('determinant = {}\n'.format(determinant))

    return rot_matrix


def create_translation_matrix_from_point(point, origin=np.zeros(3)):
    """
    Create a 4x4 homogeneous translation matrix
    :param point: np_array of x,y,z wrt global cs
    :param origin:  np_array of x,y,z point. vector from local_origin to this point defines X
    :return: translation_matrix : homogeneous translation matrix
    """
    translation_matrix = np.identity(4)
    deltas = point - origin
    translation_matrix[0:3, 3] = deltas

    # print('translation_matrix = \n{}'.format(translation_matrix))

    return translation_matrix


def invert_homogeneous_transformation_matrix(transformation_matrix):
    """
    Invert a homogeneous transformation matrix
    :param transformation_matrix: homogeneous transformation matrix
    :return: inverse_transform: inverted matrix
    """

    rot_matrix = transformation_matrix[0:3, 0:3]
    translation = transformation_matrix[0:3, 3]
    # for orthogonal arrays the transpose is equal to the inverse
    rot_inverse = rot_matrix.transpose()

    trans_inv = -rot_inverse.dot(translation)

    inverse_transform = np.identity(4)
    inverse_transform[0:3, 0:3] = rot_inverse
    inverse_transform[0:3, 3] = trans_inv

    return inverse_transform


def calculate_homogeneous_transforms(local_origin, local_x_dir, local_xy_plane):

    rot_matrix_4x4 = create_rotation_matrix_from_points(local_origin, local_x_dir, local_xy_plane)
    translation_matrix = create_translation_matrix_from_point(local_origin)

    # This is the key step showing how the transformation_matrix is created
    # using matrix multiplication of the translation and rotation matrices.
    # Order is CRITICAL:
    #     translation_matrix.dot(rot_matrix_4x4) IS NOT EQUAL TO rot_matrix_4x4.dot(translation_matrix)
    #     (except it trivial cases where it provides the correct answer but is still wrong)
    transformation_matrix = translation_matrix.dot(rot_matrix_4x4)
    inverse_transform = invert_homogeneous_transformation_matrix(transformation_matrix)

    return transformation_matrix, inverse_transform


def test_pure_rotation():
    print('{}'.format('-'*80))
    print('testing test_pure_rotation')

    local_origin = np.asarray([0.0, 0.0, .0])
    local_x_dir = np.asarray([0.0,  1.0, 0.0])
    local_xy_plane = np.asarray([-1.0, 1.0, 0.0])

    print('    local_origin = {}'.format(local_origin))
    print('    local_x_dir = {}'.format(local_x_dir))
    print('    local_xy_plane = {}'.format(local_xy_plane))

    transformation_matrix, inverse_transform = calculate_homogeneous_transforms(local_origin,
                                                                                local_x_dir,
                                                                                local_xy_plane)

    tm_str = '     {}'.format(transformation_matrix)
    tm_str = tm_str.replace('\n', '\n     ')
    print('\n    transformation_matrix = \n{}'.format(tm_str))

    point = np.asarray([1.0, 1.0, 2, 1.0])

    point_local = inverse_transform.dot(point)
    print('\n    point {} in local cs = {}'.format(point, point_local))

    point_global = transformation_matrix.dot(point_local)
    print('    local point {} in global cs = {}\n'.format(point_local, point_global))

    assert np.isclose(point, point_global).all()
    assert np.isclose(point_local, [1.0, -1, 2, 1.]).all()

    print('    Successfully completed test of test_pure_rotation\n')
    return None


def test_pure_translation():
    print('{}'.format('-'*80))
    print('testing test_pure_translation')

    local_origin = np.asarray([0.0, 0.0, 1.0])
    local_x_dir = np.asarray([1.0,  0, 1.0])
    local_xy_plane = np.asarray([0.0, 1.0, 1.0])

    print('    local_origin = {}'.format(local_origin))
    print('    local_x_dir = {}'.format(local_x_dir))
    print('    local_xy_plane = {}'.format(local_xy_plane))

    transformation_matrix, inverse_transform = calculate_homogeneous_transforms(local_origin,
                                                                                local_x_dir,
                                                                                local_xy_plane)

    tm_str = '     {}'.format(transformation_matrix)
    tm_str = tm_str.replace('\n', '\n     ')
    print('\n    transformation_matrix = \n{}'.format(tm_str))
    point = np.asarray([1.0, 1.0, 0, 1.0])

    point_local = inverse_transform.dot(point)
    print('\n    point {} in local cs = {}'.format(point, point_local))

    point_global = transformation_matrix.dot(point_local)
    print('    local point {} in global cs = {}\n'.format(point_local, point_global))

    assert np.isclose(point, point_global).all()
    assert np.isclose(point_local, [1.0, 1, -1, 1.]).all()

    print('    Successfully completed test of test_pure_translation\n')
    return None


def test_rotation_and_translation():
    print('{}'.format('-'*80))
    print('testing test_rotation_and_translation')

    local_origin = np.asarray([1.0, 1.0, 1.0])
    local_x_dir = np.asarray([.0,  1, 1.0])
    local_xy_plane = np.asarray([1.0, 2.0, 1.0])

    print('    local_origin = {}'.format(local_origin))
    print('    local_x_dir = {}'.format(local_x_dir))
    print('    local_xy_plane = {}'.format(local_xy_plane))

    transformation_matrix, inverse_transform = calculate_homogeneous_transforms(local_origin,
                                                                                local_x_dir,
                                                                                local_xy_plane)

    tm_str = '     {}'.format(transformation_matrix)
    tm_str = tm_str.replace('\n', '\n     ')
    print('\n    transformation_matrix = \n{}'.format(tm_str))

    point = np.asarray([-1.0, 2.0, 2, 1.0])
    # point = np.asarray([1.0, 1.0, 1, 1.0])
    # point = np.asarray([0.0, 0.0, 1, 1.0])

    point_local = inverse_transform.dot(point)
    print('\n    point {} in local cs = {}'.format(point, point_local))

    point_global = transformation_matrix.dot(point_local)
    print('    local point {} in global cs = {}\n'.format(point_local, point_global))

    assert np.isclose(point, point_global).all()
    assert np.isclose(point_local, [2.0, 1, -1, 1.]).all()

    print('    Successfully completed test of test_rotation_and_translation\n')
    return None


if __name__ == '__main__':
    test_pure_rotation()
    test_pure_translation()
    test_rotation_and_translation()
    print('')