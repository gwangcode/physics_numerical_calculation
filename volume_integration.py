import numpy as np
from scipy.integrate import dblquad


__func_dict = {
    '1': [lambda x, y, z, x0, y0, z0: 0, lambda x, y, z, x0, y0, z0: 0, lambda x, y, z, x0, y0, z0: z], # f, g, h;

    '1/r': [ # potential
        lambda x, y, z, x0, y0, z0: 0.5*(x-x0)/((x-x0)**2 + (y-y0)**2 + (z-z0)**2)**0.5,
        lambda x, y, z, x0, y0, z0: 0.5*(y-y0)/((x-x0)**2 + (y-y0)**2 + (z-z0)**2)**0.5,
        lambda x, y, z, x0, y0, z0: 0.5*(z-z0)/((x-x0)**2 + (y-y0)**2 + (z-z0)**2)**0.5
    ],

    'x/r3': [ # field, (x-x0)/|r-r0|**3
            lambda x, y, z, x0, y0, z0: 1/((x-x0)**2 + (y-y0)**2 + (z-z0)**2)**0.5,
            lambda x, y, z, x0, y0, z0: 0,
            lambda x, y, z, x0, y0, z0: 0
        ],

    'y/r3': [ # field, (y-y0)/|r-r0|**3
                lambda x, y, z, x0, y0, z0: 0,
                lambda x, y, z, x0, y0, z0: 1/((x-x0)**2 + (y-y0)**2 + (z-z0)**2)**0.5,
                lambda x, y, z, x0, y0, z0: 0
            ],

    'z/r3': [ # field, (z-z0)/|r-r0|**3
                lambda x, y, z, x0, y0, z0: 0,
                lambda x, y, z, x0, y0, z0: 0,
                lambda x, y, z, x0, y0, z0: 1/((x-x0)**2 + (y-y0)**2 + (z-z0)**2)**0.5
            ],

    'xy': [ # moment of inertia, I_xy
        lambda x, y, z, x0, y0, z0: 0,
        lambda x, y, z, x0, y0, z0: 0.5*x*y*y,
        lambda x, y, z, x0, y0, z0: 0
    ],

    'xz': [ # moment of inertia, I_xz
        lambda x, y, z, x0, y0, z0: 0.5*x*x*z,
        lambda x, y, z, x0, y0, z0: 0,
        lambda x, y, z, x0, y0, z0: 0
    ],

    'yz': [ # moment of inertia, I_yz
        lambda x, y, z, x0, y0, z0: 0,
        lambda x, y, z, x0, y0, z0: 0.5*y*y*z,
        lambda x, y, z, x0, y0, z0: 0
    ],

    'x2': [ # moment of inertia, I_xx
        lambda x, y, z, x0, y0, z0: x**3/3,
        lambda x, y, z, x0, y0, z0: 0,
        lambda x, y, z, x0, y0, z0: 0
    ],

    'y2': [ # moment of inertia, I_yy
        lambda x, y, z, x0, y0, z0: 0,
        lambda x, y, z, x0, y0, z0: y**3/3,
        lambda x, y, z, x0, y0, z0: 0
    ],

    'z2': [ # moment of inertia, I_zz
        lambda x, y, z, x0, y0, z0: 0,
        lambda x, y, z, x0, y0, z0: 0,
        lambda x, y, z, x0, y0, z0: z**3/3
    ],

    'x': [ # dipole, x-x0
            lambda x, y, z, x0, y0, z0: 0.5*(x-x0)**2,
            lambda x, y, z, x0, y0, z0: 0,
            lambda x, y, z, x0, y0, z0: 0
        ],

    'y': [ # dipole, x-x0
                lambda x, y, z, x0, y0, z0: 0,
                lambda x, y, z, x0, y0, z0: 0.5*(y-y0)**2,
                lambda x, y, z, x0, y0, z0: 0
            ],

    'z': [ # dipole, x-x0
                lambda x, y, z, x0, y0, z0: 0,
                lambda x, y, z, x0, y0, z0: 0,
                lambda x, y, z, x0, y0, z0: 0.5*(z-z0)**2
            ]

}


# get the functions from the above dict
def __get_fgh(key, func_dict = __func_dict): return __func_dict[key]


def __unit_normals(hull):
    """
    Compute the unit normal vectors for all triangles defined by vertices and faces.
    """
    # Get the vertices for each face
    vertices = hull.vertices
    faces = hull.faces
    triangles = vertices[faces]

    # Compute the vectors for two edges of each triangle
    AB = triangles[:, 1] - triangles[:, 0]
    AC = triangles[:, 2] - triangles[:, 0]

    # Compute the cross product of the edges to get the normals
    normals = np.cross(AB, AC)

    # Normalize the vectors
    lengths = np.linalg.norm(normals, axis=1, keepdims=True)
    unit_normals = normals / lengths

    return unit_normals, lengths, AB, AC, triangles[:, 0]


def __normal_inside(mesh, unit_normals):
    """
    Check if the unit normals are directed outward and correct them if not.
    """
    # Compute the centroids of each triangle
    centroids = np.mean(mesh.vertices[mesh.faces], axis=1)

    # Compute a point slightly along the normal direction
    epsilon = 1e-6
    outside_points = centroids + epsilon * unit_normals

    # Check if these points are inside the mesh
    is_inside = mesh.contains(outside_points)

    return is_inside


def __vector_triangle_integration(normal, length, normal_inside, v1, v2, A, f, g, h, point):
    '''

    A, B and C are vertices of the triangle
    :param A: (x_A, y_A, z_A), (3,) array
    :param B: (x_B, y_B, z_B), (3,) array
    :param C: (x_C, y_C, z_C), (3,) array
    :param f: f(x, y, z) function object with inputs x, y and z
    :param g: g(x, y, z) function object with inputs x, y and z
    :param h: h(x, y, z) function object with inputs x, y and z
    :return:
    '''

    fx = lambda u, v: f(v2[0]*v + v1[0]*u + A[0], v2[1] * v + v1[1] * u + A[1], v2[2] * v + v1[2] * u + A[2], point[0], point[1], point[2])
    fy = lambda u, v: g(v2[0]*v + v1[0]*u + A[0], v2[1] * v + v1[1] * u + A[1], v2[2] * v + v1[2] * u + A[2], point[0], point[1], point[2])
    fz = lambda u, v: h(v2[0]*v + v1[0]*u + A[0], v2[1] * v + v1[1] * u + A[1], v2[2] * v + v1[2] * u + A[2], point[0], point[1], point[2])

    integrand = lambda v, u: (normal[0] * fx(u, v) + normal[1] * fy(u, v) + normal[2] * fz(u, v))*length

    u_min, u_max = 0, 1

    v_min = 0

    v_max = lambda u: 1-u

    result, error = dblquad(integrand, u_min, u_max, v_min, v_max)
    if normal_inside: result *= -1

    return result, error


def __integrate(normals, lengths, normals_inside, v1, v2, A, f, g, h, point):
    total = 0
    err_max = 0
    for i in np.arange(normals.shape[0]):
        result, error = __vector_triangle_integration(normals[i], lengths[i], normals_inside[i], v1[i], v2[i], A[i], f, g, h, point)
        total += result
        err_max = max(err_max, error)

    return total, err_max


def volume_integration(hulls_dict, point = np.zeros((3,)), kernel = '1'):
    '''
    Perform volumetric integration to hulls in hulls_dict. It actually performs the surface integration of the divergence of the field of the quantity.
    :param hulls_dict: dict of hull_dict:
        hulls = {
            'hull_name_1': hull_dict_1,
            'hull_name_2': hull_dict_2,
            ...
            'hull_name_n': hull_dict_n
        }

        hull_dict = {
                'hull': Hull object
                'density': float: (mass/charge) density of the hull
        }

        point: 1d array: the objective point that you want to calculate

        kernel: the choice of the kernel function (fgh function)
    :return: the integration result dict
            {
                'shape_name': {'integration': integration_result, 'error': integration_error},
                ...
            }
    '''

    result = {}

    for hull_name in hulls_dict:
        hull_dict = hulls_dict[hull_name]
        hull = hull_dict['hull']
        density = hull_dict['density']

        normals, lengths, v1, v2, A = __unit_normals(hull)
        normals_inside = __normal_inside(hull, normals)

        f, g, h = __get_fgh(kernel)
        total, error = __integrate(normals, lengths, normals_inside, v1, v2, A, f, g, h, point)

        result[hull_name] = {
            'integration': density * total,
            'error': error
        }

    return result








