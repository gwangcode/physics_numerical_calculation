from volume_integration import volume_integration
from conductor_potential import charge_distribution
from magnetic_integration import A_field, B_field, inductance_sections
from numpy import pi, array, newaxis, sum as sm
from numpy.linalg import norm

eps0 = 8.85418782e-12 # permittivity of free space
Ke = 1/(4*pi*eps0)
G = 6.6743e-11
B0 = 1e-7 # magnetic permittivity of free space / (4*pi)


# mass
def mass(hulls_dict):
    '''
    Integrate the volume of the objects in hulls_dict to calculate the masses for each of them.
    :param hulls_dict: The dict of the hulls (objects) in the following form
                        hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
    :return: {'object_name': {'integration': integration result (mass or charge), 'error': calculation error}}
    '''

    return volume_integration(hulls_dict, kernel = '1')


# moment of inertia
def momentum(hulls_dict):
    '''
    Integrate the momentum of inertia tensor of the objects in hulls_dict for each of them.
    :param hulls_dict: The dict of the hulls (objects) in the following form
                        hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
    :return: {'object_name': array([[ Ixx,  Ixy, Ixz],
                                    [ Iyx,  Iyy, Iyz],
                                    [ Izx,  Izy, Izz]])}
    '''
    Ix2 = volume_integration(hulls_dict, kernel='x2')
    Iy2 = volume_integration(hulls_dict, kernel='y2')
    Iz2 = volume_integration(hulls_dict, kernel='z2')
    Ixy = volume_integration(hulls_dict, kernel='xy')
    Ixz = volume_integration(hulls_dict, kernel='xz')
    Iyz = volume_integration(hulls_dict, kernel='yz')

    r = {}
    for obj in hulls_dict:
        r[obj] = array([
        [Iy2[obj]['integration'] + Iz2[obj]['integration'], -Ixy[obj]['integration'], -Ixz[obj]['integration']],
        [-Ixy[obj]['integration'], Ix2[obj]['integration'] + Iz2[obj]['integration'], -Iyz[obj]['integration']],
        [-Ixz[obj]['integration'], -Iyz[obj]['integration'], Ix2[obj]['integration'] + Iy2[obj]['integration']]
    ])

    return r


# potential at a given point (x, y, z)
def potential(hulls_dict, point, electric = True):
    '''
    Integrate the electric potential of the objects in hulls_dict for each of them at the point of (x, y, z): tuple, list or 1darray.
    :param hulls_dict: The dict of the hulls (objects) in the following form
                        hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
    :param point: the point of (x, y, z): tuple, list or 1darray
    :param electric: if True, it calculates the electric potential, else it calculates the gravitational potential.
    :return: {'object_name': {'integration': potential at the point, 'error': calculation error}}
    '''
    for obj in hulls_dict:
        if electric: hulls_dict[obj]['density'] *= Ke
        else: hulls_dict[obj]['density'] *= G
    return volume_integration(hulls_dict, point=point, kernel='1/r')


# field at a given point (x, y, z)
def field(hulls_dict, point, electric = True):
    '''
    Integrate the electric field vector of the objects in hulls_dict for each of them at the point of (x, y, z): tuple, list or 1darray.
    :param hulls_dict: the dict of the hulls (objects) in the following form
                        hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
    :param point: the point of (x, y, z): tuple, list or 1darray
    :param electric: if True, it calculates the electric field, else it calculates the gravitational field.
    :return: {'object_name': array([Ex, Ey, Ez])}
    '''
    for obj in hulls_dict:
        if electric: hulls_dict[obj]['density'] *= Ke
        else: hulls_dict[obj]['density'] *= G

    Ex = volume_integration(hulls_dict, point=point, kernel='x/r3')
    Ey = volume_integration(hulls_dict, point=point, kernel='y/r3')
    Ez = volume_integration(hulls_dict, point=point, kernel='z/r3')

    r = {}
    for obj in hulls_dict:
        r[obj] = array([Ex[obj]['integration'], Ey[obj]['integration'], Ez[obj]['integration']])

    return r


# calculate dipole of hulls with respect to a point
def dipole(hulls_dict, point, electric = True):
    '''
    Integrate the electric dipole vector of the objects in hulls_dict for each of them with respect to the point of (x, y, z): tuple, list or 1darray.
    :param hulls_dict: the dict of the hulls (objects) in the following form
                        hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'density': e.g. 1 # mass or charge density
                                        },
                                        ...
                                    }
    :param point: the point of (x, y, z): tuple, list or 1darray
    :param electric: if True, it calculates the electric dipole, else it calculates the gravitational dipole.
    :return: {'object_name': array([dipole in the x-direction, dipole in the y-direction, dipole in the z-direction])}
    '''
    for obj in hulls_dict:
        if electric: hulls_dict[obj]['density'] *= Ke
        else: hulls_dict[obj]['density'] *= G

    dipole_x = volume_integration(hulls_dict, point=point, kernel='x')
    dipole_y = volume_integration(hulls_dict, point=point, kernel='y')
    dipole_z = volume_integration(hulls_dict, point=point, kernel='z')

    r = {}
    for obj in hulls_dict:
        r[obj] = array([dipole_x[obj]['integration'], dipole_y[obj]['integration'], dipole_z[obj]['integration']])

    return r


def conductor_charge(hulls_dict, division = 0.1):
    '''
    Calculate the charge distribution on each conductor object in hulls_dict
    :param hulls_dict: the dict of the hulls (objects) in the following form
                        hulls_dict = {
                                        '[object_name]': { # pick your object name
                                        'hull': a hull object, # e.g.sphere()
                                        'potential': e.g. 1 # potential of the object
                                        },
                                        ...
                                    }
    :param division: the spacing between adjacent points, default = 0.1.
    :return: dict of the charge distribution:
                {'object_name':
                                {'charge': {the charge on every point},
                                'points': {array[x, y, z], ...} # the coordinates (x, y, z) of every point of calculation on the surface of the conductor.
                                'total': total charge on the conductor
                                }
                }
    '''

    r = {}
    cd = charge_distribution(hulls_dict, division)
    for obj in cd:
        r[obj] = cd[obj]
        cd[obj]['charge'] /= Ke
        r[obj]['total'] = cd[obj]['charge'].sum() # total charge

    return r


def conductor_potential(charge_dist, point, epsilon = 0.1):
    '''
    Calculate potential contribution at a point (x, y, z) of a conductor from the conductor_charge result
    :param charge_dist: the charge distribution of a conductor resulted from charge_distribution() calculation
    :param point: (x, y, z) of tuple, list or 1darray, the space point NOT on the surface of the conductor
    :param epsilon: a small value to prevent singularity of very small distance, default is 0.1. It is recommended to be equal to the value of the division of the hulls.
    :return: tuple of potential contribution,
            ( { 'object_name': potential,
                ...
                }, # potential contributed by each object

              total_potential
            )
    '''
    potential = {}
    total = 0
    for obj in charge_dist:
        delta_r = norm(array(point) - charge_dist[obj]['points'], axis=1)
        delta_r[delta_r < epsilon] = epsilon
        obj_potential = Ke * charge_dist[obj]['charge'] / delta_r # this needs norm!
        potential[obj] = obj_potential.sum()
        total += potential[obj]

    return potential, total


def conductor_field(charge_dist, point, epsilon = 0.1):
    '''
    Calculate electric field contribution at a point (x, y, z) of a conductor from the conductor_charge result
    :param charge_dist: the charge distribution of a conductor resulted from charge_distribution() calculation
    :param point: (x, y, z) of tuple, list or 1darray, the space point NOT on the surface of the conductor
    :param epsilon: a small value to prevent singularity of very small distance, default is 0.1. It is recommended to be equal to the value of the division of the hulls.
    :return: tuple of potential contribution,
            ( { 'object_name': field,
                ...
                }, # electric field contributed by each object

              total_electric_field
            )
    '''
    field = {}
    total = 0
    for obj in charge_dist:
        delta_r = array(point) - charge_dist[obj]['points']
        norms = norm(delta_r, axis = 1)
        norms[norms < epsilon] = epsilon
        delta_r2_vec = delta_r / norms[:, newaxis] ** 3
        obj_field = Ke * charge_dist[obj]['charge'][:, newaxis] * delta_r2_vec
        field[obj] = sm(obj_field, axis = 0)
        total += field[obj]

    return field, total


def mfield(dict_shape3d, point):
    '''
    Calculate the magnetic field B of a dict of 3d curves at the point (x, y, z)
    :param dict_shape3d: { 'curve_name': {'curve': shap3d object,
                           'current': electric current through the curve
                            },
                            ...
                        }
    :param point: (x, y, z) of tuple, list or 1darray, the space point.
    :return: ({'object_name': magnetic field at point,
                ...
              },
              total magnetic field
              )
    '''

    r = {}
    B_total = 0
    for obj in dict_shape3d:
        curve = dict_shape3d[obj]['curve'].vertices
        coeff = B0 * dict_shape3d[obj]['current']
        B = B_field(point, curve[:-1], curve[1:])

        B = coeff * B

        B_total += B
        r[obj] = B

    return r, B_total


def mpotential(dict_shape3d, point):
    '''
    Calculate the magnetic vector potential A of a dict of 3d curves at the point (x, y, z)
    :param dict_shape3d: { 'curve_name': {'curve': shap3d object,
                           'current': electric current through the curve
                            },
                            ...
                        }
    :param point: (x, y, z) of tuple, list or 1darray, the space point.
    :return: ({'object_name': magnetic potential at point,
                ...
              },
              total magnetic potential
              )
    '''

    r = {}
    A_total = 0
    for obj in dict_shape3d:
        curve = dict_shape3d[obj]['curve'].vertices
        coeff = B0 * dict_shape3d[obj]['current']
        A = A_field(point, curve[:-1], curve[1:])

        A = coeff * A

        A_total += A
        r[obj] = A

    return r, A_total


def mutual_inductance(shape3d_1, shape3d_2):
    '''
    Calculate mutual inductance between two wires
    :param shape3d_1: Wire 1 in Shape3d
    :param shape3d_2: Wire 2 in Shape3d
    :return: mutual inductance in float
    '''
    return B0 * inductance_sections(shape3d_1.vertices, shape3d_2.vertices)


def self_inductance(shape3d):
    '''
    Calculate mutual inductance of wire
    :param shape3d: Wire in Shape3d
    :return: self inductance in float
    '''
    return B0 * inductance_sections(shape3d.vertices, shape3d.vertices)








