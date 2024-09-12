import numpy as np


def __inverse_distance_matrix(points):
    """
    Construct an N x N matrix where the (i,j)-th entry is 1/r_ij,
    with r_ij being the Euclidean distance between points r_i and r_j.

    :param points: A (N, 3) ndarray representing N points in 3D space.
    :return: A (N, N) ndarray of inverse distances.
    """
    # Calculate the pairwise Euclidean distances
    diffs = points[:, np.newaxis, :] - points[np.newaxis, :, :]  # Shape: (N, N, 3)
    dists = np.linalg.norm(diffs, axis=2)  # Shape: (N, N)

    # Compute the inverse of the distances
    with np.errstate(divide='ignore'):  # Suppress divide by zero warnings
        inv_dists = 1.0 / dists

    # Set the diagonal to 0
    np.fill_diagonal(inv_dists, 0)

    return inv_dists


def __triangle_area(vertices):
    """
    Calculate the area of a triangle given its vertices.

    :param vertices: A (3, 3) ndarray representing the vertices of the triangle.
    :return: The area of the triangle.
    """
    A, B, C = vertices
    return 0.5 * np.linalg.norm(np.cross(B - A, C - A))


def __generate_points_in_triangle(vertices, n):
    """
    Generate n points evenly distributed inside a triangle in 3D space.

    :param vertices: A (3, 3) ndarray representing the vertices of the triangle.
    :param n: Number of points to generate.
    :return: A (n, 3) ndarray of points inside the triangle.
    """
    A, B, C = vertices

    # Estimate the number of points per edge
    num_points_per_edge = int(np.ceil(np.sqrt(n)))
    #total_points = num_points_per_edge * (num_points_per_edge + 1) // 2

    # Generate barycentric coordinates
    u = np.linspace(0, 1, num_points_per_edge)
    v = np.linspace(0, 1, num_points_per_edge)
    U, V = np.meshgrid(u, v)
    U = U.flatten()
    V = V.flatten()

    # Filter points inside the triangle
    mask = U + V <= 1
    U = U[mask]
    V = V[mask]

    # Calculate the corresponding Cartesian coordinates
    W = 1 - U - V
    points = U[:, np.newaxis] * A + V[:, np.newaxis] * B + W[:, np.newaxis] * C

    # If more points are generated, sample from them
    if points.shape[0] > n:
        points = points[:n]

    return points


def __generate_points_in_triangles(triangles, division):
    """
    Generate points evenly distributed inside a collection of triangles using vectorization.

    :param triangles: A (N, 3, 3) ndarray representing N triangles with their vertices.
    :param division: A factor to determine the number of points based on the area of the triangle.
    :return: A (M, 3) ndarray of points inside all triangles.
    """
    N = triangles.shape[0]

    # Calculate the area of each triangle
    areas = np.array([__triangle_area(triangles[i]) for i in range(N)])

    # Determine the number of points per triangle based on the area
    n_points_per_triangle = (areas / division**2).astype(int)

    # Calculate the total number of points
    total_points = np.sum(n_points_per_triangle)

    # Initialize an empty list to hold the points
    all_points = []

    # Generate points for each triangle
    for i in range(N):
        vertices = triangles[i]
        points = __generate_points_in_triangle(vertices, n_points_per_triangle[i])
        all_points.append(points)

    # Combine all points into a single array
    all_points = np.vstack(all_points)

    return all_points


def charge_distribution(hulls_dict, division = 0.1):
    '''
    Calculate charge on the surfaces of hulls in hull dict
    :param hulls_dict: dict of hull_dict:
        hulls = {
            'hull_name_1': hull_dict_1,
            'hull_name_2': hull_dict_2,
            ...
            'hull_name_n': hull_dict_n
        }

        hull_dict = {
                'hull': Hull object
                'potential': float: potential of the object
        }
    :return: a 1d array of charges for each point on the surfaces of hulls. ???
    '''

    hull_points = None
    potentials = []
    lengths = []
    hull_names = []

    # iterate each hull
    for hull_name in hulls_dict:
        hull_dict = hulls_dict[hull_name]
        hull = hull_dict['hull']
        potential = hull_dict['potential']

        vertices = hull.vertices
        faces = hull.faces
        triangle_vertices = vertices[faces]

        hpoints = __generate_points_in_triangles(triangle_vertices, division)
        hpoints = np.unique(hpoints, axis = 0)

        if hull_points is None: hull_points = hpoints
        else: hull_points = np.concatenate([hull_points, hpoints], axis = 0)

        number_of_points = hpoints.shape[0]
        lengths.append(number_of_points)
        hull_names.append(hull_name)
        potentials += [potential]*number_of_points

    matrix = __inverse_distance_matrix(hull_points)
    q = np.linalg.inv(matrix) @ np.array(potentials) # charge on every point

    if q.size > 0:
        result = {}
        begin_index = 0
        for hull_name, length in zip(hull_names, lengths):
            result[hull_name] = {
                'charge': q[begin_index:begin_index+length],
                'points': hull_points[begin_index:begin_index+length]
            }
            begin_index += length

        return result





