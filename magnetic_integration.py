import numpy as np
from scipy.integrate import dblquad


def B_field(r0, r1_array, r2_array):
    """
    Calculate the magnetic field at a point r0 due to multiple straight wire segments
    each defined by corresponding points in r1_array and r2_array carrying current I,
    using the Biot-Savart law.

    Parameters:
    I        : float
               Current in the wire (in amperes)
    r0       : numpy array of shape (3,)
               Position vector of the point where the magnetic field is calculated (x0, y0, z0)
    r1_array : numpy array of shape (n, 3)
               Array of position vectors for the first endpoints of the wire segments (x1, y1, z1)
    r2_array : numpy array of shape (n, 3)
               Array of position vectors for the second endpoints of the wire segments (x2, y2, z2)

    Returns:
    B        : numpy array of shape (3,)
               Magnetic field vector at point r0 (Bx, By, Bz)
    """

    # Initialize the magnetic field vector
    B_total = np.array([0.0, 0.0, 0.0])

    # Loop through each segment
    for r1, r2 in zip(r1_array, r2_array):
        # Differential element dl
        dl = r2 - r1
        L = np.linalg.norm(dl)
        dl_unit = dl / L

        # Midpoint of the segment
        midpoint = (r1 + r2) / 2.0

        # Displacement vector r' from the midpoint to the observation point r0
        r_prime = r0 - midpoint

        # Magnitude of the displacement vector
        r_prime_mag = np.linalg.norm(r_prime)

        # Cross product of dl and r'
        cross_product = np.cross(dl_unit, r_prime)

        # Integration result (analytically calculated for a straight segment)
        factor = np.dot(dl_unit, r_prime)
        int_result = L / r_prime_mag / (r_prime_mag - factor)

        # Magnetic field due to this segment
        B_segment = int_result * cross_product / r_prime_mag ** 2

        # Add the contribution of this segment to the total magnetic field
        B_total += B_segment

    return B_total


def A_field(r0, r1_array, r2_array):
    """
    Calculate the vector potential at a point r0 due to multiple straight wire segments
    each defined by corresponding points in r1_array and r2_array carrying current I,
    using the Biot-Savart law for vector potential.

    Parameters:
    I        : float
               Current in the wire (in amperes)
    r0       : numpy array of shape (3,)
               Position vector of the point where the vector potential is calculated (x0, y0, z0)
    r1_array : numpy array of shape (n, 3)
               Array of position vectors for the first endpoints of the wire segments (x1, y1, z1)
    r2_array : numpy array of shape (n, 3)
               Array of position vectors for the second endpoints of the wire segments (x2, y2, z2)

    Returns:
    A        : numpy array of shape (3,)
               Vector potential at point r0 (Ax, Ay, Az)
    """

    # Calculate differential elements dl
    dl = r2_array - r1_array  # Shape (n, 3)
    L = np.linalg.norm(dl, axis=1)  # Shape (n,)
    dl_unit = dl / L[:, np.newaxis]  # Shape (n, 3)

    # Calculate midpoints of the segments
    midpoints = (r1_array + r2_array) / 2.0  # Shape (n, 3)

    # Displacement vectors r' from the midpoints to the observation point r0
    r_prime = r0 - midpoints  # Shape (n, 3)
    r_prime_mag = np.linalg.norm(r_prime, axis=1)  # Shape (n,)

    # Dot product of dl_unit and r_prime
    dot_product = np.einsum('ij,ij->i', dl_unit, r_prime)  # Shape (n,)

    # Integration result (vectorized)
    int_result = np.log((r_prime_mag + dot_product + L / 2) /
                        (r_prime_mag + dot_product - L / 2))  # Shape (n,)

    # Vector potential due to each segment
    A_segments = int_result[:, np.newaxis] * dl_unit  # Shape (n, 3)

    # Sum all contributions to get the total vector potential
    A_total = np.sum(A_segments, axis=0)  # Shape (3,)

    return A_total


def __integrand(t1, t2, dl1, dl2, segment1, segment2):
    r1 = segment1[0] + t1 * dl1
    r2 = segment2[0] + t2 * dl2
    r_diff = r1 - r2
    r_mag = np.linalg.norm(r_diff)
    return np.dot(dl1, dl2) / r_mag


def mutual_inductance_segment(segment1, segment2):
    """Calculate the mutual inductance between two segments."""
    dl1 = segment1[1] - segment1[0]
    dl2 = segment2[1] - segment2[0]

    # Perform the double integration
    M, _ = dblquad(lambda t1, t2: __integrand(t1, t2, dl1, dl2, segment1, segment2), 0, 1, lambda t1: 0, lambda t1: 1)

    return M # * (mu_0 / (4 * np.pi))


def inductance_sections(r1_array, r2_array):
    """
    Calculate the mutual inductance between two sections, each composed of multiple straight wire segments.
    If two sections are the same (r1_array = r2_array), calculate the self induction.

    Parameters:
    r1_array : numpy array of shape (n, 3)
               Array of position vectors for the first endpoints of the wire segments in section 1.
    r2_array : numpy array of shape (m, 3)
               Array of position vectors for the second endpoints of the wire segments in section 2.

    Returns:
    M_total : float
              Total mutual inductance between the two sections (in henries).
    """

    # Calculate the total mutual inductance by summing over all pairs of segments
    M_total = 0
    for i in range(len(r1_array) - 1):
        for j in range(len(r2_array) - 1):
            segment1 = r1_array[i:i + 2]
            segment2 = r2_array[j:j + 2]
            if not np.array_equal(segment1, segment2):
                M_total += mutual_inductance_segment(segment1, segment2)

    return M_total