import numpy as np

hkl_plane = [3, 1, 0]

normalise_G_vectors_to_unit_sphere = True

gaussian_scatter_width = 2.0/2.355

sample_normal = [0.0, 0.0, 1.0]

k_inc = np.array([-1.0, 0.0, -1.0])

wavelength = 1.378e-10

a_lattice = 3.3e-10

num_scattering_vectors = 200

num_points_per_vector = 500

num_bins = 8E2

plot_stereographs = True

plot_quivers = True


def find_family_of_vectors(hkl_plane, normalise):

    import numpy as np

    h, k, l = hkl_plane

    degenerate_family_of_vectors = [hkl_plane]

    # The following includes all permutations of the original hkl to the empty list above.

    new_plane = [k, l, h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [h, l, k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [l, k, h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [k, h, l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [l, h, k]

    degenerate_family_of_vectors.append(new_plane)

    # The following adds single negative values to the list.

    new_plane = [-h, k, l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-k, l, h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-h, l, k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-l, k, h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-k, h, l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-l, h, k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [h, -k, l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [k, -l, h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [h, -l, k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [l, -k, h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [k, -h, l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [l, -h, k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [h, k, -l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [k, l, -h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [h, l, -k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [l, k, -h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [k, h, -l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [l, h, -k]

    degenerate_family_of_vectors.append(new_plane)

    # The following adds double negatives to the list.

    new_plane = [-h, -k, l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-k, -l, h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-h, -l, k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-l, -k, h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-k, -h, l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-l, -h, k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [h, -k, -l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [k, -l, -h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [h, -l, -k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [l, -k, -h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [k, -h, -l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [l, -h, -k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-h, k, -l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-k, l, -h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-h, l, -k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-l, k, -h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-k, h, -l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-l, h, -k]

    degenerate_family_of_vectors.append(new_plane)

    # The following adds triple negatives to the list.

    new_plane = [-h, -k, -l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-k, -l, -h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-h, -l, -k]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-l, -k, -h]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-k, -h, -l]

    degenerate_family_of_vectors.append(new_plane)

    new_plane = [-l, -h, -k]

    degenerate_family_of_vectors.append(new_plane)

    if normalise is True:

        # The following removes any degeneracy and normalises the vectors onto the unit circle.

        accepted_family_of_vectors = [list(degenerate_family_of_vectors[0]/ np.linalg.norm(degenerate_family_of_vectors[0]))]

        for vector in degenerate_family_of_vectors:

            normalised_vector = list(vector/np.linalg.norm(vector))

            if normalised_vector not in accepted_family_of_vectors:

                accepted_family_of_vectors.append(normalised_vector)

    else:

        # The following removes any degeneracy without normalising.

        accepted_family_of_vectors = [degenerate_family_of_vectors[0]]

        for vector in degenerate_family_of_vectors:

            if vector not in accepted_family_of_vectors:

                accepted_family_of_vectors.append(vector)

    return accepted_family_of_vectors


def create_powder_euler_angles(num_points):

    import numpy as np

    rotation_array = np.random.normal(27.5, 20.0, num_points).tolist()

    all_euler_angles = []

    for psi in rotation_array:

        single_euler_angles = [0.0, psi, 45.0]

        all_euler_angles.append(single_euler_angles)

    return all_euler_angles


def create_alpha_fibre_euler_angles(num_points):

    import numpy as np

    rotation_array = np.random.uniform(0.0, 55.0, num_points).tolist()

    all_euler_angles = []

    for psi in rotation_array:

        single_euler_angles = [0.0, psi, 45.0]

        all_euler_angles.append(single_euler_angles)

    return all_euler_angles


def rotate_vector_family_thru_euler_angles(vector_family, euler_angles, gaussian_scatter_width, num_points):

    import numpy as np
    from numpy import sin, cos

    s_1 = np.deg2rad(np.random.normal(0.0, gaussian_scatter_width, num_points))
    s_psi = np.deg2rad(np.random.normal(0.0, gaussian_scatter_width, num_points))
    s_2 = np.deg2rad(np.random.normal(0.0, gaussian_scatter_width, num_points))

    # The following rotates the crystallite such that it is oriented along the fibre axis.
    # This formula is given in "Tables for Texture Analysis of Cubic Crystals" pg 9, by Hansen, Pospiech & Lucke

    rotated_vectors = []

    for j, euler_set in enumerate(euler_angles):

        phi_1, psi, phi_2 = np.deg2rad(euler_set[0]), np.deg2rad(euler_set[1]), np.deg2rad(euler_set[2])

        scattered_phi_1 = phi_1 + s_1[j]
        scattered_psi = psi + s_psi[j]
        scattered_phi_2 = phi_2 + s_2[j]

        euler_angle_rotation_matrix = np.array([[(cos(scattered_phi_1) * cos(scattered_phi_2)) - (sin(scattered_phi_1) * sin(scattered_phi_2) * cos(scattered_psi)),
                                                   (sin(scattered_phi_1) * cos(scattered_phi_2)) + (cos(scattered_phi_1) * sin(scattered_phi_2) * cos(scattered_psi)),
                                                   sin(scattered_phi_2) * sin(scattered_psi)],
                                                  [(- cos(scattered_phi_1) * sin(scattered_phi_2)) - (sin(scattered_phi_1) * cos(scattered_phi_2) * cos(scattered_psi)),
                                                   (-sin(scattered_phi_1) * sin(scattered_phi_2)) + (cos(scattered_phi_1) * cos(scattered_phi_2) * cos(scattered_psi)),
                                                   cos(scattered_phi_2) * sin(scattered_psi)],
                                                  [sin(scattered_phi_1) * sin(scattered_psi),
                                                   - cos(scattered_phi_1) * sin(scattered_psi),
                                                   cos(scattered_psi)]])

        transpose_rotation_matrix = np.transpose(euler_angle_rotation_matrix)

        for i, vector in enumerate(vector_family):

            vector = np.dot(transpose_rotation_matrix, vector)

            rotated_vectors.append(vector)

    return rotated_vectors


def find_diffraction_vectors(k_inc, wavelength, a_lattice, hkl_plane, num_scattering_vectors):

    import numpy as np

    h, k, l = hkl_plane

    d = a_lattice / (np.sqrt((h ** 2) + (k ** 2) + (l ** 2)))

    bragg_angle = np.arcsin((wavelength / (2.0 * d)))

    def rotate_vector_about_another_vector(v_1, axis_vector, rotation_angle):

        v_1 = v_1 / np.linalg.norm(k_inc)
        axis_vector = axis_vector / np.linalg.norm(axis_vector)
        parallel_component = (np.dot(v_1, axis_vector)) * axis_vector
        perpendicular_component = v_1 - parallel_component
        w = np.cross(axis_vector, perpendicular_component)
        w = w / np.linalg.norm(w)
        x1 = np.cos(rotation_angle)
        x2 = np.sin(rotation_angle)

        rotation_component = (x1 * perpendicular_component) + (x2 * w)

        normalised_rotated_vector = rotation_component + parallel_component

        return normalised_rotated_vector

    normalised_k_dif = rotate_vector_about_another_vector(k_inc, np.array([0.0, 1.0, 0.0]), 2 * bragg_angle)

    powder_angles = np.linspace(0.0, 2 * np.pi, num_scattering_vectors, endpoint=False)

    all_k_dif = []
    all_G = []

    for phi in powder_angles:

        normalised_rotated_k_dif = rotate_vector_about_another_vector(normalised_k_dif, k_inc, phi)

        G = k_inc - normalised_rotated_k_dif
        G = G/np.linalg.norm(G)

        all_k_dif.append(normalised_rotated_k_dif)
        all_G.append(G)

    return all_k_dif, all_G


def unnormalise_vector_quiver(vector_quiver, hkl_plane, a_lattice, wavelength):

    h, k, l = hkl_plane

    d = a_lattice / (np.sqrt((h ** 2) + (k ** 2) + (l ** 2)))

    bragg_angle = np.arcsin((wavelength / (2.0 * d)))

    length = np.sin((0.5 * np.pi) - bragg_angle)/(d * np.sin(2 * bragg_angle))

    vector_quiver = np.array(vector_quiver)

    unnormalised_vector_quiver = length * vector_quiver

    return unnormalised_vector_quiver, length


def map_vectors_using_Johns_method(fibre_vectors, diffraction_vectors, num_bins, plot_stereographs):

    import matplotlib.pyplot as plt

    def project_vectors_johns_method(vectors):

        all_theta_pole = []
        all_phi_pole = []
        all_x_num = []
        all_y_num = []

        for vector in vectors:


            mag = np.linalg.norm(vector)
            theta_pole = np.arccos(vector[2]/mag)
            phi_pole = np.arctan2(vector[1]/mag, vector[0]/mag)
            x_num = np.tan(theta_pole * 0.5) * np.cos(phi_pole)
            y_num = np.tan(theta_pole * 0.5) * np.sin(phi_pole)

            all_theta_pole.append(theta_pole)
            all_phi_pole.append(phi_pole)
            all_x_num.append(x_num)
            all_y_num.append(y_num)

        all_theta_pole = np.array(all_theta_pole)
        all_phi_pole = np.array(all_phi_pole)

        return all_theta_pole, all_phi_pole, all_x_num, all_y_num

    all_theta_pole, all_phi_pole, all_x_num, all_y_num = project_vectors_johns_method(fibre_vectors[0])

    bin_population, x_bins, y_bins = np.histogram2d(all_x_num, all_y_num, bins=num_bins, range=((-1.0, 1.0), (-1.0, 1.0)))
    plt.imshow(bin_population, cmap='jet', interpolation='None', origin='Low')
    plt.show()

    return


def map_vectors_stereographically_onto_plane(fibre_vectors, diffraction_vectors, num_bins, plot_stereographs):

    import numpy as np
    import matplotlib.pyplot as plt

    def project_vectors_stereographically(vector_quiver):

        all_X = []
        all_Y = []

        for vector in list(vector_quiver):
            x, y, z = vector

            X = 1.0 * (x / (1.0 - z))
            Y = 1.0 * (y / (1.0 - z))

            all_X.append(X)

            all_Y.append(Y)

        all_X = np.nan_to_num(all_X)
        all_Y = np.nan_to_num(all_Y)
        """
        for i, X in enumerate(all_X):

            if (X**2) + (all_Y[i]**2) > 1.0:

                all_X[i] = np.nan
                all_Y[i] = np.nan
        """

        bin_population, x_bins, y_bins = np.histogram2d(all_X, all_Y, bins = num_bins, range=((-1.0, 1.0), (-1.0, 1.0)))

        return bin_population, x_bins, y_bins

    fibre_vector_bins, fibre_x, fibre_y = project_vectors_stereographically(fibre_vectors)

    diffraction_vector_bins, diffraction_x, diffraction_y = project_vectors_stereographically(diffraction_vectors)

    plt.imshow(fibre_vector_bins, cmap='jet', interpolation='None', origin='Low')#, extent=[fibre_x[0], fibre_x[1], fibre_y[0], fibre_y[1]])
    if plot_stereographs is True:
        plt.show()
    plt.imshow(diffraction_vector_bins, cmap='jet', interpolation='None', origin='Low')#, extent=[diffraction_x[0], diffraction_x[1], diffraction_y[0], diffraction_y[1]])
    if plot_stereographs is True:
        plt.show()
    plt.close()

    return


def plot_hedgehog(vector_list):

    import numpy as np
    from mayavi import mlab
    from mayavi.mlab import quiver3d

    vector_array = np.array(vector_list)

    x = [0.0] * len(vector_list)
    y = [0.0] * len(vector_list)
    z = [0.0] * len(vector_list)
    u, v, w = vector_array[:, 0], vector_array[:, 1], vector_array[:, 2]

    quiver3d(x, y, z, u, v, w, scale_factor=1, color=(1, 1, 1), mode='2ddash')
    mlab.axes()
    mlab.show()

    return

initial_vectors = find_family_of_vectors(hkl_plane, normalise_G_vectors_to_unit_sphere)

num_points = len(initial_vectors) * num_points_per_vector

print "Total number of points will be " + str(num_points)

all_euler_sets = create_alpha_fibre_euler_angles(num_points)

#all_euler_sets = [[0.0, 0.0, 0.0]] * num_points

normalised_diffraction_vector_quiver, normalised_G_vector_quiver = find_diffraction_vectors(k_inc, wavelength, a_lattice, hkl_plane, num_scattering_vectors)

normalised_fibre_vector_quiver = rotate_vector_family_thru_euler_angles(initial_vectors, all_euler_sets, gaussian_scatter_width, num_points)

#diffraction_vector_quiver = unnormalise_vector_quiver(normalised_diffraction_vector_quiver, hkl_plane, a_lattice, wavelength)

fibre_vector_quiver = unnormalise_vector_quiver(normalised_fibre_vector_quiver, hkl_plane, a_lattice, wavelength)

map_vectors_using_Johns_method(fibre_vector_quiver, normalised_diffraction_vector_quiver, num_bins, plot_stereographs)

map_vectors_stereographically_onto_plane(normalised_fibre_vector_quiver, normalised_G_vector_quiver, num_bins, plot_stereographs)

if plot_quivers is True:

    plot_hedgehog(initial_vectors)

    plot_hedgehog(normalised_fibre_vector_quiver)

    plot_hedgehog(normalised_G_vector_quiver)
