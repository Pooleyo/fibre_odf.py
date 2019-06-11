import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

wavelength = [1.398E-10, 1.367E-10]  # He-alpha emission has low= 1.398 A; high = 1.367 A
a_lattice = 3.3E-10

hkl_plane = [4, 1, 1]
gaussian_width = 2.0/2.355  # Only applies to fibre texture.
fibre_counter = 1000
powder_counter = 1000000
num_bins = 961
num_circles = 10

polefig_fibre_texture = True
polefig_powder_texture = False
fibre_profile = False

alpha = 40.0
twotheta = [2 * np.rad2deg(np.arcsin(0.5 * wavelength[0] * np.linalg.norm(hkl_plane) / a_lattice)),
            2 * np.rad2deg(np.arcsin(0.5 * wavelength[1] * np.linalg.norm(hkl_plane) / a_lattice))]

diffraction_width_multipliers = [0.97, 1.0]

num_phi_bins = 361
phi_limits = [-113.62725450901803, 114.34869739478961]  # In degrees. This limiter is applied AFTER the bins have been setup.

show_polefig = False
show_diffractionfig = False
show_intersectionfig = False

def MakeVectorFamily(hkl_plane):

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

    # The following removes any degeneracy without normalising.

    vector_family = [degenerate_family_of_vectors[0]]

    for vector in degenerate_family_of_vectors:

        if vector not in vector_family:

            vector_family.append(vector)

    return vector_family


def RotMatBunge(G, phi1, PHI, phi2):

    import numpy as np

    phi1 = np.deg2rad(phi1)
    PHI = np.deg2rad(PHI)
    phi2 = np.deg2rad(phi2)

    cp1 = np.cos(phi1)
    sp1 = np.sin(phi1)
    ct = np.cos(PHI)
    st = np.sin(PHI)
    cp2 = np.cos(phi2)
    sp2 = np.sin(phi2)

    RotationMatrix = np.empty((3,3))

    RotationMatrix[0, 0] = cp1*cp2 - sp1*sp2*ct
    RotationMatrix[1, 0] = - cp1*sp2 - sp1*cp2*ct
    RotationMatrix[2, 0] = sp1 * st
    RotationMatrix[0, 1] = sp1*cp2 + cp1*sp2*ct
    RotationMatrix[1, 1] = - sp1*sp2 + cp1*cp2*ct
    RotationMatrix[2, 1] = - cp1*st
    RotationMatrix[0, 2] = sp2*st
    RotationMatrix[1, 2] = cp2*st
    RotationMatrix[2, 2] = ct

    RotationMatrix = np.transpose(RotationMatrix)

    G_new = np.dot(RotationMatrix, G)

    return G_new


def RotateXtal2(G, theta, phi, angle):

    import numpy as np

    theta = np.deg2rad(theta)
    phi = np.deg2rad(phi)
    angle = np.deg2rad(angle)

    u = np.empty(3)

    u[0] = np.sin(theta)*np.cos(phi)
    u[1] = np.sin(theta)*np.sin(phi)
    u[2] = np.cos(theta)

    RotationMatrix = np.empty((3, 3))

    RotationMatrix[0][0] = np.cos(angle) + (u[0] * u[0] * (1 - np.cos(angle)))
    RotationMatrix[1][0] = (u[1] * u[0] * (1 - np.cos(angle))) + (u[2] * np.sin(angle))
    RotationMatrix[2][0] = (u[2] * u[0] * (1 - np.cos(angle))) + (u[1] * np.sin(angle))
    RotationMatrix[0][1] = (u[0] * u[1] * (1 - np.cos(angle))) - (u[2] * np.sin(angle))
    RotationMatrix[1][1] = np.cos(angle) + (u[1] * u[1] * (1 - np.cos(angle)))
    RotationMatrix[2][1] = (u[2] * u[1] * (1 - np.cos(angle))) + (u[0] * np.sin(angle))
    RotationMatrix[0][2] = (u[0] * u[2] * (1 - np.cos(angle))) + (u[1] * np.sin(angle))
    RotationMatrix[1][2] = (u[1] * u[2] * (1 - np.cos(angle))) - (u[0] * np.sin(angle))
    RotationMatrix[2][2] = np.cos(angle) + (u[2] * u[2] * (1 - np.cos(angle)))

    G_new = np.dot(RotationMatrix, G)

    return G_new


def RotMatFibonacci(samples):

    import numpy as np

    offset = 2.0/samples

    increment = np.pi * (3.0 * np.sqrt(5.0))

    all_points = np.empty((samples,3))

    for s in range(samples):

        y = ((s * offset) - 1.0) + (offset/2.0)

        r = np.sqrt(1 - (y * y))

        phi = s % samples * increment

        x = np.cos(phi) * r
        z = np.sin(phi) * r

        all_points[s] = np.array([x, y, z])

    return all_points


def ReflectPolefig(polefig):

    master_polefig = polefig + np.flip(polefig, 0) + np.flip(polefig, 1) + np.flip(polefig, (0, 1))
    plt.imsave(output_folder + '/0_polefig.tif', master_polefig, origin='lower', cmap='viridis')
    plt.imshow(master_polefig, origin="lower", cmap='viridis')
    if show_polefig is True:
        plt.show()

    plt.close()

    return master_polefig


def PlotHedgehog(vectors):

    import numpy as np
    from mayavi import mlab
    from mayavi.mlab import quiver3d

    shape = np.shape(vectors)
    num_vectors = shape[0]

    x = [0.0] * num_vectors
    y = [0.0] * num_vectors
    z = [0.0] * num_vectors
    u, v, w = vectors[:, 0], vectors[:, 1], vectors[:, 2]

    quiver3d(x, y, z, u, v, w, scale_factor=1, color=(1, 1, 1), mode='2ddash')
    mlab.axes()
    mlab.show()

    return


def SingleFibreTexture(h, k, l, counter, gaussian_width, polefig, polefig_bins):

    import numpy as np

    all_G = np.zeros((counter, 3))

    for i in range(counter):

        a = np.array([1.0, 0.0, 0.0])
        b = np.array([0.0, 1.0, 0.0])
        c = np.array([0.0, 0.0, 1.0])

        G = h*a + k*b + l*c

        G = RotMatBunge(G, 0.0, np.random.uniform(0.0, 55.0), 45.0)

        G_mag = np.linalg.norm(G)

        if G[2] > 0.0:  # TODO why? What happens if we remove the if statement? ANS: This selects for hemispheres in the
            # TODO hedgehog. One hemisphere will project inside the unit circle, the other will project outside.

            G = RotateXtal2(G, 0.0, 0.0, np.random.normal(scale=gaussian_width))
            G = RotateXtal2(G, 90.0, 0.0, np.random.normal(scale=gaussian_width))
            G = RotateXtal2(G, 90.0, 90.0, np.random.normal(scale=gaussian_width))

            theta_pole = np.arccos(G[2]/G_mag)
            phi_pole = np.arctan2(G[1]/G_mag, G[0]/G_mag)

            x_num = (np.tan(theta_pole*0.5))*np.cos(phi_pole)
            y_num = (np.tan(theta_pole*0.5))*np.sin(phi_pole)

            x_bin = np.argmin(abs(polefig_bins - x_num))
            y_bin = np.argmin(abs(polefig_bins - y_num))

            polefig[x_bin, y_bin] += 1

        all_G[i] = G/G_mag

    return all_G


def PowderTexture(counter, polefig, polefig_bins):

    all_G = RotMatFibonacci(counter)

    for G in all_G:

        if G[2] > 0.0:  # TODO why? What happens if we remove the if statement? ANS: This selects for hemispheres in the
            # TODO hedgehog. One hemisphere will project inside the unit circle, the other will project outside.

            G_mag = np.linalg.norm(G)

            theta_pole = np.arccos(G[2] / G_mag)
            phi_pole = np.arctan2(G[1] / G_mag, G[0] / G_mag)

            x_num = (np.tan(theta_pole * 0.5)) * np.cos(phi_pole)  # TODO look into DimOffset
            y_num = (np.tan(theta_pole * 0.5)) * np.sin(phi_pole)  # TODO look into DimOffset

            x_bin = np.argmin(abs(polefig_bins - x_num))
            y_bin = np.argmin(abs(polefig_bins - y_num))

            polefig[x_bin, y_bin] += 1

    return


def PoleFigureXY(alpha, twotheta, phi):

    from numpy import sin, cos, tan, deg2rad, empty, arccos, arctan2
    from numpy.linalg import norm

    alpha = deg2rad(alpha)
    twotheta = deg2rad(twotheta)
    phi = deg2rad(phi)

    k_inc = empty(3)
    k_inc[0] = sin(alpha)
    k_inc[1] = 0.0
    k_inc[2] = - cos(alpha)

    k_dif_prime = empty(3)
    k_dif_prime[0] = sin(twotheta)*cos(phi)
    k_dif_prime[1] = sin(twotheta)*sin(phi)
    k_dif_prime[2] = - cos(twotheta)

    k_dif = empty(3)
    k_dif[0] = k_dif_prime[0]*cos(alpha) - k_dif_prime[2]*sin(alpha)
    k_dif[1] = k_dif_prime[1]
    k_dif[2] = k_dif_prime[0]*sin(alpha) + k_dif_prime[2]*cos(alpha)

    G = k_dif - k_inc

    G_mag = norm(G)
    G = G/G_mag

    theta_prime = arccos(G[2])
    phi_prime = arctan2(G[1], G[0])

    x_polefig = tan(theta_prime * 0.5)*cos(phi_prime)
    y_polefig = tan(theta_prime * 0.5)*sin(phi_prime)

    return x_polefig, y_polefig


def FindCircleFeatures(alpha, twotheta):

    x0, y0 = PoleFigureXY(alpha, twotheta, 0.0)
    x1, y1 = PoleFigureXY(alpha, twotheta, 180.0)

    circle_origin = [(y0 + y1) * 0.5, (x0 + x1) * 0.5]  # y value goes first since numpy has format [row, col]

    circle_radius = 0.5 * np.linalg.norm(np.array([y1 - y0, x1 - x0]))

    return circle_origin, circle_radius


def CalcPhiDiffractionRing(circle, origin, num_bins):

    import numpy as np

    x_ind, y_ind = np.nonzero(circle)

    xx, yy = np.mgrid[-1.2:1.2:(num_bins * 1j), -1.2:1.2:(num_bins * 1j)]

    xx = xx[x_ind, y_ind]
    yy = yy[x_ind, y_ind]

    centred_xx = xx - origin[0]
    centred_yy = yy - origin[1]

    phi = np.rad2deg(np.arctan2(centred_xx, centred_yy))

    for p in phi:

        if p < 0.0:

            p = 360.0 + p

    return phi


def BinDiffractionRing(texture_profile, phi, num_bins):

    import numpy as np
    import matplotlib.pyplot as plt

    bins = np.linspace(-180.0, 180.0, num_bins)

    indices = np.digitize(phi, bins)

    summed_texture_profile_across_phi = np.zeros(np.size(bins))
    population = np.zeros(np.size(bins))

    for i, ind in enumerate(indices):

        summed_texture_profile_across_phi[ind] += texture_profile[i]

        population[ind] += 1

    summed_texture_profile_across_phi[summed_texture_profile_across_phi < 1] = np.nan

    integrated_texture_profile_across_phi = summed_texture_profile_across_phi / population

    integrated_texture_profile_across_phi = np.nan_to_num(integrated_texture_profile_across_phi)

    np.savetxt(output_folder + '/texture_profile_vs_phi.csv', zip(bins, integrated_texture_profile_across_phi), delimiter=',')

    plt.plot(bins, integrated_texture_profile_across_phi)
    plt.xlabel('phi (degrees)')
    plt.ylabel('intensity (arbitrary)')
    plt.savefig(output_folder + '/3_texture_profile.png')
    plt.close()

    return integrated_texture_profile_across_phi, bins


def MakeDiffractionAnnulusArray(array_shape, origin_high_energy, origin_low_energy, radius_high_energy,
                                radius_low_energy, num_circles, twotheta):

    import numpy as np

    def Circle(origin, lower_width, upper_width, plot_name):

        xx, yy = np.mgrid[-1.2:1.2:(num_bins * 1j), -1.2:1.2:(num_bins * 1j)]

        r_sqr = (xx - origin[0]) ** 2 + (yy - origin[1]) ** 2

        circle = np.logical_and(r_sqr > lower_width**2, r_sqr < upper_width**2)

        plt.imsave(output_folder + '/' + plot_name, circle, origin='lower', cmap='viridis')
        if show_diffractionfig is True:
            plt.imshow(circle, origin="lower", cmap='viridis')
            plt.show()
        plt.close()

        return circle

    circle_high_energy = Circle(origin_high_energy, radius_high_energy, 1.02 * radius_high_energy,
                                "1_circle_high_energy")

    circle_low_energy = Circle(origin_low_energy, 0.98 * radius_low_energy, radius_low_energy, "2_circle_low_energy")

    edge_circles = circle_low_energy + circle_high_energy

    plt.imsave(output_folder + '/edge_circles', edge_circles, origin='lower', cmap='viridis')
    plt.close()

    all_circles = edge_circles.copy()

    twotheta_scan_values = np.linspace(twotheta[0], twotheta[1], num_circles, endpoint=False)

    dr = 0.02

    for i, twotheta in enumerate(twotheta_scan_values):

        o, r = FindCircleFeatures(alpha, twotheta)

        circ = Circle(o, r, (1.0 + dr) * r, "circ_" + str(i))

        all_circles += circ

    plt.imsave(output_folder + '/all_circles', all_circles, origin='lower', cmap='viridis')
    if show_diffractionfig is True:
        plt.imshow(all_circles, origin="lower", cmap='viridis')
        plt.show()
    plt.close()

    return all_circles


def FindArrayIntersection(polefig_arr, dif_arr):

    rows, cols = np.nonzero(dif_arr)

    texture_profile = []

    ring_value = 1.0 * np.max(polefig_arr)

    intersection_arr = np.copy(polefig_arr)

    for row, col in zip(rows, cols):

        texture_profile.append(polefig_arr[row, col])

        intersection_arr[row, col] = ring_value

    plt.imsave(output_folder + '/4_intersection_fig.tif', intersection_arr, origin='lower', cmap='viridis')
    if show_intersectionfig is True:
        plt.imshow(intersection_arr, origin="lower", cmap='viridis')
        plt.show()
    plt.close()

    texture_profile = np.array(texture_profile)

    return texture_profile


def CalcIntegratedIntensity(profile, phi_bins, phi_limits):


    start_ind = np.argmin(np.abs(phi_bins - phi_limits[0]))
    finish_ind = np.argmin(np.abs(phi_bins - phi_limits[1]))

    restricted_profile = profile[start_ind:finish_ind]

    integrated_intensity = np.sum(restricted_profile)

    return integrated_intensity


def CalcCorrectionFactor(fibre_intensity, powder_intensity, fibre_counter, powder_counter):

    # This accounts for the mismatch in the number of G vectors present in the powder relative to the fibre. Note that
    # half of the fibre vectors are projected out of the unit circle projection, hence the factor of 2.

    G_mismatch_factor = 2 * (powder_counter / fibre_counter)

    corrected_fibre_intensity = G_mismatch_factor * fibre_intensity

    # The intensity correction factor is the multiplier which will correct the intensity as if it were a powder.

    intensity_correction_factor = powder_intensity / corrected_fibre_intensity

    np.savetxt(output_folder + '/correction_factor.dat', np.array([intensity_correction_factor]))

    print "The correction factor for the "+str(hkl_plane[0])+str(hkl_plane[1])+str(hkl_plane[2])+' reflection = ' \
          + str(intensity_correction_factor)

    return intensity_correction_factor


def runPolefigFibre(hkl_plane, counter, alpha, twotheta, diffraction_width_mutipliers, num_phi_bins):

    # The following runs the SingleFibreTexture function for every accepted plane.

    polefig = np.zeros((num_bins, num_bins))
    polefig_bins = np.linspace(-1.2, 1.2, num_bins, endpoint=True)

    vector_family = MakeVectorFamily(hkl_plane)

    for plane in vector_family:

        h, k, l = plane

        SingleFibreTexture(h, k, l, counter, gaussian_width, polefig, polefig_bins)

    master_polefig = ReflectPolefig(polefig)

    origin_low_energy, radius_low_energy = FindCircleFeatures(alpha, twotheta[0])

    origin_high_energy, radius_high_energy = FindCircleFeatures(alpha, twotheta[1])

    dif_array = MakeDiffractionAnnulusArray(polefig.shape, origin_high_energy, origin_low_energy, radius_high_energy,
                                            radius_low_energy, num_circles, twotheta)

    texture_profile = FindArrayIntersection(master_polefig, dif_array)

    phi = CalcPhiDiffractionRing(dif_array, origin_low_energy, num_bins)

    integrated_texture_profile, bins = BinDiffractionRing(texture_profile, phi, num_phi_bins)

    fibre_intensity = CalcIntegratedIntensity(integrated_texture_profile, bins, phi_limits)

    return fibre_intensity


def runPolefigPowder(counter, alpha, twotheta, diffraction_width_multipliers):

    polefig = np.zeros((num_bins, num_bins))
    polefig_bins = np.linspace(-1.2, 1.2, num_bins, endpoint=True)

    PowderTexture(counter, polefig, polefig_bins)

    master_polefig = ReflectPolefig(polefig)
    # TODO clean up this code
    x0, y0 = PoleFigureXY(alpha, twotheta, 0.0)
    x1, y1 = PoleFigureXY(alpha, twotheta, 180.0)

    circle_origin = [(y0 + y1) * 0.5, (x0 + x1) * 0.5]  # y value goes first since numpy has format [row, col]

    circle_radius = 0.5 * np.linalg.norm(np.array([y1-y0, x1-x0]))

    inner_radius = diffraction_width_multipliers[0] * circle_radius
    outer_radius = diffraction_width_multipliers[1] * circle_radius

    dif_array = MakeDiffractionAnnulusArray(polefig.shape, circle_origin, inner_radius, outer_radius)

    texture_profile = FindArrayIntersection(master_polefig, dif_array)

    phi = CalcPhiDiffractionRing(dif_array, circle_origin)

    integrated_profile, bins = BinDiffractionRing(texture_profile, phi, num_phi_bins)

    powder_intensity = CalcIntegratedIntensity(integrated_profile, bins, [-180.0, 180.0])

    return powder_intensity


def runProfileFibre(hkl_plane, counter, alpha, twotheta, diffraction_width_mutipliers, num_phi_bins):

    # The following runs the SingleFibreTexture function for every accepted plane.

    polefig = np.zeros((num_bins, num_bins))
    polefig_bins = np.linspace(-1.2, 1.2, num_bins, endpoint=True)

    vector_family = MakeVectorFamily(hkl_plane)

    for plane in vector_family:

        h, k, l = plane

        all_G = SingleFibreTexture(h, k, l, counter, gaussian_width, polefig, polefig_bins)

    PlotHedgehog(all_G)

    """
    master_polefig = ReflectPolefig(polefig)

    origin_low_energy, radius_low_energy = FindCircleFeatures(alpha, twotheta[0])

    origin_high_energy, radius_high_energy = FindCircleFeatures(alpha, twotheta[1])


    dif_array = MakeDiffractionRingArray(polefig.shape, origin_high_energy, origin_low_energy, radius_high_energy,
                                         radius_high_energy)
    exit()
    texture_profile = FindArrayIntersection(master_polefig, dif_array)

    phi = CalcPhiDiffractionRing(dif_array, circle_origin)

    integrated_texture_profile, bins = BinDiffractionRing(texture_profile, phi, num_phi_bins)

    fibre_intensity = CalcIntegratedIntensity(integrated_texture_profile, bins, phi_limits)
    """
    return


if polefig_fibre_texture is True:

    output_folder = 'output_polefig_fibre_' + str(hkl_plane[0]) + str(hkl_plane[1]) + str(hkl_plane[2])

    if os.path.exists(output_folder):

        shutil.rmtree(output_folder)

    os.mkdir(output_folder)

    fibre_intensity = runPolefigFibre(hkl_plane, fibre_counter, alpha, twotheta, diffraction_width_multipliers, num_phi_bins)

if polefig_powder_texture is True:

    output_folder = 'output_polefig_powder_' + str(hkl_plane[0]) + str(hkl_plane[1]) + str(hkl_plane[2])

    if os.path.exists(output_folder):

        shutil.rmtree(output_folder)

    os.mkdir(output_folder)

    powder_intensity = runPolefigPowder(powder_counter, alpha, twotheta, diffraction_width_multipliers)

if polefig_powder_texture and polefig_fibre_texture is True:

    correction_factor = CalcCorrectionFactor(fibre_intensity, powder_intensity, fibre_counter, powder_counter)

if fibre_profile is True:

    output_folder = 'output_profile_fibre_' + str(hkl_plane[0]) + str(hkl_plane[1]) + str(hkl_plane[2])

    if os.path.exists(output_folder):

        shutil.rmtree(output_folder)

    os.mkdir(output_folder)

    fibre_intensity = runProfileFibre(hkl_plane, fibre_counter, alpha, twotheta, diffraction_width_multipliers, num_phi_bins)
