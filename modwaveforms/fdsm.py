import numpy as np
from bilby.gw.source import lal_binary_black_hole
from . import geomoptics
from . import waveoptics

def one_image_BBH(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, Delta_phase, **kwargs):
    """
    Generate the frequency domain waveform for lensed image of a binary black hole
    
    Delta_phase is the phase difference w.r.t. the unlensed waveform
    """
    # Actually generate the waveform by calling the generator
    wf = lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )
    F_one_images = geomoptics.one_image_BBH(frequency_array, Delta_phase)

    wf["plus"] = wf["plus"]*F_one_images
    wf["cross"] = wf["cross"]*F_one_images
    return wf

def two_images_BBH(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, mu_rel,Delta_t,Delta_phase, **kwargs):
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    mu_rel is the relative amplitude of the second image
    Delta_t is the time delay between the two images
    Delta_phase is the phase difference between the two images
    """
    # Actually generate the waveform by calling the generator
    wf = lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )
    F_two_images = geomoptics.two_images_BBH(frequency_array, mu_rel, Delta_t, Delta_phase)

    wf["plus"] = wf["plus"]*F_two_images
    wf["cross"] = wf["cross"]*F_two_images
    return wf

def fold_caustic(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, Delta_t, positive_phase, **kwargs):
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    Delta_t is the time delay between the two images
    positive_phase is the phase of the positive parity image
    """
    # Actually generate the waveform by calling the generator
    wf = lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )
    F_fold_two_images = geomoptics.F_fold_caustic(frequency_array, Delta_t, positive_phase)

    wf["plus"] = wf["plus"]*F_fold_two_images
    wf["cross"] = wf["cross"]*F_fold_two_images
    return wf

def cusp_caustic(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, Delta_t_10, Delta_t_20, mu_rel, positive_phase, **kwargs):
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    Delta_t is the time delay between the two images
    positive_phase is the phase of the positive parity image
    """
    # Actually generate the waveform by calling the generator
    wf = lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )
    F_cusp_three_images = geomoptics.F_cusp_three_images(frequency_array, Delta_t_10, Delta_t_20, mu_rel, positive_phase)
 
    wf["plus"] = wf["plus"]*F_cusp_three_images
    wf["cross"] = wf["cross"]*F_cusp_three_images
    return wf

def pointlens(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
    phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, ML, y, **kwargs
    ):
    """
    Generate the frequency domain waveform for a binary black hole merger lensed by a point mass lens

    Parameters
    ----------
    frequency_array : array
        Array of frequencies in Hz
    mass_1 : float
        Mass of the primary black hole in solar masses
    mass_2 : float
        Mass of the secondary black hole in solar masses
    luminosity_distance : float
        Luminosity distance to the source in Mpc
    a_1 : float
        Dimensionless spin magnitude of the primary black hole
    tilt_1 : float
        Tilt angle of the primary black hole spin in radians
    phi_12 : float
        Azimuthal angle between the two black hole spins in radians
    a_2 : float
        Dimensionless spin magnitude of the secondary black hole
    tilt_2 : float
        Tilt angle of the secondary black hole spin in radians
    phi_jl : float
        Azimuthal angle between the total angular momentum and the line of sight in radians
    theta_jn : float
        Angle between the total angular momentum and the line of sight in radians
    phase : float
        Reference phase of the binary in radians
    ML : float
        Redshifted mass of the lens in solar masses
    y : float
        Dimensionless impact parameter of the lens
    **kwargs : dict
        Additional keyword arguments to pass to the waveform generator

    Returns
    -------
    wf : dict
        Dictionary containing the plus and cross polarizations of the lensed waveform
    """
    # Actually generate the waveform by calling the generator
    wf = lal_binary_black_hole(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs
    )
    F = waveoptics.F_pointlens(frequency_array, ML, y)

    wf["plus"] = wf["plus"]*F
    wf["cross"] = wf["cross"]*F
    return wf
