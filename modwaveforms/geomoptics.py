import numpy as np

def one_image_BBH(frequency_array, Delta_phase):
    """
    Generate the frequency domain waveform for lensed image of a binary black hole
    
    Delta_phase is the phase difference w.r.t. the unlensed waveform
    """
    F_one_images = np.exp(-1j*np.sign(frequency_array)*Delta_phase)

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    F_one_images = np.conjugate(F_one_images)

    return F_one_images


def two_images_BBH(frequency_array, mu_rel, Delta_t, Delta_phase):
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    mu_rel is the relative amplitude of the second image
    Delta_t is the time delay between the two images
    Delta_phase is the phase difference between the two images
    """
    F_two_images = (1.+np.sqrt(mu_rel)*np.exp(1j*(2*np.pi*frequency_array*Delta_t-np.sign(frequency_array)*Delta_phase)))

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    F_two_images = np.conjugate(F_two_images)

    return F_two_images

def fold_caustic(frequency_array, Delta_t, positive_phase):
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    Delta_t is the time delay between the two images
    positive_phase is the phase of the positive parity image
    """
    F_fold_two_images = positive_phase*(1. + np.exp(1j*(2*np.pi*frequency_array*Delta_t-positive_phase*np.sign(frequency_array)*np.pi/2)))

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    F_fold_two_images = np.conjugate(F_fold_two_images)

    return F_fold_two_images

def cusp_caustic(frequency_array, Delta_t_10, Delta_t_20, mu_rel, positive_phase):
    """
    Generate the frequency domain waveform for a binary black hole system with two images

    Delta_t is the time delay between the two images
    positive_phase is the phase of the positive parity image
    """
    F_cusp_three_images = positive_phase*(1. + np.sqrt(abs(mu_rel))*np.exp(1j*(2*np.pi*frequency_array*Delta_t_10-positive_phase*np.sign(frequency_array)*np.pi/2))
                                          + np.sqrt(1-abs(mu_rel))*np.exp(1j*(2*np.pi*frequency_array*Delta_t_20-positive_phase*np.sign(frequency_array)*np.pi/2)))

    # complex conjugate because of different sign convention in the Fourier transform of 
    # LALSuite (engineering) compared to the one used in wave optics (physics)
    F_cusp_three_images = np.conjugate(F_cusp_three_images)
    
    return F_cusp_three_images
