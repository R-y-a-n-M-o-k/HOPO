import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

def normalised_profile(pgopher_input, band_center, output_file_name):
    """This function takes the pgopher input (broadened spectrum) and returns a 
       normalised profile with x-axis shifted to 0 cm-1

    Input: pgopher_input: .out file, first column is wavenumbers, second column is intensity (string)
           band_center: center of the band in cm-1 (float)
           output_file_name: name of the output file to save the normalised profile (string)
    Output: .prof file of two columns containing centered wavenumbers, normalised distribution values  
    """

    wavenumber_axis, intensity = np.loadtxt(pgopher_input, unpack=True)
    # Shift the wavenumber axis to 0 cm-1
    wavenumber_axis -= band_center

    # Normalize intensity values
    area = simpson(intensity, wavenumber_axis)
    normalised_intensity = intensity / area

    # Save the normalised profile to a .prof file
    distribution = np.column_stack((wavenumber_axis, normalised_intensity))
    np.savetxt(output_file_name, distribution, fmt="%.6f %.5e")

    return wavenumber_axis, normalised_intensity

def plot_profiles(profile_name):
    """This function plots the normalised profile from a .prof file

    Input: profile_name: name of the .prof file containing the normalised profile (string)
    Output: Plot of the normalised profile
    """

    wavenumber_axis, normalised_intensity = np.loadtxt(profile_name, unpack=True)

    plt.figure(figsize=(8, 6))
    plt.plot(wavenumber_axis, normalised_intensity, label='Normalised Profile', color='blue')
    plt.yscale('log')
    plt.title('Normalised Profile')
    plt.xlabel('Wavenumber (cm$^{-1}$)')
    plt.ylabel('Normalised Intensity')
    plt.legend()
    #plt.grid()
    plt.show()


test_wave, test_norm = normalised_profile('test.out', 'test.prof')
plot_profiles('test.prof')


