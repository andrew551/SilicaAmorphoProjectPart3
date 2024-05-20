import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import statsmodels.api as sm
lowess_sm = sm.nonparametric.lowess

hplanck = 6.6261e-34
kb = 1.380649e-23

eJ=1.602176487e-19
a_to_m=1.0e-10
amu=1.660599e-27
freq_to_wthz=np.sqrt(eJ/amu/a_to_m**2)/1.0e12
hbar = 1.054571e-34

quartic_unit = eJ / amu**2 / a_to_m**4 * hbar

c_cyan_nature=np.array([68/255, 128/255, 244/255,1.])
c_red=np.array([206./255., 30./255., 65./255.,1.0])
dark_red='#9e0b00'
c_green_nature='#63cf8c'#np.array([96./255., 172./255., 63./255.,1.0])
c_blue_nature=np.array([54./255., 79./255., 156./255.,1.0])
ec=[c_blue_nature,c_green_nature,c_red, 'yellow', 'purple']

def get_occupation(frequency, T):
    arg = frequency * hbar / (kb * T)
    arg = np.maximum(np.minimum(arg, 100), -100) # clip large values
    mode_occupation = 1/(np.exp(arg) - 1)
    return mode_occupation

'''
solve the SCP equation iteratively
'''
def compute_correction(phi4, frequency_guess, T):
    frequency_guess0 = frequency_guess
    for i in range(1):
        occupations = get_occupation(frequency_guess, T)
        inv_freqs = 1 / frequency_guess
        correction = 1 / (2) * quartic_unit * np.einsum('ij, j, j -> i', phi4, occupations + 0.5, inv_freqs)
        frequency_guess = np.sqrt(frequency_guess**2 + correction)
    #frequency_guess[:100] = frequency_guess0[:100]
    eigenvalues_corrected = (frequency_guess / freq_to_wthz / 1e12) ** 2

    return frequency_guess, correction, eigenvalues_corrected

'''
compute smoothed frequency shifts (in Hz)

'''

def compute_frequency_shifts_smoothed(phi4, frequencies0, T):
    new_freqs, correction, eigenvalues = compute_correction(phi4, frequencies0, T)
    delta = new_freqs - frequencies0
    trend = lowess_sm(delta, frequencies0, frac =1/10, it=5, return_sorted=False)
    return trend

if __name__ == '__main__':

    phi4 = np.loadtxt('10_quartic/pet.dat')
    n = phi4.shape[0]
    eigenvalues = np.loadtxt('3_diagonalise/eigenvalues.dat')
    frequencies0 = np.sqrt(eigenvalues) * freq_to_wthz * 1e12

    ### remove translational modes

    n = n-3
    phi4 = phi4[3:, 3:]
    eigenvalues = eigenvalues[3:]
    frequencies0 = frequencies0[3:]
    ###




    '''
    new_freqs, correction, _ = compute_correction(frequencies0, 500)
    plt.plot(frequencies0)
    plt.plot(new_freqs)
    plt.savefig("10_quartic/freq_corrected_plot.png")

    plt.clf()

    plt.plot(correction)
    plt.savefig("10_quartic/correction_img.png")

    np.savetxt('10_quartic/correction_values.dat', correction)
    '''

    ### make some plots

    Ts = [1, 100, 300, 600]
    new_freqs_arr = []
    plt.clf()
    for i,T in enumerate(Ts):
        new_freqs, correction, eigenvalues = compute_correction(phi4, frequencies0, T)
        print(np.any(np.isnan(new_freqs)))
        print(new_freqs)
        new_freqs_arr.append(new_freqs)
        delta = new_freqs - frequencies0
        plt.scatter(frequencies0, delta, s=2, label = str(T), color=ec[i])
        trend = lowess_sm(delta, frequencies0, frac =1/10, it=5, return_sorted=False)
        plt.plot(frequencies0, trend, color=ec[i])
    plt.legend()
    plt.savefig("10_quartic/shifts_T.png")

    plt.clf()
    plt.xlabel('$\hbar\omega$ $(cm^{-1})$')
    plt.ylabel('DOS')
    for T in Ts:
        new_freqs, correction, eigenvalues = compute_correction(phi4, frequencies0, T)
        v = np.sqrt(np.abs(eigenvalues)*9.648e27)/18.8e10
        x = np.linspace(-np.max(v)*1.02, np.max(v)*1.02, 4000)


        sigma = 10
        conv = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-x**2/(2*sigma**2))

        histo, _ = np.histogram(v, bins=np.append(x, x[-1]*2 - x[-2])) # bins have one more element than x as they represent edges
        dos = np.convolve(conv, histo, mode = 'same')

        plt.plot(x[x > 0], dos[x > 0], label = f'T={T}K')
    plt.legend()    
    plt.savefig('10_quartic/DOS_plots_T.png', dpi=400)

    plt.clf()
    plt.xlabel('$\hbar\omega$ $(cm^{-1})$')
    plt.ylabel('$g(\omega) / \omega^2$')
    for T in Ts:
        new_freqs, correction, eigenvalues = compute_correction(phi4, frequencies0, T)
        v = np.sqrt(np.abs(eigenvalues)*9.648e27)/18.8e10
        x = np.linspace(-np.max(v)*1.02, np.max(v)*1.02, 4000)
        sigma = 10
        conv = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-x**2/(2*sigma**2))

        histo, _ = np.histogram(v, bins=np.append(x, x[-1]*2 - x[-2])) # bins have one more element than x as they represent edges
        dos = np.convolve(conv, histo, mode = 'same')
        sub = np.logical_and(10 < x, x < 200)
        plt.plot(x[sub], dos[sub] / x[sub]**2, label = f'T={T}K')
    plt.legend()
    plt.savefig('10_quartic/DOS_over_omega_sq_T.png', dpi=400)
    '''
        fig, ax = plt.subplots()

        
        
        plt.savefig(outfile2, dpi=400)
        plt.show()
    '''


    #####



