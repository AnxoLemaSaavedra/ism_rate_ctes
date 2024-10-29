# Standar libraries
#--------------------------------------#
import numpy           as np

# Modules from common library, Dr. David Ferro-Costas
#--------------------------------------#
import common.physcons as pc

# Modules from Anxo Lema-Saavedra
#--------------------------------------#
import rotation        as rot


#----------------------------------------------------#
def rrkm_classical_rotational_density_of_states(molecule, energy, deltaE=1, rot_symm=1):
    '''
    Calculate classical rotational density of states
    assuming top symmetry, but only for the non-degenerate
    axis of inertia. The degenerate one will not be an active
    mode in the RRKM theory:
                  Robert G. Gilbert and Sean C. Smith
                  Theory of Unimolecular and Recombinations Reactions
                  Pages: 54-55, from 70 onwards and 145

    Input: Molecule: molecule
           float:    energy (cm-1), deltaE
           integer:  rot_symm
    
    Output: list<float>: sum_of_states (cm), energies (cm-1)
    '''
    # Check that Deltae value is not so bad
    freqs     = np.array(molecule._icfreqs) * pc.H2CM
    if deltaE > (min(freqs) / 10.): deltaE = min(freqs) / 10.

    # Approximate frequencies to the grid space
    freqs     = np.round(freqs / deltaE).astype(int)
    nE_points = int(round(energy / deltaE)) + 1
    energies  = np.arange(nE_points) * deltaE

    # Calculate rotational constants in the symmetric top approximation
    # i.e. one non degenerate and one double degenerate
    A, B = rot.prolate_rotational_constants(molecule)
    A    = A * pc.JOULE / pc.H_SI / pc.C0_SI / 100 # From Hartrees to cm-1
    B    = B * pc.JOULE / pc.H_SI / pc.C0_SI / 100 # From Hartrees to cm-1

    density_of_states = np.zeros(nE_points)
    for i in range(nE_points):
        # 2-D degenerate rotor does not couple with vibrations
        # therefore should not be taken into accout for 
        # sum of states in RRKM theory

        # Calculate rotational Sum of States for the single-degenerate axis of inertia
        # under RRKM approximation of active/inactive rotational modes
        if i == 0: density_of_states[i] = 1
        else:      density_of_states[i] = (1 / (energies[i] * A))**0.5 / rot_symm

        # Ensure quantum consistecy at very low energies
        if density_of_states[i] < 1 and density_of_states[i] > 0:
            density_of_states[i] = 1.

    return density_of_states, energies
#----------------------------------------------------#

#----------------------------------------------------#
def rrkm_classical_rotational_sum_of_states(molecule, energy, deltaE=1, rot_symm=1):
    '''
    Calculate classical rotational sum of states
    assuming top symmetry, but only for the non-degenerate
    axis of inertia. The degenerate one will not be an active
    mode in the RRKM theory:
                  Robert G. Gilbert and Sean C. Smith
                  Theory of Unimolecular and Recombinations Reactions
                  Pages: 54-55, from 70 onwards and integral over equation 4.3.31 in 145

    Input: Molecule: molecule
           float:    energy, deltaE
           integer:  rot_symm
    
    Output: list<float>: sum_of_states (cm), energies (cm-1)
    '''
    # Check that Deltae value is not so bad
    freqs     = np.array(molecule._icfreqs) * pc.H2CM
    if deltaE > (min(freqs) / 10.): deltaE = min(freqs) / 10.

    # Approximate frequencies to the grid space
    freqs     = np.round(freqs / deltaE).astype(int)
    nE_points = int(round(energy / deltaE)) + 1
    energies  = np.arange(nE_points) * deltaE

    # Calculate rotational constants in the symmetric top approximation
    # i.e. one non degenerate and one double degenerate
    A, B = rot.prolate_rotational_constants(molecule)
    A    = A * pc.JOULE / pc.H_SI / pc.C0_SI / 100 # From Hartrees to cm-1
    B    = B * pc.JOULE / pc.H_SI / pc.C0_SI / 100 # From Hartrees to cm-1

    ## Check if exists symmetry

    sum_of_states = np.zeros(nE_points)
    for i in range(nE_points):
        # Como seria par aun sigma /= 1 no rotor 2D?
        # O rotor 2-dexenerado non acopla coas vibracions
        # polo que non debemos telo en conta para a suma de estados

        # Calculate rotational Sum of States for the single-degenerate axis of inertia
        # under RRKM approximation of active/inactive rotational modes
        sum_of_states[i] = 2 * (energies[i] / A)**0.5 / rot_symm

        # Ensure quantum consistecy at very low energies
        if sum_of_states[i] < 1: sum_of_states[i] = 1.

    return sum_of_states, energies
#----------------------------------------------------#

#----------------------------------------------------#
def quantum_vibrational_density_of_states(molecule, energy, deltaE = 1):
    '''
    Calculate the density of states for each energy
    in the grid of deltaE between zero and energy
    for given Harmonic-Oscillator frequencies through
    Beyer-Swinehart direct-count algorithm

    Algorithm described in:
                  Robert G. Gilbert and Sean C. Smith
                  Theory of Unimolecular and Recombinations Reactions
                  Page 153

    Input: Molecula: molecule
           float:    energy (cm-1), deltaE (cm-1)

    Output: list<float>: density_of_states (cm), energies (cm-1)
    '''
    # Check that deltaE value is not so bad
    freqs     = np.array(molecule._icfreqs) * pc.H2CM
    if deltaE > (min(freqs) / 10.): deltaE = min(freqs) / 10.

    # Approximate frequencies to the grid space
    freqs     = np.round(freqs / deltaE).astype(int)
    nE_points = int(round(energy / deltaE)) + 1
    energies  = np.arange(nE_points) * deltaE

    # Initialize Density of States
    density_of_states = np.zeros(nE_points)
    density_of_states[0] = 1

    # Beyer-Swinehart direct-count algorithm
    for j, freq in enumerate(freqs):
        for i in range(freq, nE_points):
            density_of_states[i] += density_of_states[i - freq]

    return density_of_states, energies
#----------------------------------------------------#

#----------------------------------------------------#
def quantum_vibrational_sum_of_states(molecule, energy, deltaE=1):
    '''
    Calculate the sum of states for each energy
    in the grid of deltaE between zero and energy
    for given Harmonic-Oscillator frequencies through
    the modified Beyer-Swinehart direct-count algorithm

    Algorithm described in:
                  Robert G. Gilbert and Sean C. Smith
                  Theory of Unimolecular and Recombinations Reactions
                  Page 153

    Input: Molecula: molecule
           float:    energy (cm-1), deltaE (cm-1)

    Output: list<float>: density_of_states (cm), energies (cm-1)
    '''
    # Check that Deltae value is not so bad
    freqs     = np.array(molecule._icfreqs) * pc.H2CM
    if deltaE > (min(freqs) / 10.): deltaE = min(freqs) / 10.

    # Approximate frequencies to the grid space
    freqs     = np.round(freqs / deltaE).astype(int)
    nE_points = int(round(energy / deltaE)) + 1
    energies  = np.arange(nE_points) * deltaE

    # Initialize Density of States
    sum_of_states = np.ones(nE_points)

    # Modified Beyer-Swinehart direct-count algorithm
    for j, freq in enumerate(freqs):
        for i in range(freq, nE_points):
            sum_of_states[i] += sum_of_states[i - freq]

    return sum_of_states, energies
#----------------------------------------------------#

#----------------------------------------------------#
def rrkm_semiclassical_rovib_density_of_states(molecule, energy, deltaE=1):
    '''
    Calculate the density of states for each energy
    in the grid of deltaE between zero and energy
    for given Harmonic-Oscillator frequencies through
    modified Beyer-Swinehart employing direc-count for
    vibrations and and classical expresion for rotations.

    Algorithm described in:
                  Robert G. Gilbert and Sean C. Smith
                  Theory of Unimolecular and Recombinations Reactions
                  Page 157

    Input: Molecula: molecule
           float:    energy (cm-1), deltaE (cm-1)

    Output: list<float>: density_of_states (cm), energies (cm-1)
    '''
    # Check that Deltae value is not so bad
    freqs     = np.array(molecule._icfreqs) * pc.H2CM
    if deltaE > (min(freqs) / 10.): deltaE = min(freqs) / 10.

    # Approximate frequencies to the grid space
    freqs     = np.round(freqs / deltaE).astype(int)
    nE_points = int(round(energy / deltaE)) + 1
    energies  = np.arange(nE_points) * deltaE

    # Intialize density with rotational density
    # evaluated classicaly
    density_of_states, _ = rrkm_classical_rotational_density_of_states(molecule, energy, deltaE=deltaE)

    # Convolute over vibrational levels
    for j, freq in enumerate(freqs):
        for i in range(freq, nE_points):
            density_of_states[i] += density_of_states[i - freq]

    return density_of_states, energies
#----------------------------------------------------#

#----------------------------------------------------#
def rrkm_semiclassical_rovib_sum_of_states(molecule, energy, deltaE=1):
    '''
    Calculate the sum of states for each energy
    in the grid of deltaE between zero and energy
    for given Harmonic-Oscillator frequencies through
    modified Beyer-Swinehart employing direc-count for
    vibrations and and classical expresion for rotations.

    Algorithm described in:
                  Robert G. Gilbert and Sean C. Smith
                  Theory of Unimolecular and Recombinations Reactions
                  Page 157

    Input: Molecula: molecule
           float:    energy (cm-1), deltaE (cm-1)

    Output: list<float>: sum_of_states (cm), energies (cm-1)
    '''
    # Check that Deltae value is not so bad
    freqs     = np.array(molecule._icfreqs) * pc.H2CM
    if deltaE > (min(freqs) / 10.): deltaE = min(freqs) / 10.

    # Approximate frequencies to the grid space
    freqs     = np.round(freqs / deltaE).astype(int)
    nE_points = int(round(energy / deltaE)) + 1
    energies  = np.arange(nE_points) * deltaE

    # Intialize density with rotational density
    # evaluated classicaly
    sum_of_states, _ = rrkm_classical_rotational_sum_of_states(molecule, energy, deltaE=1)

    # Convolute over vibrational levels
    for j, freq in enumerate(freqs):
        for i in range(freq, nE_points):
            sum_of_states[i] += sum_of_states[i - freq]

    return sum_of_states, energies
#----------------------------------------------------#

#----------------------------------------------------#
def vibrational_probability(molecule, i_mode, energy=None, deltaE=1, temperature=None, ensemble="microcanonical", nmax=100):
    '''
    probability_i(n,E) = density_{n-1,i}(E-n*h*freq_i) / density(E)
    where density_{n-1,i} is the density of states for the molecule
    in the absence of the ith degree of freedom

    Employs the Beyer-Swinehart direct-count algorithm
    Described in:
                  Robert G. Gilbert and Sean C. Smith
                  Theory of Unimolecular and Recombinations Reactions
                  Page 156

    Calculate the probability of vibrational mode i
    being in level n for all n \ (n*freq[i_mode] < energy or n < nmax)
    !!!! PURE VIBRATIONAL. ROTATION IS CONSIDERED AS J = 0
    
    Input: list<float>: frequencies(cm-1)
           float:       energy (cm-1), deltaE (cm-1), temperature (K)
           string:      ensemble
           integer:     i_mode, nmax
    
    Output: list<float>: probability
    '''

    if ensemble == "microcanonical":
        if energy is None:
            print("Vibrational probability can NOT be called for microcanonical ensemble without provide an energy")
            raise Exception

        # Check that Deltae value is not so bad
        freqs     = np.array(molecule._icfreqs) * pc.H2CM
        if deltaE > (min(freqs) / 10.): deltaE = min(freqs) / 10.

        #CALCULATE DENSISTY OF STATES FOR ALL MODES
        # Approximate frequencies to the grid space
        freqs     = np.round(freqs / deltaE).astype(int)
        # Add 1 to number of points to include the zero and the provided energy
        nE_points = int(round(energy / deltaE)) + 1
        energies  = np.arange(nE_points) * deltaE

        # Initialize Density of States
        density_of_states = np.zeros(nE_points)
        density_of_states[0] = 1

        # Beyer-Swinehart direct-count algorithm
        for j, freq in enumerate(freqs):
            for i in range(freq, nE_points):
                density_of_states[i] += density_of_states[i - freq]


        #CALCULATE DENSISTY OF STATES WITHOUT MODE i
        # Remove the frequency of the i mode
        freqs_minus_i = [freq for i,freq in enumerate(freqs) if i != i_mode]

        # Initialize Density of States without i mode
        density_minus_i = np.zeros(nE_points)
        density_minus_i[0] = 1

        # Beyer-Swinehart direct-count algorithm without i mode
        for j, freq in enumerate(freqs_minus_i):
            for i in range(freq, nE_points):
                density_minus_i[i] += density_minus_i[i - freq]


        # Initialize probabilty as 0
        probability = np.zeros(nmax)
        #pos = int(round(energy / deltaE))

        # If doesnt exist a level at desired energy, will not be level
        # with energy minus n*h*v => forbiden level returns 0 for all n
        if density_of_states[nE_points-1] < 1.e-15:
            return probability

        for n in range(nmax):
            pos_minus_i = nE_points - n*freqs[i_mode] - 1

            # If energy was exceeded by n
            if pos_minus_i < 0:
                return probability

            else:
                # Get probability
                probability[n] = density_minus_i[pos_minus_i] / density_of_states[nE_points-1]


    if ensemble == "canonical":
        if temperature is None:
            print("Vibrational probability can NOT be called for canonical ensemble without provide a temperature")
            raise Exception

        beta = 1 / (pc.KB * temperature)
        # Frequency of normal mode i
        freq_i = molecule._icfreqs[i_mode]
        # Get probability
        probability = [np.exp(-n * pc.HBAR * freq_i * beta) *\
                       ( 1 - np.exp(-pc.HBAR * freq_i * beta)) for n in range(nmax)]


    return probability
#----------------------------------------------------#




