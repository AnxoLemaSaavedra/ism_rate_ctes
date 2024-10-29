# Standard libraries
#--------------------------------------#
import os
import numpy as np
from   sys   import argv

# Modules from common library, Dr. David Ferro-Costas
#--------------------------------------#
import common.physcons as pc
from   common.internal import ics_from_geom
from   common          import fncs

# Modules from Anxo Lema-Saavedra
#--------------------------------------#
import rotation              as rot
import energyStates          as enst
from   MoleculeExtended      import MoleculeExtended as Molecule

# Constants
#--------------------------------------#
DEBYE   = 1. / 0.393430307 # Bohr electron to Debye


#------------------------------------------------------------------------#
def association_rrc(reactant_A, reactant_B, temperature):
    '''
    Calculate the association rate constant
    employing Klippenstein for dipole-dipole systems

    Input: Molecule: reactant_A, reactant_B
           float:    temperature

    output: float: reaction rate constant (cm3 molecule-1 s-1)
    '''
    # Get dipole moments in Debyes
    dip_mom_A = reactant_A._dipoleMoment * DEBYE
    dip_mom_B = reactant_B._dipoleMoment * DEBYE

    # Calculate reduced mass in AMU
    mass_A       = np.sum(reactant_A._masses)
    mass_B       = np.sum(reactant_B._masses)
    reduced_mass = mass_A * mass_B / (mass_A + mass_B) * pc.AMU

    # Get reaction rate constant
    return 1.83E-9 * reduced_mass**(-0.5) * (dip_mom_A * dip_mom_B)**(2./3.) * temperature**(-1./6.) # 1.83e-9 conversion factor to (cm3 molecule-1 s-1)
#------------------------------------------------------------------------#

#------------------------------------------------------------------------#
def radiative_rrc(complx, temperature=None, energy=None, angular_moment=0, ensemble="microcanonical", nmax=100, deltaE=1):
    '''
    Calculate the radiative rate constant
    for vibration emission
    If angular momentum is not zero rotational energy
    are substracted from available energy

    Input: Molecule: complx
           float:    energy (Hartree), angular_moment (a.u.), temperature (Kelvin), deltaE (cm-1)
           string:   ensemble (microcanonical or canonical)
           
    Output: float: rate_constant (s-1)
    '''
    if ensemble=="microcanonical":
        if energy is None:
            print("Vibrational probability can NOT be called for microcanonical ensemble without provide an energy")
            raise Exception
        # Get prolate inertia moments
        imomnt_s, imomnt_d = rot.prolate_inertia_moments(complx)
        # Get energy available for vibrational excitations
        E  = energy - complx._icV1 - angular_moment / 2. / imomnt_s
        E *= pc.H2CM # to cm-1
    
    else:
        if temperature is None:
            print("Vibrational probability can NOT be called for canonical ensemble without provide a temperature")
            raise Exception

    # Get frequencies
    freqs = complx._icfreqs
    # Transform to cm-1
    freqs = [freq / 2 / np.pi / pc.C0 / pc.METER / 100 for freq in freqs]
    # Get IR intensities (km mol-1)
    ir_intensities = complx._irInten

    rate_constant = 0.
    for i_mode, (freq,ir_inten) in enumerate(zip(freqs, ir_intensities)):
        # Get probabiliy of vibrational mode i_mode
        # to be in energy level n
        prob_n = enst.vibrational_probability(complx, i_mode, energy=E,\
                 temperature=temperature, ensemble=ensemble, nmax=nmax, deltaE=deltaE)

        # Sumation part
        rate_constant += np.sum([1.25E-7 * n * prob_n[n] * ir_inten * freq**2 for n in range(nmax)])
    
    return rate_constant
#------------------------------------------------------------------------#

#------------------------------------------------------------------------#
def variational_sum_of_states(path_points, energy, angular_moment=0, rotation=False, deltaE=1):
    '''
    Calculates the position of the variational transition state
    for the given path, and the sum of states for that structure

    Input: list<Molecule>: path_points
           float:          energy (Hartrees), angular_moment (a.u.)

    Output: float:    sum_of_states
            Molecule: ts
    '''


    SoS_points = []
    # Loop over all points in the path
    for path_point in path_points:

        # If rotation energy levels are wanted to be counted
        if rotation:
            # Get inertia moments under prolate approximation
            # One 1D inertia moment and one 2D intertia moment
            imomnt_s, imomnt_d = rot.prolate_inertia_moments(path_point)
            # Get available energy in cm-1
            E  = energy - path_point._icV1 - angular_moment / 2. / imomnt_s
            E *= pc.H2CM
            # Get sum of states for the point in the path
            sum_of_states, energies = enst.rrkm_semiclassical_rovib_sum_of_states(complx, E, deltaE=deltaE)

        # If only vibrational energy levels are wanted to be counted
        else:
            # Get available energy in cm-1
            E  = energy - path_point._icV1
            E *= pc.H2CM
            # Get sum of states for the point in the path
            sum_of_states, energies = enst.quantum_vibrational_sum_of_states(path_point, E, deltaE=deltaE)

        SoS_points.append(sum_of_states[-1])

    # Get minimum along path
    SoS_ts = min(SoS_points)
    # Get structure of the transition state
    ts     = path_points[SoS_points.index(SoS_ts)]

    return SoS_ts, ts
#------------------------------------------------------------------------#

#------------------------------------------------------------------------#
def rrkm_rrc(reactant, path_points, energy=None, angular_moment=0, deltaE=1, temperature=None, ensemble="microcanonical", rotation=False):
    '''
    Calculate variational rrkm reaction rate constant for the path
    given in path_points

    Input: Molecule:       reactant
           list<Molecule>: path_points
           float:          energy(Hartrees), angular_moment(a.u.), temperature(K)
           string:         ensemble (canonical or microcanonical)
           boolean:        rotation
    '''

    # Get variational transition state
    # and sum of states of the variational ts
    SoS_ts, ts = variational_sum_of_states(path_points, energy, angular_moment)
    
    # Calculate available energy for energy levels excitation
    imomnt_s, imomnt_d = rot.prolate_inertia_moments(reactant)
    E  = energy - reactant._icV1 - angular_moment / 2. / imomnt_s
    E *= pc.H2CM

    # Calculate density of states for the reactant
    if rotation:
        DoS_complex, energies = enst.rrkm_semiclassical_rovib_density_of_states(reactant, E, deltaE=deltaE)
    else:
        DoS_complex, energies = enst.quantum_vibrational_density_of_states(reactant, E, deltaE=deltaE)

    # Check that the desired energy is not a forbidden level
    if DoS_complex[-1] == 0:
        return 0.

    # Get rate constant
    else:
        return SoS_ts / (DoS_complex[-1] / 100 / pc.C0_SI)
#------------------------------------------------------------------------#
