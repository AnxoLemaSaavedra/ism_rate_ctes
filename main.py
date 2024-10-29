'''

*************************************
** Main author: Anxo Lema-Saavedra **
** Date:        06/06/2024         **
*************************************

This program was developed by Anxo Lema-Saavedra
at the Universidade de Santiago de Compostela (Galiza).

It employs the common library created by Dr. David Ferro-Costas
from the Universidade de Santiago de Compostela and available
from the Cathedral Package, code for chemical kinetics calculations:
     https://github.com/cathedralpkg

'''

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
import reactionRateConstants as rrc
import rotation              as rot
from   MoleculeExtended      import MoleculeExtended as Molecule

# Constants
#--------------------------------------#
DEBYE   = 1. / 0.393430307 # Bohr electron to Debye




#----------------------------------------------------------------------------#
def HELP():
    help_string = '''
This program calculates the following reaction rate constants:

    + Association (k_a) (employing the Klippenstein equations: ref)

    + Radiation (k_r) (employing equation:

    + RRKM for a given path k_rrkm

    + Radiative association as: k_a * k_r / (k_r + k_d)
                where k_d is the dissociation rate constant
                calculated through RRKM


#------------------------------------------------------------------#
Use:
      python main.py [Options] input_file
#------------------------------------------------------------------#
Options:

    --help/-h: prints help message

    --assoc  : Calculate association rate constant for
               the reactants in the input file
               Need the reactants names in the input file

    --rrkm   : Calculate RRKM rate constant for
               the path in the input file
               Need the complex and path names in the input file

    --rad    : Calculate IR radiative rate constant for
               the complex in the input file
               Need the complex name in the input file
#------------------------------------------------------------------#

If no option is included the three rate constants are calculated
as well as the radiative association rate constant

The input must have the following structure:
*****************************************************************************
**    # This is a comment
**    # Calculation keywords
**    ensemble         microcanonical   # Can be canonical or microcanonical
**                                      # Default: microcanonical

**    rotation         4                # set rotation=True
**                                      # Default: False
**
**    temps            10 15 20 25 30   # If energy are desired, this should be removed
**                                      # Units are Kelvin
**                                      # If microcanonical and temperatures, the energy for each
**                                      # molecule is assumed to be average thermal energy: (E = k_b * T)
**
**    energies         30000 35000 cm-1 # If temperatures are desired, this should be removed
**                                      # Last term should be units: cm-1 kcalmol hartrees
**                                      # If energy and canonical error
**                                      # If energy, temperature for radiation rate constant
**                                      # are asumed from average thermal energy: (T = E / k_b)
**
**    scale_factor_zpe 0.975            # Scale factor for the ZPE energies
**                                      # Default: 1
**
**    reaction_coord   3 4 reverse      # optional: number of atoms defining
**                                      # distance, angle or dihedral angle (start in 1)
**                                      # if present ignores numbering of path points
**                                      # and order them by reaction coordinate increasing.
**                                      # If reverse order is decreasing
**
**    # Names
**    reactants        AcNT OH          # Names of reactant files
**    complex          GLN              # Name of Complex/Product file
**    path             DISSOCIATION     # Name of directory with path points
*****************************************************************************
'''
    return help_string
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
def handle_options(args):

    allowed_options = ["--input", "-i", "--help", "-h", "--assoc", "--rad", "--rrkm"]
    for option in args:
        if option.startswith("-") and option not in allowed_options:
            print(f"Error! The option {option} is not valid!\n\n")
            raise Exception

    help_opt    = False
    association = False
    radiative   = False
    rrkm        = False
    input_file  = None

    if len(args) == 0:
        help_opt = True
        return input_file, help_opt, association, radiative, rrkm

    if "--help"      in args or "-h" in args: help_opt = True
    if "--assoc"     in args                : association = True
    if "--rad"       in args                : radiative = True
    if "--rrkm"      in args                : rrkm = True

    if "--input" in args:
        idx = args.index("--input")
        try: input_file = args[idx+1]
        except IndexError:
            print("Error! Input file was not provided!")
            raise Exception

    elif "-i" in args:
        idx = args.index("-i")

        try: input_file = args[idx+1]
        except IndexError:
            print("Error! Input file was not provided!")
            raise Exception

    else:
        input_file = args[-1]

    if input_file.startswith("-"):
        print("Error! Input file was not provided!")
        raise Exception


    if not (association or radiative or rrkm):
        association = True
        radiative   = True
        rrkm        = True

    return input_file, help_opt, association, radiative, rrkm
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
def read_input_file(inp_file):
    '''
    Read input file for the program

    Returns: string:      ensemble, reac_A_name, reac_name_B, path_name, complx_name
             boolean:     rotation
             list<float>: energies, temps
             float:       dipolar_moment_A, dipolar_moment_B, scale_factor_zpe
    '''

    possible_ensembles = ("canonical", "microcanonical",)
    ensemble           = "microcanonical"
    reac_name_A        = None
    reac_name_B        = None
    complx_name        = None
    path_name          = None
    energies           = None
    temps              = None
    rotation           = False
    scale_factor_zpe   = 1.
    reaction_coord     = None
    reverse_path       = False
    J                  = 0
    energy_units       = "hartree"

    with open(inp_file,"r") as asdf: lines = asdf.readlines()
    lines = [line.partition("#")[0].strip() for line in lines]

    for line in lines:
        if "ensemble"         in line:
            try:               ensemble = line.split()[1]
            except IndexError: warning += "ensemble keyword was used, but no option was provided. Employing Microcanonical as default"

        if "reactants"        in line:
            reac_name_A = line.split()[1]
            try              : reac_name_B = line.split()[2]
            except IndexError: reac_name_B = None
        if "complex"          in line: complx_name = line.split()[1]
        if "path"             in line: path_name = line.split()[1]
        if "energies"         in line:
            if line.split()[-1].isnumeric():
                energies = [float(energy) for energy in line.split()[1:]]
            else:
                energies = [float(energy) for energy in line.split()[1:-1]]
                energy_units = line.split()[-1]

        if "temps"            in line: temps    = [float(temp) for temp in line.split()[1:]]
        if "rotation"         in line:
            rotation = True
            try:    J = float(line.split()[1])
            except IndexError: J = 0

        if "scale_factor_zpe" in line: scale_factor_zpe = float(line.split()[1])
        if "reaction_coord"   in line:
            reaction_coord = []
            for x in line.split()[1:]:
                try              : reaction_coord.append(int(x)-1)
                except ValueError:
                    if x == "reverse": reverse_path = True
                    else: raise Exception

    if ensemble not in possible_ensembles: raise Exception

    if energy_units   == "cm-1":     energies = [E * pc.CM2H    for E in energies]
    elif energy_units == "kcal/mol": energies = [E / pc.KCALMOL for E in energies]
    elif energy_units == "kj/mol":   energies = [E / pc.KJMOL   for E in energies]
    elif energy_units == "cal/mol":  energies = [E / pc.CALMOL  for E in energies]
    elif energy_units == "j/mol":    energies = [E / pc.JMOL    for E in energies]
    elif energy_units == "j":        energies = [E / pc.JOULE   for E in energies]
    elif energy_units == "cal":      energies = [E / pc.CAL     for E in energies]

    return ensemble, reac_name_A, reac_name_B, complx_name, path_name,\
           energies, temps, rotation, scale_factor_zpe, reaction_coord, reverse_path, J
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
def read_molecule_from_file(mol_name, project_gradient=False, scale_factor_zpe=1.):
    '''
    Create molecule from the mol_name.log file

    Return: Molecule molecule
    '''
    # Run Molecule constructor
    molecule = Molecule()
    # Set molecule data from gaussian output or formatted checkpoint file
    if   mol_name.endswith(".fchk"): molecule.set_from_fchk(mol_name)
    elif mol_name.endswith(".log"):  molecule.set_from_gauout(mol_name)
    else:
        if os.path.exists(f"{mol_name}.fchk"):
            molecule.set_from_fchk(f"{mol_name}.fchk")

        elif os.path.exists(f"{mol_name}.log"):
            molecule.set_from_gauout(f"{mol_name}.log")

        else:
            print(f"There is NO Gaussian .log file nor Gaussian .fchk file for {mol_name} molecule.")
            raise Exception

    molecule.setup(projgrad=project_gradient)
    # Create set of internal coordinates
    internal_coordinates = ics_from_geom(molecule._xcc, molecule._symbols)
    # Set scale factor for the frequencies
    molecule._fscal = scale_factor_zpe
    # Get frequencies in internal coordinates
    molecule.icfreqs(internal_coordinates, bool_pg=project_gradient)
    molecule.ana_freqs("ic")
    return molecule
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
def create_path(path_name, reaction_coord=None, scale_factor_zpe=1., reverse_path=False):
    '''
    Create the path of the reaction as a list of Molecules
    if Reaction coordinate was provided in the input,
    the list will be ordered by increasing this coordinate
    (unless reverse order was set in input file)
    elif, numerical ordered will be search in the file names
    else Exception will be raised

    Input: list<string>:   path_point_files
           reaction_coord: list<int>

    Ouput: list<Molecule>: path_points
    '''
    # Get ts files names
    if os.path.isdir(path_name):
        path_points_files   = [f"{path_name}/{path_point}" for path_point in os.listdir(path_name)]
    else:
        path_points_files   = [path_point for path_point in os.listdir() if path_point.startswith(path_name)]
    #path_points_files   = sorted(path_points_files, key = lambda x: x.name, reverse = True)

    path_points = []
    for path_point_file in path_points_files:
        # Create molecule object
        path_point = read_molecule_from_file(path_point_file, project_gradient=True, scale_factor_zpe=scale_factor_zpe)
        # Check it exist any imaginary frequency after project out the gradient
        # If imaginary frequency exists does not use that structure
        imaginary_freqs = False
        for icfreq in path_point._icfreqs:
            if icfreq < 0.: imaginary_freqs = True
        if imaginary_freqs: continue
        path_points.append(path_point)

    if reaction_coord is not None:
        if len(reaction_coord) == 2:
            at1, at2 = reaction_coord
            distances = [np.linalg.norm(np.array(point._xcc[3*at2:3*at2+3])\
                                      - np.array(point._xcc[3*at1:3*at1+3]))\
                                                     for point in path_points]
        elif len(reation_coord) == 3:
            at1, at2, at3 = reaction_coord
            distances = [fncs.angle(point._xcc[3*at3:3*at3+3],\
                                    point._xcc[3*at2:3*at2+3],\
                                    point._xcc[3*at1:3*at1+3])\
                                       for point in path_points]

        elif len(reaction_coord == 4):
            at1, at2, at3, at4 = reaction_coord
            distances = [fncs.dihedral(point._xcc[3*at4:3*at4+3],\
                                       point._xcc[3*at3:3*at3+3],\
                                       point._xcc[3*at2:3*at2+3],\
                                       point._xcc[3*at1:3*at1+3])\
                                          for point in path_points]

        else:
            print("The reaction coordinate specified in the input file contains more than 4 atoms or less than 2")
            raise Exception

        if len(distances) != len(set(distances)):
            print("There are two points along path with the same distance in the reaction coordinate.")
            raise Exception
        # Order by distances
        decrease_path = not reverse_path # To decrease direction sorted ordering should be reverse
        path_points   = sorted(path_points, key = lambda i:\
                               distances[path_points.index(i)],\
                               reverse=decrease_path)

    else:
        for j,frag in enumerate(path_point_files[0].split(".")):
            if frag.isnumeric(): break
        else: raise Exception

        path_points = sorted(path_points, key = lambda i: path_point_files[path_points.index(i)].split(".")[j])

    return path_points
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
def get_system_total_energy(reactant_A, reactant_B=None,\
                temperature=None, energy=None, J=None):
    '''
    Calculate the total energy of the system as:
        energy = E_{el,A} + E_{el,B} + E_{zpe,A} + E_{zpe,B}

    In the infinitely separation the inertial moment becomes huge
    and we asume the rotational energy to be zero.
    We are going to asume later the angular momentum from the
    first point calculated along the path

    '''
    if reactant_B is None: E0 = reactant_A._icV1
    else                 : E0 = reactant_A._icV1 + reactant_B._icV1

    if temperature is not None:
        # Add thermal energy (Hartrees)
        # We are including 1/2 KbT for reaction coordinate and 1/2 KbT for rotation
        E0 += pc.KB * temperature

    elif energy is not None:
        E0 += energy

        if reactant_B is None:
            Bs, Bd = rot.prolate_rotational_constants(reactant_A)
            E0 += Bs * J * (J+1)

    return E0
#----------------------------------------------------------------------------#

##----------------------------------------------------------------------------#
#def get_system_total_energy(reactant_A, reactant_B, temperature=None, energy=None):
#    '''
#    Calculate the total energy of the system as:
#        energy = E_{el,A} + E_{el,B} + E_{zpe,A} + E_{zpe,B}
#
#    In the infinitely separation the inertial moment becomes huge
#    and we asume the rotational energy to be zero.
#    We are going to asume later the angular momentum from the
#    first point calculated along the path
#
#    '''
#    # Electronic energy of reactants + ZPE
#    E0 = reactant_A._icV1 + reactant_B._icV1
#    if temperature is not None:
#        # Add thermal energy (Hartrees)
#        # We are including 1/2 KbT for reaction coordinate and 1/2 KbT for rotation
#        E0 += pc.KB * temperature
#
#    elif energy is not None:
#        E0 += energy
#
#    return E0
##----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
def get_average_angular_momentum(molecule, temperature=None, J=0):
    '''
    Calculate the average angular momentum at a given
    temperature assuming rotational energy from
    the equipartition theorem (E_rot = k_b * T)

    Input: Molecule:    molecule
           list<float>: temperature (Kelvin), energy (Hartrees)
           boolean:     rotation

    Ouput: float: L2 (a.u.)
    '''

    # Average angular momentum at temperature
    if temperature is not None:
        Is, Id = rot.prolate_inertia_moments(molecule)
        L2 = pc.KB * temperature * Id # hartrees

    # Calculate angular momentum from quantum number J
    else:
        L2 = pc.HBAR**2 * J * (J + 1)

    return L2
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
def get_system_conditions(reactant_A, reactant_B=None, path_points=None, reverse_path=False,\
                               rotation=False, temperatures=None, energies=None, J=0):
    '''
    Calculate the total energy of the system
    and the angular momentum

    Input: Molecule:       reactant_A, reactant_B
           list<Molecule>: path_points
           boolean:        reverse_path
           list<float>:    temperatures, energies
           integer:        J

    Ouput: list<float>: energies, angular_moments
    '''

    # If tempratures were given
    if temperatures is not None:
        # Calculate total energies
        energies = [get_system_total_energy(reactant_A, reactant_B, temperature=temperature)\
                    for temperature in temperatures]

        # If rotation it true calculate angular moments
        if rotation:
            if reverse_path: angular_moments = [get_average_angular_momentum(path_points[-1], temperature=temperature)\
                                           for temperature in temperatures]
            else           : angular_moments = [get_average_angular_momentum(path_points[0],  temperature=temperature)\
                                           for temperature in temperatures]

        # If rotation is false angular moments = 0
        else       : angular_moments = [0. for temperature in temperatures]


    # If energies were given
    elif energies is not None:

        temperatures = [energy / pc.KB for energy in energies]

        energies = [get_system_total_energy(reactant_A, reactant_B, energy=energy, J=J)\
                           for energy in energies]

        # Calculate angular momentum
        if reverse_path: angular_moments = [get_average_angular_momentum(path_points[-1], J=J)\
                           for energy in energies]
        else           : angular_moments = [get_average_angular_momentum(path_points[0],  J=J)\
                           for energy in energies]
    else:
        print("No energies nor temperatures were given in the input file. Exiting!")
        raise Exception

    return temperatures, energies, angular_moments
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
def print_results(kassoc, krad, krrkm, keff, temperatures=None, energies=None, reactant_A=None, reactant_B=None):

    assoc = False
    rad   = False
    rrkm  = False
    eff   = False

    if len(kassoc) > 0: assoc = True
    if len(krad)   > 0: rad   = True
    if len(krrkm)  > 0: rrkm  = True
    if len(keff)   > 0: eff   = True

    # Header
    print_data = ""
    if temperatures is not None: print_data += "Temp(K)   "
    print_data += "Energy(cm-1)"
    if assoc: print_data += "   association(cm3 molecule s-1)"
    if rad  : print_data += "   radiative(s-1)"
    if rrkm : print_data += "   rrkm(s-1)"
    if eff  : print_data += "   effective(s-1)"
    print_data += "\n"

    # Data
    for i in range(len(energies)):
        data = []
        if temperatures is not None:
            data.append(f"{temperatures[i]:>7.2f}")

        energy = energies[i]
        if reactant_A is not None: energy -= reactant_A._icV1
        if reactant_B is not None: energy -= reactant_B._icV1
        energy *= pc.H2CM

        data.append(f"{energy:>12.2f}")

        if assoc : data.append(f"{kassoc[i]:>12.6e}")
        if rad   : data.append(f"{krad[i]:>12.6e}")
        if rrkm  : data.append(f"{krrkm[i]:>12.6e}")
        if eff   : data.append(f"{keff[i]:>12.6e}")

        print_data += "   ".join(data) + "\n"


    return print_data
#----------------------------------------------------------------------------#


#==========================#
#==========================#
# MAIN PROGRAM STARTS HERE #
#==========================#
#==========================#
def main(args):

    try:
        input_file, help_opt, association, radiative, rrkm = handle_options(args)
    except:
        print(HELP())
        return

    if help_opt:
        print(HELP())
        return

    # Read input file
    # Return energies in Hartrees
    ensemble, reac_name_A, reac_name_B, complx_name, path_name,\
            energies, temperatures, rotation, scale_factor_zpe, reaction_coord,\
            reverse_path, J = read_input_file(input_file)

    print(temperatures)
    print(energies)
    if association and reac_name_B is None:
        print("The association rate constant calculation was requested")
        print("but only one reactant was given in the input file. Exiting!")
        raise Exception

    print("ensemble, ", ensemble, "reac_name_A, ", reac_name_A, "reac_name_B, ", reac_name_B, "complx_name, ", complx_name, "path_name, ", path_name,\
            "energies, ", energies, "temps, ", temperatures, "rotation, ", rotation, "scale_factor_zpe, ", scale_factor_zpe, "reaction_coord, " ,reaction_coord,\
            "reverse_path ",reverse_path, "J ", J)
    print()
    print()
    print()
    # Create ractant A object
    print(reac_name_A, scale_factor_zpe)
    print(type(scale_factor_zpe))
    reactant_A = read_molecule_from_file(reac_name_A, scale_factor_zpe=scale_factor_zpe)
    reactant_A.set_dipole_moment_from_gauout(f"{reac_name_A}.log")
    # Create ractant B object
    if reac_name_B is not None:
        reactant_B = read_molecule_from_file(reac_name_B, scale_factor_zpe=scale_factor_zpe)
        reactant_B.set_dipole_moment_from_gauout(f"{reac_name_B}.log")
    else: reactant_B = None
    # Create complex object
    complx     = read_molecule_from_file(complx_name, scale_factor_zpe=scale_factor_zpe)
    complx.ir_intensities_from_gauout(f"{complx_name}.log")

    # Create list of points along path
    path_points = create_path(path_name, reaction_coord, scale_factor_zpe, reverse_path)
    # OLLO PORQUE PATH-POINTS NON CONTEN O COMPLEXO

    # Calculate total energies and angular momentum of the system
    temperatures, energies, angular_moments = get_system_conditions(reactant_A, reactant_B,\
                                 path_points, reverse_path=reverse_path, rotation=rotation,\
                                           temperatures=temperatures, energies=energies, J=J)



    # Loop over temperatures/energies
    # and calculate rate constants
    k_a, k_r, k_b, k_eff = [], [], [], []
    for temperature,energy,angular_moment in zip(temperatures,energies,angular_moments):

        print((energy - complx._icV1) * pc.H2CM)
        print("\n\nTemperature: ", temperature)
        # Calculate association constant
        if association:
            k_a.append(rrc.association_rrc(reactant_A, reactant_B, temperature))
            print("k_association: ", k_a[-1])

        # Calculate radiative constant
        if radiative:
            k_r.append(rrc.radiative_rrc(complx, energy=energy, temperature=temperature, ensemble=ensemble,nmax=100))
            print("k_radiative: ", k_r[-1])

        # Calculate dissociation constant
        if rrkm:
            k_b.append(rrc.rrkm_rrc(complx, path_points, energy, angular_moment, ensemble=ensemble, rotation=rotation))
            print("k_rrkm: ", k_b[-1])

        # Calculate effective rate constant
        if association and radiative and rrkm:
            k_eff.append( k_a[-1] * k_r[-1] / (k_b[-1] + k_r[-1]))
            print("k_effective: ", k_eff[-1])


    print_data = print_results(k_a, k_r, k_b, k_eff, temperatures, energies, reactant_A=reactant_A, reactant_B=reactant_B)
    with open("rate_constants.out","w") as asdf: asdf.write(print_data)




if __name__ == "__main__":
    args = argv[1:]
    main(args)


