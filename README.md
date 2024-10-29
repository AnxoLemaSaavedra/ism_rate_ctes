# ism_rate_ctes

This is a Beta version, the program is still in development. 

It works well (I think), but the standard output is so ugly, please don't be angry.

<br>

This program calculates rate constants of interest for ISM processes.

This program calculates three types of rate constants:

  - RRKM. Calculates the unimolecular reaction rate constant
    under the RRKM approximation for a given path.
    If several points along the path are provided
    the variational calculation is done. Else, the
    standard RRKM is calculated.

  - Association. Calculates the association rate constant for
    two molecules employing the Klippenstein equation for
    dipole-dipole systems. Other kind of systems can be easily
    got by taking the correct Klippenstein equation[[1]](#longRangeInteractions).
    
  - Radiative. Calculates the radiative rate constant
    for the vibration emission (IR). The energy of the
    emitted phton is not explicitily know, but is
    typically dominated by 1E-3 cm-1 region[[2]](#radiation).

  - The effective rate constant for reactions where the
    three processes are present.
    

# How to run
```python main.py [ Options ] input_file```

### Options:
  - ```--help/-h```: prints help message.
  
  - ```--assoc```  : Calculate the association rate constant.
  
  - ```--rrkm```   : Calculate the RRKM rate constant for a given path.
  
  - ```--rad```    : Calculate the radiative rate constant.

If NO option is provided, the three processes are calculated,
and also the effective rate constant, asuming the path for RRKM is
the dissociation reaction (competitive with association).

# Input example
```
# This is a comment
# Calculation keywords
ensemble         microcanonical   # Can be canonical or microcanonical
                                  # Default: microcanonical

rotation         4                # set rotation=True
                                  # Default: False 

temps            10 15 20 25 30   # If energy are desired, this should be removed
                                  # Units are Kelvin
                                  # If microcanonical and temperatures, the energy for each
                                  # molecule is assumed to be average thermal energy: (E = k_b * T)
 
energies         30000 35000 cm-1 # If temperatures are desired, this should be removed
                                  # Last term should be units: cm-1 kcalmol hartrees
                                  # If energy and canonical error
                                  # If energy, temperature for radiation rate constant
                                  # are asumed from average thermal energy: (T = E / k_b)

scale_factor_zpe 0.975            # Scale factor for the ZPE energies
                                  # Default: 1

reaction_coord   3 4 reverse      # optional: number of atoms defining
                                  # distance, angle or dihedral angle (start in 1)
                                  # if present ignores numbering of path points
                                  # and order them by reaction coordinate increasing. 
                                  # If reverse order is decreasing

# Names                                 
reactants        AcNT OH          # Names of reactant files
complex          GLN              # Name of Complex/Product file
path             DISSOCIATION     # Name of directory with path points
                                  # If only one point in the directory
                                  # standard RRKM instead of variational
```

#### Considerations
Rotation can be included through the prolate approximation,
that treats the molecule rotations as a one 2-dimensions degenerated
rotation and one 1-dimension non degenerated one. The non-degenerated
external rotation is considered to be an active mode and is allowed to 
change energy with reaction coordinate. The degenerate rotation is considered
to be a good approximation for the total angular momentum, which should be considered.
Due to this effect, the centrifugal barrier appears.

In general, including the rotation does not affect significantly to the reaction rate constant.
If you are not sure about its importance, just not include it.


## References

<a id="longRangeInteractions">[1]</a>
J. Chem. Phys. 122, 194103 (2005); doi: [10.1063/1.1899603](https://doi.org/10.1063/1.1899603)

<a id="radiation">[2]</a> 
J. Chem. Phys. 104, 4502 (1996); doi: [10.1063/1.471201](https://doi.org/10.1063/1.471201)
