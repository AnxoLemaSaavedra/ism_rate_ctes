import numpy           as np
import common.physcons as pc



#----------------------------------------------------#
def center_of_mass(masses, xcc):
    '''
    Calculate center of masses of a give molecule
    
    Input: [masses], [[x_i, y_i, z_i]]
    
    Output: [x, y, z] center of masses
    '''
    M = np.sum(masses)
    cm = np.sum((masses*xcc.T).T, axis = 0)
    cm = cm / M
    return cm
#----------------------------------------------------#

#----------------------------------------------------#
def rotational_constants(molecule):
    '''
    Get rotational constants from molecule

    Input: Molecule: molecule

    Output: list<float>: rotational_constants (a.u.)
    '''
    Ia, Ib, Ic = np.array(molecule._imoms)

    A, B, C    = [pc.HBAR / (2 * I) for I in (Ia, Ib, Ic)]
    A, B, C    = [pc.HBAR / (4 * np.pi * pc.C0 * I) for I in (Ia, Ib, Ic)]
    return A, B, C
#----------------------------------------------------#

#----------------------------------------------------#
def prolate_rotational_constants(molecule):
    '''
    Get rotational constants from molecule
    under prolate approximation (one 1-D and one 2-D rotors)

    Input: Molecule: molecule

    Output: list<float>: rotational_constants (a.u.)
    '''
    # Calculate moments of inertia in the symmetric top approximation
    # i.e. one non degenerate and one double degenerate
    Is, Id = prolate_inertia_moments(molecule)

    Bs, Bd = [pc.HBAR**2 / (2 * I) for I in (Is,Id)]
    return Bs, Bd
#----------------------------------------------------#

#----------------------------------------------------#
def rotational_constants_SI(inertia_moment):
    '''
    Get rotational constants form moments of inertia

    Input:  list<float>: inertia_moment (kg * m**2)

    Output: list<float>: Rotational constants (m-1)
    '''
    # If inertia_moment is an iterable
    try:
        return [pc.H_SI / (8 * np.pi**2 * pc.C0_SI * I_i) for I_i in inertia_moment]

    # If inertia_moment is a float
    except TypeError:
        return pc.H_SI / (8 * np.pi**2 * pc.C0_SI * inertia_moment)
    #return np.array([h**2 / (8 * np.pi**2 * I_i) / 100 for I_i in inertia_moment])
#----------------------------------------------------#

#----------------------------------------------------#
def prolate_inertia_moments(molecule):
    '''
    Calculate prolate top rotational moments of inertia in SI units

    Input: Molecule: molecule

    Output: float: Is(1-D), Id(2-D) (a.u.)
    '''
    Ia, Ib, Ic = np.array(molecule._imoms) # a.u.

    # Approximate to symmetric top
    if abs(Ia - Ib) < abs(Ia - Ic):
        if abs(Ia - Ib) < abs(Ib - Ic):
            Id = (Ia * Ib)**0.5
            Is = Ic
        else:
            Id = (Ib * Ic)**0.5
            Is = Ia
    else:
        if abs(Ia - Ic) < abs(Ib - Ic):
            Id = (Ia * Ic)**0.5
            Is = Ib
        else:
            Id = (Ib * Ic)**0.5
            Is = Ia
    return Is, Id
#----------------------------------------------------#

#----------------------------------------------------#
def moments_of_inertia(masses, xcc):
    '''
    Calculate principal moments of inertia

    Input: list<float>:       masses
           list<list<float>>: xcc

    Ouput: float: Ia, Ib, Ic (mass * dist**2)
    '''
    # Get center of mass
    cm = center_of_mass(masses, xcc)
    # Move cartesians to center of mass
    xcc = xcc - cm
    #Initialize Inertia tensor
    I = np.zeros((3,3))

    # Build inertia tensor
    for ii in range(3):
        # Diagonal terms
        I[ii][ii] = np.sum([masses[kk] * (np.linalg.norm(xcc[kk])**2 - xcc[kk][ii]**2)\
                            for kk in range(len(masses))])
        # Non-diagonal terms
        for jj in range(ii+1,3):
            I[ii][jj] = np.sum([-masses[kk] * xcc[kk][ii] * xcc[kk][jj]\
                                for kk in range(len(masses))])
            I[jj][ii] = I[ii][jj]

    # Diagonalize inertia tensor
    (Ia, Ib, Ic), vect = np.linalg.eigh(I)
    return Ia, Ib, Ic
#----------------------------------------------------#
