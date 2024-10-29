#------------------------------------------------------------------------#
def zmat_from_zmatfile(input_name):
    '''
    Reads Z-matrix from .zmat file

    Output ex.: [
                 ["at1"],
                 ["at2", (1,), ("dist12")],
                 ["at3", (2, 1,), (dist23, ang123,)], ...
                ],
                {"dist12": distancia12, "ang123": angle123, ...}
    '''
    with open(input_name,"r") as asdf: lines = asdf.readlines()
    zmatvals = {}
    zmat, rings, mult_bonds = [], [], []
    bool_variables, bool_rings, bool_multbonds = False, False, False
    bool_zmat = True

    for i,line in enumerate(lines):
        if "Variables:" in line:
            bool_zmat = False
            bool_variables = True
            continue
        if "Rings:" in line:
            bool_variables = False
            bool_rings = True
            continue
        if "MultBonds:" in line:
            bool_variables, bool_rings = False
            bool_multbonds = True
            continue
        if line == "" or line == "\n":
            if bool_zmat:
                bool_zmat = False
                bool_variables = True
                continue
            if bool_variables:
                bool_variables = False
                bool_rings = True
                continue
            if bool_rings:
                bool_variables, bool_rings = False, False
                bool_multbonds = True
                continue

        # Z-matrix and variables names
        if bool_zmat:
            if i == 0:
                zmat.append([line.split()[0]])
            elif i == 1:
                zmat.append([line.split()[0], (int(line.split()[1])-1,), (line.split()[2],)])
            elif i == 2:
                zmat.append([line.split()[0], (int(line.split()[1])-1, int(line.split()[3])-1,), (line.split()[2], line.split()[4],)])
            else:
                zmat.append([line.split()[0], (int(line.split()[1])-1, int(line.split()[3])-1, int(line.split()[5])-1), (line.split()[2], line.split()[4], line.split()[6])])

        # Read variables values
        if bool_variables:
            zmatvals[f"{line.split()[0]}"] = float(line.split()[-1])

        # Read ring(s) specification
        if bool_rings:
            rings.append([])
            for atom in line.split(","):
                if "-" in atom:
                    for j in range(int(atom.split("-")[0])-1, int(atom.split("-")[1])): rings[-1].append(j)
                else:
                    rings[-1].append(int(atom)-1)

        # Read multiple bonds if present
        if bool_multbonds:
            if   "-" in line: mult_bonds.append([int(line.split("-")[0])-1, int(line.split("-")[1])-1])
            elif "," in line: mult_bonds.append([int(line.split(",")[0])-1, int(line.split(",")[1])-1])
            elif " " in line: mult_bonds.append([int(line.split()[0])-1   , int(line.split()[1])-1])

    return zmat, zmatvals, rings, mult_bonds
#------------------------------------------------------------------------#

