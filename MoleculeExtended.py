import numpy as np
from   common.Molecule import Molecule


class MoleculeExtended(Molecule):
    def __init__(self,label=None):
        super().__init__(label)
        self._irInten       = None
        self._dipoleMoment = None


          #---------------------#


    def ir_intensities_from_gauout(self, gauout):
        '''
        Read IR intensities from Gaussian output file
        Store intensities in km mol-1
        '''
        with open(gauout,"r") as asdf: lines = asdf.readlines()
        lines = [line.split() for line in lines]
        
        # Get IR intensities from gaussian output file
        # units in km / mol
        ir_intensities = []
        for line in lines:
            try:
                if line[:2] == ["IR", "Inten"]:
                    for inten in line[3:]: ir_intensities.append(float(inten))

            except: pass

        self._irInten = ir_intensities


          #---------------------#


    def set_dipole_moment(self, dipole_moment):
        self._dipoleMoment = dipole_moment


          #---------------------#


    def set_dipole_moment_from_gauout(self, gauout):
        '''
        Read dipole moment from Gaussian output file
        Store dipole moment in atomic units
        '''
        with open(gauout,"r") as asdf: lines = asdf.readlines()
        data = "".join([line.strip() for line in lines])
        data = data.split("Dipole=")[-1].split("\\")[0].split(",")

        dipole_moment = [float(x) for x in data]
        self._dipoleMoment = np.linalg.norm(dipole_moment)


          #---------------------#


    def set_dipole_moment_from_fchk(self, fchk):
        '''
        Read dipole moment from Gaussian output file
        Store dipole moment in atomic units
        '''
        with open(fchk,"r") as asdf: lines = asdf.readlines()

        # Loop over lines in Gaussian Output File
        for i,line in enumerate(lines):
            # If dipole moment was found do NOT break loop
            if "Dipole Moment" in line:
                dipole_moment = [float(x) for x in lines[i+1].split()]
                break
        self._dipoleMoment = np.linalg.norm(dipole_moment)
        self._dipoleMoment = dipole_moment
#----------------------------------------------------#

