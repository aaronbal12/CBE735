#!/.conda/envs/ase/bin/python

import ase
from ase import io
from ase.io import read,write
from ase.neb import NEB
from ase.constraints import FixScaled, FixAtoms
import glob, os
import shutil
import sys
from getopt import getopt



class NEBSkeleton(object):
    def __init__(self, initial = [], final = [], nImages = 7, interpolation = 'idpp', nebtype = 'nebnormal', fmax=0.0001):
        if not initial or not final:
            nImages = 0
        images = self.setupNEBImages(initial, final, nImages)
        self.Images = images
        self.Interpolation = interpolation
        self.NEBType = nebtype
        self.Fmax = fmax
        self.Silent = False

    @property
    def nImages(self):
        return len(self.Images)

    def setupNEBImages(self, initial, final, nImages):
        images = [initial]
        images += [initial.copy() for i in range(int(nImages))]
        images += [final]
        return images

    def printInfo(self):
        print("ase version is: %s" % ase.__version__)
        print("creating neb with %s images, %s interpolation" % (str(self.nImages - 2), self.Interpolation))
        print("type of neb is %s\n" % self.NEBType )
    
    def interpolate(self):
        constraints = self.Images[0].constraints
        if self.NEBType == 'nebtr':  # only remove constraints for nebtr, as it goes wrong if we don't do so. for idpp, removing constraints leads to wild movement of fixed atoms
            for image in self.Images:
                image.constraints = []

        if self.NEBType == 'nebtr':
            try:
                neb = NEB(self.Images, remove_rotation_and_translation=True)
            except:
                neb = NEB(self.Images)
        else:
                neb = NEB(self.Images)


        if self.Interpolation == 'linear':
            neb.interpolate('linear', mic=True)
        elif self.Interpolation == 'idpp':
            neb.interpolate('linear', mic=True)
            neb.idpp_interpolate(mic=True, fmax=self.Fmax, steps=250)
        elif self.Interpolation == 'idpp-direct':
            neb.idpp_interpolate(mic=True, fmax=self.Fmax, steps=250)
        elif self.Interpolation == 'none':
            pass
        else:
            print('Error, interpolation method %s not supported!' % self.Interpolation)
            exit(1)

        # Set back original constraints
        if self.NEBType == 'nebtr':
            self.appendConstraints(neb, constraints)
        return neb

    def appendConstraints(self, neb, constraints):
        # Blank fixed scaled constraint so that ase.io.write prints selective dynamics and all T T T
        blankConstraint = FixScaled(self.Images[0].cell,0, mask=(0,0,0))
        for image in neb.images:
            newConstraints = [constraint for constraint in image.constraints]
            newConstraints.append(constraints)
            newConstraints = [constraint for constraint in newConstraints if constraint]

            if not newConstraints:
                newConstraints.append(blankConstraint)
            image.set_constraint(newConstraints)



class NEBWriter(object):
    def __init__(self, nebdir, idximagestowrite='all'):
        self.Dir = nebdir
        self.IdxImagesToWrite = idximagestowrite



    def write(self, neb):
        if self.IdxImagesToWrite == 'all':
            images = [neb.images[i] for i in range(len(neb.images))]
        else:
            images = [neb.images[i] for i in self.IdxImagesToWrite]

        executingdir = os.getcwd()
        self.prepareNEBFolder(executingdir)
        os.chdir(self.Dir)
        for [imagecounter, image] in enumerate(images):
            folder = '%.2i' % imagecounter
            os.system("mkdir %s" % folder)
            image.write('POSCAR_tmp', format='vasp',direct=True,vasp5=True)
            shutil.copy('POSCAR_tmp', folder + '/POSCAR')
        
        os.chdir(executingdir)
    
    def prepareNEBFolder(self, executingdir):
        if not os.path.abspath(self.Dir) == executingdir:
            if not os.path.exists(self.Dir):
                os.makedirs(self.Dir)
            
            self.copyFilesToNEBFolder()
    
    def copyFilesToNEBFolder(self):
        auxilaryFiles = ['INCAR','KPOINTS','POTCAR','OUTCAR_start','OUTCAR_end','POSCAR_start','POSCAR_end']
        auxilaryFiles = [os.path.expanduser(file) for file in auxilaryFiles]
        auxilaryFiles = [os.path.join(os.getcwd(), filename) if not filename.startswith("/") else filename for filename in auxilaryFiles ]

        for file in auxilaryFiles:
            try:
                shutil.copy(file, self.Dir)
            except:
                pass
    
    def copyOutcarsToStartAndEndFolders(self, OUTCAR_START="OUTCAR_start", OUTCAR_END="OUTCAR_end"):
        outcarfiles = [OUTCAR_START, OUTCAR_END]
        outcarfiles = [os.path.expanduser(file) for file in outcarfiles]
        outcarfiles = [os.path.join(os.getcwd(), filename) if not filename.startswith("/") else filename for filename in outcarfiles ]

        folders = glob.glob('[0-9][0-9]')
        folders.sort()

        try:
            shutil.copy(OUTCAR_START, '%s/OUTCAR' % folders[0])
        except:
            pass
        
        try:
            shutil.copy(OUTCAR_END, '%s/OUTCAR' % folders[-1])
        except:
            pass

class NullNEBWriter(NEBWriter):
    def __init__(self):
        NEBWriter.__init__(self, nebdir=None)

    def write(self, neb):
        return

    def copyOutcarsToStartAndEndFolders(self, OUTCAR_START="OUTCAR_start", OUTCAR_END="OUTCAR_end"):
        return

###############################################################################
def parseOptions():
    options = getDefaultOptions()
    opts, args = getopt(sys.argv[1:], "fn:s:r")
    for opt, arg in opts:
        if opt == '-n':
            options["nImages"] = int(arg)
        elif opt == '-s':
            options["ScanningAtomIdx"] = int(arg) # index must be in python indexing
        elif opt == '-r':
            options["NEBType"] = 'nebtr'
        else:	
            print("unknown option!\n")
            print("options available: o - optimization")
            sys.exit()

    return options


def getDefaultOptions(options = {}):
    keys = ["nImages", "ScanningAtomIdx", "NEBType", "Interpolation"]
    defaultValues = [7, -1, "nebnormal", "idpp"]
    
    for (key, defaultValue) in zip(keys, defaultValues):
        if not key in options:
            options[key] = defaultValue 

    return options


def readNEBInitialAndFinalPoscars(INITIAL_FILENAME = "POSCAR_start", FINAL_FILENAME = "POSCAR_end"):
    def checkIfStartAndEndPoscarsExist(startFile, endFile):
        if not os.path.isfile(startFile):
            print("Starting configuration file %s does not exist!" % startFile)
        if not os.path.isfile(endFile):
            print("Ending configuration file %s does not exist!" % endFile)
    
    def checkStartAndEndPoscars(startFile, endFile):
        checkIfStartAndEndPoscarsExist(startFile, endFile)

    
    checkStartAndEndPoscars(INITIAL_FILENAME, FINAL_FILENAME)

    initialConfig = io.read(INITIAL_FILENAME)
    finalConfig = io.read(FINAL_FILENAME)
    
    return (initialConfig, finalConfig)


def copyStartAndEndFoldersPoscarToContcar():
    slicingFolders = glob.glob('[0-9][0-9]')
    slicingFolders.sort()
    
    startFolderPoscar = os.path.join(slicingFolders[0], 'POSCAR')
    endFolderPoscar = os.path.join(slicingFolders[-1], 'POSCAR')
    
    startFolderContcar = os.path.join(slicingFolders[0], 'CONTCAR')
    endFolderContcar = os.path.join(slicingFolders[-1], 'CONTCAR')
    
    try:
        shutil.copy(startFolderPoscar, startFolderContcar)
    except:
        pass
    
    try:
        shutil.copy(endFolderPoscar, endFolderContcar)
    except:
        pass


def main(nebskeleton, nebwriter, silent=False):
    if not silent:
        nebskeleton.printInfo()

    neb = nebskeleton.interpolate()
    nebwriter.write(neb)
    nebwriter.copyOutcarsToStartAndEndFolders()
    
    return neb

#################################################################################

if __name__ == "__main__":  # Runs when called as a script, but not when imported
    options = parseOptions()

    (initialConfig, finalConfig) = readNEBInitialAndFinalPoscars()


    nebskele = NEBSkeleton(initial=initialConfig, final=finalConfig, nImages=options["nImages"],
                           nebtype=options["NEBType"])

    nebwriter = NEBWriter(nebdir='./')
    main(nebskele, nebwriter)
