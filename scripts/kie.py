#!/home/vempatti/.conda/envs/ase/bin/python

# ---------------------- README -----------------------------------------------
#
# Module for calculating vibrational frequencies
# and Infrared intensity from VASP vibrational calculation results
#
# Files needed: OUTCAR, POSCAR, OSZICAR
#
# Need to read mass-weighted Hessian and dynamatical dipolemoment from OUTCAR
#
# by Guowen Peng 03/28/2009
#
# ------------------------------------------------------------------------------
# Version 02082019
# Refactored and adapted for python3
# IR module not tested.
# Added symmetry detection using pymatgen modules
#
# by Benjamin Chen 02/08/2019
# ------------------------------------------------------------------------------
# Acknowledgements: BC, GP, LG, RN, FE, YX, and modules in ASE3.
# ------------------------------------------------------------------------------

import logging
import os
import sys
import warnings
from contextlib import contextmanager
from abc import abstractmethod, ABC
from copy import deepcopy
from functools import lru_cache
from getopt import getopt
from math import sqrt, pi, log, exp

import ase.io
import numpy as np
from ase.constraints import FixAtoms, FixScaled
from scipy import optimize
from ase import Atoms
try:
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
    from pymatgen.core.structure import Molecule, Structure
    from pymatgen.io.xyz import XYZ
    from pymatgen.io.ase import AseAtomsAdaptor
except ImportError:
    print('Could not load pymatgen modules. Symmetry detection for gas phase molecules will not work.\n')

warnings.simplefilter("ignore", np.exceptions.ComplexWarning)
logger = logging.getLogger(__name__)

# Constants
_c = 299792458.  # speed of light
_mu0 = 4.e-7 * pi  # permeability of vacuum
_eps0 = 1 / _mu0 / _c ** 2  # permittivity of vacuum
_Grav = 6.67259e-11  # gravitational constant
_hplanck = 6.6260755e-34  # Planck constant
_hbar = _hplanck / (2 * pi)  # Planck constant / 2pi
_e = 1.60217733e-19  # elementary charge
_me = 9.1093897e-31  # electron mass
_mp = 1.6726231e-27  # proton mass
_Nav = 6.0221367e23  # Avogadro number
_k = 1.380658e-23  # Boltzmann constant
_amu = 1.6605402e-27  # amu
_aa = 1E-10  # angstrom
Debye = 1e11 * _e * _c  # Debye factor:  eA/Debye
eps = 1E-8
NA = 6.0221367e23  # Avogadro's number
C = 299792458. * 100.  # speed of light in cm/s
R = 8.3144598  # Gas Constant J/mol/K
h = 6.62606755E-34  # Planck's Constant m^2*kg/s
kb = R / NA  # Boltzmann's Constant J/K
eV2kJ = _e * _Nav / 1000
eV2J = eV2kJ * 1E3

# Settings
FREQ_CUTOFF = 10  # Frequency cutoff in cm-1 for the Entropy, Cp, ZPE calculation
POSCAR = 'POSCAR'
OUTCAR = 'OUTCAR'



def main(options):
    if not os.path.exists('POSCAR'):
        raise ValueError('POSCAR does not exist.')

    if not os.path.exists('OUTCAR'):
        if os.path.exists('OUTCAR.gz'):
            print('Unzipping OUTCAR.gz')
            os.system('gzip -d OUTCAR.gz')
        else:
            raise ValueError('OUTCAR does not exist.')

    T = 298.15
    P = 101325

    freq_calculator = _create_frequency_calculator(poscar=POSCAR, outcar=OUTCAR, options=options)
    freq_calculator.print_results()
    energy = extract_U0()

    if options['calc_thermo_properties']:
        _calc_thermo(T, P, energy, freq_calculator)

    if options['calc_shomate_parameters']:
        _calc_shomate(T, P, energy, freq_calculator)

    if options['animate']:
        _animate(freq_calculator)

    if options['calc_ir_spectrum']:
        freqs, nonvibs, modes, modetypes = freq_calculator.calc()
        _calculate_ir_spectrum(dof=freq_calculator.dof, frequencies=freqs, outcarfile=OUTCAR,
                               masses=freq_calculator.masses, modes=modes)

def _animate(freq_calculator):
    try:
        from visualize import MovieWriter
    except ImportError:
        raise ImportError('Animation failed as we could not import MovieWriter from visualize.py.')

    print(f'\n-----------------------Animation-----------------------')
    freqs, _, modes, modetypes = freq_calculator.calc()
    logger.info(f'Modes: {modes}')
    atoms = freq_calculator.atoms
    positions = atoms.get_positions(wrap=False)
    symbols = atoms.get_chemical_symbols()
    relaxed_atoms_mask = freq_calculator.mask
    for idx, (freq, mode, modetype) in enumerate(zip(freqs, modes, modetypes)):
        if np.isreal(freq):
            freq_string = f'{np.real(freq):3.1f}'
        else:
            freq_string = f'{np.imag(freq):3.1f}j'


        print(f'Animating {modetype} mode with frequency {freq_string} cm-1...')
        images = []
        for i in [-1, 0, 1, 0, -1]:
            positions_new = positions.copy()
            positions_new[relaxed_atoms_mask] += mode*i/4
            images.append(Atoms(positions=positions_new, symbols=symbols))

        logger.info(f'Modetype {modetype}')
        writer = MovieWriter(images=images, views=('xtopview','xfrontview'))

        filename = f'{modetype}{idx}_{freq_string}.gif'

        logger.info(f'Filename {filename}')
        writer.write(filename)
        os.system(f"convert -fill blue -pointsize 18 -draw \"text 15,20 \'Freq. = {freq_string} [1/cm]\'\" {filename} {filename}")

    print(f'Animations complete.')
    print(f'-------------------------------------------------------\n')

def _calc_shomate(T, P, energy, freq_calculator):
    entropy_calculator = EntropyCalculator.from_freq_calculator(freq_calculator=freq_calculator)
    enthalpy_calculator = EnthalpyCalculator.from_freq_calculator_and_energy(energy=energy,
                                                                             freq_calculator=freq_calculator)
    shomate_calculator = ShomateCalculator(entropy_calculator=entropy_calculator,
                                           enthalpy_calculator=enthalpy_calculator)
    Trange1 = [10, 300]
    Trange2 = [300, 1200]
    shomate_calculator.print_results(Trange1, P)
    shomate_calculator.print_results(Trange2, P)


def _calc_thermo(T, P, energy, freq_calculator):
    # T = float(get_input("Temperature [298.15K] : ", 298.15)) # K
    # P = float(get_input("Pressure [1 atm] : ", 1)) * 101325 # Pa
    entropy_calculator = EntropyCalculator.from_freq_calculator(freq_calculator=freq_calculator)
    enthalpy_calculator = EnthalpyCalculator.from_freq_calculator_and_energy(energy=energy,
                                                                             freq_calculator=freq_calculator)
    printer = ThermoPropertiesPrinter(atoms=freq_calculator.relaxed_atoms,
                                      entropy_calculator=entropy_calculator,
                                      enthalpy_calculator=enthalpy_calculator)
    printer.print(T, P)


def _create_frequency_calculator(poscar, outcar, options):

    if options['internal_coords']:
        freq_calculator = InternalCoordsFrequencyCalculator.from_outcar_and_poscar(poscarfile=poscar, outcarfile=outcar,
                                                                                   gaussian_hessian=options[
                                                                                       'gaussian_hessian'])
    else:
        freq_calculator = MassWeightedCoordsFrequencyCalculator.from_outcar_and_poscar(poscarfile=poscar,
                                                                                       outcarfile=outcar,
                                                                                       gaussian_hessian=options[
                                                                                           'gaussian_hessian'])
    return freq_calculator


def get_input(msg, default):
    user_input = input(msg)
    if user_input:
        return user_input
    else:
        return default


def _calculate_ir_spectrum(dof, frequencies, outcarfile, masses, modes):
    with open(outcarfile) as f:
        lines = f.readlines()
    # read displacement
    for line in lines:
        if line.rfind('Finite differences POTIM') > -1:
            delta = float(line.split()[3])

    # Find # of SCF for each ionic steps [(DOF +1) in total]
    scf = []
    with open('OSZICAR', 'r') as f:
        for i, line in enumerate(f):
            if 'F=' in line:
                scf.append(i)
    nscf = np.array(scf)
    nm_scf = np.zeros(2 * dof + 1, dtype=np.int)
    for j in range(2 * dof + 1):
        nm_scf[j] = nscf[j] - 2 * (j + 1)

        # display nm_scf

        # Read dynamical dipolement for +- displacement
        #
        # LDIPOL = .TURE. and .FALSE. have differnt stories
        #
    dip = []
    tflag = 0
    for line in lines:
        if line.rfind('LDIPOL =      T') > -1:
            tflag = 1
            print("LDIPOL = .TRUE.")

        if line.find('dipolmoment') > -1:
            dip.append(np.array([float(f) for f in line.split()[1:4]]))

    # LDIPOL = .FALSE., one dipolemoment per each ionic step
    #   elif line.find('dipolmoment') > -1:
    #        dip.append(np.array([float(f) for f in line.split()[1:4]]))
    #   else:
    #     display "Cannot find dipole momenet from OUTCAR. I don't know how to calculate IR intensity. Quit. "
    #     sys.exit(1)
    #
    # tmpd: temporary array for dipole moment
    # tflag: 1 for LDIPOLE = .TRUE.; 0 for .FALSE.
    if tflag:
        tmpd = np.zeros([2 * dof + 1, 3])
        j = 0
        for i in nm_scf:
            tmpd[j] = dip[i - 1]
            j += 1
    else:
        tmpd = dip

    # display 'Static dipole moment (Debye):', tmpd[0]* Debye
    dipole = np.array(tmpd)[1:2 * dof + 1]
    # calculate dp/dx=(dipole(+)-dipole(-))/2*delta
    dpdx = np.zeros([dof, 3])
    for j in range(0, dof):
        dpdx[j] = (dipole[2 * j + 1] - dipole[2 * j]) / (2 * delta)

    # mass-weighted dp/dx
    dpdq = np.array([dpdx[j] / sqrt(masses[int(j / 3)]) for j in range(dof)])
    # --- calculate dpdq * modes----
    #
    # X_ki are the elecments of the i-th eigenvector
    #
    # dp/dQ_i=sum(k=1,3N) (dp/dq_k * X_ki)
    modes = modes.T.copy()
    dpdQ = np.dot(modes, dpdq)
    # display dpdQ
    # calculate IR intensities
    # IR_i ~ |dp/dQ_i|**2
    intensities = np.array([sum(dpdQ[j] ** 2) for j in range(dof)])
    intensities = intensities * (Debye) ** 2
    print('-------------------------------------')
    print(' Mode    Frequency      IR Intensity')
    print('  #       cm^-1       (D/A)^2 amu^-1')
    print('-------------------------------------')
    for n, e in enumerate(frequencies):
        if e.imag != 0:
            c = 'i'
            e = e.imag
        else:
            c = ' '
        print('%3d   %7.1f%s  %10.6f' % (n, e, c, intensities[n]))
    print('-------------------------------------')
    print('')
    fd = open('ir.dat', 'w')
    fd.write('-------------------------------------\n')
    fd.write(' Mode    Frequency      IR Intensity\n')
    fd.write('  #       cm^-1       (D/A)^2 amu^-1\n')
    fd.write('-------------------------------------\n')
    for n, e in enumerate(frequencies):
        if e.imag != 0:
            c = 'i'
            e = e.imag
        else:
            c = ' '
        fd.write('%3d   %7.1f%s  %10.6f \n' % (n, e, c, intensities[n]))
    fd.close()
    write_spectra(frequencies, intensities, 'ir-spectra.dat', 1000, 4000, 3001, 20, 'Lorentzian')


def get_spectrum(frequencies, intensities, start=1000, end=4000, npts=None, width=20, type='Gaussian'):
    """Get infrared spectrum.
       The method returns wavenumbers in cm^-1 with corresonding absolute infrared intensity.
       Start and end point, and width of the Gaussian/Lorentzian should be given in cm^-1."""

    type = type.lower()
    assert type in ['gaussian', 'lorentzian']
    if not npts:
        npts = (end - start) / width * 10 + 1
    if type == 'lorentzian':
        intensities = intensities * width * pi / 2.
    else:
        sigma = width / 2. / sqrt(2. * log(2.))
    # Make array with spectrum data
    spectrum = np.empty(npts, np.float)
    ediff = (end - start) / float(npts - 1)
    energies = np.arange(start, end + ediff, ediff)
    for i, energy in enumerate(energies):
        energies[i] = energy
        if type == 'lorentzian':
            spectrum[i] = (intensities * 0.5 * width / pi / ((frequencies - energy) ** 2 + 0.25 * width ** 2)).sum()
        else:
            spectrum[i] = (intensities * np.exp(-(frequencies - energy) ** 2 / 2. / sigma ** 2)).sum()
    return [energies, spectrum]


def write_spectra(frequencies, intensities, out='ir-spectra.dat', start=1000, end=4000, npts=None, width=20,
                  type='Gaussian'):
    """Write out infrared spectrum to file.

    First column is the wavenumber in cm^-1, the second column the absolute infrared intensities, and
    the third column the absorbance scaled so that data runs from 1 to 0. Start and end
    point, and width of the Gaussian/Lorentzian should be given in cm^-1."""
    energies, spectrum = get_spectrum(frequencies, intensities, start, end, npts, width, type)

    # Write out spectrum in file. First column is absolute intensities.
    # Second column is absorbance scaled so that data runs from 1 to 0
    # spectrum2 = 1. - spectrum/spectrum.max()
    spectrum2 = spectrum / spectrum.max()
    outdata = np.empty([len(energies), 3])
    outdata.T[0] = energies
    outdata.T[1] = spectrum
    outdata.T[2] = spectrum2
    fd = open(out, 'w')
    for row in outdata:
        fd.write('%.3f  %15.5e  %15.5e \n' % (row[0], row[1], row[2]))
    fd.close()


def extract_U0():
    file = open('OSZICAR', 'r')
    ulines = file.readlines()
    file.close()
    for line in ulines:
        if line.rfind(' 1 F=') > -1:
            U0 = float(line.split()[4])
    return U0


def parse_args():
    optlist, args = getopt(sys.argv[1:], "invsatgp")
    options = {'internal_coords': False,
               'animate': False,
               'kie': False,
               'verbose': False,
               'gaussian_hessian': False,
               'calc_ir_spectrum': False,
               'calc_shomate_parameters': False,
               'calc_thermo_properties': False}

    for opt, arg in optlist:
        if opt == '-i':
            options['internal_coords'] = True
        elif opt == '-v':
            options['verbose'] = True
        elif opt == '-a':
            options['animate'] = True
        elif opt == '-g':
            options['gaussian_hessian'] = True
        elif opt == '-s':
            options['calc_ir_spectrum'] = True
        elif opt == '-p':
            options['calc_shomate_parameters'] = True
            options['calc_thermo_properties'] = True
        elif opt == '-t':
            options['calc_thermo_properties'] = True
        elif opt == '-k':
            options['kie'] = True

    if options['internal_coords'] and options['calc_ir_spectrum']:
        raise NotImplementedError('Cannot use internal coords to calculate ir spectrum!')

    if args:
        print("Usage: kie.py")
        print("-v verbose")
        print('-s IR intensity')
        print('-t isotopic effect')
        print('-i use hessian in internal coordinates (translation/rotation is projected out')
        print('-a animation')
        print('-g gaussian method for filling in hessian')
        print('-p calc shomate parameters')
        sys.exit(1)

    if options['animate'] and options['internal_coords']:
        warnings.warn('Internal coordinates may give animations that do not correspond to the correct vibrational mode. '
                      'Use with caution.')

    return options

def array2string(array):
    freqs_string = ','.join([f'{np.real(f):.1f}'  if np.isreal(f) else f'{np.imag(f):.1f}j' for f in array])
    freqs_string = freqs_string.replace(' ', '')
    return freqs_string

# -------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------
class FrequencyCalculator(ABC):
    def __init__(self, atoms, dof, raw_hessian, freq_cutoff=FREQ_CUTOFF, gaussian_hessian=False):
        atoms.center()
        self.gaussian_hessian = gaussian_hessian
        self.factor = sqrt(_e / (_aa) ** 2 / _amu)
        self.atoms = atoms
        self.freq_cutoff = freq_cutoff
        self.dof = dof
        self.raw_hessian = raw_hessian

    @abstractmethod
    def calc(self):
        raise NotImplementedError

    def _replace_low_freqs_with_cutoff(self, freqs):
        freqsnew = freqs.copy()
        for i, freq in enumerate(freqs):
            if freq <= self.freq_cutoff and freq.imag == 0 and abs(freq) > eps:
                logger.info('***')
                logger.info('WARNING: Frequency %4.3f is below threshhold of %4.3f!' % (freq, self.freq_cutoff))
                logger.info('Replacing frequency by threshold value.')
                logger.info('If this is not a TS calculations, consider relaxing the endpoints further.')
                logger.info('***')
                freqsnew[i] = self.freq_cutoff
            elif abs(freq) <= eps:
                logger.info('Frequency #%d is very close to zero, probably from internal coords freq calculation' % i)

        return freqsnew

    def is_gas_species(self):
        freqs, nonvibs, modes, modetypes = self.calc()
        return nonvibs >= 5

    def is_gas_species_linear(self):
        freqs, nonvibs, modes, modetypes = self.calc()
        return nonvibs == 5

    @property
    def relaxed_atoms(self):
        return self.atoms[self.relaxed_atoms_idx]

    @property
    def relaxed_atoms_idx(self):
        constraints = self.atoms.constraints
        fixed_atoms_constraint = FixAtoms(indices=[])
        for constraint in constraints:
            if isinstance(constraint, FixScaled):
                raise NotImplementedError('Cannot calculate for atoms that are partially fixed!')
            else:
                fixed_atoms_constraint = constraint
        return [i for i, atom in enumerate(self.atoms) if i not in fixed_atoms_constraint.index]

    @property
    def hessian(self):
        H = self.raw_hessian.copy()
        # dynmass: mass repeated in x, y, z for one atom, mass matrix
        dynmass = []
        for i in range(self.dof):
            massIdx = self.relaxed_atoms_idx[int(np.floor(i / 3))]
            dynmass.append(self.masses[massIdx])
        # SYMMETRIZE H -- adding Guowen's edit to symmetrize hessian for vasp5.
        for i in range(self.dof):
            H[:, i] = 0.5 * (H[:, i] + H[i, :])
            H[i, :] = H[:, i]

        for i in range(self.dof):
            H[:, i] = H[:, i] / sqrt(dynmass[i])
            H[i, :] = H[i, :] / sqrt(dynmass[i])

        return -H  # in VASP, H is the first derevative of forces, H should times -1

    @property
    def mask(self):
        # Which coordinates are relaxed (1's) and which aren't (0's)
        num_atoms = len(self.atoms)
        mask = np.zeros([num_atoms, 3])
        for i in range(self.dof):
            idx = self.relaxed_atoms_idx[int(np.floor(i / 3))]
            mask[idx, :] = 1
        return mask.astype('bool')

    @property
    def masses(self):
        masses = []
        for i, atom in enumerate(self.atoms):
            if atom.symbol == 'X':
                warnings.warn(f'Atom {i} with symbol X is treated as deuterium with mass of 2.014')
                masses.append(2.014)
            else:
                masses.append(atom.mass)
        return masses

    @property
    def mass_matrix(self):
        Mg = np.zeros([self.dof, self.dof])
        counter3 = 0
        for m in range(len(self.atoms)):
            for n in range(3):
                if self.mask[m][n] == 1:
                    Mg[counter3][counter3] = 1. / sqrt((float(self.masses[m])))
                    counter3 = counter3 + 1
        return Mg

    def print_results(self):
        freqs, nonvibs, modes_ic, modetypes = self.calc()
        print('----------------------------------------')
        print(' %-5s  %-10s  %15s' % ('Mode#', 'ModeTypes', 'Frequency(1/cm)'))
        print('----------------------------------------')
        for n, freq in enumerate(freqs):
            if freq.imag != 0:
                c = 'i'
                freq = freq.imag
            else:
                c = ' '
                freq = freq.real

            print(' %-5s  %-10s  %15.1f%s ' % (n, modetypes[n], freq, c))
        print('----------------------------------------')
        print('')

        string = ''
        for i in np.real(freqs):
            if i == 0:
                continue
            string = string + ' %3.1f' % i
        print('All non-imaginary freqs: %s' % string)
        print('All freqs: %s' % array2string(freqs))

    @classmethod
    def from_outcar_and_poscar(cls, outcarfile, poscarfile, **kwargs):
        outcar = Outcar.from_file(outcarfile)
        atoms = ase.io.read(poscarfile)
        return cls(dof=outcar.dof, atoms=atoms, raw_hessian=outcar.raw_hessian, **kwargs)

    @classmethod
    def from_folder(cls, folder, **kwargs):
        outcarfile = os.path.join(folder, 'OUTCAR')
        poscarfile = os.path.join(folder, 'POSCAR')
        return cls.from_outcar_and_poscar(outcarfile=outcarfile, poscarfile=poscarfile, **kwargs)


class MassWeightedCoordsFrequencyCalculator(FrequencyCalculator):
    def __init__(self, atoms, dof, raw_hessian, gaussian_hessian=False):
        FrequencyCalculator.__init__(self, atoms=atoms, dof=dof, raw_hessian=raw_hessian,
                                     gaussian_hessian=gaussian_hessian)

    @lru_cache()
    def calc(self):
        ModeTypes = self.dof * ['vib']
        nonvibs = 0
        omega2, modes = np.linalg.eigh(self.hessian)
        modes = modes.T

        frequencies = self.factor * omega2.astype(complex) ** 0.5 / (_c * 100 * pi * 2)
        frequencies = self._replace_low_freqs_with_cutoff(frequencies)
        return frequencies, nonvibs, modes, ModeTypes


class InternalCoordsFrequencyCalculator(FrequencyCalculator):
    def __init__(self, atoms, dof, raw_hessian, gaussian_hessian=False):
        FrequencyCalculator.__init__(self, atoms=atoms, dof=dof, raw_hessian=raw_hessian,
                                     gaussian_hessian=gaussian_hessian)

    @lru_cache()
    def calc(self):
        logger.info("Using internal coordinates instead of mass weighted coordinates.")
        Hmwc = deepcopy(self.hessian)
        H, D, nonvibs, mol_mass = self._convert_mass_weighted_hessian_to_internal_coords()

        # omega2_ic, modes_ic = scipy.linalg.eigh(H)
        omega2_ic, modes_ic = np.linalg.eigh(H)
        logger.info(f'{self.factor * omega2_ic.astype(complex) ** 0.5 / (_c * 100 * pi * 2)}')
        logger.info(f'{modes_ic}')
        evectors_ic = modes_ic.T
        omega2 = np.zeros([self.dof])
        evectors = np.zeros([self.dof, self.dof])

        ModeTypes = ['vib'] * self.dof
        ModeTypes[0] = ModeTypes[1] = 'trans'
        if nonvibs >= 5:
            ModeTypes[2] = 'trans'
            ModeTypes[3] = ModeTypes[4] = 'rot'

        if nonvibs == 6:
            ModeTypes[5] = 'rot'

        for i in range(nonvibs):
            evectors[i][i] = 1

        for i in range(nonvibs, self.dof):
            omega2[i] = omega2_ic[i - nonvibs]
            for j in range(nonvibs, self.dof):
                evectors[i][j] = evectors_ic[i - nonvibs][j - nonvibs]

        if nonvibs == 2:
            # surface species, calculate the translational frequencies separately
            # translational modes along x,y plane are in D1 and D2
            Dtrans = np.zeros([self.dof, 2])
            Dtrans[:, 0] = D[:, 0]
            Dtrans[:, 1] = D[:, 1]
            Htrans_int = np.dot(Dtrans.T, np.dot(Hmwc, Dtrans))
            omega2_trans, modes_trans = np.linalg.eigh(Htrans_int)
            evectors_trans = modes_trans.T

            for i in range(nonvibs):
                omega2[i] = omega2_trans[i]
                for j in range(nonvibs):
                    evectors[i][j] = evectors_trans[i][j]

        # Transfrom the eigenvectors back into cartesian space
        evectors_cart = (np.dot(self.mass_matrix, np.dot(D, evectors.T))).T

        # normalize the evectors
        for v in range(len(evectors_cart)):
            vec = evectors_cart[v]
            l = sqrt(np.dot(vec, vec))
            if abs(l) < eps:
                continue
            evectors_cart[v] = vec / l

        frequencies = self.factor * omega2.astype(complex) ** 0.5 / (_c * 100 * pi * 2)
        frequencies = self._replace_low_freqs_with_cutoff(frequencies)
        return frequencies, nonvibs, evectors_cart, ModeTypes
        #return frequencies, nonvibs, modes_ic.T, ModeTypes

    def _convert_mass_weighted_hessian_to_internal_coords(self):
        atoms = self.atoms.copy()
        H = self.hessian
        # Calculate com (center of mass)
        com = np.zeros([3])
        sum_m = 0.0
        nonvibs = 6
        num_atoms = len(atoms)
        for i in range(len(self.atoms)):
            if self.mask[i][0] == 0:
                # atoms must be fully relaxed or fixed!!!
                nonvibs = 2
                continue
            for j in range(3):
                com[j] += atoms[i].position[j] * self.masses[i]
            sum_m += self.masses[i]

        for j in range(3):
            com[j] = com[j] / sum_m
        com = atoms.get_center_of_mass()
        logger.info(f'Center of Mass: {com}')

        # now we shift the origin to center of mass
        atoms.set_positions(atoms.get_positions() - com)

        # Calculate the moment of inertia tensor
        (principal_I, X) = self._calc_inertia_tensor()

        logger.info(f'The principal moments of inertia are: {principal_I}.')
        logger.info(f'The principal axis of inertia are: {X}')

        # Create rotating and translating frame, only fully relaxed or static atoms work.
        D = np.zeros([self.dof, self.dof])
        P = np.zeros([int(self.dof / 3), 3])
        counter4 = 0
        for i in range(num_atoms):
            mass = self.masses[i]
            if self.mask[i][0] == 0:
                continue

            # Generate the P dot products (P_x,y,z)_i for atom i
            for j in range(3):
                P[counter4, j] = np.dot(mass * atoms[i].position, X[j])

            D[3 * counter4][0] = sqrt(mass)
            D[3 * counter4 + 1][1] = sqrt(mass)
            if not (nonvibs == 2):  # not a surface species
                D[3 * counter4 + 2][2] = sqrt(mass)
                for j in range(3):
                    D[3 * counter4 + j][3] = (P[counter4, 1] * X[j][2] - P[counter4, 2] * X[j][1]) / sqrt(mass)
                    D[3 * counter4 + j][4] = (P[counter4, 2] * X[j][0] - P[counter4, 0] * X[j][2]) / sqrt(mass)
                    D[3 * counter4 + j][5] = (P[counter4, 0] * X[j][1] - P[counter4, 1] * X[j][0]) / sqrt(mass)
            counter4 += 1

        # Normalize the D vectors
        linear = False
        for i in range(nonvibs):
            # take scalar product of each vector with itself
            sc_prod = np.dot(D[:, i], D[:, i])
            if sc_prod <= eps:  # this is not actually a normal mode
                logger.info("\nThis must be a linear molecule. If not, something is screwed up.")
                linear = True
                if self.gaussian_hessian:
                    warnings.warn("Do NOT use the -G option for linear molecules.\n")
            else:
                D[:, i] /= sqrt(sc_prod)
        if linear:
            nonvibs = 5

        # fill D matrix using null space
        Dvib = self._null(D.T)
        for i in range(nonvibs, self.dof):
            D[:, i] = Dvib[:, i - nonvibs]

        logger.info(f'D Transformation Matrix: {D}')

        # Calculate Hessian in internal coordinates Hint
        if self.gaussian_hessian:  # do the selections the Gaussian way using a submatrix
            print("Using the NvibxNvib submatrix approach used in Gaussian.\n")
            Hsubint = np.zeros([self.dof - nonvibs, self.dof - nonvibs])
            Hint = np.dot(D.T, np.dot(H, D))
            for i in range(self.dof - nonvibs):
                for j in range(self.dof - nonvibs):
                    Hsubint[i][j] = Hint[i + nonvibs][j + nonvibs]
        else:
            Hsubint = np.dot(Dvib.T, np.dot(H, Dvib))

        return (Hsubint, D, nonvibs, sum_m)

    def _calc_inertia_tensor(self):
        principal_I, X = self.relaxed_atoms.get_moments_of_inertia(vectors=True)
        return principal_I, X

    @staticmethod
    def _null(A, eps=1e-4):
        u, s, v = np.linalg.svd(A)
        null_mask = (s <= eps)
        null_space = np.compress(null_mask, v, axis=0)
        return null_space.T


class EnthalpyCalculator(ABC):
    def __init__(self, energy, freqs):
        self.freqs = freqs
        self.energy = energy

    def calc(self, T):
        return self.energy \
               + self._calc_vibrational(T) \
               + self._calc_rotational(T) \
               + self._calc_translational(T) \
               + self._calc_pv_term(T)

    def _calc_vibrational(self, T):
        U = 0
        for freq in self.freqs:
            if freq.imag > 0 or abs(freq) < eps:
                continue

            x = h * C * freq / (kb * T)
            U += R * (0.5 * x * T + x * T / (exp(x) - 1))
        U = U.real
        h_vib = U / eV2J  # Return in eV
        logger.info(f'H_vib: {h_vib:4.3f}')
        return h_vib

    @abstractmethod
    def _calc_translational(self, T):
        raise NotImplementedError

    @abstractmethod
    def _calc_rotational(self, T):
        raise NotImplementedError

    @abstractmethod
    def _calc_pv_term(self, T):
        raise NotImplementedError

    @property
    def ZPE(self):
        ZPE = 0
        for freq in self.freqs:
            if freq.imag > 0 or abs(freq) < eps:
                continue

            ZPE += 0.5 * h * _c * 100 * freq.real
        return ZPE / _e  # Return in eV

    @staticmethod
    def from_freq_calculator_and_energy(freq_calculator, energy):
        freqs, _, _, _ = freq_calculator.calc()
        if freq_calculator.is_gas_species():
            is_species_linear = freq_calculator.is_gas_species_linear()
            return GasSpeciesEnthalpyCalculator(energy=energy, freqs=freqs, is_species_linear=is_species_linear)
        else:
            return SurfSpeciesEnthalpyCalculator(energy=energy, freqs=freqs)


class GasSpeciesEnthalpyCalculator(EnthalpyCalculator):
    def __init__(self, energy, freqs, is_species_linear):
        super().__init__(freqs=freqs, energy=energy)
        self.is_species_linear = is_species_linear

    def _calc_translational(self, T):
        h_trans = (3 / 2 * R * T) / eV2J
        logger.info(f'H_trans: {h_trans:4.3f}')
        return h_trans

    def _calc_rotational(self, T):
        if self.is_species_linear:
            # Cp = R
            h_rot = (R * T) / eV2J
        else:
            # Cp = 3/2R
            h_rot = (3 / 2 * R * T) / eV2J
        logger.info(f'H_rot: {h_rot:4.3f}')
        return h_rot

    def _calc_pv_term(self, T):
        h_pv = (R * T) / eV2J  # ideal gas, pv=nRT
        logger.info(f'H_pv: {h_pv:4.3f}')
        return h_pv


class SurfSpeciesEnthalpyCalculator(EnthalpyCalculator):
    def __init__(self, energy, freqs):
        super().__init__(freqs=freqs, energy=energy)

    def _calc_translational(self, T):
        return 0

    def _calc_rotational(self, T):
        return 0

    def _calc_pv_term(self, T):
        return 0


class EntropyCalculator(ABC):
    RRHO_WEIGHTING_CONSTANT = 4

    def __init__(self, freqs, atoms, modetypes, rrhocutoff):
        # Atoms are only unfixed atoms
        self.rrhocutoff = rrhocutoff
        self.freqs = freqs
        self.atoms = atoms
        self.modetypes = modetypes

    @property
    def mass(self):
        return sum([atom.mass for atom in self.atoms])

    def calc(self, T, P):
        Strans = self._calc_translational(T, P)
        Srot = self._calc_rotational(T)
        Svib = self._calc_vibrational(T)
        logger.info('Strans (%3.2f): %3.2f' % (T, Strans))
        logger.info('Srot (%3.2f): %3.2f' % (T, Srot))
        logger.info('Svib (%3.2f): %3.2f' % (T, Svib))
        return Strans + Srot + Svib

    def _calc_vibrational(self, T):
        if self.rrhocutoff == 0:
            return self._calc_vibrational_harmonic_oscillator(T, weights=None)
        else:
            freqs = np.array([f for f in self.freqs.real if f > 0])
            weights = 1 / (1 + (self.rrhocutoff / freqs) ** self.RRHO_WEIGHTING_CONSTANT)
            harmonic_oscillator = self._calc_vibrational_harmonic_oscillator(T, weights=weights)
            rotational_osciallator = self._calc_vibrational_rrho(T, weights=1 - weights)

            return harmonic_oscillator + rotational_osciallator

    def _calc_vibrational_harmonic_oscillator(self, T, weights=None):
        freqs = np.array([f for f in self.freqs.real if f > 0])
        x = h * C * freqs / (kb * T)
        S = R * (x / (np.exp(x) - 1) - np.log(1 - np.exp(-x)))
        if weights is not None:
            S = S * weights
        return np.sum(S)

    def _calc_vibrational_rrho(self, T, weights=None):
        freqs = np.array([f for f in self.freqs.real if f > 0])
        v = freqs * C
        u = h / (8 * pi ** 2 * v)

        average_u = 1E-44  # constant by Grimme
        uprime = (u * average_u) / (u + average_u)
        logger.info('uprime: {uprime}')
        S = R * (0.5 + 0.5 * np.log(8 * pi ** 3 * uprime * kb * T / h ** 2))
        if weights is not None:
            S = S * weights
        return np.sum(S)

    @abstractmethod
    def _calc_translational(self, T, P):
        raise NotImplementedError

    @abstractmethod
    def _calc_rotational(self, T):
        raise NotImplementedError

    @staticmethod
    def from_freq_calculator(freq_calculator, **kwargs):
        freqs, nonvibs, modes, modetypes = freq_calculator.calc()
        factory = EntropyCalculatorFactory(freqs=freqs, nonvibs=nonvibs, modetypes=modetypes,
                                           atoms=freq_calculator.atoms,
                                           is_gas_species=freq_calculator.is_gas_species(),
                                           **kwargs)
        return factory.build()


class GasSpeciesEntropyCalculator(EntropyCalculator):
    def __init__(self, freqs, atoms, modetypes, sno, rrhocutoff=0):
        EntropyCalculator.__init__(self, freqs, atoms, modetypes, rrhocutoff)
        self.sno = sno

    @property
    def _is_species_linear(self):
        principal_I, X = self._get_inertia_tensor()
        principal_I = [p / 1000 / NA * 1E-20 for p in principal_I if p > eps]
        if len(principal_I) == 2:
            is_linear = True
        elif len(principal_I) == 3:
            is_linear = False
        else:
            raise NotImplementedError('Error in determining linearity of molecule. Too many small inertia moments!')

        return is_linear

    def _calc_translational(self, T, P):
        mass = self.mass / 1000 / NA
        return R * ((log(
            (2 * pi * mass * kb * T / (h ** 2)) ** 1.5 * kb * T / P)) + 2.5)  # 3D Translational contribution

    def _calc_rotational(self, T):
        # T is in K, P is in Pa
        principal_I, X = self._get_inertia_tensor()
        principal_I = [p / 1000 / NA * 1E-20 for p in principal_I if p > eps]
        ProdI = np.prod(principal_I)

        if self._is_species_linear:
            Imax = max(principal_I)
            S = R * (log(8 * pi ** 2 * Imax * kb * T / (self.sno * h ** 2)) + 1)
            logger.info("Using rotational entropy for linear molecule with I " + str(Imax) + ".\n")
        else:
            S = R * (log(8 * pi ** 2 * sqrt(8 * pi ** 3 * ProdI) * (kb * T) ** 1.5 / (self.sno * h ** 3)) + 1.5)
            logger.info("Using rotational entropy for nonlinear molecule.")

        # Cp = (1 + nonvibs / 2.0) * R
        return S

    def _get_inertia_tensor(self):
        principal_I, X = self.atoms.get_moments_of_inertia(vectors=True)
        return principal_I, X


class SurfSpeciesEntropyCalculator(EntropyCalculator):
    def __init__(self, freqs, atoms, modetypes, rrhocutoff=0):
        EntropyCalculator.__init__(self, freqs, atoms, modetypes, rrhocutoff)

    def _calc_translational(self, T, P):
        return 0

    def _calc_rotational(self, T):
        return 0


class EntropyCalculatorFactory(object):
    def __init__(self, freqs, modetypes, nonvibs, atoms, is_gas_species, rrhocutoff=0):
        self.freqsnew = freqs
        self.modetypesnew = modetypes
        self.nonvibs = nonvibs
        self.atoms = atoms
        self.is_gas_species = is_gas_species
        self.rrhocutoff = rrhocutoff

    def build(self):
        if self.is_gas_species:
            sno = self._get_symmetry_number()
            entropy_calculator = GasSpeciesEntropyCalculator(atoms=self.atoms,
                                                             freqs=self.freqsnew,
                                                             sno=sno,
                                                             modetypes=self.modetypesnew,
                                                             rrhocutoff=self.rrhocutoff)
        else:
            entropy_calculator = SurfSpeciesEntropyCalculator(atoms=self.atoms,
                                                              freqs=self.freqsnew,
                                                              modetypes=self.modetypesnew,
                                                              rrhocutoff=self.rrhocutoff)

        return entropy_calculator

    def _get_symmetry_number(self):
        try:
            sno = self._get_symmetry_number_via_symmetry_detection()
        except:
            sno = self._get_symmetry_no_via_user_input()
        return sno

    def _get_symmetry_number_via_symmetry_detection(self):
        print(f'\n-----------------------Symmetry Detection-----------------------')
        atoms = self.atoms
        dist_vectors = atoms.get_distances(0, range(len(atoms)), mic=True, vector=True)

        pos = atoms.get_positions()
        pos = pos[0] + dist_vectors

        atoms.set_positions(pos - np.mean(pos, 0) + atoms.cell @ np.array([0.5, 0.5, 0.5]))
        structure = AseAtomsAdaptor.get_structure(atoms)
        m = Molecule.from_str(str(XYZ(structure)), fmt='xyz')
        a = PointGroupAnalyzer(m, tolerance=0.3)

        symbol = a.sch_symbol
        if a.rot_sym:
            _, sno = max(a.rot_sym, key=lambda v: v[1])
        else:
            if symbol in ['D*h']:
                sno = 2
            elif symbol in ['C*v', 'Cs', 'C1']:
                sno = 1
            else:
                raise ValueError(f'Symmetry number for {symbol} unknown.')
        print(f'Point group     :  {symbol}')
        print(f'Symmetry number :  {sno}')
        print(f'----------------------------------------------------------------')
        return sno

    def _get_symmetry_no_via_user_input(self):
        print("The rotational symmetry number s is the number of indistinguishable orientations of the molecule"
              "generated by rotation. If the molecule has no center of symmetry (e.g. HCN) s = 1, whereas if the"
              "molecule has a center of symmetry (e.g. CO2) then s = 2.\n")
        sno = int(get_input("Symmetry Number for this molecule [1] : ", 1))
        return sno


class ShomateCalculator(ABC):
    def __init__(self, entropy_calculator, enthalpy_calculator):
        self.entropy_calculator = entropy_calculator
        self.enthalpy_calculator = enthalpy_calculator

    def calc(self, temp_range, P):
        params = self._calc_shomate(temp_range, P)
        return params

    def _calc_shomate(self, temp_range, P):
        ndata = 100
        dT = (temp_range[-1] - temp_range[0]) / ndata
        T = temp_range[0]
        S_data = []
        T_data = []
        iniparam = (1, 1, 1, 1, 1, 1)
        for i in range(ndata + 1):
            S_data.append(self.entropy_calculator.calc(T, P))
            T_data.append(T)
            T += dT

        T_data = np.array(T_data)
        output = optimize.least_squares(self._Shomate_S_fit, iniparam,
                                        args=(T_data, S_data), xtol=1E-8, gtol=1E-8,
                                        jac='3-point', max_nfev=500)
        param = output.x

        H298 = self.enthalpy_calculator.calc(298) * eV2kJ
        shomate = [param[0], param[1], param[2], param[3], param[4], 0, param[5], H298]
        # t0 = Trange[0] / 1000.0
        t0 = 298 / 1000

        shomate[5] = H298 \
                     - param[0] * t0 \
                     - param[1] * t0 ** 2 / 2 \
                     - param[2] * t0 ** 3 / 3 \
                     - param[3] * t0 ** 4 / 4 \
                     + param[4] / t0

        return [s.real for s in shomate]

    def print_results(self, Trange, P):
        params = self.calc(Trange, P)
        print('')
        print('------Least squares fitting of Shomate equation for Entropy-----')
        print('Using DFT Energies as Enthalpy reference')
        print('Shomate Parameters A, B, C, D, E, F, G, H, F for T range %s' % Trange)
        for i, pname in enumerate(('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')):
            print("%s:  %15.7e" % (pname, params[i]))
        print('All shomate: %7.7e %7.7e %7.7e %7.7e %7.7e %7.7e %7.7e %7.7e' % tuple(params))
        print('-----------------------------------------------------------------')
        print('')

    @staticmethod
    def _Shomate_H(param, T):
        t = T / 1000.0
        H = param[0] * t + param[1] * t ** 2 / 2 + param[2] * t ** 3 / 3 + param[3] * t ** 4 / 4 - param[4] / t + param[
            5]
        return H

    @staticmethod
    def _Shomate_S_fit(param, T, S_data):
        t = T / 1000.0
        S = param[0] * np.log(t) + param[1] * t + param[2] * t * t / 2 + param[3] * t * t * t / 3 - param[4] / (
                2 * t * t) + param[5]
        Sdiff = S - S_data
        return Sdiff


class ThermoPropertiesPrinter(object):
    def __init__(self, atoms, enthalpy_calculator, entropy_calculator):
        self.atoms = atoms
        self.enthalpy_calculator = enthalpy_calculator
        self.entropy_calculator = entropy_calculator

    def print(self, T, P):
        enthalpy = self.enthalpy_calculator.calc(T)
        entropy = self.entropy_calculator.calc(T, P)

        print('')
        print('---------------------Thermodynamic Properties-------------------')
        print("Molecular Weight                 : %10.3f g" % (sum([atom.mass for atom in self.atoms])))
        print("Electronic Energy                : %10.3f eV" % self.enthalpy_calculator.energy)
        print("Zero-point Energy                : %10.3f eV" % self.enthalpy_calculator.ZPE)
        print("Enthalpy                         : %10.3f eV at %.2f Kelvin" % (enthalpy, T))
        print("Entropy                          : %10.3f J/mol/K at %.2f Kelvin, %.3f Pa" % (entropy, T, P))
        print('----------------------------------------------------------------')
        print('')


class Outcar(object):
    def __init__(self, dof, raw_hessian):
        self.dof = dof
        self.raw_hessian = raw_hessian

    @classmethod
    def from_file(cls, outcar):
        def _get_degrees_of_freedom(lines):
            for line in lines:
                if line.rfind('Degrees of freedom') > -1:
                    dof = int(line.split()[5])
            assert dof % 3 == 0
            return dof

        def _get_raw_hessian(lines, dof):
            n = 0
            hes = []
            for line in lines:
                if line.find('SECOND DERIVATIVES (NOT SYMMETRIZED)') > -1:
                    for i in range(dof):
                        hes.append(np.array([float(f) for f in lines[n + 3 + i].split()[1:dof + 1]]))
                n += 1

            if not hes:
                raise ValueError('Could not read Hessian from OUTCAR.')

            return np.array(hes)
           
        with zopen(outcar, 'rt') as f:
            lines = f.readlines()
        
        dof = _get_degrees_of_freedom(lines)
        raw_hessian = _get_raw_hessian(lines, dof)
        return cls(dof=dof, raw_hessian=raw_hessian)

@contextmanager
def zopen(file, mode):
    if os.path.exists(file):
        f = open(file, mode)
        yield f
    elif os.path.exists(file + '.gz'):
        from monty.io import zopen as monty_zopen
        f = monty_zopen(file + '.gz', mode)
        yield f
    else:
        raise FileNotFoundError('Could not find file {file}.')

    f.close()
##############################################
if __name__ == '__main__':
    options = parse_args()
    logger = logging.getLogger()
    logger.handlers = []
    if options['verbose']:
        formatter = logging.Formatter('%(asctime)s - %(levelname)-6s - [%(processName)s %(name)s] %(message)s')
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)
        logger.setLevel(logging.DEBUG)
    main(options)
