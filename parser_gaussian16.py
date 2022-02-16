#!/usr/bin/env python
import re
import os

from .exceptions import *
from .tools import *
from .parser_nbo7 import NaturalBondOrbital7
from .tools.tools import atom_from_atomic_number

class GaussianOutput(NaturalBondOrbital7):
    def __init__(self,output_file):
        self.output_file = output_file
    
    def get_number_of_atoms(self) -> int:
        '''Fetches number of atoms from Gaussian16 output file. 

        Args:
            None

        Returns:
            The number of atoms in the molecule as int. 

        Raises:
            PropertyNotFoundError: When 'NAtoms=' is not found on the output file 
        '''
        with open(self.output_file) as output:
            for line in output:
                if re.search('NAtoms=', line):
                    number_of_atoms = int(line.split()[1])
                    return number_of_atoms
        raise PropertyNotFoundError("Output does not contain number of atoms")
    
    def get_scf_energies(self) -> list:
        '''Fetches SCF Energies from Gaussian16 output file. 

        Args:
            None

        Returns:
            All SCF Energies from the output in a list. 

        Raises:
            PropertyNotFoundError: When 'SCF Done' is not found on the output file 
        '''
        scf_energies = list()
        with open(self.output_file) as output:
            for line in output:
                if re.search('SCF Done', line):
                    scf_energy = float(line.split()[4])
                    scf_energies.append(scf_energy)
        if scf_energies:
            return scf_energies
        else:
            raise PropertyNotFoundError("Output does not contain SCF Energy")

    def get_dipole(self) -> float:
        '''Fetches the dipole magnitude from Gaussian16 output file. 

        Args:
            None

        Returns:
            Dipole vector norm in Debye (float) 

        Raises:
            PropertyNotFoundError: When 'Electric dipole moment (input orientation)' is not found on the output file 
        '''
        with open(self.output_file) as output:
            for line in output:
                if re.search('Electric dipole moment \(input orientation\)', line):
                        for _ in range(2):
                            output.readline()
                        dipole = float(output.readline().split()[1].replace('D','E'))
                        return dipole
        raise PropertyNotFoundError("Output does not contain dipole information")

    def get_polarizability(self) -> tuple:
        '''Fetches the dipole polarizability from Gaussian16 output file. 

        Args:
            None

        Returns:
            Tuple containing (isotropic_pol, anisotropic_pol) in Debye. 

        Raises:
            PropertyNotFoundError: When 'Dipole polarizability, Alpha (input orientation)' is not found on the output file. 
        '''
        with open(self.output_file) as output:
            for line in output:
                if re.search('Dipole polarizability, Alpha \(input orientation\)', line):
                        for _ in range(3):
                            output.readline()
                        isotropic_polarizability = float(output.readline().split()[1].replace('D','E'))
                        anisotropic_polarizability = float(output.readline().split()[1].replace('D','E'))
                        return(isotropic_polarizability, anisotropic_polarizability)
        raise PropertyNotFoundError("Output does not contain dipole polarizability information")

    def get_orbitals_energies(self) -> tuple:
        '''Fetches the molecular orbitals energies from Gaussian16 output file. 

        Args:
            None

        Returns:
            Tuple containing ([occupied_orbs], [empty_orbs]) energies in Eh. The lists are 0-indexed  

        Raises:
            PropertyNotFoundError: When population analysis using SCF Density is not present on the output file. 
        '''
        with open(self.output_file) as output:
            for line in output:
                if re.search('Population analysis using the SCF Density', line):
                        for _ in range(3):
                            output.readline()
                        occupied_orbitals_eigenvalues = list()
                        empty_orbitals_eigenvalues = list()
                        while True:
                            current_line = output.readline()
                            if re.search('Alpha  occ. eigenvalues', current_line):
                                energies_from_line = current_line.split()[4:]
                                for energy in energies_from_line:
                                    occupied_orbitals_eigenvalues.append(float(energy))
                            elif re.search('Alpha virt. eigenvalues', current_line):
                                energies_from_line = current_line.split()[4:]
                                for energy in energies_from_line:
                                    empty_orbitals_eigenvalues.append(float(energy))
                            else:
                                break
                        return(occupied_orbitals_eigenvalues, empty_orbitals_eigenvalues)
        raise PropertyNotFoundError("Output does not contain orbitals energies from SCF Density")

    def extract_nbo7_output(self) -> list:
        '''Fetches data from NBO7 output within the Gaussian16 output. 

        Args:
            None

        Returns:
            List of all lines of the NBO7 output. 

        Raises:
            PropertyNotFoundError: When NBO analysis is not present in Gaussian16 output file. 
        '''
        nbo_output = list()
        with open(self.output_file) as output:
            for line in output:
                if re.search(' NBO 7.0 ', line):
                    nbo_output.append(line)
                    while True:
                        current_line = output.readline()
                        if re.search('NBO analysis completed', current_line):
                            return nbo_output
                        else:
                            nbo_output.append(current_line)
            raise PropertyNotFoundError("Output does not contain NBO7 output")

    def get_thermochemistry(self) -> dict:
        pass

    def get_geometries(self):
        pass

    def write_geometry(self):
        pass

    def _standard_method(self):
        with open(self.output_file) as output:
            for line in output:
                pass
