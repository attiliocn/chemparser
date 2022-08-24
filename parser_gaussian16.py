#!/usr/bin/env python
import re
import os

import numpy as np

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

    def get_hirshfeld_charges(self):
        '''Fetches the charges from Hirshfeld population analysis from Gaussian16 output file. 

        Args:
            None

        Returns:
            Dict containing charges for all atoms (currently not available) and charges with 
            hydrogens summed into heavy atoms. The dict_keys are 0-indexed atom numbers 

        Raises:
            PropertyNotFoundError: When population analysis using Hirshfeld is not present on the output file. 
        '''
        hirshfeld_charges = {
            'all_atoms': {},
            'without_H': {}
        }
        with open(self.output_file) as output:
            for line in output:
                if re.search('Hirshfeld charges with hydrogens summed into heavy atoms:', line):
                    output.readline()
                    while True:
                        hirshfeld_output = output.readline()
                        if re.search("\s+[0-9]+\s+[A-Za-z]+\s+-?[0-9]+\.[0-9]+", hirshfeld_output):
                            hirshfeld_output_split = hirshfeld_output.split()
                            atom_number = int(hirshfeld_output_split[0])
                            hirshfeld_charges['without_H'][atom_number-1] = {
                                'atom_number': atom_number,
                                'element': hirshfeld_output_split[1],
                                'hirshfeld_charge': float(hirshfeld_output_split[2]),
                                'cm5_charge': float(hirshfeld_output_split[3])
                            }
                        else:
                            break
        if hirshfeld_charges['without_H']:
            return hirshfeld_charges
        else: 
            raise PropertyNotFoundError("Output does not contain Hirshfeld charges")
    
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
                        regex_energy = re.compile('-?[0-9]+.[0-9]+')
                        while True:
                            current_line = output.readline()
                            if re.search('Alpha  occ. eigenvalues', current_line):
                                energies_from_line = regex_energy.findall(current_line)
                                for energy in energies_from_line:
                                    occupied_orbitals_eigenvalues.append(float(energy))
                            elif re.search('Alpha virt. eigenvalues', current_line):
                                energies_from_line = regex_energy.findall(current_line)
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

    def get_nmr_tensors(self):
        '''Fetches NMR Magnetic shielding tensors from Gaussian Output. 

        Args:
            None

        Returns:
            Dict of NMR Shielding Tensors

        Raises:
            PropertyNotFoundError: When NMR calculation is not present in Gaussian16 output file. 
        '''
        regex_tensors = re.compile('\s*([A-Z]{2})=\s*(-?[0-9]+.[0-9]+)')
        tensors = dict()
        with open(self.output_file) as output:
            for line in output:
                if re.search('SCF GIAO Magnetic shielding tensor', line):
                    while True:
                        current_line = output.readline()
                        if re.search('\s+[0-9]+\s+[A-Za-z]{1,2}\s+Isotropic', current_line):
                            atom_number = int(current_line.split()[0])
                            atom_tensors = dict()
                            for _ in range(4):
                                current_line = output.readline()
                                if not re.search('Eigenvalues', current_line):
                                    all_tensors = regex_tensors.findall(current_line)
                                    for tensor in all_tensors:
                                        atom_tensors[tensor[0]] = float(tensor[1])
                                    tensors[atom_number] = atom_tensors
                        else:
                            break
        if tensors:
            return tensors
        else:
            raise PropertyNotFoundError("Output does not contain NMR Shielding Tensors")

                            

    def get_thermochemistry(self) -> dict:
        pass

    def get_geometries(self):
        '''Fetches al XYZ Coordinates from Gaussian Output. 

        Args:
            None

        Returns:
            List of numpy arrays containing the molecular geometries

        Raises:
            PropertyNotFoundError: When no geometry is found in the output file 
        '''
        geometries = list()
        with open(self.output_file) as output:
            for line in output:
                if re.search('Coordinates \(Angstroms\)', line):
                    current_geometry = list()
                    output.readline()
                    output.readline()
                    while True:
                        current_line = output.readline()
                        if re.search('\s-{69}', current_line):
                            break
                        else:
                            current_geometry.append(current_line.strip().split())
                    current_geometry = np.array(current_geometry)[:,[1,3,4,5]] # Select only the columns: atomic number and x, y and z coordinates
                    current_geometry[:,0] = [atom_from_atomic_number(int(i)) for i in current_geometry[:,0]] # Replace atomic number for element
                    geometries.append(current_geometry)
        if geometries:
            return geometries
        else:
            raise PropertyNotFoundError("Output does not contain any geometry information")



    def write_geometry(self):
        pass

    def _standard_method(self):
        with open(self.output_file) as output:
            for line in output:
                pass
