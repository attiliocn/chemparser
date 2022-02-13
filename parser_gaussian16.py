#!/usr/bin/env python
import re
import os

from .tools import *
from .parser_nbo7 import NaturalBondOrbital7

class GaussianOutput(NaturalBondOrbital7):
    def __init__(self,output_file):
        self.output_file = output_file

    def parse(self):
        with open(self.output_file) as output:
            if not re.search('Entering Gaussian System', output.readline()):
                print('Not Gaussian Output!') #TODO: This should return error, not a print.
            else:
                gaussian_output_data = {
                    'number_of_atoms': None,
                    'scf_energies': list(),
                    'dipole': None,
                    'isotropic_polarizability': None,
                    'anisotropic_polarizability': None,
                    'occupied_orbitals_eigenvalues': None,
                    'virtual_orbitals_eigenvalues': None,
                }
                for line in output:
                    # Set number_of_atoms attribute from first match on the output
                    if not hasattr(self, 'number_of_atoms'):
                        if re.search('NAtoms=', line):
                            number_of_atoms = int(line.split()[1])
                            gaussian_output_data['number_of_atoms'] = number_of_atoms
                    # Obtain SCF Energies
                    if re.search('SCF Done', line):
                        scf_energy = float(line.split()[4])
                        gaussian_output_data['scf_energies'].append(scf_energy)
                    # Obtain the molecular dipole
                    if re.search('Electric dipole moment \(input orientation\)', line):
                        for _ in range(2):
                            output.readline()
                        dipole = float(output.readline().split()[1].replace('D','E'))
                        gaussian_output_data['dipole'] = dipole
                    # Obtain the polarizability
                    if re.search('Dipole polarizability, Alpha \(input orientation\)', line):
                        for _ in range(3):
                            output.readline()
                        isotropic_polarizability = float(output.readline().split()[1].replace('D','E'))
                        anisotropic_polarizability = float(output.readline().split()[1].replace('D','E'))
                        gaussian_output_data['isotropic_polarizability'] = isotropic_polarizability
                        gaussian_output_data['anisotropic_polarizability'] = anisotropic_polarizability
                    # Obtain the orbital energies from SCF Density
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
                        gaussian_output_data['occupied_orbitals_eigenvalues'] = occupied_orbitals_eigenvalues
                        gaussian_output_data['virtual_orbitals_eigenvalues'] = empty_orbitals_eigenvalues
                    # Parse NBO7 output
                    if re.search(' NBO 7.0 ', line):
                        with open('nbo_temp.out', mode='w') as nbo_output:
                            nbo_output.write(line)
                            while True:
                                current_line = output.readline()
                                if re.search('NBO analysis completed', current_line):
                                    break
                                else:
                                    nbo_output.write(current_line)
                        nbo_object = NaturalBondOrbital7('nbo_temp.out')
                        nbo_object.check_this_out()
                        os.remove('./nbo_temp.out')

    def get_geometries(self):
        pass

    def write_geometry(self):
        pass

    def print_filename(self):
        print(self.output_file)
