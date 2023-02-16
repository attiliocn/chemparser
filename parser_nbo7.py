import os
import re

from .exceptions import *

class NaturalBondOrbital7():
    def __init__(self,output_file):
        self.output_file = output_file
    
    def get_natural_population_analysis(self) -> list:
        '''Fetches natural population analysis table from NBO7 output 

        Args:
            None

        Returns:
            List of natural population analysis for each atom in the molecule. This
            output is 0-indexed.  

        Raises:
            PropertyNotFoundError: When NPA analysis is not present in the output. 
        '''
        with open(self.output_file) as output:
            natural_population_analysis = list()
            for line in output:
                if re.search('Summary of Natural Population Analysis', line):
                    regex_npa = re.compile('([A-Z][a-z]?)\s{0,2}([0-9]{1,3})\s*(-?[0-9]*\.[0-9]*)\s*(-?[0-9]*\.[0-9]*)\s*(-?[0-9]*\.[0-9]*)\s*(-?[0-9]*\.[0-9]*)\s*(-?[0-9]*\.[0-9]*)')
                    for _ in range(5):
                        output.readline()
                    while True:
                        current_line = output.readline()
                        if re.search('={68}', current_line):
                            break
                        else:
                            current_atom_npa = regex_npa.findall(current_line)[0]
                            atom_natural_population = {
                                'atom': str(current_atom_npa[0]),
                                'atom_number': int(current_atom_npa[1]),
                                'natural_charge': float(current_atom_npa[2]),
                                'core_population': float(current_atom_npa[3]),
                                'valence_population': float(current_atom_npa[4]),
                                'rydberg_population': float(current_atom_npa[5]),
                                'total_population': float(current_atom_npa[6])
                            }
                            natural_population_analysis.append(atom_natural_population)
                            # TODO: Make shure that this list is populated on atom-order. Current implementation
                            # relies on NBO output table order, which I assume to be already sorted from 1 to n atoms
        if natural_population_analysis:
            return natural_population_analysis
        else:
            raise PropertyNotFoundError("Output does not contain NPA analysis")

    def get_natural_bond_orbitals(self) -> list:
        '''Fetches natural orbitals from NBO7 output. 

        Args:
            None

        Returns:
            0-indexed list containing dicts of parsed natural bond orbitals.  

        Raises:
            PropertyNotFoundError: When NBO orbitals are not present in the output. 
        '''
        with open(self.output_file) as output:
            natural_bond_orbitals_raw = list()
            for line in output:
                if re.search('NATURAL BOND ORBITALS \(Summary\):', line):
                    for _ in range(6):
                        output.readline()
                    while True:
                        current_line = output.readline()
                        if re.search('NBO analysis completed', current_line):# re.search('^ +-{31}', current_line):
                            break
                        else:
                            natural_bond_orbitals_raw.append(current_line.strip())
            
        if not natural_bond_orbitals_raw:
            raise PropertyNotFoundError("Output does not contain NBO Orbitals Summary")
        
        # regular expressions 
        regex_nbo_parser = re.compile('^([0-9]+)\. +(CR|LP|BD|BD\*|RY) +([0-9]+) +([a-zA-Z]+.*[0-9]+).*(-?[0-9]+\.[0-9]+) +(-?[0-9]+\.[0-9]+)')
        #regex_nbo_identifier = re.compile('^[0-9]+\. ')
        #regex_nbo_float = re.compile('-?[0-9]+\.[0-9]+')
        #regex_nbo_participants = re.compile('[A-Za-z]+\s*([0-9]+)')
        regex_delocalizations = re.compile('[0-9]+\([a-z]\)')

        all_nbo_parsed = list()               
        for nbo_output in natural_bond_orbitals_raw:
            nbo_string = re.sub('[()-]', ' ', nbo_output)
            nbo_string_parsed = regex_nbo_parser.findall(nbo_string)

            if nbo_string_parsed:
                nbo_string_parsed = nbo_string_parsed[0]
                nbo_parsed = {
                    'nbo_number': int(nbo_string_parsed[0]),
                    'nbo_type': nbo_string_parsed[1],
                    'nbo_bond_order': int(nbo_string_parsed[2]),
                    'nbo_occupancy': float(nbo_string_parsed[4]),
                    'nbo_energy': float(nbo_string_parsed[5]),
                    'nbo_participants': [int(i) for i in re.sub('[A-Za-z]','',nbo_string_parsed[3]).split()],
                    'nbo_delocalizations': regex_delocalizations.findall(nbo_output)
                }
                all_nbo_parsed.append(nbo_parsed)
            else:
                nbo_delocalizations = regex_delocalizations.findall(nbo_output)
                for delocalization in nbo_delocalizations:
                    all_nbo_parsed[-1]['nbo_delocalizations'].append(delocalization)
        return all_nbo_parsed

    def get_perturbation_analysis(self):
        '''
        '''
        with open(self.output_file) as output:
            perturbation_raw = list()
            for line in output:
                if re.search('SECOND ORDER PERTURBATION THEORY', line):
                    for _ in range(7):
                        output.readline()
                    while True:
                        perturbation_output = output.readline().strip()
                        if re.search('NATURAL BOND ORBITALS', perturbation_output):
                            break
                        elif re.search('^[0-9]+\. ', perturbation_output):
                            #print(perturbation_output)
                            perturbation_raw.append(perturbation_output)
        
        perturbation_regex = re.compile('([0-9]+\.) (.*) ([0-9]+\.) (.*) ([0-9]+\.[0-9]+) +([0-9]+\.[0-9]+) +([0-9]+\.[0-9]+)')
        perturbations_parsed = list()
        for entry in perturbation_raw:
            #print(entry)
            perturbation_data = perturbation_regex.findall(entry)[0]
            perturbation_data = list(perturbation_data)
            perturbation_data = [i.strip() for i in perturbation_data]
            #print(perturbation_data)

            donor_orbital_number = int(perturbation_data[0].replace('.',''))
            donor_orbital_raw = perturbation_data[1]
            #print(donor_orbital_raw)
            donor_orbital = self.parse_nbo_participants(self,donor_orbital_raw)

            acceptor_orbital_number = int(perturbation_data[2].replace('.',''))
            acceptor_orbital_raw = perturbation_data[3]
            acceptor_orbital = self.parse_nbo_participants(self,acceptor_orbital_raw)

            energy = float(perturbation_data[4])
            E_term = float(perturbation_data[5])
            F_term = float(perturbation_data[6])

            perturbation_parsed = {
                'donor orbital': donor_orbital_number,
                'donor type': donor_orbital['nbo_type'],
                'donor bond order': donor_orbital['nbo_bond_order'],
                'donor participants': donor_orbital['nbo_participants'],
                'acceptor orbital': acceptor_orbital_number,
                'acceptor type': acceptor_orbital['nbo_type'],
                'acceptor bond order': acceptor_orbital['nbo_bond_order'],
                'acceptor participants': acceptor_orbital['nbo_participants'],
                'energy': energy,
                'E_term': E_term,
                'F_term': F_term,
            }
            perturbations_parsed.append(perturbation_parsed)
        return perturbations_parsed



    @staticmethod
    def parse_nbo_participants(self, nbo_string):
        # regular expressions 
        regex_nbo_float = re.compile('-?[0-9]+\.[0-9]+')
        regex_nbo_participants = re.compile('[A-Za-z]+\s*([0-9]+)')
        regex_delocalizations = re.compile('[0-9]+\([a-z]\)')

        split_nbo_content = re.sub('[()]', ' ', nbo_string).split()
        #nbo_number = split_nbo_content[0].replace('.', '')
        nbo_type = split_nbo_content[0]
        nbo_bond_order = split_nbo_content[1]
        #nbo_occupancy, nbo_energy = regex_nbo_float.findall(nbo_string)
        nbo_participants = regex_nbo_participants.findall(nbo_string)
        nbo_delocalizations = regex_delocalizations.findall(nbo_string)

        nbo_parsed = {
            #'nbo_number': int(nbo_number),
            'nbo_type': str(nbo_type),
            'nbo_bond_order': int(nbo_bond_order),
            #'nbo_occupancy': float(nbo_occupancy),
            #'nbo_energy': float(nbo_energy),
            'nbo_participants': [int(i) for i in nbo_participants],
            'nbo_delocalizations': nbo_delocalizations
            }
        #print(nbo_parsed)
        return nbo_parsed
                    
