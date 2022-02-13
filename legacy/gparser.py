import re

## utils module ##
##################
def search_lines(lines_list, pattern):
    results = list()
    for i,content in enumerate(lines_list):
        if re.search(pattern,content):
            results.append([i,content])
    return results

def parse_energies_list(energies_list, pattern):
    for i in range(len(energies_list)):
        energies_list[i][1] = float(re.search(pattern,energies_list[i][1]).group(0))
    energies = [energies_list[i][1] for i in range(len(energies_list))]
    return energies

PERIODIC_TABLE = ["Bq","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo","X"]
 

#######################################################################################

class GaussianOutput():
    def __init__(self, file_path):
        
        self.file_path = file_path
        self.name = file_path.split('/')[-1]
        with open(self.file_path) as file: self.output_content = file.readlines()
            
        self.n_atoms = self.number_of_atoms()
        self.coordinates = self.xyz_coordinates()
        self.last_xyz = self.coordinates[-1]
            
        self.electronic_energies = self.scf_energies()
        self.electronic_energy = self.electronic_energies[-1]
        
        self.thermo = self.thermochemistry()
        self.gibbs_energy = None
        if self.thermo:
            self.gibbs_energy = self.thermo['gibbs_energy'][-1]
        
        self.dipole_moment = self.dipole()
        
        self.polar = self.polarizability()
        self.isotropic_pol = None
        if self.polar:
            self.isotropic_pol = self.polar['isotropic']
            
        self.nbo = self.nbo_analysis()
        self.npa = None
        self.nbo_orbitals = None
        if self.nbo:
            self.npa = self.nbo['npa']
            self.nbo_orbitals = self.nbo['nbo_orbitals']
            
        self.orb_energies = self.orbitals_energies()
        self.homo_energy = None
        self.lumo_energy = None
        if self.orb_energies:
            self.homo_energy = self.orb_energies['occ'][-1]
            self.lumo_energy = self.orb_energies['virt'][0]

    def number_of_atoms(self):
        search_results = search_lines(self.output_content, 'NAtoms=')
        n_atoms = int(search_results[0][1].split()[1])
        return n_atoms
    
    def xyz_coordinates(self):
        search_results = search_lines(self.output_content, 'Coordinates \(Angstroms\)')
        coordinates = list()
        
        for result in search_results:
            entry_coordinates = list()
            first_line = result[0]+3
            last_line = first_line + self.number_of_atoms()
            
            for line in range(first_line,last_line):
                data = self.output_content[line].split()
                element_number = int(data[1])
                element = PERIODIC_TABLE[element_number]
                float_format = "{:.5f}"
                x_coordinate = float(float_format.format(float(data[3])))
                y_coordinate = float(float_format.format(float(data[4])))
                z_coordinate = float(float_format.format(float(data[5])))
                entry_coordinates.append([element,x_coordinate,y_coordinate,z_coordinate])     
            coordinates.append(entry_coordinates)
        return coordinates
    
    def write_xyz_file(self, coordinate_number=-1):
        filename = f"{self.name.split('.')[0]}.xyz"
        with open(filename,mode='w') as xyz_file:
            xyz_file.write(f"{str(self.n_atoms)}\n")
            xyz_file.write(f"{filename}\n")
            for coordinate in self.coordinates[coordinate_number]:
                coordinates_line = ' '.join(str(x) for x in coordinate)
                xyz_file.write(f"{coordinates_line}\n")
    
    def split_elements_coords(self, coordinate_number=-1):
        elements = list()
        coords = list()

        for coordinate in self.coordinates[coordinate_number]:
            elements.append(coordinate[0])
            coords.append(coordinate[1:])
        
        return elements,coords
    
    def scf_energies(self):
        search_results = search_lines(self.output_content, 'SCF Done')
        electronic_energies = parse_energies_list(search_results, '-[0-9]+.[0-9]+')
        
        return electronic_energies
        
    def thermochemistry(self):
        if search_lines(self.output_content,'- Thermochemistry -'):
            thermo = dict()
            # Vibrational Spectrum
            
            # Gibbs Free Energy
            search_results = search_lines(self.output_content,'Sum of electronic and thermal Free Energies')
            gibbs_energies = parse_energies_list(search_results, '-[0-9]+.[0-9]+')
            thermo['gibbs_energy'] = gibbs_energies
            
            return thermo
        else:
            return None
            
    def dipole(self):
        search_results = search_lines(self.output_content, 'Electric dipole moment \(input')
        if search_results:
            dipole_line = search_results[0][0]+3
            dipole_moment = float(self.output_content[dipole_line].split()[1].replace('D', 'E'))
            return dipole_moment     
    
    def polarizability(self):
        polar = dict()
        search_results = search_lines(self.output_content, 'Dipole polarizability, Alpha \(input')
        if search_results:
            isotropic_pol_line = search_results[0][0]+4
            isotropic_pol = float(self.output_content[isotropic_pol_line].split()[1].replace('D', 'E'))
            polar['isotropic'] = isotropic_pol

            anisotropic_pol_line = search_results[0][0]+5
            anisotropic_pol = None
            polar['anisotropic'] = anisotropic_pol

            return polar
        
    def nbo_analysis(self):
        if search_lines(self.output_content,'N A T U R A L   B O N D   O R B I T A L   A N A L Y S I S'):
            nbo_data = dict()
            
            # natural population analysis (NPA)
            search_results = search_lines(self.output_content, 'Summary of Natural Population Analysis')
            first_atom_line = search_results[0][0]+6
            last_atom_line = first_atom_line + self.n_atoms
            
            nbo_pop = dict()
            for i in range(first_atom_line,last_atom_line):
                natural_data = self.output_content[i].split()
                nbo_pop[int(natural_data[1])] = {
                    'atom_type': natural_data[0],
                    'natural_charge': float(natural_data[2]),
                }
                nbo_data['npa'] = nbo_pop

            # natural bond orbitals (NBO)
            search_results = search_lines(self.output_content, 'NATURAL BOND ORBITALS \(Summary\):')
            first_nbo_line = search_results[0][0] + 7

            search_results = search_lines(self.output_content, '          -------------------------------')
            last_nbo_line = search_results[0][0]

            nbo_orbitals = dict()
            for i in range(first_nbo_line, last_nbo_line):
                content_line = self.output_content[i]
                if re.search('^\s{1,4}[0-9]',content_line):
                    content_line = content_line.strip()
                    if re.search(r'[0-9]{1,3}\([rvg]\)', content_line): 
                        contain_delocalization = True
                    else:
                        contain_delocalization = False
                    #print(content_line) #print for debug purposes

                    content_line = re.sub(r'(\))([A-Z])',r'\1 \2', content_line) #correction for 2-letters atoms (Br, Ca, etc...)
                    content_line = re.sub(r'([0-9]-)([A-Z])',r'\1 \2', content_line) #correction for 2-letters atoms (Br, Ca, etc...)

                    content_line = content_line.replace('(','')
                    content_line = content_line.replace(')','')
                    content_line = re.sub(' +',' ', content_line)
                    content_line = content_line.split()
                    if contain_delocalization: content_line = content_line[:-1]
                    #print(content_line) #print for debug purposes

                    nbo_number = content_line[0]
                    nbo_type = content_line[1]
                    nbo_bond_order = content_line[2]
                    nbo_participants = content_line[3:-2]
                    nbo_occ = float(content_line[-2])
                    nbo_energy = float(content_line[-1])

                    nbo_number = int(nbo_number.replace('.',''))
                    nbo_bond_order = int(nbo_bond_order)

                    nbo_participants = [i.replace('-', '') for i in nbo_participants]
                    nbo_participants_numbers = list()
                    for participant in nbo_participants:
                        try:
                            nbo_participants_numbers.append(int(participant ))
                        except ValueError:
                            pass

                    #print(content_line)
                    #print(f'''
                    #number {nbo_number}
                    #type {nbo_type}
                    #bond order {nbo_bond_order}
                    #participants {nbo_participants_numbers}
                    #occ {nbo_occ}
                    #energy {nbo_energy}
                    #''')

                    #print(nbo_participants)

                    nbo_orbitals[nbo_number] = {
                        'number':nbo_number,
                        'type':nbo_type,
                        'bond order': nbo_bond_order,
                        'participants': nbo_participants_numbers,
                        'occ':nbo_occ,
                        'energy':nbo_energy
                    }
            nbo_data['nbo_orbitals'] = nbo_orbitals
  
            return nbo_data
        else:
            return None
    
    def orbitals_energies(self):
        if search_lines(self.output_content,'Population analysis using the SCF Density'):
            orb_energies = dict()
            
            search_occupied_orbitals = search_lines(self.output_content, 'Alpha  occ. eigenvalues')
            search_virt_orbitals = search_lines(self.output_content, 'Alpha virt. eigenvalues')
            
            occ_energies = list()
            for search_entry in search_occupied_orbitals:
                orbital_energies = search_entry[1].split()[4:]
                for energy in orbital_energies:
                    occ_energies.append(float(energy))
            
            virt_energies = list()
            for search_entry in search_virt_orbitals:
                orbital_energies = search_entry[1].split()[4:]
                for energy in orbital_energies:
                    virt_energies.append(float(energy))
            
            orb_energies['occ'] = occ_energies
            orb_energies['virt'] = virt_energies
            
            return orb_energies
        else:
            return None
