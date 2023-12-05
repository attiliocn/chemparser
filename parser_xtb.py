import re
import numpy as np

class xtbOutput():
    def __init__(self, output_file):
        self.output_file = output_file

    def get_fukui_indexes(self):
        with open(self.output_file) as output:
            for line in output:
                if re.search('(#)\s+(f\(\+\))\s+(f\(\-\))\s+(f\(0\))', line):
                    fukui_content = []
                    while True:
                        current_line = output.readline()
                        if re.search('Property Printout', current_line):
                            fukui_content = fukui_content[:-1]
                            break
                        else:
                            parsed_line = current_line.strip()
                            parsed_line = re.sub(' +',' ',parsed_line)
                            fukui_content.append(parsed_line.split(' '))
            parsed_fukui = {}
            for content in fukui_content:
                atom_number, element = re.search('([0-9]+)([a-zA-Z]{1,2})', content[0]).groups()
                parsed_fukui[atom_number] = {
                    'element': element,
                    'f(+)': float(content[1]),
                    'f(-)': float(content[2]),
                    'f(0)': float(content[3])
                }
        return parsed_fukui