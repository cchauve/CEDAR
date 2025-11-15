"""
CEDAR: manipulating phylogenetic rooted trees representations as vectors
Generic read/write functions
"""

__author__ = "Cedric Chauve"
__credits__ = ["Cedric Chauve", "Louxin Zhang"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Cedric Chauve"
__email__ = "cedric.chauve@sfu.ca"
__status__ = "Release"

def __read_file(input_file):
    """
    Read a file and returns a list of strings, one per line
    """
    lines = []
    with open(input_file) as in_file:
        for line in in_file.readlines():
            if line.strip():
                lines.append(line.rstrip())
    return lines

def __write_file(in_lines, output_file):
    """
    Write a list of strings into a file
    """
    with open(output_file, "w") as out_file:
        for line in in_lines:
            out_file.write(f"{line}\n")

            
