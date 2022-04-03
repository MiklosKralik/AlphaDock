# all utilitiy and helper functions
# it is assumed that vina & obabel are within the environment
import os
from subprocess import run
from AlphaDock import dir_path, cache_folder_path, get_full_path


def convert_to_pdbqt(outFile_name, file_path):
    file_path = get_full_path(file_path)
    obabel = 'obabel'
    outfile = '/'.join([cache_folder_path, outFile_name + '.pdbqt'])
    command = [obabel, file_path, '-o', 'pdbqt', '-O', outfile]
    process = run(command, capture_output=True, text=True)
    print(process.stdout)
    print(process.stderr)
    return outfile

if __name__ == '__main__':
    print(dir_path)
    print(cache_folder_path)
