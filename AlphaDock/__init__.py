# initialize python package
import os
import glob

def get_full_path(path):
    # use to only use absolute paths
    return os.path.abspath(path)

dir_path = get_full_path(os.path.dirname(os.path.realpath(__file__)))
cache_folder_path = get_full_path(dir_path + '/cache_folder')

def clear_cache_folder():
    files = glob.glob(cache_folder_path + '/*')
    print('removing files -------')
    for f in files:
        if f.split('.')[-1] in ['txt', 'pdb', 'pdbqt', 'csv', 'json', 'jpeg']:
            print(f)
            os.remove(f)
