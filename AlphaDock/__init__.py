# initialize python package
import os

def get_full_path(path):
    # use to only use absolute paths
    return os.path.abspath(path)

dir_path = get_full_path(os.path.dirname(os.path.realpath(__file__)))
cache_folder_path = get_full_path(dir_path + '/cache_folder')
