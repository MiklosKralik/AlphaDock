# all utilitiy and helper functions
# it is assumed that vina & obabel are within the environment
import os
from subprocess import run
from AlphaDock import dir_path, cache_folder_path, get_full_path
from pymol import cmd
from pymol import stored
import pandas as pd

def convert_to_pdbqt(outFile_name, file_path, obabel):
    file_path = get_full_path(file_path)
    outfile = '/'.join([cache_folder_path, outFile_name + '.pdbqt'])
    command = [obabel, file_path, '-o', 'pdbqt', '-O', outfile]
    process = run(command, capture_output=True, text=True)
    print(process.stdout)
    print(process.stderr)
    return outfile

def make_inflexible(file_path):
    print(f'making {file_path} inflexible...')

    with open(file_path) as f:
        file = f.readlines()
    print('original_length: {}'.format(len(file)))

    strip_file = []
    for line in file:
        if line[:4] == 'ATOM':
            strip_file.append(line)
    strip_file_path = '_inflexible.'.join(file_path.split('.'))
    print('new_length: {}'.format(len(strip_file)))

    with open(strip_file_path, 'w') as f:
        f.writelines(strip_file)

    return strip_file_path

def find_bounding_box(molecule_path, padding=1.0):
    #load protein into pymol and retrieve center cooridinates & edge lengths
    cmd.delete('all')
    cmd.load(molecule_path)
    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent('all')

    box_attributes = {}
    # box edges
    box_attributes['lenX'] = str((maxX - minX) + padding)
    box_attributes['lenY'] = str((maxY - minY) + padding)
    box_attributes['lenZ'] = str((maxZ - minZ) + padding)
    # box centre
    box_attributes['cenX'] = str(minX + (maxX - minX) / 2)
    box_attributes['cenY'] = str(minY + (maxY - minY) / 2)
    box_attributes['cenZ'] = str(minZ + (maxZ - minZ) / 2)

    return box_attributes

def run_dock(cache, vina, protein_name, ligand_name, job_name, num_modes):
    # use vina to dock protein and define docking parameters
    dump_data = 'outputs'
    box_attributes = cache['box_attributes']
    protein_path = cache['cached_outfile_protein_stripped']
    ligand_path = cache['cached_outfile_ligand_pdbqt']
    out_path = dump_data + '/' + job_name + protein_name + ligand_name + '.pdbqt'
    log_path = dump_data + '/' + job_name + protein_name + ligand_name + '.txt'

    command = [vina, '--receptor', protein_path,
    '--ligand', ligand_path,
    '--center_x', box_attributes['cenX'],
    '--center_y', box_attributes['cenY'],
    '--center_z', box_attributes['cenZ'],
    '--size_x', box_attributes['lenX'],
    '--size_y', box_attributes['lenY'],
    '--size_z', box_attributes['lenZ'],
    '--out', out_path,
    '--log', log_path,
    '--num_modes', str(num_modes)]

    print(f'Docking {protein_name} & {ligand_name} Job: {job_name} with vina....')
    process = run(command, capture_output=True, text=True)
    print(process.stdout)
    print(process.stderr)
    return get_full_path(out_path), get_full_path(log_path)



def collect_surroundings(radius, protein_path, cache, job_name):
    # collect the surroundings of each protein in pymol
    docked_ligand_path = cache['dock_out_path']
    docked_ligand_name = docked_ligand_path.split('/')[-1].split('.')[0]
    protein_name = protein_path.split('/')[-1].split('.')[0]
    radius = str(radius)
    print('Collecting surroundings ' + docked_ligand_name + protein_name + job_name + str(radius) + '......')
    # make sure environment is clear
    cmd.delete('all')

    #load molecules
    cmd.load(docked_ligand_path)
    cmd.load(protein_path)
    stored.data = []

    for state in range(1, 10):
        state_selection_name = f'surroundings{state}'
        cmd.select(state_selection_name, docked_ligand_name, state=state)
        string_selection = ' '.join(['byres', protein_name, 'within', radius, 'of', state_selection_name])
        cmd.select('temp', string_selection)
        stored.state = state
        cmd.iterate('temp', 'stored.data.append([resn, resi, chain, name, stored.state])')
    #print(stored.data)
    columns = ['resn', 'resi', 'chain', 'name', 'ligand_state']
    df_surroundings = pd.DataFrame(data=stored.data, columns=columns).drop_duplicates()

    # parse log file for affinity vaues
    with open(cache['dock_log_path']) as f:
        log_file = f.readlines()
    #print(log_file[-1:])
    return df_surroundings

if __name__ == '__main__':
    print(dir_path)
    print(cache_folder_path)
