# all utilitiy and helper functions
# it is assumed that vina & obabel are within the environment
import os
from subprocess import run
from AlphaDock import dir_path, cache_folder_path, get_full_path
from pymol import cmd
from pymol import stored
import pandas as pd
import json

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
    dump_data = cache_folder_path
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

    # adding all information containing lines to file
    info_lines = []

    for line in log_file[-2::-1]:
        info_lines.append(line)
        if line.split(' ')[3] == '1':
            break
    # invert as is in the wrong order
    info_lines = info_lines[::-1]
    #print(info_lines[0].split(' '))

    affinity = {}
    rmsd_lb = {}
    rmsd_ub = {}

    for line in info_lines:
        line_split = []

        for i in line.split(' '):
            if i != '':
                line_split.append(i)


        affinity[int(line_split[0])] = float(line_split[1])
        rmsd_lb[int(line_split[0])] = float(line_split[2])
        rmsd_ub[int(line_split[0])] = float(line_split[3][:-1])

    for state in df_surroundings['ligand_state'].unique():
        df_surroundings.loc[df_surroundings['ligand_state']==state, 'affinity'] = affinity[state]
        df_surroundings.loc[df_surroundings['ligand_state']==state, 'rmsd_lb'] = rmsd_lb[state]
        df_surroundings.loc[df_surroundings['ligand_state']==state, 'rmsd_ub'] = rmsd_ub[state]

    print(df_surroundings.info())
    return df_surroundings

def cache_cache(cache, job_name, action):
    # create full cache so different jobs are not mixed up
    cache_path = cache_folder_path + '/' + f'cache{job_name}.json'

    if os.path.exists(cache_path):
        if action == 'read':
            with open(cache_path, 'r+') as file:
                cache = json.load(file)
        if action == 'clear':
            os.remove(cache_path)

    if action == 'write':
        with open(cache_path, 'w+') as file:
            cache_obj = json.dumps(cache)
            file.write(cache_obj)

    return cache


if __name__ == '__main__':
    print(dir_path)
    print(cache_folder_path)
