# all utilitiy and helper functions
# it is assumed that vina & obabel are within the environment
import os
from subprocess import run
from AlphaDock import dir_path, cache_folder_path, get_full_path
from pymol import cmd
from pymol import stored
import pandas as pd
import json
import random
import matplotlib.pyplot as plt

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

def run_dock(cache, vina, protein_name, ligand_name, job_name, num_modes=9,
            exhaustiveness=8, seed=random.randint(1, 1000000), energy_range=3):
    # use vina to dock protein and define docking parameters
    dump_data = cache_folder_path
    box_attributes = cache['box_attributes']
    protein_path = cache['cached_outfile_protein_stripped']
    ligand_path = cache['cached_outfile_ligand_pdbqt']
    out_path_full = dump_data + '/' + '_'.join([job_name, protein_name, ligand_name]) + '.pdb'
    log_path_full = dump_data + '/' + '_'.join([job_name, protein_name, ligand_name]) + '.txt'

    command = [vina, '--receptor', protein_path,
    '--ligand', ligand_path,
    '--center_x', box_attributes['cenX'],
    '--center_y', box_attributes['cenY'],
    '--center_z', box_attributes['cenZ'],
    '--size_x', box_attributes['lenX'],
    '--size_y', box_attributes['lenY'],
    '--size_z', box_attributes['lenZ'],
    '--out', out_path_full,
    '--log', log_path_full,
    '--num_modes', str(num_modes),
    '--exhaustiveness', str(exhaustiveness),
    '--seed', str(seed),
    '--energy_range', str(energy_range)]

    print(f'Docking {protein_name} & {ligand_name} Job: {job_name} with vina....')
    process = run(command, capture_output=True, text=True)
    print(process.stdout)
    print(process.stderr)
    return get_full_path(out_path_full), get_full_path(log_path_full)



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

    for state in range(1, cache['num_modes']+1):
        try:
            state_selection_name = f'surroundings{state}'
            cmd.select(state_selection_name, docked_ligand_name, state=state)
            # removed byres because that constricts information about the proximity of the lipid
            string_selection = ' '.join([protein_name, 'within', radius, 'of', state_selection_name])
            cmd.select('temp', string_selection)
            stored.state = state
            cmd.iterate('temp', 'stored.data.append([resn, resi, chain, name, stored.state])')
        except:
            break
            print('NUM MODES IS LARGER THAN NUM STATES')
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

# function that takes labels as input and outputs a list of colors accordingly
# input can contain chains from A-G and all residues

def consistent_colors(labels):

    # dict to store corresponding labels and colors
    color_dict = {'A':'blue', 'B':'green', 'C':'red', 'D':'cyan',
                  'E':'magenta', 'F':'yellow', 'G':'white', 'H':'silver',
                 'ARG':'brown', 'HIS':'tomato', 'LYS':'peru', 'ASP':'darkorange',
                  'GLU':'tan', 'SER':'gold', 'THR':'yellow', 'ASN':'limegreen',
                  'GLN':'mintcream','CYS':'teal', 'GLY':'steelblue', 'PRO':'navy',
                  'ALA':'indigo','VAL':'plum', 'ILE':'hotpink', 'LEU':'crimson',
                  'MET':'moccasin', 'PHE':'gray', 'TYR':'maroon', 'TRP':'peachpuff'}
    #print([mcolors.CSS4_COLORS[i] for i in color_dict.values()])
    colors = [color_dict[i] for i in labels]
    return colors

def calculate_percentage(values):
    total = sum(values)
    values = [(i/total)*100 for i in values]
    return values

def plot_chains(df, ax):
    chains = df.groupby('chain').size().sort_values(ascending=False)
    labels = chains.keys()
    values = chains.values
    ax.set_title('Chains')
    ax.set_ylabel('Percent of docking states')
    ax.set_xlabel('Chains')
    ax.axes.xaxis.set_visible(False)
    colors = consistent_colors(labels)
    values = calculate_percentage(values)

    return ax.bar_label(ax.bar([i for i in range(len(values))], values, color=colors, edgecolor='dimgrey'), labels=labels)

def plot_residues(df, ax):
    residues = df.groupby('resn').size().sort_values()
    labels = residues.keys()
    values = residues.values
    ax.set_xlim(0, 100)
    ax.set_xlabel('Percent of docking states')
    ax.set_title('Residues')
    ax.axes.yaxis.set_visible(False)
    colors = consistent_colors(labels)
    values = calculate_percentage(values)

    return ax.bar_label(ax.barh([i for i in range(len(values))], values, color=colors, edgecolor='dimgrey'), labels=labels)

def plot_specific_residues_chain_affinity(df, ax):
    df['resn_resi_chain_aff'] = df['resn'] + df['resi'].astype(str) + df['chain']  + df['affinity'].astype(str)
    residues = df.groupby('resn_resi_chain_aff').size().sort_values()
    labels = residues.keys()[-20:]
    values = residues.values[-20:]
    ax.set_xlim(0, 100)
    ax.set_xlabel('Percent of docking states')
    ax.set_title('Residues & Chains & Affinities')
    ax.axes.yaxis.set_visible(False)
    # labels are in different format for coloring
    c_labels = [i[:3] for i in labels]
    colors = consistent_colors(c_labels)
    values = calculate_percentage(values)

    return ax.bar_label(ax.barh([i for i in range(len(values))], values, color=colors, edgecolor='dimgrey'), labels=labels)


def plot_function(cache, job_name, protein_name, ligand_name, show):
    df = pd.read_csv(cache['df_path'])
    # unionized colors
    # general chains
    # general residue types
    # binding affinity boxplots
    ############################
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(1,3, figsize=(20, 10))
    plt.tight_layout()

    plot_chains(df, ax[0])
    plot_residues(df, ax[1])
    plot_specific_residues_chain_affinity(df, ax[2])
    fig_path = cache_folder_path + '/' + '_'.join([job_name, protein_name, ligand_name, 'charts']) + '.jpeg'

    plt.savefig(fig_path, dpi=200)
    if show:
        plt.show()

    return fig_path

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
