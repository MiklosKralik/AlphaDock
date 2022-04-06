# class that encompasses all features for a docking job and has the ability to cache
# version 1 should have be a simple whole procedure implementation
from AlphaDock.utils import *
from AlphaDock import dir_path, cache_folder_path, get_full_path
import pandas as pd
import json
import random
from shutil import copy2
#from __init__ import dir_path, cache_folder_path, get_full_path


class DJ:

	def __init__(self, protein_path, ligand_path, job_name,
					obabel='obabel', vina='vina'):
		'''version_1:
		func convert to pdbqt,
		func strip
		dock'''
		self.ProteinPath = protein_path
		self.ProteinName = protein_path.split('/')[-1].split('.')[0]
		self.LigandPath = ligand_path
		self.LigandName = ligand_path.split('/')[-1].split('.')[0]
		self.JobName = job_name
		self.Cache = {}
		self.Obabel = obabel
		self.Vina = get_full_path(vina)


	def to_pdbqt(self):
		self.Cache = cache_cache(self.Cache, self.JobName, action='read')
		print('Converting protein & ligand to .pdbqt with obabel...')
		# protein
		outFile_name_protein = self.JobName + '_' + self.ProteinName
		cached_outfile_protein_pdbqt = convert_to_pdbqt(outFile_name_protein, self.ProteinPath, self.Obabel)
		self.Cache['cached_outfile_protein_pdbqt'] = cached_outfile_protein_pdbqt
		# ligand
		outFile_name_ligand = self.JobName + '_' + self.LigandName
		cached_outfile_ligand_pdbqt = convert_to_pdbqt(outFile_name_ligand, self.LigandPath, self.Obabel)
		self.Cache['cached_outfile_ligand_pdbqt'] = cached_outfile_ligand_pdbqt
		self.Cache = cache_cache(self.Cache, self.JobName, action='write')

	def strip_protein(self, source='cache', file_path=None):
		self.Cache = cache_cache(self.Cache, self.JobName, action='read')
		# strips protein of all lines except ATOM
		if source =='cache':
			file_path = self.Cache['cached_outfile_protein_pdbqt']
		elif source =='original':
			file_path = self.ProteinPath

		strip_file_path = make_inflexible(file_path)
		self.Cache['cached_outfile_protein_stripped'] = strip_file_path
		self.Cache = cache_cache(self.Cache, self.JobName, action='write')

	def strip_ligand(self, source='cache', file_path=None):
		self.Cache = cache_cache(self.Cache, self.JobName, action='read')
		# strips protein of all lines except ATOM
		if source =='cache':
			file_path = self.Cache['cached_outfile_protein_pdbqt']
		elif source =='original':
			file_path = self.ProteinPath

		strip_file_path = make_inflexible(file_path)
		self.Cache['cached_outfile_ligand_stripped'] = strip_file_path
		self.Cache = cache_cache(self.Cache, self.JobName, action='write')

	def box(self, padding=1.0):
		self.Cache = cache_cache(self.Cache, self.JobName, action='read')
		# gets box attributes and stores them
		box_attributes = find_bounding_box(self.ProteinPath)
		self.Cache['box_attributes'] = box_attributes
		print(box_attributes)
		self.Cache = cache_cache(self.Cache, self.JobName, action='write')

	def dock(self, num_modes=9, exhaustiveness=8,
	 seed=random.randint(1, 1000000), energy_range=3):
		self.Cache = cache_cache(self.Cache, self.JobName, action='read')

		# uses a set of parameters as a dictionary for docking
		out_path, log_path = run_dock(self.Cache, self.Vina,
		 					self.ProteinName, self.LigandName, self.JobName, num_modes=num_modes,
							exhaustiveness=exhaustiveness, seed=seed, energy_range=energy_range)
		self.Cache['dock_out_path'] = out_path
		self.Cache['dock_log_path'] = log_path
		self.Cache['num_modes'] = num_modes
		self.Cache = cache_cache(self.Cache, self.JobName, action='write')


	def surroundings(self, radius=3):
		self.Cache = cache_cache(self.Cache, self.JobName, action='read')
		df_surroundings = collect_surroundings(radius, self.ProteinPath, self.Cache, self.JobName)
		print(df_surroundings)
		df_path = cache_folder_path + '/' + '_'.join([self.JobName, self.ProteinName, self.LigandName]) + '.csv'
		df_surroundings.to_csv(df_path, index=False)
		self.Cache['df_path'] = df_path
		self.Cache = cache_cache(self.Cache, self.JobName, action='write')

	def plot(self, show=False):
		self.Cache = cache_cache(self.Cache, self.JobName, action='read')
		fig_path = plot_function(self.Cache, self.JobName, self.ProteinName, self.LigandName, show=show)
		# add path to cache and extract_cache
		print(fig_path)
		self.Cache['fig_path'] = fig_path
		self.Cache = cache_cache(self.Cache, self.JobName, action='write')

	def clear_cache(self):
		self.Cache = cache_cache(self.Cache, self.JobName, action='clear')

	def extract_cache(self, destination):
		copy2(self.Cache['dock_out_path'], destination)
		copy2(self.Cache['df_path'], destination)
		copy2(self.Cache['fig_path'], destination)
