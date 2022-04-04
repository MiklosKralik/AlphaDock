# class that encompasses all features for a docking job and has the ability to cache
# version 1 should have be a simple whole procedure implementation
from AlphaDock.utils import *
from AlphaDock import dir_path, cache_folder_path, get_full_path
import pandas as pd
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
		print('Converting protein & ligand to .pdbqt with obabel...')
		# protein
		outFile_name_protein = self.JobName + '_' + self.ProteinName
		cached_outfile_protein_pdbqt = convert_to_pdbqt(outFile_name_protein, self.ProteinPath, self.Obabel)
		self.Cache['cached_outfile_protein_pdbqt'] = cached_outfile_protein_pdbqt
		# ligand
		outFile_name_ligand = self.JobName + '_' + self.LigandName
		cached_outfile_ligand_pdbqt = convert_to_pdbqt(outFile_name_ligand, self.LigandPath, self.Obabel)
		self.Cache['cached_outfile_ligand_pdbqt'] = cached_outfile_ligand_pdbqt

	def strip_protein(self, source='cache', file_path=None):
		# strips protein of all lines except ATOM
		if source =='cache':
			file_path = self.Cache['cached_outfile_protein_pdbqt']
		elif source =='original':
			file_path = self.ProteinPath

		strip_file_path = make_inflexible(file_path)
		self.Cache['cached_outfile_protein_stripped'] = strip_file_path

	def strip_ligand(self, source='cache', file_path=None):
		# strips protein of all lines except ATOM
		if source =='cache':
			file_path = self.Cache['cached_outfile_protein_pdbqt']
		elif source =='original':
			file_path = self.ProteinPath

		strip_file_path = make_inflexible(file_path)
		self.Cache['cached_outfile_ligand_stripped'] = strip_file_path

	def box(self, padding=1.0):
		# gets box attributes and stores them
		box_attributes = find_bounding_box(self.ProteinPath)
		self.Cache['box_attributes'] = box_attributes
		print(box_attributes)

	def dock(self, num_modes=9):
		# uses a set of parameters as a dictionary for docking
		out_path, log_path = run_dock(self.Cache, self.Vina,
		 					self.ProteinName, self.LigandName, self.JobName, num_modes)
		self.Cache['dock_out_path'] = out_path
		self.Cache['dock_log_path'] = log_path

	def surroundings(self, radius=3):
		df_surroundings = collect_surroundings(radius, self.ProteinPath, self.Cache, self.JobName)
		print(df_surroundings)
