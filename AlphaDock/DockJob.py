# class that encompasses all features for a docking job and has the ability to cache
# version 1 should have be a simple whole procedure implementation
from AlphaDock.utils import *
from AlphaDock import dir_path, cache_folder_path, get_full_path
#from __init__ import dir_path, cache_folder_path, get_full_path


class DJ:

	def __init__(self, protein_path, ligand_path, job_name):
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

	def to_pdbqt(self):
		print('Converting protein & ligand to .pdbqt with obabel...')
		# protein
		outFile_name_protein = self.JobName + '_' + self.ProteinName
		cached_outfile_protein_pdbqt = convert_to_pdbqt(outFile_name_protein, self.ProteinPath)
		self.Cache['cached_outfile_protein_pdbqt'] = cached_outfile_protein_pdbqt
		# ligand
		outFile_name_ligand = self.JobName + '_' + self.LigandName
		cached_outfile_ligand_pdbqt = convert_to_pdbqt(outFile_name_ligand, self.LigandPath)
		self.Cache['cached_outfile_ligand_pdbqt'] = cached_outfile_ligand_pdbqt
		
