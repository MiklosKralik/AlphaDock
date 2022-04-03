# general troubleshooting
from AlphaDock.DockJob import DJ
from AlphaDock import get_full_path
JobTest = DJ('test_structures/MTB4_MCE.pdb', 'test_structures/cholesterol.sdf', 'test')
JobTest.to_pdbqt()
print(JobTest.Cache)
