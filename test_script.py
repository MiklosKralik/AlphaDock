# general troubleshooting
from AlphaDock.DockJob import DJ
from AlphaDock import get_full_path
JobTest = DJ('test_structures/MTB4_MCE.pdb', 'test_structures/cholesterol.sdf', 'test',
                vina='AlphaDock/vina')
JobTest.to_pdbqt()
JobTest.strip_protein()
JobTest.box(padding=1.0)
JobTest.dock()
JobTest.surroundings()
