# general troubleshooting
from AlphaDock.DockJob import DJ
from AlphaDock import get_full_path, clear_cache_folder
JobTest = DJ('test_structures/MTB4_MCE.pdb', 'test_structures/cholesterol.sdf', 'test_',
                vina='AlphaDock/vina')
clear_cache_folder()
#JobTest.clear_cache()
#JobTest.to_pdbqt()
#JobTest.strip_protein()
#JobTest.box(padding=1.0)
#JobTest.dock()
#JobTest.surroundings()
