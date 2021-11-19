from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

m1 = Chem.MolFromSmiles('Cc1ccccc1')
cyclohexane = Chem.MolFromSmiles('C1CCCCC1')
methyl_3_bromo_cyclohexene_1 = Chem.MolFromSmiles('CC1=CC(CCC1)Br')
dicyclohexyl = Chem.MolFromSmiles('C1CCCCC1C1CCCCC1') # 2 rings connected by a separate edge
biphenyl = Chem.MolFromSmiles('c1ccccc1-c2ccccc2') # 2 aromatic rings that are connected by a separate edge
spiro_undecane = Chem.MolFromSmiles('C12(CCCCC1)CCCCC2') # 2 rings that share a vertex
perhydroisoquinoline = Chem.MolFromSmiles('N1CC2CCCC2CC1') # 2 rings that share an edge
# indane = Chem.MolFromSmiles('c1ccc2CCCc2c1C1=CC=CC(CCC2)=C12') # NOTE: this one prints 4 rings instead of 2, but it is also aromatic

cubane = Chem.MolFromSmiles('C12C3C4C1C5C4C3C25')
print("cubane num rings: ", Chem.rdmolops.GetSSSR(cubane))
symmSSSRcubane = Chem.GetSymmSSSR(cubane)
for i in range(0, len(symmSSSRcubane)):
    print(list(symmSSSRcubane[i]))

# SMILES -> OUR ENCODING
# 1. find the rings (create a collection C of ring edges) 
# 2. find the branches (s.t. none of the edges in the branch are from C) using neighbors of atoms in the rings
# 3. temporarily assign branches to rings (each ring gets all the branches attached to it)
# 4. for each ring, figure out the ring traversal
#    a. for each vertex that has the longest branch incident on it, try 2 different traversals: clockwise and counterclockwise
#    b. pick the alphabetically the 1st one (if it is symmetric, it does not matter -- an arbitrary top one works)
# 5. Encode rings as strings (along with their branches) and sort them (R1, R2, ....)
# 6. Perform branch trimming (remove duplicate branches -- e.g. if branch A appears in ring R1 and ring R2, then assign it to R1)
# 7. ?? Could remove unnecessary rings

ri = perhydroisoquinoline.GetRingInfo()
symmSSr = Chem.GetSymmSSSR(perhydroisoquinoline)
print("rings info:", ri)
print("symm sssr:", symmSSr[0])
# print(ri.bondRings())
print(MolWt(perhydroisoquinoline))
print(Chem.rdmolops.GetSSSR(perhydroisoquinoline))

# getNeighbors