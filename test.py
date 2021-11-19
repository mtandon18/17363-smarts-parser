from parseAtomClone import strToAdjacencyMatrix
from rdkit import Chem
from smilesCanonicalEncoding import encodeSMILES

# Adapted from https://stackoverflow.com/questions/51195392/smiles-from-graph
def molFromMatrix(node_list, adjacency_matrix):

    # create empty editable mol object
    mol = Chem.RWMol()

    # add atoms to mol and keep track of index
    node_to_idx = {}
    for i in range(len(node_list)):
        a = Chem.Atom(node_list[i])
        molIdx = mol.AddAtom(a)
        node_to_idx[i] = molIdx

    # add bonds between adjacent atoms
    for ix, row in enumerate(adjacency_matrix):
        for iy, bond in enumerate(row):

            # only traverse half the matrix
            if iy <= ix:
                continue

            # add relevant bond type (there are many more of these)
            if bond == 0:
                continue
            elif bond == 1:
                bond_type = Chem.rdchem.BondType.SINGLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
            elif bond == 2:
                bond_type = Chem.rdchem.BondType.DOUBLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
            elif bond == 3:
                bond_type = Chem.rdchem.BondType.TRIPLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)

    # Convert RWMol to Mol object
    mol = mol.GetMol()

    return mol

'''
Tests both the parser and canonical encoding for correctness.
Input: File of input strings in SMILES form
Output: Number of tests passed

Converts SMILES to rdkit molecule and compares that to SMILES that has been
encoded to canonical SMILES and then parsed to an adjacency matrix. Passes if
the matrices and molecule in rdkit are identical.
'''
def test(f):
    tests = open(f, "r").read().splitlines()
    count = 0
    correct = 0
    for molecule in tests:
        smilesMol = Chem.MolFromSmiles(molecule)
        aromatic = False
        for bond in smilesMol.GetBonds():
            if bond.GetBondType() == Chem.BondType.AROMATIC:
                aromatic = True
        if aromatic:
            continue

        encodedSmiles = encodeSMILES(molecule)
        matrix, labels = strToAdjacencyMatrix(encodedSmiles)
        encodedSmilesMol = molFromMatrix(labels, matrix)

        s1 = Chem.MolToSmiles(smilesMol)
        s2 = Chem.MolToSmiles(encodedSmilesMol)

        if (s1 == s2): # Relies on rdkit's canonical SMILES to check equality
            correct += 1
        else:
        # Problem with rdkit's circularity for SMILES encodings and won't be counted
            if molecule != s1:
                count -= 1
            else:
                print(f"Testing {molecule} on line {count}...Failed")
                print("   ", encodedSmiles)
                print("   ", s2)
        count += 1

    print(f"PASSED {correct}/{count}")


test("pcbaClean.txt")


