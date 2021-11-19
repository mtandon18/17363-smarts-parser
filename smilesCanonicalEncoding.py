from rdkit import Chem
import copy
from collections import defaultdict


'''
Goal
    Classify all edges (bonds) as ring edges or non-ring, branch edges.

Input 
    m- rdkit Mol object
    ringAtoms- set of atoms in rings only

Output:
    ringEdges- set of edges contained in rings
    branchEdgesAll- dictionary mapping atom to all its incident branch edges 
'''
def initializeRingData(m, ringAtoms):
    branchEdgesAll = defaultdict(list)
    ringEdges = set()

    for bond in m.GetBonds():
        # This encoding does not support aromatic rings so request a new molecule
        if str(bond.GetBondType()) == "AROMATIC":
            x = input("This molecule is aromatic. Please try again.\n")
            return encodeSMILES(x)

        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()

        if bond.IsInRing():
            ringEdges.add((start,end))
        else:
            # Exclude mapping branch edges from ring atoms
            # If included, may combine branches with same root into one branch
            if end not in ringAtoms or\
                    (end in ringAtoms and start in ringAtoms):
                branchEdgesAll[end].append(start)
            if start not in ringAtoms or\
                    (end in ringAtoms and start in ringAtoms):
                branchEdgesAll[start].append(end)

    return ringEdges, branchEdgesAll


'''
Goal
    Find all connected components in the branch edges to get trees.

Input 
    branchEdgesAll- dictionary mapping atom to all its incident branch edges 
    ringAtoms- set of atoms in rings only
    m- rdkit Mol object

Output:
    allTrees- list of disconnected branch trees as tuple of atoms and bonds
              (labeled by index ID set by rdkit)
'''
def getAllTreeBranches(branchEdgesAll, ringAtoms, m):
    # while there are still non-ring edges we have not used
    allTrees = []
    while len(branchEdgesAll) > 0:
        thisTreeNodeAtoms = []
        thisTreeNodeBonds = []
        queue = set()
        # add the starting item for this tree:
        (s,neighbors) = branchEdgesAll.popitem()
        thisTreeNodeAtoms.append(s)
        thisTreeNodeAtoms.extend(neighbors)
        thisTreeNodeBonds.extend([m.GetBondBetweenAtoms(s,x).GetIdx() for x in neighbors])
        for i in neighbors:
            if i in branchEdgesAll.keys() and s in branchEdgesAll[i]:
                branchEdgesAll[i].remove(s)
                if branchEdgesAll[i] == []:
                    del branchEdgesAll[i]
            if i not in ringAtoms:
                queue.add(i) # keep the elements in the queue that we wish to explore
        # keep adding all the edges that belong to this tree
        while (len(queue) > 0):
            elem = queue.pop()
            neibs = branchEdgesAll[elem]
            thisTreeNodeAtoms.append(elem)
            thisTreeNodeBonds.extend([m.GetBondBetweenAtoms(elem,x).GetIdx() for x in neibs])
            # remove these edges from branchEdgesAll
            del branchEdgesAll[elem]
            for i in neibs:
                if i not in thisTreeNodeAtoms and i not in ringAtoms:
                    queue.add(i)
            # Moved after check to add to queue
            thisTreeNodeAtoms.extend(neibs)
        tree = (set(thisTreeNodeAtoms),set(thisTreeNodeBonds))
        if tree not in allTrees:
            allTrees.append(tree)
    return allTrees


'''
Goal
    Maps ring edges into edges which correspond to 2 opposing traversals.
    Uses ring vectors which have directionality built in.

Input 
    m- rdkit Mol object
    ringVectors- set of atoms in rings only
    ringEdges- dictionary mapping atom to all its incident branch edges 

Output:
    ringEdgesByStart- maps atom to next atom in ring in one direction
    ringEdgesByEnd- maps atom to next atom in ring in opposite direction
'''
def getRingData(m, ringVectors, ringEdges):
    ringEdgesByStart = defaultdict(dict)
    ringEdgesByEnd = defaultdict(dict)

    # All ring edges implied by rings in ringVectors
    for ring in ringVectors:
        for i in range(len(ring)):
            s = ring[i]
            e = ring[(i+1) % len(ring)]
            ringEdgesByStart[ring][s] = e
            ringEdgesByEnd[ring][e] = s
            ringEdges.difference_update([(s,e),(e,s)])

    # All edges should be found by processing ringVectors
    assert(len(ringEdges) == 0)

    return ringEdgesByStart, ringEdgesByEnd


'''
Goal
    Converts tree branches into their SMILES form for every possible root.

Input 
    allTrees- 
    m- rdkit Mol object
    ringAtoms- set of atoms in rings only

Output:
    allTreeMolecules- dictionary mapping atom index to branch list containing 
                      all branches rooted at that index where each branch is 
                      represented by a 3-tuple (SMILES with atom index tags, 
                      set of atoms in branch, SMILES without atom index tags)
    branchAtoms- set of all atoms contained in branches only
'''
def moleculesFromBranches(allTrees,molecule,ringAtoms):
    allTreeMolecules = defaultdict(list)
    branchAtoms = []
    idxSyms = []
    for i in range(molecule.GetNumAtoms()):
        idxSyms.append(molecule.GetAtomWithIdx(i).GetSymbol() + str(i))
    for tree in allTrees:
        atoms, bonds = tree
        for root in atoms:
            a = list(atoms)
            b = list(bonds)
            branchAtoms.extend(a)
            m = Chem.MolFragmentToSmiles(molecule, atomsToUse=a, bondsToUse=b,
                    rootedAtAtom=root, atomSymbols=idxSyms)
            mstr = Chem.MolFragmentToSmiles(molecule, atomsToUse=a, bondsToUse=b,
                    rootedAtAtom=root)
            i = 1
            while m[i].isdigit():
                i += 1
            allTreeMolecules[root].append((m[i:],atoms, mstr[1:]))
    return allTreeMolecules, set(branchAtoms)


'''
Goal
    Traverse the ring given the implicit direction and starting point while 
    embedding branches.

Input 
    m- rdkit Mol object
    startingIndex- integer starting point for current traversal
    ring- ring vector for ring being traversed
    ringEdges- set of all directed ring edges all in the same direction
    branchesByRoot- dictionary mapping atom index to branch list containing 
                    all branches rooted at that index where each branch is 
                    represented by a 3-tuple (SMILES with atom index tags, 
                    set of atoms in branch, SMILES without atom index tags)

Output:
    traversal- string traversal of ring with all branches embedded
    order- traversal of ring in terms of atom indices
'''
def getTraversalForStartingPoint(m, startingIndex, ring, ringEdges, branchesByRoot):
    order = [startingIndex]
    currIndex = startingIndex
    traversal = ":"
    traversal = traversal + str(m.GetAtomWithIdx(currIndex).GetSymbol())
    if currIndex in branchesByRoot.keys():
        branchStrs = []
        for _,_,b in branchesByRoot[currIndex]:
            branchStrs.append(b)
        branchStrs.sort()
        for b in branchStrs:
            traversal = traversal + "(" + b + ")"
    currIndex = ringEdges[currIndex]
  
    prevIndex = currIndex
    while (currIndex != startingIndex):
        bond = ""
        if prevIndex != currIndex:
            bondType = m.GetBondBetweenAtoms(prevIndex, currIndex).GetBondType()
            if bondType == Chem.BondType.DOUBLE:
                bond = "="
            elif bondType == Chem.BondType.TRIPLE:
                bond = "#"
        order.append(currIndex)
        traversal = traversal + bond + str(m.GetAtomWithIdx(currIndex).GetSymbol())
        if currIndex in branchesByRoot.keys():
            branchStrs = []
            for _,_,b in branchesByRoot[currIndex]:
                branchStrs.append(b)
            branchStrs.sort()
            for b in branchStrs:
                traversal = traversal + "(" + b + ")"
        prevIndex = currIndex
        currIndex = ringEdges[currIndex]
    traversal = traversal + (":")
    bond = ""
    bondType = m.GetBondBetweenAtoms(prevIndex, startingIndex).GetBondType()
    if bondType == Chem.BondType.DOUBLE:
        bond = "="
    elif bondType == Chem.BondType.TRIPLE:
        bond = "#"
    traversal = traversal + bond

    return traversal, order

'''
Goal
    Merges branches and rings into strings by exploring traversals for the rings
    and an overall ordering for each ring-branch string. Branch trimming and
    tagging are performed to reduce redundancy down to the key atoms.

Input 
    m- rdkit Mol object
    branchesByRoot- dictionary mapping atom index to branch list containing 
                    all branches rooted at that index where each branch is 
                    represented by a 3-tuple (SMILES with atom index tags, 
                    set of atoms in branch, SMILES without atom index tags)
    ringEdgesByStart- maps atom to next atom in ring in one direction
    ringEdgesByStart- maps atom to next atom in ring in opposite direction

Output:
    minTagResult- string with encoded SMILES
'''
def getRingTraversals(m, branchesByRoot, ringEdgesByStart, ringEdgesByEnd):
    completeMolecule = []
    stringToRingOrder = dict()
    ringVectors = ringEdgesByStart.keys()

    # get all traversals for each ring with all branches embedded to determine total ordering
    for ring in ringVectors:
        allTraversals = []
        traversalToOrder = dict()
        for startingIndex in ring:
            # try forward direction:
            ringTraversalForward,atomOrderingF = getTraversalForStartingPoint(m, startingIndex, ring, ringEdgesByStart[ring], branchesByRoot)
            # try reverse direction:
            ringTraversalReversed,atomOrderingR = getTraversalForStartingPoint(m, startingIndex, ring, ringEdgesByEnd[ring], branchesByRoot)

            allTraversals.append((ringTraversalForward,ring))
            allTraversals.append((ringTraversalReversed, ring))
            traversalToOrder[(ringTraversalForward, ring)] = atomOrderingF
            traversalToOrder[(ringTraversalReversed, ring)] = atomOrderingR

        allTraversals.sort(key=lambda t: t[0])
        minTraversal = allTraversals[0]
        stringToRingOrder[minTraversal] = traversalToOrder[minTraversal]
        completeMolecule.append(allTraversals[0])

    completeMolecule.sort(key=lambda t: t[0])
    moleculeOrderings = [stringToRingOrder[mol] for mol in completeMolecule]

    # initial string, with atom indices as tags for everything
    result = ""    
    freqDict = defaultdict(lambda: 0)
    idxOrder = []

    for j in range(len(moleculeOrderings)):
        traversal = moleculeOrderings[j]
        result = result + ":"
        prevIdx = traversal[0]
        for i in range(len(traversal)):
            atomIdx = traversal[i]
            idxOrder.append(atomIdx)
            freqDict[atomIdx] += 1
            bond = ""
            if prevIdx != atomIdx:
                bondType = m.GetBondBetweenAtoms(prevIdx, atomIdx).GetBondType()
                if bondType == Chem.BondType.DOUBLE:
                    bond = "="
                elif bondType == Chem.BondType.TRIPLE:
                    bond = "#"

            atom = m.GetAtomWithIdx(atomIdx).GetSymbol()
            result = result + bond + atom + str(atomIdx)

            if atomIdx in branchesByRoot.keys():
                branchStrs = []
                branchAlists = []
                for b,alist,_ in branchesByRoot[atomIdx]:
                    branchStrs.append(b)
                    branchAlists.append(alist)
                branchStrs.sort()
                for b in branchStrs:
                    result = result + "(" + b + ")"
                    i = 0
                    idxStr = ""
                    while i < len(b):
                        if b[i].isalpha() and len(idxStr) > 0:
                            idx = int(idxStr)
                            idxOrder.append(idx)
                            freqDict[idx] += 1
                            idxStr = ""
                        elif b[i].isdigit():
                            idxStr += b[i]

                        i += 1
                    if len(idxStr) > 0:
                        idx = int(idxStr)
                        idxOrder.append(idx)
                        freqDict[idx] += 1
                        idxStr = ""

                for alist in branchAlists:
                    for a in alist:
                        if a in branchesByRoot:
                            for branch in branchesByRoot[a]:
                                if branch[1] == alist:
                                    branchesByRoot[a].remove(branch)
                            if branchesByRoot[a] == []:       
                                del branchesByRoot[a]
            prevIdx = atomIdx
        bond = ""
        bondType = m.GetBondBetweenAtoms(prevIdx, traversal[0]).GetBondType()
        if bondType == Chem.BondType.DOUBLE:
            bond = "="
        elif bondType == Chem.BondType.TRIPLE:
            bond = "#"

        result = result + ":" + bond
        if j != len(moleculeOrderings)-1:
            result = result + "-"

    tagDict = dict() 
    curTag = 0
    for idx in idxOrder:
        if freqDict[idx] > 1 and idx not in tagDict:
            tagDict[idx] = curTag
            curTag += 1
    minTagResult = ""
    idxStr = ""
    i = 0
    while i < len(result):
        if result[i].isdigit():
            idxStr += result[i]
        else:
            if len(idxStr) > 0:
                idx = int(idxStr)
                if idx in tagDict:
                    minTagResult += str(tagDict[idx])
                idxStr = ""
            minTagResult += result[i]
        i += 1
    if len(idxStr) > 0:
        idx = int(idxStr)
        if idx in tagDict:
            minTagResult += str(tagDict[idx])
        idxStr = ""

    return minTagResult

'''
Goal
    Encode SMILES into our format.

Input 
    s- SMILES string of valid, non-aromatic molecule

Output:
    Encoded SMILES that is canonical and parser-friendly
'''
def encodeSMILES(s):
    m = Chem.MolFromSmiles(s) # gets the molecule

    # Asks for input again if molecule not valid SMILES
    if m == None:
        x = input("Please try again.\n")
        return encodeSMILES(x)

    ri = m.GetRingInfo()
    n = ri.NumRings()

    # If there are any rings in structure, we need to perform encoding so
    # molecule can be parsed using context-free grammar
    if n > 0:
        ringVectors = ri.AtomRings() # atoms on the rings only
        ringAtoms = set()
        for ring in ringVectors:
            for r in ring:
                ringAtoms.add(r)

        # Separate ring edges from branch edges 
        ringEdges, branchEdgesAll = initializeRingData(m, ringAtoms)

        # Extend tree branches, getting distinct connected components
        allTreeBranches = getAllTreeBranches(branchEdgesAll, ringAtoms, m)

        # Create 2 ring edge dictionaries the 2 traversal directions 
        # (arbitrarily clockwise and counter-clockwise)
        ringEdgesByStart, ringEdgesByEnd = getRingData(m, ringVectors, ringEdges)

        # Find SMILES encodings for all branches by root
        branchesByRoot, branchAtoms = moleculesFromBranches(allTreeBranches, m, ringAtoms)

        # Get final SMILES encoding in canonical order
        result = getRingTraversals(m, branchesByRoot, ringEdgesByStart, ringEdgesByEnd)
        return result
    else:
        return s
        # Output SMILES

# List of simple test molecules used to debug code
if __name__ == "__main__":
    cubane = 'C12C3C4C1C5C4C3C25'
    bicyclohexyl = 'C1CCCCC1C1CCCCC1'
    cyclohexane = 'C1CCCCC1'
    _3_propyl_4_isopropyl_1_heptene = 'C1=CC(CCC)C(C(=C)C)CCC1'
    testMultBranchSameRoot1 = 'C1=CC(CCC)(C(=C)C)CCC1'
    testMultBonds1 = 'C1=C=C=C=C=C=C=1'
    _3_ringed_multiple_branches = 'C1(C2=CC(F)=C(C=C2N(C2CC2)C=C1C(=O)O)N1CCNCC1)=O'
    _1_methyl_3_bromo_cyclohexene_1 = 'CC1=CC(CCC1)Br'
    testSymHack1 = 'C1OSN1CCCCCC1PSN1'
    testSymHack2 = 'C1CCC1C1CCC1'
    #print(encodeSMILES("O=C1NC(=O)C(C2CCCS2)(C2CCCS2)N1"))
