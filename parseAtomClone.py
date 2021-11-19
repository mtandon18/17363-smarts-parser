from arpeggio import Optional, ZeroOrMore, OneOrMore, EOF, StrMatch
from arpeggio import RegExMatch as _
from arpeggio import ParserPython
from arpeggio import PTNodeVisitor, visit_parse_tree
import sys
import copy

# From Arpeggio tutorial, prevents string from appearing in parse tree
class SuppressStr(StrMatch):
    suppress = True

# Parser Functions
def atomUnit():
    # Second option lists accepted atoms, can add more if necesary
    return [(SuppressStr('('), atomList, SuppressStr(')')),
            (Optional(['=','#']), _(r'Br|Cl|[BCFHINOPS]'), Optional(_(r'\d*')))]

def atomList():
    return ZeroOrMore(atomUnit)

def ring():
    return SuppressStr(':'), atomList, \
        SuppressStr(':'), Optional(['=','#'])

def fragment():
    return [ring, atomList]

def fragmentList():
    return fragment, ZeroOrMore(SuppressStr('-'),fragment), EOF

# Visitor class that walks parse tree
# Returns list of fragments,
# Each fragment is tuple of 2 lists,
# First list is starting branch,
# Second list is ring
# Any atomList in parenthesis is its own list
# Any atomUnit is tuple of Atom string and int tag
class AtomCloneVistor(PTNodeVisitor):
    def visit_fragmentList(self, node, children):
        return list(children)

    def visit_fragment(self, node, children):
        d = dict()
        for i in range(len(node)):
            if node[i].rule_name == "atomList":
                d[node[i].rule_name] = children[i]
            elif node[i].rule_name == "ring":
                d[node[i].rule_name], bond = children[i]
                d["ringBond"] = bond

        return d

    def visit_ring(self, node, children):
        if len(children) == 0 or children[0] in ['=','#']:
            return None
        bond = 1
        if len(children) == 2 and children[1] == '=':
            bond = 2
        elif len(children) == 2 and children[1] == '#':
            bond = 3

        return (children[0], bond)

    def visit_atomList(self, node, children):
        return list(children)

    def visit_atomUnit(self, node, children):
        # Parenthesized branch
        if node[0].rule_name == "atomList":
            return children[0]

        noBondChildren = children
        bond = 1
        tag = None
        if children[0] == '=':
            bond = 2
            noBondChildren = children[1:]
        elif children[0] == '#':
            bond = 3
            noBondChildren = children[1:]

        # Atom and tag
        if len(noBondChildren) == 2:
            d= dict()
            d["atom"] = noBondChildren[0]
            d["tag"] = int(noBondChildren[1])
            d["bond"] = bond
            return d
        # Atom without tag
        else:
            d= dict()
            d["atom"] = noBondChildren[0]
            d["bond"] = bond
            return d


def initParser():
    """
    Set parser global variable to arpeggio parser for our format
    """
    global parser
    parser = ParserPython(fragmentList)

def fragmentAdjacencyMatrix(fragment):
    """
    Takes tuple representing a fragment and returns adjacency matrix,
    label list, and tag dictionary.

    arguments:
    fragment - tuple representing fragment as returned from visitor class

    returns:
    adjacency matrix as list of lists, label matrix as list of strings,
    dictionary mapping tags to index
    """
    numAtoms = 0
    ringStart = None
    labels = []
    bonds = []
    tags = dict()
    fidxStack = [0]
    atomStack = [-1]
    struct = "atomList"
    if struct not in fragment:
        struct = "ring"
        ringStart = 0
    # Traverse fragment in DFS order with iterative stack
    while len(fidxStack) > 0:
        # Use struct and fidxStack to get current object
        obj = fragment[struct]
        for i in fidxStack:
            if i < len(obj):
                obj = obj[i]
            else:
                # When index at one level goes out of bounds,
                # return to previous level
                fidxStack.pop()
                # Check if finished with ring, add ring closing bond
                if len(fidxStack) == 0 and struct == "ring":
                    bonds.append((ringStart, atomStack[0], fragment["ringBond"]))
                atomStack.pop()
                if len(fidxStack) > 0:
                    fidxStack[-1] += 1
                obj = None
                break
        # Handle individual atom by adding tag, label, and implicit bond
        if isinstance(obj, dict):
            labels.append(obj["atom"])
            if "tag" in obj:
                # assert (obj["tag"] not in tags), "Duplicate tag in same fragment"
                tags[obj["tag"]] = numAtoms

            bonds.append((atomStack[-1],numAtoms,obj["bond"]))
            fidxStack[-1] += 1
            atomStack[-1] = numAtoms
            numAtoms += 1
        # Handle branch (atomList) by setting up parent atom for next bond
        elif isinstance(obj, list):
            fidxStack.append(0)
            # atomStack.append(numAtoms-1)
            atomStack.append(atomStack[-1])

    # Remove dummy starter bond
    bonds = bonds[1:]

    # Make adjacency matrix with all zeros
    mat = []
    for i in range(numAtoms):
        mat.append([0]*numAtoms)

    for (i,j, b) in bonds:
        mat[i][j] = b
        mat[j][i] = b

    return mat, labels, tags

def mergeFragments(a, b):
    """
    Merges two fragment adjacency matrices and labels using tags.

    arguments:
    a - Fragment tuple (matrix, label, tag)
    b - Second fragment tuple (matrix, label, tag)

    returns:
    Merged fragment tuple (matrix, label, tag)
    """
    matA = a[0]
    labA = a[1]
    tagA = a[2]

    matB = b[0]
    labB = b[1]
    tagB = b[2]

    # Compute final amount of atoms in mergeed fragment
    newNumAtom = len(labA) + len(labB)
    intersect = [k for k in tagB if k in tagA]
    newNumAtom -= len(intersect)

    newLab = copy.copy(labA)
    # Construct list mapping index in B to combined index (accounting for tags)
    idxB = [None] * len(labB)
    for k in intersect:
        assert labA[tagA[k]] == labB[tagB[k]], "Same tag on different atom types"
        idxB[tagB[k]] = tagA[k]

    currIdx = len(labA)
    for i in range(len(idxB)):
        if idxB[i] is None:
            idxB[i] = currIdx
            newLab.append(labB[i])
            currIdx += 1
    for k in tagB.keys():
        tagB[k] = idxB[tagB[k]]
    newTag = {**tagB, **tagA}

    newMat = []
    for i in range(newNumAtom):
        newMat.append([0]*newNumAtom)

    # Copy A matrix
    for i in range(len(matA)):
        for j in range(len(matA)):
            newMat[i][j] = matA[i][j]

    # Add bonds from B matrix
    for i in range(len(matB)):
        for j in range(len(matB)):
            newMat[idxB[i]][idxB[j]] = max(matB[i][j],
                    newMat[idxB[i]][idxB[j]])

    return (newMat, newLab, newTag)

def strToAdjacencyMatrix(s):
    """
    Parses string and returns adjacency matrix and label array

    arguments:
    s - string to be parsed

    returns:
    adjacency matrix as list of lists, label matrix as list of strings
    """
    initParser()
    tree = parser.parse(s)
    listForm = visit_parse_tree(tree, AtomCloneVistor())
    currFragment = fragmentAdjacencyMatrix(listForm[0])
    for nextList in listForm[1:]:
        nextFragment = fragmentAdjacencyMatrix(nextList)
        currFragment = mergeFragments(currFragment, nextFragment)

    return currFragment[0:2]

def printMatrix(mat):
    """
    Prints matrix (list of lists) in nice format
    (assumes all entries take same number of characters)

    arguments:
    mat - list of lists representing matrix
    """
    for row in mat:
        print(' '.join([str(x) for x in row]))

def printLabels(lab):
    """
    Prints label list in nice format

    arguments:
    mat - list of lists representing matrix
    """
    for i,c in enumerate(lab):
        print(str(i)+" "+c)

# initParser()
# test = ":C(C)C(C)C0C1S:-:C(CNC2)CCCO:-:C2C0C1NCN:"
# mat,lab = strToAdjacencyMatrix(test)
# print(test)
# printMatrix(mat)
# printLabels(lab)
#assert(len(sys.argv) == 2)

# initParser()

# with open(sys.argv[1]) as f:
#     for line in f:
#         cleanLine = line.strip()
#         print("String: ", cleanLine)
#         matrix, lab = strToAdjacencyMatrix(cleanLine)
#         printMatrix(matrix)
#         printLabels(lab)
