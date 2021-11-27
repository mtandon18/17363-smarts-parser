from arpeggio import Optional, ZeroOrMore, EOF, StrMatch
from arpeggio import RegExMatch as _
from arpeggio import ParserPython
from arpeggio import PTNodeVisitor, visit_parse_tree

# From Arpeggio tutorial, prevents string from appearing in parse tree
class SuppressStr(StrMatch):
    suppress = True

def atom():
    return [(_(r'Br|Cl|[BCFHINOPS]')), (_(r'br|cl|[bcfhinops]'))]

def charge():
    return [(_(r'\+\d*')), (_(r'-\d*'))]

def atomWithOptCharge():
    return [(SuppressStr('['), atomUnit, Optional(charge), SuppressStr(']')), (atom)]

def bond():
    return ['-', '/', '\\', '/?', '\?', '=', '#', ':', '~', '@']

def wildCard():
    return '*'

def aromatic():
    return 'A'

def degree():
    return (_(r'D\d*'))

def totalHCount():
    return (_(r'H\d*'))

def implicitHCount():
    return (_(r'h\d*'))

def valence():
    return (_(r'v\d*'))

def connectivity():
    return (_(r'X\d*'))

def ringConnectivity():
    return (_(r'x\d*'))

def ringMembership():
    return (_(r'R\d*'))

def ringSize():
    return (_(r'r\d*'))

def atomicNumber():
    return (_(r'#\d*'))

def atomicPrimitive():
    return [charge, wildCard, aromatic, degree, totalHCount, implicitHCount,
    valence, connectivity, ringConnectivity, ringMembership, ringSize, 
    atomicNumber]

# Parser Functions
def atomUnit():
    # Second option lists accepted atoms, our subset will only support these atoms
    # Because of PEG, we match with a branch first and only if this doesn't work
    # do we attempt to match with a ring and then just an atom
    return [(SuppressStr('('), ZeroOrMore(atomUnit), SuppressStr(')')),
            (Optional(bond), ZeroOrMore(atomWithOptCharge), Optional(ringClosure))]

def fragment():
    return atom, ZeroOrMore(atomUnit), EOF

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
    parser = ParserPython(fragment)

def parseSMARTS(s):
    initParser()
    tree = parser.parse(s)
    #listForm = visit_parse_tree(tree, AtomCloneVistor())


# f = open("pcbaClean.txt").read()  
# initParser()
# for line in f.splitlines():
#     tree = parser.parse(line)
#     print(tree)