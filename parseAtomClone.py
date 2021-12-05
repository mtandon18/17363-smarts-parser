from arpeggio import Optional, ZeroOrMore, EOF, StrMatch
from arpeggio import RegExMatch as _
from arpeggio import ParserPython
from arpeggio import PTNodeVisitor, visit_parse_tree
from atomClass import AtomWithContext, Atom, Bond, Context

# From Arpeggio tutorial, prevents string from appearing in parse tree
class SuppressStr(StrMatch):
    suppress = True

def specific_aromatic_atom():
    return (_(r'br|cl|[bcfhinops]'))

def specific_aliphatic_atom():
    return (_(r'Br|Cl|[BCFHINOPS]'))

def atom():
    return [specific_aliphatic_atom, 
            specific_aromatic_atom, 
            wildCard, 
            aromatic, 
            aliphatic]

def charge():
    return [(_(r'\+\d*')), (_(r'-\d*'))]

def bond():
    return ['-', '=', '#', ':', '~']

def wildCard():
    return '*'

def aromatic():
    return 'a'

def aliphatic():
    return 'A'

def degree():
    return (SuppressStr("D"), _(r'\d*'))

def totalHCount():
    return _(r'H\d*')

def implicitHCount():
    return (SuppressStr("h"), _(r'\d*'))

def valence():
    return (SuppressStr("v"), _(r'\d*'))

def connectivity():
    return (SuppressStr("X"), _(r'\d*'))

def ringConnectivity():
    return (SuppressStr("x"), _(r'\d*'))

def ringMembership():
    return (SuppressStr("R"), _(r'\d*'))

def ringSize():
    return (SuppressStr("r"), _(r'\d*'))

def atomContext():
    return [charge, degree, totalHCount, implicitHCount,
    valence, connectivity, ringConnectivity, ringMembership, ringSize]

def ringClosure():
    return (_(r'\d*'))

def logicalOperator():
    return [_(r'&'), _(r','), _(r';')]

def not_op():
    return (_(r'!'))

def atomWithContext():
    return (Optional(not_op), atom, ZeroOrMore(atomContext))

def bracketedAtom():
    return [(SuppressStr("["), atomWithContext, ZeroOrMore(logicalOperator, [atomWithContext, ZeroOrMore(atomContext)]), SuppressStr("]")), atom]

def smilesBranch():
    return (SuppressStr('('), ZeroOrMore(smilesAtomUnit), SuppressStr(')'))

def smartsBranch():
    return (SuppressStr('('), ZeroOrMore(smartsAtomUnit), SuppressStr(')'))

def smilesAtomUnit():
    return [smilesBranch,
            (Optional(bond), 
            ZeroOrMore(bracketedAtom), 
            Optional(ringClosure))]

def smartsAtomUnit():
    return [smartsBranch,
            Optional(bond), 
            bracketedAtom, 
            Optional(ringClosure)]

def smarts():
    return ZeroOrMore(smartsAtomUnit), EOF

def smiles():
    return atom, ZeroOrMore(smilesAtomUnit), EOF

# Visitor class that walks parse tree
# Returns list of fragments,
# Each fragment is tuple of 2 lists,
# First list is starting branch,
# Second list is ring
# Any atomList in parenthesis is its own list
# Any atomUnit is tuple of Atom string and int tag
class TreeVisitor(PTNodeVisitor):
    def visit_smarts(self, node, children):
        return list(children)

    def smartsAtomUnit(self, node, children):
        return children

    def visit_bracketedAtom(self, node, children):
        return children
    
    def visit_wildcard(self, node, children):
        return Atom({"is_wildcard": True})

    def visit_aromatic(self, node, children):
        return Atom(**{"is_wildcard": True, "is_aromatic": True})

    def visit_aliphatic(self, node, children):
        return Atom(**{"is_wildcard": True, "is_aromatic": False})

    def visit_specific_aliphatic_atom(self, node, children):
        return Atom(**{"is_wildcard": False, "is_aromatic": False, "atomic_symbol": node})

    def visit_connectivity(self, node, children):
        return int(children[0])

    def visit_atom(self, node, children):
        return children[0]

    def visit_atomWithContext(self, node, children):
        atom = children[0]
        context = Context()
        for i in range(1, len(children)):
            child = children[i]
            rule = node[i]
            if rule == "atomContext":
                contextRule = child.rule
                if contextRule == "connectivity":
                    context.set_connectivity(child)
        return AtomWithContext(atom, context)

def initParser():
    """
    Set parser global variable to arpeggio parser for our format
    """
    global parser
    #smilesParser = ParserPython(smiles)
    parser = ParserPython(smarts)

def parseSMARTS(s):
    initParser()
    tree = parser.parse(s)
    #listForm = visit_parse_tree(tree, AtomCloneVistor())

def print_list_units(result):
    for elem in result:
        if type(elem) != AtomWithContext:
            print_list_units(elem)
        else:
            print(elem)

f = open("smartsTests.txt").read()  
initParser()
for line in f.splitlines():
    print(f"\n{line}")
    tree = parser.parse(line)
    print(tree.tree_str(), "\n")
    result = visit_parse_tree(tree, TreeVisitor())
    print_list_units(result)

