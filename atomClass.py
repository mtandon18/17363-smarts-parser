class AtomWithContext():
    def __init__(self, atoms, contexts):
        self.atoms = atoms  # List of Atom objects, implicit that these are combined with the OR operator
        self.contexts = contexts  # List of corresponding Context objects for the atoms list above

class Atom():
    def __init__(self, **kwargs):
        self.is_aromatic = kwargs.get("is_aromatic", False)
        self.is_wildcard = kwargs.get("is_wildcard", False)
        self.is_anti_pattern = kwargs.get("is_anti_pattern", False) # Used to describe that the substructure should NOT be present
        self.atomic_number = kwargs.get("atomic_number", None)
        self.atomic_symbol = kwargs.get("atomic_symbol", None)

        if self.is_wildcard and (self.atomic_number or self.atomic_symbol):
            raise "A wildcard should not have known atomic symbol or number"
        
class Context():
    def __init__(self, **kwargs):
        self.degree = kwargs.get("degree", None)
        self.explicit_h_count = kwargs.get("explicit_h_count", None)
        self.implicit_h_count = kwargs.get("implicit_h_count", None)
        self.ring_membership = kwargs.get("ring_membership", None)
        self.ring_size = kwargs.get("ring_size", None)
        self.valence = kwargs.get("valence", None)
        self.connectivity = kwargs.get("connectivity", None)
        self.ring_connectivity = kwargs.get("ring_connectivity", None)
        self.charge = kwargs.get("charge", None)

    # Need this for every parameter, will help with typechecking
    def set_degree(self, degree):
        if self.degree:
            raise "Degree for this atomic unit has already been set"
        else:
            self.degree = degree     

class Bond():
    def __init__(self, bond_type):
        self.bond_type = bond_type