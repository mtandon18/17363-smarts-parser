class AtomWithContext():
    def __init__(self, atoms, contexts):
        if (type(atoms) != list):
            atoms = [atoms]
            contexts = [contexts]

        self.atoms = atoms  # List of Atom objects, implicit that these are combined with the OR operator
        self.contexts = contexts  # List of corresponding Context objects for the atoms list above

    def __str__(self):
        result = ""
        for i in range(len(self.atoms)):
            atom = self.atoms[i]
            ctx = self.contexts[i]
            result = result + f"ATOM\n{str(atom)}\n\nCONTEXT\n{str(ctx)}\n"
        return result
            

class Atom():
    def __init__(self, **kwargs):
        self.is_aromatic = kwargs.get("is_aromatic", False)
        self.is_wildcard = kwargs.get("is_wildcard", False)
        self.is_anti_pattern = kwargs.get("is_anti_pattern", False) # Used to describe that the substructure should NOT be present
        self.atomic_symbol = kwargs.get("atomic_symbol", None)

        if self.is_wildcard and self.atomic_symbol:
            raise "A wildcard should not have known atomic symbol or number"

    def __str__(self):
        return f"  Aromatic: {self.is_aromatic}\n  Wildcard: {self.is_wildcard}\n  Negated: {self.is_anti_pattern}\n  Atomic Symbol: {self.atomic_symbol}"   

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

    def set_explicit_h_count(self, explicit_h_count):
        if self.explicit_h_count:
            raise "Explicit H count for this atomic unit has already been set"
        else:
            self.explicit_h_count = explicit_h_count

    def set_implicit_h_count(self, implicit_h_count):
        if self.implicit_h_count:
            raise "Implicit H count for this atomic unit has already been set"
        else:
            self.implicit_h_count = implicit_h_count
    
    def set_ring_membership(self, ring_membership):
        if self.ring_membership:
            raise "Ring membership for this atomic unit has already been set"
        else:
            self.ring_membership = ring_membership

    def set_ring_size(self, ring_size):
        if self.ring_size:
            raise "Ring size for this atomic unit has already been set"
        else:
            self.ring_size = ring_size

    def set_valence(self, valence):
        if self.valence:
            raise "Valence for this atomic unit has already been set"
        else:
            self.valence = valence

    def set_connectivity(self, connectivity):
        if self.connectivity:
            raise "Connectivity for this atomic unit has already been set"
        else:
            self.connectivity = connectivity

    def set_ring_connectivity(self, ring_connectivity):
        if self.ring_connectivity:
            raise "Ring connectivity for this atomic unit has already been set"
        else:
            self.ring_connectivity = ring_connectivity

    def set_charge(self, charge):
        if self.charge:
            raise "Charge for this atomic unit has already been set"
        else:
            self.charge = charge

    def __str__(self):
        degree = f"  Degree: {self.degree}"
        exp_h_count = f"  Explicit H Count: {self.explicit_h_count}"
        imp_h_count = f"  Implicit H Count: {self.implicit_h_count}"
        ring_member = f"  Ring Membership: {self.ring_membership}"
        ring_size = f"  Ring Size: {self.ring_size}"
        valence = f"  Valence: {self.valence}"
        connectivity = f"  Ring Connectivity: {self.connectivity}"
        ring_connectivity = f"  Ring Connectivity: {self.ring_connectivity}"
        charge = f"  Charge: {self.charge}"
        return f"{degree}\n{exp_h_count}\n{imp_h_count}\n{ring_member}\n{ring_size}\n{valence}\n{connectivity}\n{ring_connectivity}\n{charge}"
            
class Bond():
    def __init__(self, bond_type):
        self.bond_type = bond_type