from CoolProp.CoolProp import PropsSI

class InputError(Exception):
    pass

class Chemical:
    
    def __init__(self, name):
        self.name = name
        self.mol_weight = None
        self.cp = None
        self.cv = None
        self.R = None
        self.gamma = None

        self.T = None
        self.rho = None
        self.P = None
        self.s = None
        self.h = None

    def set_const(self,mol_weight, cp, cv):
        self.mol_weight = mol_weight
        self.cp = cp
        self.cv = cv
        self.R = cp - cv
        self.gamma = cp / cv

    def set_mol_weight(self,mol_weight):
        self.mol_weight = mol_weight

    def set_cp(self,cp):
        self.cp = cp

    def get_mol_weight(self):
        if (self.mol_weight == None): raise ImportError("Molecular Weight has not been defined yet!")
        return self.mol_weight
    
    def get_cp(self):
        if (self.cp == None): raise ImportError("Isobaric specific heat capacity has not been defined yet!")
        return self.cp


class CoolProp_Chemical(Chemical):

    def __init__(self, name, T, rho):
        self.name = name # must be a valid coolprop chemical name!
        self.mol_weight = PropsSI(self.name,"molemass")/1000
        self.cp = PropsSI("C", "T", T, "D", rho, self.name)
        self.cv = PropsSI("O", "T", T, "D", rho, self.name)
        self.R = self.cp - self.cv

        self.T = T
        self.rho = rho
        self.P = PropsSI("P", "T", T, "D", rho, self.name)
        self.s = PropsSI("S", "T", T, "D", rho, self.name)
        self.h = PropsSI("H", "T", T, "D", rho, self.name)