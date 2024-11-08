from algorithm import *

class Chemical: # class for holding the properties of a chemical

    def __init__(self, name:str):
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
        self.rho_solid = None

    def set_const(self,mol_weight:float, cp:float, cv:float):
        global N_A, k_B
        self.mol_weight = mol_weight
        self.cp = cp
        self.cv = cv
        self.R = (N_A * k_B) / mol_weight
        self.gamma = cp / cv

    def set_mol_weight(self,mol_weight:float):
        global N_A, k_B
        self.mol_weight = mol_weight
        self.R = (N_A * k_B) / mol_weight

    def set_cp(self,cp:float):
        self.cp = cp

    def set_solid_density(self,rho_solid:float):
        self.rho_solid = rho_solid

    def get_name(self):
        return self.name

    def get_mol_weight(self):
        if (self.mol_weight == None): raise ImportError("Molecular Weight has not been defined yet!")
        return self.mol_weight
    
    def get_cp(self):
        if (self.cp == None): raise ImportError("Isobaric specific heat capacity has not been defined yet!")
        return self.cp
    
    def get_solid_density(self):
        if (self.rho_solid == None): raise ImportError("Density of solid phase has not been defined yet!")
        return self.rho_solid   

    def get_R(self):
        if (self.R == None): raise ImportError("Gas constant has not been defined yet!")
        return self.R  

class ChemicalSet: # class for holding the amount (mass) and rate of change of amount (massflow) of many chemicals

    def __init__(self):
        self.species = {}
        self.m_species = {}
        self.m_dot_species = {}

    def append_chemical(self, aChemical:Chemical): # append a new chemical to the set
        self.species[aChemical.get_name()] = aChemical

    def append_chemical_list(self, aChemicalList:list[Chemical]): # append multiple new chemicals to the set
        for k in range(len(aChemicalList)):
            self.species[aChemicalList[k].get_name()] = aChemicalList[k] 
    
    def set_chemical_mass(self, aChemical:Chemical, new_mass): # set the mass of a chemical
        self.m_species[aChemical.get_name()] = new_mass

    def set_multiple_chemical_mass(self, aChemicalList:list[Chemical], new_mass_list:list[float]): # set the masses of multiple chemicals at once
        for k in range(len(aChemicalList)):
            self.m_species[aChemicalList[k].get_name()] = new_mass_list[k]

    def set_chemical_massflow(self, aChemical:Chemical, new_massflow): # set the mass flow rate of a chemical
        self.m_dot_species[aChemical.get_name()] = new_massflow
    
    def set_multiple_chemical_massflow(self, aChemicalList:list[Chemical], new_massflow_list:list[float]): # set the mass flow rates of multiple chemicals at once
        for k in range(len(aChemicalList)):
            self.m_species[aChemicalList[k].get_name()] = new_massflow_list[k]

    def get_chemical_dict(self): # return the dictionary of chemicals 
        return self.species
    
    def get_chemical_mass(self,aChemical:Chemical): # return the mass of this chemical
        return self.m_species[aChemical:Chemical.get_name()]
    
    def get_chemical_mass_by_index(self,k:int): # return the mass of the 'k_th' chemical
        return list(self.m_species.values)[k]

    def get_chemical_massflow(self,aChemical:Chemical): # return the massflow of this chemical
        return self.m_dot_species[aChemical:Chemical.get_name()]

    def get_chemical_massflow_by_index(self,k:int): # return the mass flow rate of the 'k_th' chemical
        return list(self.m_species.values)[k]

    def len(self): # return the length of the dict
        return len(self.species)

    # def set_atmospheric(self,Vol:float):
    #     self.m_N2 = (101325 * Vol) / (N2.get_R() * 300) * 0.8
    #     self.m_O2 = (101325 * Vol) / (O2.get_R() * 300) * 0.2

# class CoolProp_Chemical(Chemical):

#     def __init__(self, name, T, rho):
#         self.name = name # must be a valid coolprop chemical name!
#         self.mol_weight = PropsSI(self.name,"molemass")/1000
#         self.cp = PropsSI("C", "T", T, "D", rho, self.name)
#         self.cv = PropsSI("O", "T", T, "D", rho, self.name)
#         self.R = self.cp - self.cv

#         self.T = T
#         self.rho = rho
#         self.P = PropsSI("P", "T", T, "D", rho, self.name)
#         self.s = PropsSI("S", "T", T, "D", rho, self.name)
#         self.h = PropsSI("H", "T", T, "D", rho, self.name)