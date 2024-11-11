from algorithm import *

class Chemical: 
    # class for holding the properties of a chemical

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

    def set_mol_weight(self, mol_weight:float):
        global N_A, k_B
        self.mol_weight = mol_weight
        self.R = (N_A * k_B) / mol_weight * 1000

    def set_cp(self, cp:float):
        self.cp = cp
        if (self.cv != None):
            self.gamma = cp / self.cv
    
    def set_cv(self, cv:float):
        self.cv = cv
        if (self.cp != None):
            self.gamma = self.cp / cv

    def set_solid_density(self,rho_solid:float):
        self.rho_solid = rho_solid

    def get_name(self) -> str:
        return self.name

    def get_mol_weight(self) -> float:
        if (self.mol_weight == None): raise ImportError("Molecular Weight has not been defined yet!")
        return self.mol_weight
    
    def get_cp(self) -> float:
        if (self.cp == None): raise ImportError("Isobaric specific heat capacity has not been defined yet!")
        return self.cp
    
    def get_solid_density(self) -> float:
        if (self.rho_solid == None): raise ImportError("Density of solid phase has not been defined yet!")
        return self.rho_solid   

    def get_R(self) -> float:
        if (self.R == None): raise ImportError("Gas constant has not been defined yet!")
        return self.R  
    
    def get_gamma(self) -> float:
        if (self.gamma == None): raise ImportError("Heat capacity ratio has not been defined yet!")
        return self.gamma

class ChemicalSet: 
    # class for holding the amount (mass) {mols are a dumb unit invented by chemists because they are cowards} and rate of change of amount (massflow) of many chemicals

    def __init__(self):
        self.species: dict[str, Chemical] = {}
        self.m_species: dict[str, float] = {}
        self.m_dot_species: dict[str, float] = {}

    def append_chemical(self, aChemical:Chemical): 
        # append a new chemical to the set
        self.species[aChemical.get_name()] = aChemical

    def append_chemical_list(self, aChemicalList:list[Chemical]): 
        # append multiple new chemicals to the set
        for k in range(len(aChemicalList)):
            self.species[aChemicalList[k].name] = aChemicalList[k] 
    
    def set_chemical_mass(self, aChemical:Chemical, new_mass:float): 
        # set the mass of a chemical
        self.m_species[aChemical.name] = new_mass

    def set_multiple_chemical_mass(self, aChemicalList:list[Chemical], new_mass_list:list[float]): 
        # set the masses of multiple chemicals at once
        for k in range(len(aChemicalList)):
            self.m_species[aChemicalList[k].name] = new_mass_list[k]

    def set_chemical_massflow(self, aChemical:Chemical, new_massflow:float): 
        # set the mass flow rate of a chemical
        self.m_dot_species[aChemical.name] = new_massflow
    
    def set_multiple_chemical_massflow(self, aChemicalList:list[Chemical], new_massflow_list:list[float]): 
        # set the mass flow rates of multiple chemicals at once
        for k in range(len(aChemicalList)):
            self.m_dot_species[aChemicalList[k].name] = new_massflow_list[k]

    def print_chemical_mass(self):
        for key, value in self.m_species.items():
            print(f"Key: {key}, Value: {value}, Value Type: {type(value)}")

    def get_chemical_dict(self) -> dict: 
        # return the dictionary of chemicals 
        return self.species

    def get_chemical_name_by_index(self,k:int) -> str: 
        # return the name of the 'k_th' chemical
        return list(self.species.keys())[k]
    
    def get_chemical_by_index(self,k:int) -> Chemical:
        # return the 'k_th' chemical
        return list(self.species.values())[k]

    def get_chemical_mass(self,aChemical:Chemical) -> float: 
        # return the mass of this chemical
        return self.m_species[aChemical:Chemical.get_name()]
    
    def get_chemical_mass_by_index(self,k:int) -> float: 
        # return the mass of the 'k_th' chemical
        return list(self.m_species.values())[k]

    def get_chemical_massflow(self,aChemical:Chemical) -> float: 
        # return the massflow of this chemical
        return self.m_dot_species[aChemical:Chemical.get_name()]

    def get_chemical_massflow_by_index(self,k:int) -> float: 
        # return the mass flow rate of the 'k_th' chemical
        return list(self.m_species.values())[k]

    def len(self) -> int: 
        # return the number of chemicals in the chemical set
        return len(self.species)
    
    def get_total_mass(self) -> float:
        # return the total mass of all the chemicals in the chemical set
        if (len(self.species) != len(self.m_species)):
            warnings.warn("Warning: some chemical species have undefined masses and have not contributed to the total mass.", UserWarning)
            print("len(self.species)", len(self.species), "len(self.m_species)", len(self.m_species))
        m_total = 0
        for k in range(len(self.m_species)):
            m_total += self.get_chemical_mass_by_index(k)
        return m_total

    def get_ave_cp(self) -> float:
        # return mass-averaged isobaric specific heat capacity
        if (len(self.species) != len(self.m_species)):
            warnings.warn("Warning: some chemical species have undefined masses and have not contributed to the total mass.", UserWarning)
            print("len(self.species)", len(self.species), "len(self.m_species)", len(self.m_species))
        m_total = 0
        cp_total = 0
        for k in range(len(self.species)):
            m_total += self.get_chemical_mass_by_index(k)
            cp_total += self.get_chemical_mass_by_index(k) * self.get_chemical_by_index(k).get_cp()
        return cp_total / m_total
    
    def get_ave_R(self) -> float:
        # return mass-averaged gas constant
        if (len(self.species) != len(self.m_species)):
            warnings.warn("Warning: some chemical species have undefined masses and have not contributed to the total mass.", UserWarning)
            print("len(self.species)", len(self.species), "len(self.m_species)", len(self.m_species))
        m_total = 0
        R_total = 0
        for k in range(len(self.species)):
            m_total += self.get_chemical_mass_by_index(k)
            R_total += self.get_chemical_mass_by_index(k) * self.get_chemical_by_index(k).get_R()
        return R_total / m_total           

    def get_ave_gamma(self) -> float:
        # return mass-averaged specific heat capacity ratio
        if (len(self.species) != len(self.m_species)):
            warnings.warn("Warning: some chemical species have undefined masses and have not contributed to the total mass.", UserWarning)
            print("len(self.species)", len(self.species), "len(self.m_species)", len(self.m_species))
        m_total = 0
        gamma_total = 0
        for k in range(len(self.species)):
            m_total += self.get_chemical_mass_by_index(k)
            gamma_total += self.get_chemical_mass_by_index(k) * self.get_chemical_by_index(k).get_gamma()
        return gamma_total / m_total

    def set_atmospheric_conditions(self, pressure:float, temperature:float, volume:float):
        # set the mass of N2 and O2 assuimg Earth's chemical abundance
        if ("N2" in self.species and "O2" in self.species):
            for k in range(len(self.species)):
                if (self.get_chemical_name_by_index(k) == "N2"):
                    self.set_chemical_mass(self.get_chemical_by_index(k), (0.79 * pressure * volume) / (self.get_chemical_by_index(k).get_R() * temperature))
                elif (self.get_chemical_name_by_index(k) == "O2"):
                    self.set_chemical_mass(self.get_chemical_by_index(k), (0.21 * pressure * volume) / (self.get_chemical_by_index(k).get_R() * temperature))
                else:
                    self.set_chemical_mass(self.get_chemical_by_index(k), 0.0)
        else:
            raise KeyError("No instances of \"N2\" or \"O2\" exist in your ChemicalSet")

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