from chemical import *

class ChemicalReactionError(Exception):
    pass

class Chemical_Reaction:

    def __init__(self, reactants, products, stoich_coeff, Q_reaction):
        if ((len(reactants) + len(products)) != len(stoich_coeff)): raise ChemicalReactionError("Number of chemical reaction coefficient must match number of reactants and products!")
        self.reactants = reactants # array-like of valid coolprop chemical names
        self.products = products # array-like of valid coolprop chemical names
        self.reactant_stoich_coeff = stoich_coeff[:len(reactants)]
        self.product_stoich_coeff = stoich_coeff[len(reactants):]
        self.Q_reaction= Q_reaction

    def complete_reaction(self, reactant_flowrates):
        if (len(self.reactants) != len(reactant_flowrates)): raise ChemicalReactionError("Must specify a reactant flow rate for each reactant!")
        
        reactant_molflowrates = [0] * len(self.reactants)
        stoichmetric_molflowrates = [0] * len(self.reactants)
        for k in range(len(self.reactants)):
            reactant_molflowrates[k] = reactant_flowrates[k] / (self.reactants[k].get_mol_weight())
            stoichmetric_molflowrates[k] = reactant_molflowrates[k] / self.reactant_stoich_coeff[k]

        molflowrate = np.min(stoichmetric_molflowrates) 

        for k in range(len(self.reactants)):
            reactant_molflowrates[k] = reactant_molflowrates[k] - molflowrate * self.reactant_stoich_coeff[k]
            reactant_flowrates[k] = (reactant_molflowrates[k] * self.reactants[k].get_mol_weight())

        product_molflowrates = [0] * len(self.products)
        product_flowrates = [0] * len(self.products)
        for k in range(len(self.products)):
            product_molflowrates[k] = molflowrate * self.product_stoich_coeff[k]
            product_flowrates[k] = (product_molflowrates[k] * self.products[k].get_mol_weight())

        P_total = (self.Q_reaction * molflowrate)
        return reactant_flowrates, product_flowrates, P_total