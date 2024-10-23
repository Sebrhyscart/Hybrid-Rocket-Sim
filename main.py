from chemical import *
from chemical_reaction import *
from injector import *
from combustion import *
from nozzle import *
from printout import *

def main():
    print_title()

    N2O = Chemical("N2O")
    N2O.set_mol_weight(34)

    N2 = Chemical("N2")
    N2.set_mol_weight(28)

    O2 = Chemical("O2")
    O2.set_mol_weight(32)

    paraffin = Chemical("paraffin")
    paraffin.set_mol_weight(12*32+66)

    H2O = Chemical("H2O")
    H2O.set_mol_weight(18)

    CO2 = Chemical("CO2")
    CO2.set_mol_weight(44)


    decomposition = Chemical_Reaction([N2O], [N2, O2], [1,1,1/2], 1)

    combustion = Chemical_Reaction([paraffin, O2], [H2O, CO2], [1,97/2,33,32], 1)

    aNozzle = Nozzle(0.01,0.06)
    aNozzle.set(1.35,400,1e5,1e6,3500)
    print(aNozzle.get_thrust_Isp())

main()