from Run import *
import matplotlib.pyplot as plt

## Base SI units are used for all quantities throughout
## Time - [s]
## Mass - [kg]
## Length - [m]
## Temperature - [k]
## Force - [N]
## Pressure - [Pa]
## Density - [kg/m3]
## Molecular weight - [g/mol]
## Entropy or heat capcity - [J/kg/k]

n2 = ChemicalSet()

def main():
    print_title()

    # ============================================================================================================
    # MATERIALS
    # ============================================================================================================
    N2O = Chemical("N2O")
    N2O.set_mol_weight(34)
    N2O.set_cp(0.88e3)
    N2O.set_cv(0.69e3)

    N2 = Chemical("N2")
    N2.set_mol_weight(28)
    N2.set_cp(1.04e3)
    N2.set_cv(0.743e3)

    O2 = Chemical("O2")
    O2.set_mol_weight(32)
    O2.set_cp(0.919e3)
    O2.set_cv(0.659e3)    

    paraffin = Chemical("paraffin")
    paraffin.set_mol_weight(12*32+66)
    paraffin.set_cp(1.67e3) # value for propane not paraffin vapor lol, it probably doesnt matter that much
    paraffin.set_cv(1.48e3) # value for propane not paraffin vapor lol, it probably doesnt matter that much 
    paraffin.set_solid_density(900)

    H2O = Chemical("H2O")
    H2O.set_mol_weight(18)
    H2O.set_cp(2.26e3)
    H2O.set_cv(1.76e3) 

    CO2 = Chemical("CO2")
    CO2.set_mol_weight(44)
    CO2.set_cp(0.844e3)
    CO2.set_cv(0.655e3) 

    decomposition = ChemicalReaction("decomposition", [N2O], [N2, O2], [1,1,1/2], 0)
    vaporization = ChemicalReaction("vaporization", [paraffin],[paraffin],[1,1], (-1.46e5-1.6e6)*paraffin.mol_weight/1000)
    combustion = ChemicalReaction("combustion", [paraffin, O2], [H2O, CO2], [1,97/2,33,32], 19.87e6)

    aChemicalSet = ChemicalSet()
    aChemicalSet.append_chemical_list([N2O, N2, O2, paraffin, H2O, CO2])

    aChemicalReactionSet = ChemicalReactionSet()
    aChemicalReactionSet.append_chemical_reaction_list([decomposition,vaporization,combustion])

    # ============================================================================================================
    # GEOMETRY
    # ============================================================================================================
    combustion_chamber_outer_wall_radius = 0.06
    combustion_chamber_port_radius = 0.02
    combustion_chamber_pre_fuel_length = 0.05
    combustion_chamber_fuel_length = 0.2
    combustion_chamber_post_fuel_length = 0.05
    nozzle_throat_radius = 0.01
    nozzle_exit_radius = 0.02
    m_oxidizer = 20
    m_dot_injector = 0.2

    aTank = Tank(m_oxidizer)
    aInjector = Injector(m_dot_injector)
    aCombustionChamber = CombustionChamber(combustion_chamber_outer_wall_radius, combustion_chamber_port_radius, combustion_chamber_fuel_length,
                                           combustion_chamber_pre_fuel_length, combustion_chamber_post_fuel_length)
    aNozzle = Nozzle((np.pi*nozzle_throat_radius**2), (np.pi*nozzle_exit_radius**2))

    # ============================================================================================================
    # RUN
    # ============================================================================================================

    aRun = Run(aChemicalSet, aChemicalReactionSet, aTank, aInjector, aCombustionChamber, aNozzle)

    aRun.run(VERBOSE=True, dt=1e-3, endtime=0.2)

    # r_throat_list = [0.0025,0.005,0.0075,0.01,0.0125,0.015,0.0175,0.02]
    # for r_throat in r_throat_list:
    #     filename = "rt=" + str(r_throat)

    #     aRun.aNozzle.A_throat = np.pi * r_throat**2
    #     aRun.run(output_file_name=filename, dt=1e-3, endtime=5)
    #     print("Done run ", filename)

main()