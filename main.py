from Run import *

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

def main():
    print_title()

    # ============================================================================================================
    # MATERIALS
    # ============================================================================================================
    # NOX - oxidizer
    N2O = Chemical("N2O")
    N2O.set_mol_weight(34)
    N2O.set_cp(0.88e3)
    N2O.set_cv(0.69e3)

    # Nitrogen gas - oxidizer by-product
    N2 = Chemical("N2")
    N2.set_mol_weight(28)
    N2.set_cp(1.04e3)
    N2.set_cv(0.743e3)

    # Oxygen gas - Oxidizer after decomp
    O2 = Chemical("O2")
    O2.set_mol_weight(32)
    O2.set_cp(0.919e3)
    O2.set_cv(0.659e3)    

    # paraffin wax - fuel
    paraffin = Chemical("fuel")
    paraffin.set_mol_weight(12*32+66)
    paraffin.set_cp(1.67e3) # value for propane not paraffin vapor lol, it probably doesnt matter that much
    paraffin.set_cv(1.48e3) # value for propane not paraffin vapor lol, it probably doesnt matter that much 
    paraffin.set_solid_density(900)

    # Hydroxyl-terminated Poly Butane - fuel
    HTPB = Chemical("fuel")
    HTPB.set_mol_weight(2700)
    HTPB.set_cp(1.67e3) # value for propane not paraffin vapor lol, it probably doesnt matter that much
    HTPB.set_cv(1.48e3) # value for propane not paraffin vapor lol, it probably doesnt matter that much 
    HTPB.set_solid_density(900)

    # Water vapor - combustion product
    H2O = Chemical("H2O")
    H2O.set_mol_weight(18)
    H2O.set_cp(2.26e3)
    H2O.set_cv(1.76e3) 

    # CO2 - combustion product
    CO2 = Chemical("CO2")
    CO2.set_mol_weight(44)
    CO2.set_cp(0.844e3)
    CO2.set_cv(0.655e3) 

    # NOX Decomposition reaction
    decomposition_nox = ChemicalReaction("decomposition", [N2O], [N2, O2], [1,1,1/2], 0)
    
    # Paraffin vaporization and combustion reactions
    vaporization_paraffin = ChemicalReaction("vaporization", [paraffin],[paraffin],[1,1], (-1.46e6)*paraffin.mol_weight/1000)
    combustion_paraffin = ChemicalReaction("combustion", [paraffin, O2], [H2O, CO2], [1,19,26,25], 42e6)
    
    # HTPB vaporization and combustion reactions
    vaporization_HTPB = ChemicalReaction("vaporization", [HTPB],[HTPB],[1,1], (-1.0e6)*HTPB.mol_weight/1000)
    combustion_HTPB = ChemicalReaction("combustion", [HTPB, O2], [H2O, CO2], [1,2.75,3,4], 46e6)

    # append chemicals to chemical sets
    paraffinChemicalSet = ChemicalSet()
    paraffinChemicalSet.append_chemical_list([N2O, N2, O2, paraffin, H2O, CO2])
    HTPBChemcialSet = ChemicalSet()
    HTPBChemcialSet.append_chemical_list([N2O, N2, O2, HTPB, H2O, CO2])

    # append chemical reactions to chemcial reaction sets
    paraffinChemicalReactionSet = ChemicalReactionSet()
    paraffinChemicalReactionSet.append_chemical_reaction_list([decomposition_nox,vaporization_paraffin,combustion_paraffin])
    HTPBChemicalReactionSet = ChemicalReactionSet()
    HTPBChemicalReactionSet.append_chemical_reaction_list([decomposition_nox,vaporization_HTPB,combustion_HTPB])

    # ============================================================================================================
    # GEOMETRY
    # ============================================================================================================
    combustion_chamber_outer_wall_radius = 0.0857/2
    combustion_chamber_port_radius = 0.0303
    combustion_chamber_pre_fuel_length = 0.05
    combustion_chamber_fuel_length = 0.2446
    combustion_chamber_post_fuel_length = 0.05
    nozzle_throat_radius = inches_to_meters(0.843/2)
    nozzle_exit_radius = inches_to_meters(1.83/2)
    nozzle_diverging_section_length = inches_to_meters(2)
    injector_hole_radius = 0.002
    m_oxidizer = 2.5
    m_dot_injector = 0.2

    macTank = Tank(m_oxidizer)
    macInjector = Injector(1, (np.pi*injector_hole_radius**2))
    macInjector.flowrate = m_dot_injector
    macCombustionChamber = CombustionChamber(combustion_chamber_outer_wall_radius, combustion_chamber_port_radius, combustion_chamber_fuel_length,
                                           combustion_chamber_pre_fuel_length, combustion_chamber_post_fuel_length)
    macNozzle = Nozzle((np.pi*nozzle_throat_radius**2), (np.pi*nozzle_exit_radius**2), nozzle_diverging_section_length)

    combustion_chamber_outer_wall_radius = 0.0857/2
    combustion_chamber_port_radius = 0.0303
    combustion_chamber_pre_fuel_length = 0.05
    combustion_chamber_fuel_length = 0.2446
    combustion_chamber_post_fuel_length = 0.05
    nozzle_throat_radius = inches_to_meters(0.869/2)
    nozzle_exit_radius = inches_to_meters(1.890/2)
    nozzle_diverging_section_length = inches_to_meters(2)
    injector_hole_radius = 0.001473
    m_oxidizer = 2.5
    m_dot_injector = 0.3

    waterlooTank = Tank(m_oxidizer)
    waterlooInjector = Injector(6, (np.pi*injector_hole_radius**2))
    waterlooInjector.flowrate = m_dot_injector
    waterlooCombustionChamber = CombustionChamber(combustion_chamber_outer_wall_radius, combustion_chamber_port_radius, combustion_chamber_fuel_length,
                                           combustion_chamber_pre_fuel_length, combustion_chamber_post_fuel_length)
    waterlooNozzle = Nozzle((np.pi*nozzle_throat_radius**2), (np.pi*nozzle_exit_radius**2), nozzle_diverging_section_length)

    # ============================================================================================================
    # RUN
    # ============================================================================================================

    mcmasterRun = Run(paraffinChemicalSet, paraffinChemicalReactionSet, macTank, macInjector, macCombustionChamber, macNozzle)
    #mcmasterRun.run(PLOT=True, timestep=1e-4, endtime=1e-5, output_name="test", INJECTOR_MODEL='complex')

    waterlooRun = Run(HTPBChemcialSet, HTPBChemicalReactionSet, waterlooTank, waterlooInjector, waterlooCombustionChamber, waterlooNozzle)
    waterlooRun.run(PLOT=True, endtime=2, output_name="waterloo", REGRESSION_MODEL='HTPB', INJECTOR_MODEL='simple')

main()