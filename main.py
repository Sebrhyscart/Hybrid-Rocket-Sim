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
    vaporization = ChemicalReaction("vaporization", [paraffin],[paraffin],[1,1], 0)
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
    nozzle_throat_radius = 0.005
    nozzle_exit_radius = 0.02
    m_oxidizer = 20
    m_dot_injector = 0.1

    aTank = Tank(m_oxidizer)
    aInjector = Injector(m_dot_injector)
    aCombustionChamber = CombustionChamber(combustion_chamber_outer_wall_radius, combustion_chamber_port_radius, combustion_chamber_fuel_length,
                                           combustion_chamber_pre_fuel_length, combustion_chamber_post_fuel_length)
    aNozzle = Nozzle((np.pi*combustion_chamber_outer_wall_radius**2), (np.pi*nozzle_throat_radius**2), (np.pi*nozzle_exit_radius**2))

    # ============================================================================================================
    # RUN
    # ============================================================================================================

    aRun = Run(aChemicalSet, aChemicalReactionSet, aTank, aInjector, aCombustionChamber, aNozzle)
    t,P,T = aRun.run(VERBOSE = True, THRUST_ISP = False, dt = 1e-3 ,endtime = 1)
    # plt.plot(t,P)

    # aNozzle.set_gamma(1.4)
    # aNozzle.set_R(287.05)
    # aNozzle.set_back_pressure(101325)
    # P_0_list = np.linspace(1.5*101325,3*101325,10)
    # T_0 = 2000

    # m_dot_list = np.zeros(len(P_0_list))
    # for i in range(len(P_0_list)):
    #     P_0 = P_0_list[i]
    #     aNozzle.calc_flowrate(P_0, T_0)
    #     m_dot_list[i] = aNozzle.flowrate

        # print("----------------------------------------------------")
        # print("Pressure Ratio P_0/P_b:",round(P_0_list[i]/101325,2))
        # print("Area Ratio A_e/A_t:",round(aNozzle.A_over_A_star,2))
        # print("Exit Pressure:",round(aNozzle.P_exit,2))
        # print("Exit Temperature:",round(aNozzle.T_exit,2))
        # print("Choked:",aNozzle.choked)
        # print("Flowrate:",round(aNozzle.flowrate,2))
        # print("Exit Mach:",round(aNozzle.M_exit,2))
        # print("Exit Velocity:",round(aNozzle.v_eff,2))

    # plt.style.use('dark_background')
    # plt.figure(figsize=(10, 5))
    # plt.plot(P_0_list, m_dot_list, label='flowrate')
    # plt.title('Flowrate vs Pressure')
    # plt.xlabel('Pressure (Pa)')
    # plt.ylabel('Flowrate')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig('flowrate_vs_pressure.png')  # Save plot as an image
    # plt.show()
        
main()