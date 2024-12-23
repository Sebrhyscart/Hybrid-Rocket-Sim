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
    nozzle_exit_radius = 0.014
    diverging_section_nozzle_length = 0.04
    m_oxidizer = 20
    m_dot_injector = 0.2

    aTank = Tank(m_oxidizer)
    aInjector = Injector(m_dot_injector)
    aCombustionChamber = CombustionChamber(combustion_chamber_outer_wall_radius, combustion_chamber_port_radius, combustion_chamber_fuel_length,
                                           combustion_chamber_pre_fuel_length, combustion_chamber_post_fuel_length)
    aNozzle = Nozzle((np.pi*nozzle_throat_radius**2), (np.pi*nozzle_exit_radius**2), diverging_section_nozzle_length)

    # ============================================================================================================
    # RUN
    # ============================================================================================================

    # aRun = Run(aChemicalSet, aChemicalReactionSet, aTank, aInjector, aCombustionChamber, aNozzle)
    # aRun.run(PLOT=True, VERBOSE=False, timestep=1e-4, endtime=4, output_name="short")

    # filename = "short"
    # # Read the header row (8th row, index 7) and data
    # with open(filename+".o", 'r') as file:
    #     headers = file.readlines()[7].strip().split()  # Get column names from the header row

    # # Load the rest of the data, skipping the first 7 rows (header and metadata)
    # data = np.loadtxt(filename+".o", skiprows=8)

    # # Extract the time column
    # time = data[:, 0]  # First column: Time (s)
    # flowrate = data[:, 9] * 6.9e6
    # pressure = data[:, 1]

    # plt.figure(figsize=(10, 5))
    # plt.plot(time,flowrate,label='flowrate')
    # plt.plot(time,pressure,label='pressure')
    # plt.title(f'P_exit vs Time')
    # plt.xlabel('P_0')
    # plt.ylabel('P_exit')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f'P_exit_vs_time.png')
    # plt.close()

    # aNozzle.set_gamma(1.4)
    # aNozzle.set_R(270)
    # aNozzle.set_back_pressure(101325)

    # P_list = np.linspace(101325,15*101325,1000)
    # T_0 = 3000

    # M_throat = []
    # M_exit = []
    # P_exit = []
    # A_star = []
    # for P_0 in P_list:
    #     aNozzle.calc_flowrate(P_0,T_0)
    #     M_throat.append(aNozzle.M_throat)
    #     M_exit.append(aNozzle.M_exit)
    #     P_exit.append(aNozzle.P_exit)
    #     A_star.append(aNozzle.A_star_temp)

    # plt.figure(figsize=(10, 5))
    # plt.plot(P_list, M_throat)
    # plt.title(f'M_throat vs Time')
    # plt.xlabel('P_0')
    # plt.ylabel('M_throat')
    # plt.tight_layout()
    # plt.savefig(f'M_throat_vs_time.png')
    # plt.close()

    # plt.figure(figsize=(10, 5))
    # plt.plot(P_list, M_exit)
    # plt.title(f'M_exit vs Time')
    # plt.xlabel('P_0')
    # plt.ylabel('M_exit')
    # plt.tight_layout()
    # plt.savefig(f'M_exit_vs_time.png')
    # plt.close()

    # plt.figure(figsize=(10, 5))
    # plt.plot(P_list, P_exit)
    # plt.title(f'P_exit vs Time')
    # plt.xlabel('P_0')
    # plt.ylabel('P_exit')
    # plt.tight_layout()
    # plt.savefig(f'P_exit_vs_time.png')
    # plt.close()

    # plt.figure(figsize=(10, 5))
    # plt.plot(P_list, A_star)
    # plt.axhline(y=aNozzle.A_throat, color='r', linestyle='--', linewidth=2)
    # plt.title(f'A_star vs Time')
    # plt.xlabel('P_0')
    # plt.ylabel('A_star')
    # plt.tight_layout()
    # plt.savefig(f'A_star_vs_time.png')
    # plt.close()

main()