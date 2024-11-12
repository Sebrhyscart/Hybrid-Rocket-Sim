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
    nozzle_throat_radius = 0.04
    nozzle_exit_radius = 0.06
    m_oxidizer = 20
    m_dot_injector = 0.1

    aTank = Tank(m_oxidizer)
    aInjector = Injector(m_dot_injector)
    aCombustionChamber = CombustionChamber(combustion_chamber_outer_wall_radius, combustion_chamber_port_radius, combustion_chamber_fuel_length,
                                           combustion_chamber_pre_fuel_length, combustion_chamber_post_fuel_length)
    aNozzle = Nozzle((np.pi*nozzle_throat_radius**2), (nozzle_exit_radius**2 / nozzle_throat_radius**2))

    # ============================================================================================================
    # RUN
    # ============================================================================================================
    
    aRun = Run(aChemicalSet, aChemicalReactionSet, aTank, aInjector, aCombustionChamber, aNozzle)

    aRun.run(VERBOSE = True)













    # eps = 1e-3 # size of residual
    # max_iter = 100 # max number of iterations in iterative methods
    # dt = 1e-3 # time step
    # i = 0

    # # Lists of values to record
    # time = 0
    # times = []
    # P_list = []
    # T_list = []
    # thrust_list = []
    # Isp_list = []

    # # Oxidizer and oxidizer flow rate
    # m_oxidizer = 20
    # injector_flow_rate = 1
    # aInjector.set_flowrate(injector_flow_rate)

    # # Geometry variables
    # r_port = aCombustionChamber.get_r_port()
    # r_wall = combustion_chamber_outer_wall_radius
    # l_fuel = combustion_chamber_fuel_length

    # m_dot_out = 0
    # V = (np.pi * r_port**2 * l_fuel) + (np.pi * r_wall**2 * (combustion_chamber_pre_fuel_length + combustion_chamber_post_fuel_length))
    # T = 300
    # aChemicalSet.set_atmospheric_conditions(V,T)

    # while((r_port < r_wall) and (m_oxidizer >= 0)):                                                         # While we're not out of fuel or oxidizer ...
    #     print("i =", i)
    #     aChemicalSet.set_chemical_massflow(N2O,aInjector.get_flowrate())                                    # set flowrate entering the combustion chamber this timestep
    #     m_dot_N2O_, (m_dot_N2, m_dot_O2), Q_dot_decomp = decomposition.complete_reaction([aInjector.get_flowrate()]) # N2O -> N2 + 1/2 O2
    #     m_dot_N2O = m_dot_N2O_
    #     aChemicalSet.set_multiple_chemical_massflow([N2O, N2, O2],[m_dot_N2O, m_dot_N2, m_dot_O2])          # add flowrates to chemical set
    #     m_oxidizer = m_oxidizer - aInjector.get_flowrate()*dt

    #     r_dot = aCombustionChamber.regression_rate(aInjector.get_flowrate() / (np.pi * r_port**2))          # regression rate
    #     m_dot_paraffin = paraffin.get_solid_density() * r_dot * 2 * np.pi * r_port * l_fuel                 # flow rate of paraffin liquid + vapor from regression
    #     aChemicalSet.set_chemical_massflow(paraffin,m_dot_paraffin)                                         # add flowrate to chemical set
    #     _, _, Q_dot_vap = vaporization.complete_reaction([m_dot_paraffin])                                    # rate latient heat is consumed by vaporizing paraffin
        
    #     V = V + 2 * np.pi * r_port * l_fuel * r_dot * dt
    #     r_port = r_port + r_dot*dt
    #     aCombustionChamber.set_r_port(r_port)

    #     (m_dot_paraffin, m_dot_O2), (m_dot_H2O, m_dot_CO2), Q_dot_comb = combustion.complete_reaction([m_dot_paraffin, m_dot_O2]) # C32H66 + 97/2 O2 -> 33 H2O + 32 CO2
    #     aChemicalSet.set_multiple_chemical_massflow([paraffin, O2, H2O, CO2],[m_dot_paraffin, m_dot_O2, m_dot_H2O, m_dot_CO2]) # add flowrate to chemical set

    #     m_total = aChemicalSet.get_total_mass()
    #     cp_ave = aChemicalSet.get_ave_cp()
    #     R_ave = aChemicalSet.get_ave_R()
    #     gamma_ave = aChemicalSet.get_ave_gamma()
    #     aNozzle.set_gamma_R(gamma_ave, R_ave)

    #     P_res = 10*eps # initilize to be > eps
    #     T_res = 10*eps # initilize to be > eps
    #     m_dot_out_res = 10*eps # initilize to be > eps
    #     j = 0
    #     P = 101325
    #     while (P_res >= eps and T_res >= eps and m_dot_out_res >= eps):
    #         T_prior = T
    #         P_prior = P
    #         m_dot_out_prior = m_dot_out

    #         T = (T_prior + (Q_dot_decomp + Q_dot_vap + Q_dot_comb)*dt / (cp_ave*m_total)) / (1 + m_dot_out*dt / m_total) # Conservation of Energy
    #         P = 0 
    #         for k in range(aChemicalSet.len()):                                                             # Conservation of Mass
    #             mass_k = aChemicalSet.get_chemical_mass_by_index(k) + aChemicalSet.get_chemical_massflow_by_index(k) * dt + (aChemicalSet.get_chemical_mass_by_index(k)/m_total) * m_dot_out * dt
    #             aChemicalSet.set_chemical_mass(aChemicalSet.get_chemical_by_index(k),mass_k) 
    #             P += mass_k * aChemicalSet.get_chemical_by_index(k).get_R() * T / V
    #         m_dot_out = aNozzle.get_flowrate_from_stagnation(P, T)
            
    #         P_res = abs((P - P_prior)/P)
    #         T_res = abs((T - T_prior)/T)
    #         m_dot_out_res = abs((m_dot_out - m_dot_out_prior)/m_dot_out)
    #         j += 1
    #         if (j > max_iter):
    #             warnings.warn("Warning: Iterator exceeded maximum number of iteration steps", UserWarning)
    #             break
        
    #     thrust, Isp = aNozzle.get_thrust_Isp()

    #     i += 1
    #     time += dt
    #     times.append(time)

    #     P_list.append(P)
    #     T_list.append(T)
    #     thrust_list.append(thrust)
    #     Isp_list.append(Isp)

    #     if i>10: break
        
    # print("Number of timesteps:", i)
main()