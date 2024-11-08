from chemical import *
from ChemicalReaction import *
from injector import *
from CombustionChamber import *
from nozzle import *
from printout import *

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
    # ============================================================================================================
    # FLAGS
    # ============================================================================================================

    # ============================================================================================================
    # MATERIALS
    # ============================================================================================================
    N2O = Chemical("N2O")
    N2O.set_mol_weight(34)

    N2 = Chemical("N2")
    N2.set_mol_weight(28)

    O2 = Chemical("O2")
    O2.set_mol_weight(32)

    paraffin = Chemical("paraffin")
    paraffin.set_mol_weight(12*32+66)
    paraffin.set_solid_density(1)

    H2O = Chemical("H2O")
    H2O.set_mol_weight(18)

    CO2 = Chemical("CO2")
    CO2.set_mol_weight(44)

    aChemicalSet = ChemicalSet()
    aChemicalSet.append_chemical_list([N2O, N2, O2, paraffin, H2O, CO2])

    decomposition = ChemicalReaction([N2O], [N2, O2], [1,1,1/2], 1)
    vaporation = ChemicalReaction([paraffin],[paraffin],[1,1], -1)
    combustion = ChemicalReaction([paraffin, O2], [H2O, CO2], [1,97/2,33,32], 1)

    # ============================================================================================================
    # INITIAL GEOMETRY
    # ============================================================================================================
    combustion_chamber_outer_wall_radius = 0.06
    combustion_chamber_port_radius = 0.02
    combustion_chamber_pre_fuel_length = 0.05
    combustion_chamber_fuel_length = 0.2
    combustion_chamber_post_fuel_length = 0.05
    nozzle_throat_radius = 0.02
    nozzle_exit_radius = 0.06

    aInjector = Injector()
    aCombustionChamber = CombustionChamber(combustion_chamber_outer_wall_radius, combustion_chamber_port_radius, combustion_chamber_fuel_length,
                                           combustion_chamber_pre_fuel_length, combustion_chamber_post_fuel_length)
    aNozzle = Nozzle((np.pi*nozzle_throat_radius**2), (nozzle_exit_radius**2 / nozzle_throat_radius**2))
    aNozzle.set_back_pressure(101325)

    # ============================================================================================================
    # TIME INTEGRATOR
    # ============================================================================================================
    # Implicit - Euler time integration scheme

    dt = 0.001
    i = 0

    # Lists of values to record
    time = []
    P_c = []
    T_c = []
    thrust = []
    Isp = []
    c_F = []
    c_star = []

    # Oxidizer and oxidizer flow rate
    m_oxidizer = 20
    injector_flow_rate = 1
    aInjector.set_flowrate(injector_flow_rate)

    # Geometry variables
    r_port = aCombustionChamber.get_r_port()
    r_wall = combustion_chamber_outer_wall_radius
    l_fuel = combustion_chamber_fuel_length

    m_dot_out = 0
    vol = (np.pi * r_port**2 * l_fuel) + (np.pi * r_wall**2 * (combustion_chamber_pre_fuel_length + combustion_chamber_post_fuel_length))
    aChemicalSet.set_atmospheric_conditions(vol)




    while((r_port < r_wall) and (m_oxidizer >= 0)):                                                         # While we're not out of fuel or oxidizer ...

        # N2O enters and decomposes
        aChemicalSet.set_chemical_massflow(N2O,aInjector.get_flowrate())                                    # set flowrate entering the combustion chamber this timestep
        m_dot_N2O, (m_dot_N2, m_dot_O2), Q_dot_decomp = decomposition.complete_reaction([m_dot_N2O])        # N2O -> N2 + 1/2 O2
        aChemicalSet.set_multiple_chemical_massflow([N2O, N2, O2],[m_dot_N2O, m_dot_N2, m_dot_O2])          # add flowrates to chemical set

        # regression rate and fuel evaporates
        r_dot = aCombustionChamber.regression_rate(aInjector.get_flowrate() / (np.pi * r_port**2))          # regression rate
        m_dot_paraffin = paraffin.get_solid_density() * r_dot * 2 * np.pi * r_port * l_fuel                 # flow rate of paraffin liquid + vapor from regression
        aChemicalSet.set_chemical_massflow(paraffin,m_dot_paraffin)                                         # add flowrate to chemical set
        _, _, Q_dot_vap = vaporation.complete_reaction([m_dot_paraffin])                                    # rate latient heat is consumed by vaporizing paraffin

        # combustion reaction
        (m_dot_paraffin, m_dot_O2), (m_dot_H2O, m_dot_CO2), Q_dot_comb = combustion.complete_reaction([m_dot_paraffin, m_dot_O2]) # C32H66 + 97/2 O2 -> 33 H2O + 32 CO2
        aChemicalSet.set_multiple_chemical_massflow([paraffin, O2, H2O, CO2],[m_dot_paraffin, m_dot_O2, m_dot_H2O, m_dot_CO2]) # add flowrate to chemical set

        # average properties and mass of chemicals
        m_total = aChemicalSet.get_total_mass()
        cp_ave = aChemicalSet.get_ave_cp()
        R_ave = aChemicalSet.get_ave_R()
        gamma_ave = aChemicalSet.get_ave_gamma()





        # Change in volume
        vol = vol - 2 * np.pi * r_port * l_fuel * r_dot * dt

        # Conservation of Energy
        Q_dot_out = m_dot_out * (cp_ave) * temp ### use of previous timestep value instead of current value - fix this later (sub into temp equation and rearrange)
        temp = temp + (Q_dot_decomp + Q_dot_vap + Q_dot_comb - Q_dot_out) * dt / (cp_ave * m_total)

        # Conservation of Mass
        press = 0
        for k in range(aChemicalSet.len()):
            mass_k = aChemicalSet.get_chemical_mass_by_index(k) + aChemicalSet.get_chemical_massflow_by_index(k) * dt + (aChemicalSet.get_chemical_mass_by_index(k)/m_total) * m_dot_out * dt
            aChemicalSet.set_chemical_mass(aChemicalSet.get_chemical_by_index(k),mass_k) 
            press += mass_k * aChemicalSet.get_chemical_by_index(k).get_R() * temp / vol

        
    print_title()
main()