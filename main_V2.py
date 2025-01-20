from Run import *

import numpy as np
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
    # DECOMPOSITION
    # ============================================================================================================
    
    # Initialize Reactants & Products
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

    # Set Up Decomposition Reaction
    decomposition = ChemicalReaction("decomposition", [N2O], [N2, O2], [1,1,1/2], 0)    

    # ============================================================================================================
    # COMBUSTION
    # ============================================================================================================

    H2O = Chemical("H2O")
    H2O.set_mol_weight(18)
    H2O.set_cp(2.26e3)
    H2O.set_cv(1.76e3) 

    CO2 = Chemical("CO2")
    CO2.set_mol_weight(44)
    CO2.set_cp(0.844e3)
    CO2.set_cv(0.655e3) 


    fuel_type=input("Enter the fuel type ('Paraffin' or 'HTPB'): ")
    # """
    # Configure the fuel based on the user input.

    # Parameters:
    #     fuel_type (str): The type of fuel ('Paraffin' or 'HTPB').

    # Returns:
    #     Chemical: Configured Chemical object for the selected fuel.
    # """
    if fuel_type.lower() == "paraffin":
        fuel = Chemical("fuel")
        fuel.set_mol_weight(345)
        fuel.set_cp(950)
        fuel.set_cv(720)
        fuel.set_solid_density(900)

        # Paraffin Reaction
        vaporization = ChemicalReaction("vaporization", [fuel],[fuel],[1,1], (-1.46e6)*fuel.mol_weight/1000)
        combustion = ChemicalReaction("combustion", [fuel, O2], [H2O, CO2], [1,38/2,25,26], 42e6)

        aChemicalSet = ChemicalSet()
        aChemicalSet.append_chemical_list([N2O, N2, O2, fuel, H2O, CO2])

        aChemicalReactionSet = ChemicalReactionSet()
        aChemicalReactionSet.append_chemical_reaction_list([decomposition, vaporization, combustion])

    elif fuel_type.lower() == "htpb":
        fuel = Chemical("fuel")
        fuel.set_mol_weight(2700)
        fuel.set_cp(1050)
        fuel.set_cv(800)
        fuel.set_solid_density(940)

        vaporization = ChemicalReaction("vaporization", [fuel],[fuel],[1,1], (-1.0e6)*fuel.mol_weight/1000)
        combustion = ChemicalReaction("combustion", [fuel, O2], [H2O, CO2], [1,2.75,4,3], 46e6)

        aChemicalSet = ChemicalSet()
        aChemicalSet.append_chemical_list([N2O, N2, O2, fuel, H2O, CO2])

        aChemicalReactionSet = ChemicalReactionSet()
        aChemicalReactionSet.append_chemical_reaction_list([decomposition, vaporization, combustion])

    else:
        raise ValueError("Invalid fuel type. Please choose 'Paraffin' or 'HTPB'.")
    


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
    oxidizer_tank_radius = inches_to_meters(5/2 - 3/16)
    oxidizer_tank_length = 0.353
    m_oxidizer = 2.5
    m_dot_injector = 0.2

    aTank = Tank((np.pi*oxidizer_tank_length*oxidizer_tank_radius**2),m_oxidizer)
    aInjector = Injector(1, (np.pi*injector_hole_radius**2))
    aInjector.flowrate = m_dot_injector
    aCombustionChamber = CombustionChamber(combustion_chamber_outer_wall_radius, combustion_chamber_port_radius, combustion_chamber_fuel_length,
                                           combustion_chamber_pre_fuel_length, combustion_chamber_post_fuel_length)
    aNozzle = Nozzle((np.pi*nozzle_throat_radius**2), (np.pi*nozzle_exit_radius**2), nozzle_diverging_section_length)



    # ============================================================================================================
    # RUN
    # ============================================================================================================
    aRun = Run(aChemicalSet, aChemicalReactionSet, aTank, aInjector, aCombustionChamber, aNozzle)
    aRun.run(PLOT=True, VERBOSE=False, timestep=1e-4, endtime=1, output_name="SimulationOutputData")

main()




# def run_simulation(fuel_type, plot=True, verbose=False):
#     """
#     Run a single simulation for the given fuel type.

#     Parameters:
#         fuel_type (str): The type of fuel ('Paraffin' or 'HTPB').
#         plot (bool): Whether to generate plots.
#         verbose (bool): Whether to print detailed run information.

#     Returns:
#         dict: Simulation results.
#     """

# def compare_fuels():
#     """
#     Compare Paraffin and HTPB simulations and generate comparison plots.
#     """
#     fuels = ["Paraffin", "HTPB"]
#     results = {}

#     for fuel in fuels:
#         results[fuel] = run_simulation(fuel, plot=False, verbose=False)

#     # Generate comparison plots for pressure, temperature, and thrust
#     time = results["Paraffin"]["time"]

#     # Pressure vs Time
#     plt.figure(figsize=(10, 5))
#     plt.plot(time, results["Paraffin"]["pressure"], label="Paraffin")
#     plt.plot(time, results["HTPB"]["pressure"], label="HTPB")
#     plt.title('Chamber Pressure vs Time for Different Fuels')
#     plt.xlabel('Time (s)')
#     plt.ylabel('Pressure (Pa)')
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig('comparison_pressure_vs_time.png')
#     plt.close()

#     # Temperature vs Time
#     plt.figure(figsize=(10, 5))
#     plt.plot(time, results["Paraffin"]["temperature"], label="Paraffin")
#     plt.plot(time, results["HTPB"]["temperature"], label="HTPB")
#     plt.title('Chamber Temperature vs Time for Different Fuels')
#     plt.xlabel('Time (s)')
#     plt.ylabel('Temperature (K)')
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig('comparison_temperature_vs_time.png')
#     plt.close()

#     # Thrust vs Time
#     plt.figure(figsize=(10, 5))
#     plt.plot(time, results["Paraffin"]["thrust"], label="Paraffin")
#     plt.plot(time, results["HTPB"]["thrust"], label="HTPB")
#     plt.title('Thrust vs Time for Different Fuels')
#     plt.xlabel('Time (s)')
#     plt.ylabel('Thrust (N)')
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig('comparison_thrust_vs_time.png')
#     plt.close()

#     print("Comparison plots generated: 'comparison_pressure_vs_time.png', 'comparison_temperature_vs_time.png', 'comparison_thrust_vs_time.png'")



# def main():
#     print("Welcome to the Hybrid Rocket Simulation!")
#     print("1. Single fuel simulation")
#     print("2. Compare Paraffin and HTPB")
#     choice = input("Enter your choice (1 or 2): ")

#     if choice == "1":
#         fuel = input("Enter the fuel type ('Paraffin' or 'HTPB'): ").capitalize()
#         results = run_simulation(fuel)
#         print(f"Simulation complete for {fuel}.")
#     elif choice == "2":
#         compare_fuels()
#         print("Fuel comparison complete.")
#     else:
#         print("Invalid choice. Exiting.")
# main()