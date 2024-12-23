import matplotlib.pyplot as plt
import os
import shutil
from Chemical import *
from ChemicalReaction import *
from Tank import *
from Injector import *
from CombustionChamber import *
from Nozzle import *
from printout import *

class Run:
    def __init__(self, aChemicalSet:ChemicalSet, aChemicalReactionSet:ChemicalReactionSet, aTank:Tank, aInjector:Injector, aCombustionChamber:CombustionChamber, aNozzle:Nozzle):
        '''
        Run class brings allows all of the other classes to talk to each other and holds the main algorithm for the code

        Attributes:
            aChemicalSet: ChemicalSet object to hold chemical data
            aChemicalReactionSet: ChemicalReactionSet object to hold chemical reaction data
            aTank: Oxidizer tank object
            aInjector: Oxidizer injector object
            aCombustionChamber: CombustionChamber object
            aNozzle: Nozzle object
        '''
        self.aChemicalSet = aChemicalSet
        self.aChemicalReactionSet = aChemicalReactionSet
        self.aTank = aTank
        self.aInjector = aInjector
        self.aCombustionChamber = aCombustionChamber
        self.aNozzle = aNozzle

    def run(self, timestep:float=1e-4, endtime:float=100, output_name:str=None, VERBOSE=False, PLOT=False, REGRESSION_MODEL = 'simple', THRUST_ISP=False, CF_CSTAR=False):
        '''
        This method is the main time-stepping loop to run the code. Time dependance is solved using a Backwards-Euler time integration scheme.

        Parameters:
            dt: (optional) time step size. Default: 1e-4 s.
            endtime: (optional) endtime of the simulation. Default: 100 s.
            output_name: (optional) the name of the output files. Default: a name will be chosen for you.
            VERBOSE: (optional) run with verbose outputs. Default: False.
            PLOT: (optional) plot time-dependent results. Default: False.
            REGRESSION_MODEL: (optional) choice of model for fuel regression rate.  Options: simple. Default: simple.
            THRUST_ISP: (optional) calculate the thrust and Isp of the engine. Default: False.
            CF_CSTAR: (optional) calcualte the thrust coefficient and characteristic velocity of the engine. Default: False.
        '''

        print(40*"-","\nRun Mode:\nREGRESSION_MODEL =", REGRESSION_MODEL, "\nVERBOSE =", VERBOSE, "\nPLOT =", PLOT, "\nTHRUST_ISP =", THRUST_ISP, "\nCF_CSTAR =", CF_CSTAR, "\n"+40*"-")
        if VERBOSE:
            new_file(output_name)
            print_file(output_name,header)
            print_file(output_name,"Run Mode:\nVERBOSE =", VERBOSE, "\nTHRUST_ISP =", THRUST_ISP, "\nCF_CSTAR =", CF_CSTAR)
        
        new_data(output_name)
        if VERBOSE or PLOT:
            print_data(output_name,"combustion_chamber_outer_wall_radius =",self.aCombustionChamber.r_wall)
            print_data(output_name,"combustion_chamber_port_radius =",self.aCombustionChamber.r_port)
            print_data(output_name,"combustion_chamber_pre_fuel_length =",self.aCombustionChamber.l_pre)
            print_data(output_name,"combustion_chamber_post_fuel_length =",self.aCombustionChamber.l_post)
            print_data(output_name,"combustion_chamber_fuel_length =",self.aCombustionChamber.l_fuel)
            print_data(output_name,"nozzle_throat_Area =",round(self.aNozzle.A_throat,5))
            print_data(output_name,"nozzle_exit_Area =",round(self.aNozzle.A_exit,5))
            print_data(output_name,
                       "time(s)",
                       "Pressure(Pa)",
                       "Temperature(k)",
                       "mdot_out(kg/s)",
                       "Thrust(N)",
                       "Isp(s)",
                       "c_F",
                       "c*(m/s)",
                       "m_ox(kg)",
                       "mdot_ox(kg/s)",
                       "regression_rate(m/s)",
                       "r_port(m)",
                       "Chamber_Volume(m^3)",
                       "gamma",
                       "R(J/kg/K)",
                       "O/F",
                       "M_exit",
                       "P_exit",
                       "x_shock(m)", 
                       "m_O2(kg)",
                       "m_fuelvapor(kg)",
                       "m_N2(kg)",
                       "m_H2O(kg)",
                       "m_CO2(kg)")
        else:
            if (THRUST_ISP == True):
                if (CF_CSTAR == True):
                    print_data(output_name,"t(s)","P(Pa)","T(k)","mdot(kg/s)","F(N)","Isp(s)","cF(ul)","c*(m/s)")
                else:
                    print_data(output_name,"t(s)","P(Pa)","T(k)","mdot(kg/s)","F(N)","Isp(s)")
            else:
                print_data(output_name,"t(s)","P(Pa)","T(k)","mdot(kg/s)")

        # Extract Chemical data
        N2O = self.aChemicalSet.species["N2O"]
        N2 = self.aChemicalSet.species["N2"]
        O2 = self.aChemicalSet.species["O2"]
        paraffin = self.aChemicalSet.species["paraffin"]
        H2O = self.aChemicalSet.species["H2O"]
        CO2 = self.aChemicalSet.species["CO2"]
        decomposition = self.aChemicalReactionSet.reactions["decomposition"]
        vaporization = self.aChemicalReactionSet.reactions["vaporization"]
        combustion = self.aChemicalReactionSet.reactions["combustion"]

        # Initial Conditions
        T = 300 # [k] initial temperature
        P = 101325 # [Pa] initial pressure
        V = self.aCombustionChamber.get_volume() # [m^3] initial engine volume
        self.aChemicalSet.set_atmospheric_conditions(P,T,V) # set the intial masses of N2 and O2 in the engine
        self.aNozzle.set_back_pressure(P)
        m_dot_out = 0

        dt = timestep
        t = 0
        i = 0

        if VERBOSE: 
            print("-----------------------------------------------------------------------------------------")
            print("timestep i =", round(i,6), "time = ", round(t,8), "[s]")
            print_file(output_name,"-----------------------------------------------------------------------------------------")
            print_file(output_name,"timestep i =", round(i,6), "time = ", round(t,8), "[s]")

        if VERBOSE or PLOT:
            print_data(output_name,
                       round(t,8),
                       round(P,4),
                       round(T,4),
                       0,
                       0,
                       0,
                       0,
                       0,
                       self.aTank.m_oxidizer,
                       0,
                       0,
                       self.aCombustionChamber.r_port,
                       round(V,8),
                       round(self.aChemicalSet.get_ave_gamma(),4),
                       round(self.aChemicalSet.get_ave_R(),4),
                       0,
                       0,
                       P,
                       0,
                       round(self.aChemicalSet.get_chemical_mass(O2),8),
                       round(self.aChemicalSet.get_chemical_mass(paraffin),8),
                       round(self.aChemicalSet.get_chemical_mass(N2),8),
                       round(self.aChemicalSet.get_chemical_mass(H2O),8),
                       round(self.aChemicalSet.get_chemical_mass(CO2),8))
        else:
            if (THRUST_ISP == True):
                if (CF_CSTAR == True):
                    print_data(output_name,round(t,8),round(P,2),round(T,2),0,0,0,0,0)
                else:
                    print_data(output_name,round(t,8),round(P,2),round(T,2),0,0,0)
            else:
                print_data(output_name,round(t,8),round(P,2),round(T,2),0)       

        T_list = [T]
        im1ChemicalSet = self.aChemicalSet.copy()

        while((self.aCombustionChamber.r_port < self.aCombustionChamber.r_wall) and (self.aTank.m_oxidizer > 0) and t < endtime): # While we're not out of fuel or oxidizer ...
            
            # start with very small timesteps
            if (i < 100 and timestep > 1e-6): dt=1e-6
            else: dt = timestep

            # Indexing update
            i += 1
            t += dt

            if VERBOSE: 
                print("-----------------------------------------------------------------------------------------")
                print("timestep i =", round(i,8), "time = ", round(t,8), "[s]\n")
                print("NOX Decomposition Reaction:")
                print("\tm_dot N2O into engine =",round(self.aInjector.flowrate,3),"[kg/s]")
                print_file(output_name,"-----------------------------------------------------------------------------------------")
                print_file(output_name,"timestep i =", round(i,8), "time = ", round(t,8), "[s]\n")
                print_file(output_name,"NOX Decomposition Reaction:")
                print_file(output_name,"\tm_dot N2O into engine =",round(self.aInjector.flowrate,3),"[kg/s]")

            # get oxidizer flow rate from the injector
            m_dot_N2O_injector = self.aInjector.flowrate
            # remove some oxidizer from the tank
            self.aTank.set_m_oxidizer(self.aTank.m_oxidizer - m_dot_N2O_injector*dt)
            # set flowrate entering the combustion chamber this timestep
            self.aChemicalSet.set_chemical_massflow(N2O,m_dot_N2O_injector)
            # N2O -> N2 + 1/2 O2
            m_dot_N2O_, (m_dot_N2, m_dot_O2), Q_dot_decomp = decomposition.complete_reaction([m_dot_N2O_injector])
            m_dot_N2O = m_dot_N2O_[0]
            # add flowrates to chemical set
            self.aChemicalSet.set_multiple_chemical_massflow([N2O, N2, O2],[m_dot_N2O, m_dot_N2, m_dot_O2])
            # Note: m_dot_N2O_injector is the NOX flowrate BEFORE decomposition, and m_dot_N2O is the NOX flowrate AFTER decomposition

            # # remove some oxidizer from the tank
            # self.aTank.set_m_oxidizer(self.aTank.m_oxidizer - self.aInjector.sinusoidal_flowrate(t)*dt)
            # # set flowrate entering the combustion chamber this timestep
            # self.aChemicalSet.set_chemical_massflow(N2O,self.aInjector.sinusoidal_flowrate(t))
            # # N2O -> N2 + 1/2 O2
            # m_dot_N2O_, (m_dot_N2, m_dot_O2), Q_dot_decomp = decomposition.complete_reaction([self.aInjector.sinusoidal_flowrate(t)])
            # m_dot_N2O = m_dot_N2O_[0]
            # # add flowrates to chemical set
            # self.aChemicalSet.set_multiple_chemical_massflow([N2O, N2, O2],[m_dot_N2O, m_dot_N2, m_dot_O2])

            if VERBOSE: 
                print("\tm_dot N2 from decomp =", round(m_dot_N2,3), "[kg/s]")
                print("\tm_dot O2 from decomp =", round(m_dot_O2,3), "[kg/s]")
                print("\tQ_dot Power from decomp =", round(Q_dot_decomp,3), "[W]\n")
                print_file(output_name,"\tm_dot N2 from decomp =", round(m_dot_N2,3), "[kg/s]")
                print_file(output_name,"\tm_dot O2 from decomp =", round(m_dot_O2,3), "[kg/s]")
                print_file(output_name,"\tQ_dot Power from decomp =", round(Q_dot_decomp,3), "[W]\n")

            # regression rate
            r_dot = self.aCombustionChamber.simple_regression_rate(m_dot_N2O_injector / (np.pi * self.aCombustionChamber.r_port**2))
            # flow rate of paraffin liquid + vapor from regression
            m_dot_paraffin_grain = paraffin.get_solid_density() * r_dot * 2 * np.pi * self.aCombustionChamber.r_port * self.aCombustionChamber.l_fuel
            # add flowrate to chemical set
            self.aChemicalSet.set_chemical_massflow(paraffin,m_dot_paraffin_grain)
            # rate latient heat is consumed by vaporizing paraffin
            _, _, Q_dot_vap = vaporization.complete_reaction([m_dot_paraffin_grain])
            # expand the port radius
            self.aCombustionChamber.set_r_port(self.aCombustionChamber.r_port + r_dot*dt)
            # calculate the new chamber volume
            V = self.aCombustionChamber.get_volume()

            if VERBOSE: 
                print("Fuel Grain Regression:")
                print("\tRegression rate =", round(r_dot*1000,2), "[mm/s]")
                print("\tm_dot paraffin from regression =", round(m_dot_paraffin_grain,3), "[kg/s]")
                print("\tNew port radius after regression =", round(self.aCombustionChamber.r_port,10), "[m]")
                print("\tNew Chamber Volume after regression =", round(V,6), "[m^3]")
                print("\tQ_dot Power from vaporization =", round(Q_dot_vap,3), "[W]\n")
                print_file(output_name,"Fuel Grain Regression:")
                print_file(output_name,"\tRegression rate =", round(r_dot*1000,2), "[mm/s]")
                print_file(output_name,"\tm_dot paraffin from regression =", round(m_dot_paraffin_grain,3), "[kg/s]")
                print_file(output_name,"\tNew port radius after regression =", round(self.aCombustionChamber.r_port,10), "[m]")
                print_file(output_name,"\tNew Chamber Volume after regression =", round(V,6), "[m^3]")
                print_file(output_name,"\tQ_dot Power from vaporization =", round(Q_dot_vap,3), "[W]\n")

            # C32H66 + 97/2 O2 -> 33 H2O + 32 CO2
            (m_dot_paraffin, m_dot_O2), (m_dot_H2O, m_dot_CO2), Q_dot_comb = combustion.complete_reaction([m_dot_paraffin_grain, m_dot_O2])
            # add flowrates of reaction products to chemical set
            self.aChemicalSet.set_multiple_chemical_massflow([paraffin, O2, H2O, CO2],[m_dot_paraffin, m_dot_O2, m_dot_H2O, m_dot_CO2])
            # Note: m_dot_paraffin_grain is the flowrate of paraffin from the fuel grain, and m_dot_paraffin is the flowrate of paraffin after combustion

            if VERBOSE: 
                print("Combustion Reaction:")
                print("\tQ_dot Power from combustion =", round(Q_dot_comb,3), "[W]")
                print("\tm_dot O2 surviving combustion =", round(m_dot_O2,3), "[kg/s]")
                print("\tm_dot paraffin surviving combustion =", round(m_dot_paraffin,3), "[kg/s]")
                print("\tm_dot H2O from combustion =", round(m_dot_H2O,3), "[kg/s]")
                print("\tm_dot CO2 from combustion =", round(m_dot_CO2,3), "[kg/s]\n")
                print_file(output_name,"Combustion Reaction:")
                print_file(output_name,"\tQ_dot Power from combustion =", round(Q_dot_comb,3), "[W]")
                print_file(output_name,"\tm_dot O2 surviving combustion =", round(m_dot_O2,3), "[kg/s]")
                print_file(output_name,"\tm_dot paraffin surviving combustion =", round(m_dot_paraffin,3), "[kg/s]")
                print_file(output_name,"\tm_dot H2O from combustion =", round(m_dot_H2O,3), "[kg/s]")
                print_file(output_name,"\tm_dot CO2 from combustion =", round(m_dot_CO2,3), "[kg/s]\n")

            # total mass of all chemical species
            m_total = self.aChemicalSet.get_total_mass()
            # average cp of mixture
            cp_ave = self.aChemicalSet.get_ave_cp()
            # average R of mixture
            R_ave = self.aChemicalSet.get_ave_R()
            # average gamma of mixture
            gamma_ave = self.aChemicalSet.get_ave_gamma()
            # set the gamma and R on the exaust in the nozzle
            self.aNozzle.set_gamma(gamma_ave)
            self.aNozzle.set_R(R_ave)

            if VERBOSE: 
                print("Total Chemical Species Info:")
                print("\tTotal mass of species in engine =", round(m_total,3), "[kg]")
                print("\tAverage R =", round(R_ave,3), "[J/kg/K]")
                print("\tAverage cp =", round(cp_ave,3), "[J/kg/K]")
                print("\tAverage gamma =", round(gamma_ave,3), "\n")
                print_file(output_name,"Total Chemical Species Info:")
                print_file(output_name,"\tTotal mass of species in engine =", round(m_total,3), "[kg]")
                print_file(output_name,"\tAverage R =", round(R_ave,3), "[J/kg/K]")
                print_file(output_name,"\tAverage cp =", round(cp_ave,3), "[J/kg/K]")
                print_file(output_name,"\tAverage gamma =", round(gamma_ave,3), "\n")

            # Solve for new temperature using previous temperature and conservation of energy
            T = (T_list[i-1] + (Q_dot_decomp + Q_dot_vap + Q_dot_comb)*dt / (cp_ave*m_total)) / (1 + m_dot_out*dt / m_total)
            # Recalculate the pressure from the partial pressures of each chemical species present
            P = 0 
            for k in range(self.aChemicalSet.len()):
                # calculate the mass of the k_th chemical species present from conservation of mass
                mass_k = im1ChemicalSet.get_chemical_mass_by_index(k) + self.aChemicalSet.get_chemical_massflow_by_index(k) * dt - (im1ChemicalSet.get_chemical_mass_by_index(k)/m_total) * m_dot_out * dt
                # no negative mass lol
                if mass_k < 0: mass_k = 0
                # save this calculated mass in the chemical set for later use (next timestep)
                self.aChemicalSet.set_chemical_mass(self.aChemicalSet.get_chemical_by_index(k),mass_k)
                # Add the partial pressure of the k_th chemical species (using ideal gas law) to a running pressure total
                P += mass_k * self.aChemicalSet.get_chemical_by_index(k).get_R() * T / V

            # calculate the flowrate out of the engine using the Temperature and Pressure in the chamber
            self.aNozzle.calc_flowrate(P, T)
            # get the flowrate out, to be used in the next timestep
            m_dot_out = self.aNozzle.flowrate
            # create a i-1 chemical set to be used throughout the next timestep
            im1ChemicalSet = self.aChemicalSet.copy()
            # append the engine temperature to a list to be used in the next timestep
            T_list.append(T)

            if VERBOSE: 
                print("Conservation Law Balance:")
                print("\tP =", round(P,1))
                print("\tT =", round(T,1))
                print("\tNozzle choked: ", self.aNozzle.is_choked())
                print("\tm_dot_out =", round(m_dot_out,5))
                print_file(output_name,"Conservation Law Balance:")
                print_file(output_name,"\tP =", round(P,1))
                print_file(output_name,"\tT =", round(T,1))
                print_file(output_name,"\tNozzle choked: ", self.aNozzle.is_choked())
                print_file(output_name,"\tm_dot_out =", round(m_dot_out,5))

            if VERBOSE or PLOT:
                thrust, Isp = self.aNozzle.get_thrust_Isp()
                c_f, c_star = self.aNozzle.get_thrust_coeff_characteristic_velocity()
                print_data(output_name,
                           round(t,8),
                           round(P,4),
                           round(T,4),
                           round(m_dot_out,8),
                           round(thrust,8),
                           round(Isp,4),
                           round(c_f,8),
                           round(c_star,4),
                           round(self.aTank.m_oxidizer,8),
                           round(m_dot_N2O_injector,8),
                           round(r_dot,8),
                           round(self.aCombustionChamber.r_port,8),
                           round(V,8),
                           round(self.aNozzle.gamma,4),
                           round(self.aNozzle.R,4),
                           round(m_dot_N2O_injector/m_dot_paraffin_grain,4),
                           round(self.aNozzle.M_exit,4),
                           round(self.aNozzle.P_exit,8),
                           round(self.aNozzle.x_shock,8),
                           round(self.aChemicalSet.get_chemical_mass(O2),8),
                           round(self.aChemicalSet.get_chemical_mass(paraffin),8),
                           round(self.aChemicalSet.get_chemical_mass(N2),8),
                           round(self.aChemicalSet.get_chemical_mass(H2O),8),
                           round(self.aChemicalSet.get_chemical_mass(CO2),8))
            else:
                if (THRUST_ISP):
                    thrust, Isp = self.aNozzle.get_thrust_Isp()
                    if (CF_CSTAR):
                        c_f, c_star = self.aNozzle.get_thrust_coeff_characteristic_velocity()
                        print_data(output_name,round(t,8),round(P,2),round(T,2),round(m_dot_out,5),round(thrust,2),round(Isp,2),round(c_f,2),round(c_star,2))
                    else:
                        print_data(output_name,round(t,8),round(P,2),round(T,2),round(m_dot_out,5),round(thrust,2),round(Isp,2))
                else:
                    print_data(output_name,round(t,8),round(P,2),round(T,2),round(m_dot_out,5))

            if (i % 1000 == 0): 
                print("\tCurrent Timestep:",i)
                   
        print(40*"-","\nCalculation Complete!")

        if VERBOSE: 
            print("=========================================================================================")
            print("Done!")
            print("Number of timesteps =", round(i,0), "Final time =", round(t,6))
            print_file(output_name,"=========================================================================================")
            print_file(output_name,"Done!")
            print_file(output_name,"Number of timesteps =", round(i,0), "Final time =", round(t,6))

        if PLOT:
            self.plot(output_name)

    def plot(self, filename):
        plt.style.use('dark_background')

        if filename == None: filename = 'data'
        directory = filename+'_plots'

        # Create a folder named 'plots'
        if os.path.exists(directory):
            # Remove all contents if the folder already exists
            shutil.rmtree(directory)
        os.makedirs(directory, exist_ok=True)

        # Read the header row (8th row, index 7) and data
        with open(filename+".o", 'r') as file:
            headers = file.readlines()[7].strip().split()  # Get column names from the header row

        # Load the rest of the data, skipping the first 7 rows (header and metadata)
        data = np.loadtxt(filename+".o", skiprows=8)

        # Extract the time column
        time = data[:, 0]  # First column: Time (s)

        # Generate a plot for each column
        for i in range(1, data.shape[1]):  # Start from 1 to skip the time column
            parameter = data[:, i]
            parameter_name = headers[i]

            # Sanitize the parameter name for saving
            sanitized_name = parameter_name.replace('/', '_over_').replace('*', '_star')

            # Plot the data
            plt.figure(figsize=(10, 5))
            plt.plot(time, parameter)
            plt.title(f'{parameter_name.split("(")[0]} vs Time')
            plt.xlabel('Time (s)')
            plt.ylabel(parameter_name)
            plt.tight_layout()
            plt.savefig(f'{directory}/{sanitized_name.split("(")[0]}_vs_time.png')
            plt.close()
                
        # Generate a combined plot for specific parameters with legend
        combined_parameters = ["m_O2(kg)", "m_fuelvapor(kg)", "m_N2(kg)", "m_H2O(kg)", "m_CO2(kg)"]
        plt.figure(figsize=(10, 5))
        for param_name in combined_parameters:
            if param_name in headers:
                param_index = headers.index(param_name)
                param_data = data[:, param_index]
                plt.plot(time, param_data, label=param_name)

        plt.title("Chemical Species Mass vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Mass of Species (kg)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{directory}/species_mass_vs_time.png')
        plt.close()
        