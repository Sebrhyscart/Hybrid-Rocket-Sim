from Chemical import *
from ChemicalReaction import *
from Tank import *
from Injector import *
from CombustionChamber import *
from Nozzle import *
from printout import *

class Run:
    def __init__(self, aChemicalSet:ChemicalSet, aChemicalReactionSet:ChemicalReactionSet, aTank:Tank, aInjector:Injector, aCombustionChamber:CombustionChamber, aNozzle:Nozzle):
        self.aChemicalSet = aChemicalSet
        self.aChemicalReactionSet = aChemicalReactionSet
        self.aTank = aTank
        self.aInjector = aInjector
        self.aCombustionChamber = aCombustionChamber
        self.aNozzle = aNozzle

    def run(self, dt:float=1e-3, endtime:float=1000, eps:float=1e-3, max_iter:int=10, VERBOSE=False, THRUST_ISP=False, CF_CSTAR=False):
        '''
        Implicit - Euler based time integration scheme
        '''

        print("Run Mode:\nVERBOSE =", VERBOSE, "\nTHRUST_ISP =", THRUST_ISP, "\nCF_CSTAR =", CF_CSTAR)
        new_file()
        print_file(header)
        print_file("Run Mode:\nVERBOSE =", VERBOSE, "\nTHRUST_ISP =", THRUST_ISP, "\nCF_CSTAR =", CF_CSTAR)

        new_data()
        if (THRUST_ISP == True):
            if (CF_CSTAR == True):
                print_data("t(s)","P(Pa)","T(k)","mdot(kg/s)","F(N)","ISP(s)","CF(ul)","CSTAR(m/s)")
            else:
                print_data("t(s)","P(Pa)","T(k)","mdot(kg/s)","F(N)","ISP(s)")
        else:
            print_data("t(s)","P(Pa)","T(k)","mdot(kg/s)")

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

        t = 0
        i = 0

        if VERBOSE: 
            print("-----------------------------------------------------------------------------------------")
            print("timestep i =", round(i,6), "time = ", round(t,8), "[s]")
            print("\tP =", round(P,1), " | T =", round(T,1), " | m_dot_out =", round(m_dot_out,2))
            print_file("-----------------------------------------------------------------------------------------")
            print_file("timestep i =", round(i,6), "time = ", round(t,8), "[s]")
            print_file("\tP =", round(P,1), " | T =", round(T,1), " | m_dot_out =", round(m_dot_out,2))

        t_list = [t]
        P_list = [P]
        T_list = [T]
        if THRUST_ISP: 
            thrust_list = [0]
            Isp_list = [0]
        if CF_CSTAR:
            c_f_list = [0]
            c_star_list = [0]

        # Requested tallies
        if (THRUST_ISP == True):
            if (CF_CSTAR == True):
                print_data(round(t,8),round(P,2),round(T,2),0,0,0,0,0)
            else:
                print_data(round(t,8),round(P,2),round(T,2),0,0,0)
        else:
            print_data(round(t,8),round(P,2),round(T,2),0)       

        im1ChemicalSet = self.aChemicalSet.copy()

        while((self.aCombustionChamber.r_port < self.aCombustionChamber.r_wall) and (self.aTank.m_oxidizer > 0) and t < endtime): # While we're not out of fuel or oxidizer ...
            
            i += 1
            t += dt

            if VERBOSE: 
                print("-----------------------------------------------------------------------------------------")
                print("timestep i =", round(i,8), "time = ", round(t,8), "[s]\n")
                print("NOX Decomposition Reaction")
                print("\tm_dot N2O into engine =",round(self.aInjector.flowrate,3),"[kg/s]")
                print_file("-----------------------------------------------------------------------------------------")
                print_file("timestep i =", round(i,8), "time = ", round(t,8), "[s]\n")
                print_file("NOX Decomposition Reaction")
                print_file("\tm_dot N2O into engine =",round(self.aInjector.flowrate,3),"[kg/s]")

            self.aTank.set_m_oxidizer(self.aTank.m_oxidizer - self.aInjector.flowrate*dt)             # remove some oxidizer from the tank
            self.aChemicalSet.set_chemical_massflow(N2O,self.aInjector.flowrate)                                    # set flowrate entering the combustion chamber this timestep
            m_dot_N2O_, (m_dot_N2, m_dot_O2), Q_dot_decomp = decomposition.complete_reaction([self.aInjector.flowrate]) # N2O -> N2 + 1/2 O2
            m_dot_N2O = m_dot_N2O_[0]
            self.aChemicalSet.set_multiple_chemical_massflow([N2O, N2, O2],[m_dot_N2O, m_dot_N2, m_dot_O2])          # add flowrates to chemical set

            if VERBOSE: 
                print("\tm_dot N2 from decomp =", round(m_dot_N2,3), "[kg/s]")
                print("\tm_dot O2 from decomp =", round(m_dot_O2,3), "[kg/s]")
                print("\tQ_dot Power from decomp =", round(Q_dot_decomp,3), "[W]\n")
                print_file("\tm_dot N2 from decomp =", round(m_dot_N2,3), "[kg/s]")
                print_file("\tm_dot O2 from decomp =", round(m_dot_O2,3), "[kg/s]")
                print_file("\tQ_dot Power from decomp =", round(Q_dot_decomp,3), "[W]\n")

            r_dot = self.aCombustionChamber.regression_rate(self.aInjector.flowrate / (np.pi * self.aCombustionChamber.r_port**2))          # regression rate
            m_dot_paraffin = paraffin.get_solid_density() * r_dot * 2 * np.pi * self.aCombustionChamber.r_port * self.aCombustionChamber.l_fuel                 # flow rate of paraffin liquid + vapor from regression
            self.aChemicalSet.set_chemical_massflow(paraffin,m_dot_paraffin)                                         # add flowrate to chemical set
            _, _, Q_dot_vap = vaporization.complete_reaction([m_dot_paraffin])                                    # rate latient heat is consumed by vaporizing paraffin
            
            self.aCombustionChamber.set_r_port(self.aCombustionChamber.r_port + r_dot*dt)
            V = self.aCombustionChamber.get_volume()

            if VERBOSE: 
                print("Fuel Grain Regression")
                print("\tRegression rate =", round(r_dot*1000,2), "[mm/s]")
                print("\tm_dot paraffin from regression =", round(m_dot_paraffin,3), "[kg/s]")
                print("\tNew port radius after regression =", round(self.aCombustionChamber.r_port,10), "[m]")
                print("\tNew Chamber Volume after regression =", round(V,6), "[m^3]")
                print("\tQ_dot Power from vaporization =", round(Q_dot_vap,3), "[W]\n")
                print_file("Fuel Grain Regression")
                print_file("\tRegression rate =", round(r_dot*1000,2), "[mm/s]")
                print_file("\tm_dot paraffin from regression =", round(m_dot_paraffin,3), "[kg/s]")
                print_file("\tNew port radius after regression =", round(self.aCombustionChamber.r_port,5), "[m]")
                print_file("\tNew Chamber Volume after regression =", round(V,6), "[m^3]")
                print_file("\tQ_dot Power from vaporization =", round(Q_dot_vap,3), "[W]\n")

            (m_dot_paraffin, m_dot_O2), (m_dot_H2O, m_dot_CO2), Q_dot_comb = combustion.complete_reaction([m_dot_paraffin, m_dot_O2]) # C32H66 + 97/2 O2 -> 33 H2O + 32 CO2
            self.aChemicalSet.set_multiple_chemical_massflow([paraffin, O2, H2O, CO2],[m_dot_paraffin, m_dot_O2, m_dot_H2O, m_dot_CO2]) # add flowrate to chemical set

            if VERBOSE: 
                print("Combustion Reaction")
                print("\tQ_dot Power from combustion =", round(Q_dot_comb,3), "[W]")
                print("\tm_dot O2 surviving combustion =", round(m_dot_O2,3), "[kg/s]")
                print("\tm_dot paraffin surviving combustion =", round(m_dot_paraffin,3), "[kg/s]")
                print("\tm_dot H2O from combustion =", round(m_dot_H2O,3), "[kg/s]")
                print("\tm_dot CO2 from combustion =", round(m_dot_CO2,3), "[kg/s]\n")
                print_file("Combustion Reaction")
                print_file("\tQ_dot Power from combustion =", round(Q_dot_comb,3), "[W]")
                print_file("\tm_dot O2 surviving combustion =", round(m_dot_O2,3), "[kg/s]")
                print_file("\tm_dot paraffin surviving combustion =", round(m_dot_paraffin,3), "[kg/s]")
                print_file("\tm_dot H2O from combustion =", round(m_dot_H2O,3), "[kg/s]")
                print_file("\tm_dot CO2 from combustion =", round(m_dot_CO2,3), "[kg/s]\n")

            m_total = self.aChemicalSet.get_total_mass()
            cp_ave = self.aChemicalSet.get_ave_cp()
            R_ave = self.aChemicalSet.get_ave_R()
            gamma_ave = self.aChemicalSet.get_ave_gamma()
            self.aNozzle.set_gamma(gamma_ave)
            self.aNozzle.set_R(R_ave)

            if VERBOSE: 
                print("Total Chemical Species Info")
                print("\tTotal mass of species in engine =", round(m_total,3), "[kg]")
                print("\tAverage R =", round(R_ave,3), "[J/kg/K]")
                print("\tAverage cp =", round(cp_ave,3), "[J/kg/K]")
                print("\tAverage gamma =", round(gamma_ave,3), "\n")
                print_file("Total Chemical Species Info")
                print_file("\tTotal mass of species in engine =", round(m_total,3), "[kg]")
                print_file("\tAverage R =", round(R_ave,3), "[J/kg/K]")
                print_file("\tAverage cp =", round(cp_ave,3), "[J/kg/K]")
                print_file("\tAverage gamma =", round(gamma_ave,3), "\n")

            T = (T_list[i-1] + (Q_dot_decomp + Q_dot_vap + Q_dot_comb)*dt / (cp_ave*m_total)) / (1 + m_dot_out*dt / m_total) # Conservation of Energy
            P = 0 
            for k in range(self.aChemicalSet.len()):                                                             # Conservation of Mass
                mass_k = im1ChemicalSet.get_chemical_mass_by_index(k) + self.aChemicalSet.get_chemical_massflow_by_index(k) * dt - (im1ChemicalSet.get_chemical_mass_by_index(k)/m_total) * m_dot_out * dt
                if mass_k < 0: mass_k = 0
                print("mass_k: ", mass_k, "[kg]")
                self.aChemicalSet.set_chemical_mass(self.aChemicalSet.get_chemical_by_index(k),mass_k) 
                P += mass_k * self.aChemicalSet.get_chemical_by_index(k).get_R() * T / V
                print("P: ", P, "[Pa]")

            self.aNozzle.calc_flowrate(P, T)
            m_dot_out = self.aNozzle.flowrate

            if VERBOSE: 
                print("Conservation Law Balance")
                print("\tP =", round(P,1))
                print("\tT =", round(T,1))
                print("\tNozzle choked: ", self.aNozzle.is_choked())
                print("\tm_dot_out =", round(m_dot_out,5))
                print_file("Conservation Law Balance")
                print_file("\tP =", round(P,1))
                print_file("\tT =", round(T,1))
                print_file("\tNozzle choked: ", self.aNozzle.is_choked())
                print_file("\tm_dot_out =", round(m_dot_out,5))
            
            im1ChemicalSet = self.aChemicalSet.copy()

            if THRUST_ISP: thrust, Isp = self.aNozzle.get_thrust_Isp()
            if CF_CSTAR: c_f, c_star = self.aNozzle.get_thrust_coeff_characteristic_velocity()

            # Requested tallies
            if (THRUST_ISP == True):
                thrust, Isp = self.aNozzle.get_thrust_Isp()
                thrust_list.append(thrust)
                Isp_list.append(Isp)
                if (CF_CSTAR == True):
                    c_f, c_star = self.aNozzle.get_thrust_coeff_characteristic_velocity()
                    c_f_list.append(c_f)
                    c_star_list.append(c_star)
                    print_data(round(t,8),round(P,2),round(T,2),round(m_dot_out,5),round(thrust,2),round(Isp,2),round(c_f,2),round(c_star,2))
                else:
                    print_data(round(t,8),round(P,2),round(T,2),round(m_dot_out,5),round(thrust,2),round(Isp,2))
            else:
                print_data(round(t,8),round(P,2),round(T,2),round(m_dot_out,5))       
                t_list.append(t)
                P_list.append(P)
                T_list.append(T)
    
        if VERBOSE: 
            print("=================================================================")
            print("Done")
            print("Number of timesteps =", round(i,0), "Final time =", round(t,6))
            print_file("=================================================================")
            print_file("Done")
            print_file("Number of timesteps =", round(i,0), "Final time =", round(t,6))

        return t_list.append(t), P_list.append(P), T_list.append(T)