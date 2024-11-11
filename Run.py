from Chemical import *
from ChemicalReaction import *
from Tank import *
from Injector import *
from CombustionChamber import *
from Nozzle import *
from printout import *
from math import isnan

class Run:
    def __init__(self, aChemicalSet:ChemicalSet, aChemicalReactionSet:ChemicalReactionSet, aTank:Tank, aInjector:Injector, aCombustionChamber:CombustionChamber, aNozzle:Nozzle):
        self.aChemicalSet = aChemicalSet
        self.aChemicalReactionSet = aChemicalReactionSet
        self.aTank = aTank
        self.aInjector = aInjector
        self.aCombustionChamber = aCombustionChamber
        self.aNozzle = aNozzle

    def run(self, dt:float=1e-3, endtime:float=1000, eps:float=1e-3, max_iter:int=10, VERBOSE=False, THRUST_ISP=False, CF_CSTAR=False):
        # Implicit - Euler time integration scheme

        print("Run Mode:\nVERBOSE =", VERBOSE, "\nTHRUST_ISP =", THRUST_ISP, "\nCF_CSTAR =", CF_CSTAR)

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

        if VERBOSE: print("-----------------------------------------------------------------------------------------")
        if VERBOSE: print("timestep i =", round(i,6), "time = ", round(t,6), "[s]")

        # Requested tallies
        t_list = [t]
        P_list = [P]
        T_list = [T]
        if THRUST_ISP: 
            thrust_list = [0]
            Isp_list = [0]
        if CF_CSTAR:
            c_f_list = [0]
            c_star_list = [0]

        while((self.aCombustionChamber.r_port < self.aCombustionChamber.r_wall) and (self.aTank.m_oxidizer > 0) and t < endtime):                                                         # While we're not out of fuel or oxidizer ...
            
            i += 1
            t += dt

            if VERBOSE: print("-----------------------------------------------------------------------------------------")
            if VERBOSE: print("timestep i =", round(i,6), "time = ", round(t,6), "[s]")

            if VERBOSE: print("m_dot N2O into engine =",round(self.aInjector.get_flowrate(),3),"[kg/s]")

            self.aTank.set_m_oxidizer(self.aTank.m_oxidizer - self.aInjector.get_flowrate()*dt)             # remove some oxidizer from the tank
            self.aChemicalSet.set_chemical_massflow(N2O,self.aInjector.get_flowrate())                                    # set flowrate entering the combustion chamber this timestep
            m_dot_N2O_, (m_dot_N2, m_dot_O2), Q_dot_decomp = decomposition.complete_reaction([self.aInjector.get_flowrate()]) # N2O -> N2 + 1/2 O2
            m_dot_N2O = m_dot_N2O_
            self.aChemicalSet.set_multiple_chemical_massflow([N2O, N2, O2],[m_dot_N2O, m_dot_N2, m_dot_O2])          # add flowrates to chemical set

            r_dot = self.aCombustionChamber.regression_rate(self.aInjector.get_flowrate() / (np.pi * self.aCombustionChamber.r_port**2))          # regression rate
            m_dot_paraffin = paraffin.get_solid_density() * r_dot * 2 * np.pi * self.aCombustionChamber.r_port * self.aCombustionChamber.l_fuel                 # flow rate of paraffin liquid + vapor from regression
            self.aChemicalSet.set_chemical_massflow(paraffin,m_dot_paraffin)                                         # add flowrate to chemical set
            _, _, Q_dot_vap = vaporization.complete_reaction([m_dot_paraffin])                                    # rate latient heat is consumed by vaporizing paraffin
            
            self.aCombustionChamber.set_r_port(self.aCombustionChamber.r_port + r_dot*dt)
            V = self.aCombustionChamber.get_volume()

            (m_dot_paraffin, m_dot_O2), (m_dot_H2O, m_dot_CO2), Q_dot_comb = combustion.complete_reaction([m_dot_paraffin, m_dot_O2]) # C32H66 + 97/2 O2 -> 33 H2O + 32 CO2
            self.aChemicalSet.set_multiple_chemical_massflow([paraffin, O2, H2O, CO2],[m_dot_paraffin, m_dot_O2, m_dot_H2O, m_dot_CO2]) # add flowrate to chemical set

            m_total = self.aChemicalSet.get_total_mass()
            cp_ave = self.aChemicalSet.get_ave_cp()
            R_ave = self.aChemicalSet.get_ave_R()
            gamma_ave = self.aChemicalSet.get_ave_gamma()
            self.aNozzle.set_gamma_R(gamma_ave, R_ave)

            P_res = 10*eps # initilize to be > eps
            T_res = 10*eps # initilize to be > eps
            m_dot_out_res = 10*eps # initilize to be > eps


            j = 0
            while (P_res >= eps and T_res >= eps and m_dot_out_res >= eps):
                if VERBOSE: print("\titeration step j =", j)
                if VERBOSE: print("\tP =", P, "T =", T, "m_dot_out =",m_dot_out)
                T_jm1 = T # T_j-1
                P_jm1 = P # P_j-1
                m_dot_out_jm1 = m_dot_out # m_dot_out_j-1

                T = (T_list[i-1] + (Q_dot_decomp + Q_dot_vap + Q_dot_comb)*dt / (cp_ave*m_total)) / (1 + m_dot_out*dt / m_total) # Conservation of Energy
                P = 0 
                for k in range(self.aChemicalSet.len()):                                                             # Conservation of Mass
                    # FIX THIS SECTION
                    # mass_k = self.aChemicalSet.get_chemical_mass_by_index(k) + self.aChemicalSet.get_chemical_massflow_by_index(k) * dt + (self.aChemicalSet.get_chemical_mass_by_index(k)/m_total) * m_dot_out * dt
                    # self.aChemicalSet.set_chemical_mass(self.aChemicalSet.get_chemical_by_index(k),mass_k) 
                    # P += mass_k * self.aChemicalSet.get_chemical_by_index(k).get_R() * T / V
                m_dot_out = self.aNozzle.get_flowrate_from_stagnation(P, T)
                
                P_res = abs((P - P_jm1)/P)
                print("P_res", P_res)
                T_res = abs((T - T_jm1)/T)
                m_dot_out_res = abs((m_dot_out - m_dot_out_jm1)/m_dot_out)
                j += 1

                if (j > max_iter):
                    warnings.warn("Warning: Iterator exceeded maximum number of iteration steps", UserWarning)
                    break
            


            if THRUST_ISP: thrust, Isp = self.aNozzle.get_thrust_Isp()
            if CF_CSTAR: c_f, c_star = self.aNozzle.get_thrust_coeff_characteristic_velocity()

            # Requested tallies
            t_list.append(t)
            P_list.append(P)
            T_list.append(T)
            if THRUST_ISP: 
                thrust, Isp = self.aNozzle.get_thrust_Isp()
                thrust_list.append(thrust)
                Isp_list.append(Isp)
            if CF_CSTAR:
                c_f, c_star = self.aNozzle.get_thrust_coeff_characteristic_velocity()
                c_f_list.append(c_f)
                c_star_list.append(c_star)

            if i>10: break # <--------------------------------------------------------------------------------- REMOVE
        
        if VERBOSE: print("Number of timesteps =", i, "Final time =", t)