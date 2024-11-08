from algorithm import *

class Nozzle: # isentropic converging-diverging (DeLaval) nozzle class

    def __init__(self, A_throat: float, A_exit: float):
        self.A_throat = A_throat                    # [m] throat area
        self.A_exit = A_exit                        # [m] exit area
        self.A_e_over_A_star = A_exit / A_throat    # expansion area ratio
        self.gamma = None                           # specific heat capacity ratio
        self.R = None                               # [J/kg/K] gas constant
        self.choked = None                          # boolean set if flow is choked or not
        self.P_b = None                             # [Pa] back pressure
        self.P_0 = None                             # [Pa] stagnation pressure (in chamber)
        self.T_0 = None                             # [k] stagnation temperature (in chamber)
        self.rho_0 = None                           # [kg/m^3] stagnation density (in chamber)
        self.flowrate = None                        # [kg/s] nozzle flow rate
        self.P_star = None                          # [Pa] critical throat pressure
        self.T_star = None                          # [k] critical throat temperature
        self.rho_star = None                        # [kg/m^3] critical throat densityow rate 
        self.P_exit = None                          # [Pa] exit pressure
        self.T_exit = None                          # [k] exit temperature
        self.rho_exit = None                        # [kg/m^3] exit denisty
        self.M_exit = None                          # exit Mach number
        self.speed_of_sound_exit = None             # [m/s] speed of sound at exit
        self.thrust = None                          # [N] engine thrust
        self.c_star = None                          # [m/s] characteristic velocity
        self.c_F = None                             # thrust coefficient
        self.Isp = None                             # [s] specific impulse (ISP)

    def set_gamma(self, gamma:float):
        self.gamma = gamma

    def set_R(self, R:float):
        self.R = R

    def set_gamma_R(self, gamma:float, R:float):
        self.gamma = gamma
        self.R = R

    def set_back_pressure(self, P_b:float):
        self.P_b = P_b

    def get_flowrate_from_stagnation(self, P_0:float, T_0:float) -> float: 
        # this method returns the flowrate through the nozzle given stagnation pressure and temperature
        if (self.gamma != None and self.R != None):
            self.P_0 = P_0
            self.T_0 = T_0
            self.rho_0 = P_0 / (self.R * T_0)

            self.P_star = self.P_0 / ((1 + (self.gamma-1)/2)**((self.gamma)/(self.gamma-1)))
            self.T_star = self.T_0 / (1 + (self.gamma-1)/2)
            self.rho_star = self.rho_0 / ((1 + (self.gamma-1)/2)**(1/(self.gamma-1)))
            chokedflowrate = self.rho_star * self.A_throat * np.sqrt(self.gamma * self.R * self.T_star)

            if (self.P_b != None):
                if (self.P_b < self.P_star):
                    self.choked = True
                    self.flowrate = chokedflowrate
                else:
                    self.choked = False
                    self.flowrate = np.sqrt((2 * self.P_b * (self.P_0 - self.P_b)) / (self.R * ((self.T_0) / (1 + (((self.gamma - 1)/(self.gamma)) * ((self.P_0 - self.P_b)/self.P_b))))))
            else:
                raise InputError("Please Define Nozzle back pressure before attempting to calculate flowrate!\nUse the method my_nozzle.set_back_pressure(P_b)")
            return self.flowrate
        else:
            raise InputError("Please Define Nozzle gamma and R before attempting to calculate flowrate!\nUse the method my_nozzle.set_gamma(gamma) and my_nozzle.set_R(R)")

    def is_choked(self) -> bool:
        if(self.choked == None): raise InputError("Define the flowrate first!")
        return self.choked
    
    def get_flowrate(self) -> float:
        if(self.choked == None): raise InputError("Define the flowrate first!")
        return self.flowrate
    
    def get_area_ratio(self) -> float:
        return self.A_e_over_A_star
    
    def __Mach_number_func_minus_Area_ratio(self, M:float) -> float:
        return (((self.gamma+1)/2)**(-(self.gamma+1)/(2*(self.gamma-1)))) * ((1 + ((self.gamma-1)/2) * M**2)**((self.gamma+1)/(2*(self.gamma-1)))) / M - self.A_e_over_A_star

    def get_exit_mach_number(self) -> float:
        self.M_exit = golden_search(self.__Mach_number_func_minus_Area_ratio, 1, 5)
        return self.M_exit
    
    def get_exit_conditions(self) -> tuple[float, float, float, float]:
        if(self.P_star == None or self.T_star == None or self.rho_star == None): raise InputError("Please set the stagnation conditions before calculating exit conditions!")
        if(self.M_exit == None):
            self.get_exit_mach_number()
        
        self.P_exit = self.P_star * (1 + ((self.gamma-1)/2)*self.M_exit**2)**(-(self.gamma)/(self.gamma-1))
        self.T_exit = self.T_star * (1 + ((self.gamma-1)/2)*self.M_exit**2)**(-1)
        self.rho_exit = self.rho_star * (1 + ((self.gamma-1)/2)*self.M_exit**2)**(-1/(self.gamma-1))
        self.speed_of_sound_exit = np.sqrt(self.gamma * self.R * self.T_exit)
        return self.P_exit, self.T_exit, self.rho_exit, self.speed_of_sound_exit

    def get_thrust_Isp(self) -> tuple[float, float]:
        if(self.P_star == None or self.T_star == None or self.rho_star == None): raise InputError("Please set the stagnation conditions before calculating exit conditions!")
        if(self.M_exit == None):
            self.get_exit_mach_number()
        if(self.speed_of_sound_exit == None or self.P_exit == None):
            self.get_exit_conditions()
        self.thrust = self.flowrate * self.M_exit * self.speed_of_sound_exit + (self.P_exit - self.P_b)*self.A_exit
        self.c_F = self.thrust / (self.P_0 * self.A_thoat)
        self.c_star = np.sqrt(((self.R * self.T_0)/self.gamma)*((self.gamma+1)/2)**((self.gamma+1)/(self.gamma-1)))
        self.Isp = (self.c_F * self.c_star)/9.81
        return self.thrust, self.Isp
    
    def get_thrust_coeff_characteristic_velocity(self) -> tuple[float, float]:
        if(self.c_F == None or self.c_star == None):
            self.get_Thrust_Isp()
        return self.c_F, self.c_star
