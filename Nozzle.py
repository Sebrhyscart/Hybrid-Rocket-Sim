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
        self.P_0 = None
        self.T_0 = None
        self.flowrate = None                        # [kg/s] nozzle flow rate
        self.P_exit = None                          # [Pa] exit pressure
        self.T_exit = None                          # [K] exit temperature
        self.M_exit = None                          # exit Mach number
        self.v_eff = None                           # effective Nozzle velocity

    def set_gamma(self, gamma:float):
        self.gamma = gamma

    def set_R(self, R:float):
        self.R = R

    def set_back_pressure(self, P_b:float):
        self.P_b = P_b

    def __Mach_number_func_minus_Area_ratio(self, M:float) -> float:
        return (((self.gamma+1)/2)**(-(self.gamma+1)/(2*(self.gamma-1)))) * ((1 + ((self.gamma-1)/2) * M**2)**((self.gamma+1)/(2*(self.gamma-1)))) / M - self.A_e_over_A_star

    def calc_flowrate(self, P_0:float, T_0:float):
        self.P_0 = P_0
        self.T_0 = T_0

        if (self.gamma != None and self.R != None and self.P_b != None):

            P_star = P_0 / ((1 + (self.gamma-1)/2)**((self.gamma)/(self.gamma-1)))
            if (self.P_b < P_star):
                self.choked = True

                self.flowrate = P_0 * np.sqrt(self.gamma / (self.R * T_0)) * self.A_throat * (1 + ((self.gamma-1)/2))**(-(self.gamma + 1)/(2 * (self.gamma - 1)))

                M_exit = golden_search(self.__Mach_number_func_minus_Area_ratio, 1, 5)
                P_exit = P_star * (1 + ((self.gamma-1)/2)*M_exit**2)**(-(self.gamma)/(self.gamma-1))
                T_exit = T_0 * ((1 + ((self.gamma-1)/2)*M_exit)**(-1))
                P_shock_exit = P_exit * (2*self.gamma*M_exit - (self.gamma - 1)) / (self.gamma + 1)
                speed_of_sound_exit = np.sqrt(self.gamma * self.R * T_exit)

                if (P_shock_exit <= self.P_b): # Shock in nozzle - aproximate as shock at exit for all shock positions
                    P_exit = P_shock_exit
                    T_exit = T_exit * (((2*self.gamma*M_exit**2 - (self.gamma - 1))) * ((self.gamma-1)*M_exit**2 + 2))/((self.gamma+1)**2 * M_exit**2)
                    M_exit = ((self.gamma-1)*M_exit**2 + 2)/(2*self.gamma*M_exit**2 - (self.gamma - 1))
                    speed_of_sound_exit = np.sqrt(self.gamma * self.R * T_exit)
            else:
                self.choked = False
                P_exit = self.P_b
                M_exit = np.sqrt((2/(self.gamma-1)) * ((P_0/self.P_b)**((self.gamma-1)/self.gamma) - 1))
                T_exit = T_0 * ((1 + ((self.gamma-1)/2)*M_exit)**(-1)) 
                speed_of_sound_exit = np.sqrt(self.gamma * self.R * T_exit)
                self.flowrate = P_0 * np.sqrt(self.gamma / (self.R * T_0)) * self.A_exit * M_exit * (1 + ((self.gamma-1)/2 * M_exit**2))**(-(self.gamma + 1)/(2 * (self.gamma - 1)))

            self.P_exit = P_exit
            self.T_exit = T_exit
            self.M_exit = M_exit
            self.v_eff = M_exit * speed_of_sound_exit

    def __Mach_number_func_minus_Area_ratio(self, M:float) -> float:
        return (((self.gamma+1)/2)**(-(self.gamma+1)/(2*(self.gamma-1)))) * ((1 + ((self.gamma-1)/2) * M**2)**((self.gamma+1)/(2*(self.gamma-1)))) / M - self.A_e_over_A_star

    def is_choked(self) -> bool:
        if(self.choked == None): raise InputError("Define the flowrate first!")
        return self.choked

    def get_thrust_Isp(self) -> tuple[float, float]:
        self.thrust = self.flowrate * self.v_eff + (self.P_exit - self.P_b)*self.A_exit
        self.c_F = self.thrust / (self.P_0 * self.A_throat)
        return self.thrust, self.Isp
    
    def get_thrust_coeff_characteristic_velocity(self) -> tuple[float, float]:
        self.thrust = self.flowrate * self.v_eff + (self.P_exit - self.P_b)*self.A_exit
        self.c_F = self.thrust / (self.P_0 * self.A_throat)
        self.c_star = np.sqrt(((self.R * self.T_0)/self.gamma)*((self.gamma+1)/2)**((self.gamma+1)/(self.gamma-1)))
        self.Isp = (self.c_F * self.c_star)/9.81
        return self.c_F, self.c_star
    
