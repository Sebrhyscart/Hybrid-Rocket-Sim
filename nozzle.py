from algorithm import *

class Nozzle: # converging-diverging (DeLaval) nozzle class

    def __init__(self, A_throat: float, A_exit: float):
        # GEOMETRY
        self.A_thoat = A_throat # [m]
        self.A_exit = A_exit # [m]
        self.A_e_over_A_star = A_exit / A_throat

        # GAS CONSTANTS
        self.gamma = None # specific heat capacity ratio
        self.R = None # [J/kg/K] gas constant

        # NOZZLE CONDITIONS
        self.choked = None
        self.P_b = None # [Pa] back pressure
        self.P_0 = None
        self.T_0 = None
        self.rho_0 = None
        self.flowrate = None

        self.P_star = None
        self.T_star = None
        self.rho_star = None
        self.chokedflowrate = None

        self.P_exit = None
        self.T_exit = None
        self.rho_exit = None
        self.M_exit = None
        self.speed_of_sound_exit = None

        self.thrust = None
        self.c_star = None
        self.c_F = None
        self.Isp = None

    def set_gamma(self, gamma):
        self.gamma = gamma

    def set_R(self, R):
        self.R = R

    def set_backpressure(self, P_b):
        self.P_b = P_b

    def set_stagnation(self, P_0, T_0):
        if (self.gamma != None and self.R != None):
            self.P_0 = P_0
            self.T_0 = T_0
            self.rho_0 = P_0 / (self.R * T_0)

            self.P_star = self.P_0 / ((1 + (self.gamma-1)/2)**((self.gamma)/(self.gamma-1)))
            self.T_star = self.T_0 / (1 + (self.gamma-1)/2)
            self.rho_star = self.rho_0 / ((1 + (self.gamma-1)/2)**(1/(self.gamma-1)))
            self.chokedflowrate = self.rho_star * self.A_thoat * np.sqrt(self.gamma * self.R * self.T_star)
        else:
            raise InputError("Please Define Nozzle gamma and R before setting the stagnation conditions!\nUse the method my_nozzle.set_gamma(gamma) and my_nozzle.set_R(R)")

    def set_flowrate_from_P0_T0_Pb(self):
        if (self.chokedflowrate != None):
            flowrate = np.sqrt((2 * self.P_b * (self.P_0 - self.P_b)) / (self.R * ((self.T_0) / (1 + (((self.gamma - 1)/(self.gamma)) * ((self.P_0 - self.P_b)/self.P_b))))))
            if(flowrate <= self.chokedflowrate):
                self.flowrate = flowrate
                self.choked = False
            else:
                self.flowrate = self.chokedflowrate
                self.choked = True
            return self.flowrate
        else:
            raise InputError("Please Define Nozzle stagnation conditions and/or back pressure before setting the flowrate using stagnation conditions!\nUse the method my_nozzle.set_stagnation(P_0, T_0) and/or my_nozzle.set_backpressure(P_b)")

    def set(self, gamma, R, P_b, P_0, T_0):
        self.gamma = gamma
        self.R = R
        self.P_b = P_b
        self.P_0 = P_0
        self.T_0 = T_0
        self.rho_0 = P_0 / (self.R * T_0)

        self.P_star = self.P_0 / ((1 + (self.gamma-1)/2)**((self.gamma)/(self.gamma-1)))
        self.T_star = self.T_0 / (1 + (self.gamma-1)/2)
        self.rho_star = self.rho_0 / ((1 + (self.gamma-1)/2)**(1/(self.gamma-1)))
        self.chokedflowrate = self.rho_star * self.A_thoat * np.sqrt(self.gamma * self.R * self.T_star)

        flowrate = np.sqrt((2 * self.P_b * (self.P_0 - self.P_b)) / (self.R * ((self.T_0) / (1 + (((self.gamma - 1)/(self.gamma)) * ((self.P_0 - self.P_b)/self.P_b))))))
        if(flowrate <= self.chokedflowrate):
            self.flowrate = flowrate
            self.choked = False
        else:
            self.flowrate = self.chokedflowrate
            self.choked = True
        return self.flowrate

    def is_choked(self):
        if(self.choked == None): raise InputError("Define the flowrate first!")
        return self.choked
    
    def get_flowrate(self):
        if(self.choked == None): raise InputError("Define the flowrate first!")
        return self.flowrate
    
    def get_area_ratio(self):
        return self.A_e_over_A_star
    
    def __Mach_number_func_minus_Area_ratio(self, M):
        return (((self.gamma+1)/2)**(-(self.gamma+1)/(2*(self.gamma-1)))) * ((1 + ((self.gamma-1)/2) * M**2)**((self.gamma+1)/(2*(self.gamma-1)))) / M - self.A_e_over_A_star

    def get_exit_mach_number(self):
        self.M_exit = golden_search(self.__Mach_number_func_minus_Area_ratio, 1, 5)
        return self.M_exit
    
    def get_exit_conditions(self):
        if(self.P_star == None or self.T_star == None or self.rho_star == None): raise InputError("Please set the stagnation conditions before calculating exit conditions!")
        if(self.M_exit == None):
            self.get_exit_mach_number()
        
        self.P_exit = self.P_star * (1 + ((self.gamma-1)/2)*self.M_exit**2)**(-(self.gamma)/(self.gamma-1))
        self.T_exit = self.T_star * (1 + ((self.gamma-1)/2)*self.M_exit**2)**(-1)
        self.rho_exit = self.rho_star * (1 + ((self.gamma-1)/2)*self.M_exit**2)**(-1/(self.gamma-1))
        self.speed_of_sound_exit = np.sqrt(self.gamma * self.R * self.T_exit)
        return self.P_exit, self.T_exit, self.rho_exit, self.speed_of_sound_exit

    def get_thrust_Isp(self):
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
    
    def get_thrust_coeff_characteristic_velocity(self):
        if(self.c_F == None or self.c_star == None):
            self.get_Thrust_Isp()
        return self.c_F, self.c_star
