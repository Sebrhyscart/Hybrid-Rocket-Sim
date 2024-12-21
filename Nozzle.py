from algorithm import *
from scipy.optimize import fsolve

class Nozzle:

    def __init__(self, A_throat:float, A_exit:float):
        '''
        Converging-Diverging (DeLaval) nozzle class. It uses isentropic compressible flow equations to model the pressure, temperature, flowrate, 
        and Mach number in the nozzle. From this, thrust and specific impulse (Isp) is calculated.
        '''
        self.A_throat = A_throat                    # [m] throat area
        self.A_exit = A_exit                        # [m] exit area
        self.gamma = None                           # specific heat capacity ratio
        self.R = None                               # [J/kg/K] gas constant
        self.P_0 = None                             # [Pa] Stagnation pressure
        self.T_0 = None                             # [K] Stagnation temperature
        self.choked = None                          # boolean set if flow is choked or not
        self.P_b = None                             # [Pa] back pressure
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

    def mach_number_solver(self, A_ratio, supersonic=True):
        """
        Solves for the Mach number given the area ratio and specific heat ratio.

        Parameters:
            A_ratio (float): Area ratio (A/A*)
            gamma (float): Specific heat ratio
            supersonic (bool): If True, find the supersonic solution. If False, find the subsonic solution.

        Returns:
            float: Mach number corresponding to the given area ratio.
        """
        gamma = self.gamma
        def func(M):
            return (((gamma + 1) / 2) ** (-(gamma + 1) / (2 * (gamma - 1)))) * \
                ((1 + ((gamma - 1) / 2) * M**2) ** ((gamma + 1) / (2 * (gamma - 1)))) / M - A_ratio

        # Initial guess for supersonic and subsonic solutions
        #   Solve by Newton's Method - supply initial guess for Mach number
        #   Solve by Golden Search - supply closed bounds containing 1 root
        if supersonic:
            # M_guess = 2.0 
            a = 1 
            b = 10
        else:
            # M_guess = 0.5
            a = 1e-6
            b = 1

        # return fsolve(func, M_guess)[0] # Newton's method -esque rootfinding algorithm
        return golden_search(func,a,b) # Golden search rootfinding algorithm

    def calc_flowrate(self, P_0, T_0):
        '''
        Calculate the mass flow rate and effective velocity at exit through this converging diverging nozzle.

        Parameters
            P_0: Stagnation Pressure, assumed to be pressure in the chamber
            T_0: Stagnation Temperture, assumed to be temperature in the chamber
        '''
        self.P_0 = P_0
        self.T_0 = T_0

        if (P_0 > self.P_b):
            # Assuming subsonic flow, calculate exit Mach number 
            M_exit = np.sqrt(2/(self.gamma-1) * ((self.P_b/P_0) ** (-(self.gamma-1)/self.gamma) - 1))

            # Find the area required to choak this flow, A*
            A_star = self.A_exit / ((((self.gamma+1)/2)**(-(self.gamma+1)/(2*(self.gamma-1)))) * ((1 + ((self.gamma-1)/2) * M_exit**2)**((self.gamma+1)/(2*(self.gamma-1)))) / M_exit)
            # Find the mach number at the throat based on the area ratio A_throat/A*
            M_throat = self.mach_number_solver(self.A_throat/A_star, supersonic=False)

            P_throat = P_0 / ((1 + (self.gamma-1)/2 * M_throat**2) ** ((self.gamma)/(self.gamma-1)))
            T_throat = T_0 / ((1 + (self.gamma-1)/2 * M_throat**2))
            self.flowrate = P_throat * np.sqrt(self.gamma/(self.R * T_throat)) * M_throat * self.A_throat

            if (abs(M_throat - 1) < 1e-6):
                M_throat = 1
                self.choked = True
                # Assuming fully expanded flow, the exit Mach number, M_exit is a function of A_e / A_t
                M_exit = self.mach_number_solver(self.A_exit/self.A_throat, supersonic=True)

            else:
                self.choked = False
                
            # Exit pressure and temperature and speed of sound are a function of M_exit
            P_exit = P_0 * (1 + ((self.gamma-1)/2)*M_exit**2)**(-(self.gamma)/(self.gamma-1))
            T_exit = T_0 * ((1 + ((self.gamma-1)/2)*M_exit)**(-1))
            speed_of_sound_exit = np.sqrt(self.gamma * self.R * T_exit)

            # Shock in the nozzle if P_shock_exit < back pressure
            P_shock_exit = P_exit * (2*self.gamma*M_exit - (self.gamma - 1)) / (self.gamma + 1)
            if(M_exit > 1 and P_shock_exit <= self.P_b): # Approximate all shocks positions in nozzle as shock at exit
                P_exit = P_exit * ((2*self.gamma*M_exit**2 - (self.gamma-1)) / (self.gamma+1))
                T_exit = T_exit * ((((self.gamma-1)*M_exit**2 + 2) * (2*self.gamma*M_exit**2 - (self.gamma-1))) / ((self.gamma-1)*M_exit**2))
                M_exit = np.sqrt(((self.gamma-1)*M_exit**2 + 2) / (2*self.gamma*M_exit**2 - (self.gamma-1)))
                speed_of_sound_exit = np.sqrt(self.gamma * self.R * T_exit)

            self.P_exit = P_exit
            self.T_exit = T_exit
            self.M_exit = M_exit
            self.v_eff = M_exit * speed_of_sound_exit
        
        else:
            self.flowrate = 0
            self.P_exit = P_0
            self.T_exit = T_0
            self.M_exit = 0
            self.v_eff = 0

    def is_choked(self) -> bool:
        if(self.choked == None): raise InputError("Define the flowrate first!")
        return self.choked

    def get_thrust_Isp(self) -> tuple[float, float]:
        '''
            Calculate the thrust and specific impulse (Isp) of the engine
        '''
        thrust = self.flowrate * self.v_eff + (self.P_exit - self.P_b)*self.A_exit
        c_F = thrust / (self.P_0 * self.A_throat)
        c_star = np.sqrt(((self.R * self.T_0)/self.gamma)*((self.gamma+1)/2)**((self.gamma+1)/(self.gamma-1)))
        Isp = (c_F * c_star)/9.81
        return thrust, Isp
    
    def get_thrust_coeff_characteristic_velocity(self) -> tuple[float, float]:
        '''
            Calculate the thrust coefficient, c_F and the characteristic velocity, c* of the engine
        '''
        thrust = self.flowrate * self.v_eff + (self.P_exit - self.P_b)*self.A_exit
        c_F = thrust / (self.P_0 * self.A_throat)
        c_star = np.sqrt(((self.R * self.T_0)/self.gamma)*((self.gamma+1)/2)**((self.gamma+1)/(self.gamma-1)))
        return c_F, c_star
