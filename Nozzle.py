from algorithm import *
from scipy.optimize import fsolve

class NozzleError(Exception):
    pass

class Nozzle:

    def __init__(self, A_throat:float, A_exit:float, length:float):
        '''
        Converging-Diverging (DeLaval) nozzle class. It uses isentropic compressible flow equations to model the pressure, temperature, flowrate, 
        and Mach number in the nozzle. From this, thrust and specific impulse (Isp) is calculated.
        '''
        self.A_throat = A_throat                    # [m] throat area
        self.A_exit = A_exit                        # [m] exit area
        self.length = length                        # [m] length of diverging section
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
        self.M_throat = None
        self.x_shock = None

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
        def func(M):
            return self.A_over_Astar(M) - A_ratio

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
        if self.P_b == None: raise NozzleError("Define a back pressure (self.P_b) before calling calc_flowrate(P_0, T_0)")

        if (P_0 > self.P_b):
            # Assuming subsonic flow, calculate exit Mach number 
            M_exit = np.sqrt(2/(self.gamma-1) * ((self.P_b/P_0) ** (-(self.gamma-1)/self.gamma) - 1))
            if M_exit > 1: M_exit = 1
            A_star = self.A_exit / self.A_over_Astar(M_exit)

            self.A_star_temp = A_star

            if A_star >= self.A_throat:
                self.choked = True
                M_throat = 1
                M_exit = self.mach_number_solver(self.A_exit/self.A_throat, supersonic=True)
            else: 
                self.choked = False 
                M_throat = self.mach_number_solver(self.A_throat/A_star, supersonic=False)
            
            self.M_throat = M_throat

            P_throat = P_0 * self.P_over_P0(M_throat)
            T_throat = T_0 * self.T_over_T0(M_throat)
            self.flowrate = P_throat * np.sqrt(self.gamma/(self.R * T_throat)) * M_throat * self.A_throat
            
            # Exit pressure and temperature and speed of sound are a function of M_exit
            P_exit = P_0 * self.P_over_P0(M_exit)
            T_exit = T_0 * self.T_over_T0(M_exit)
            speed_of_sound_exit = np.sqrt(self.gamma * self.R * T_exit)

            # Shock in the nozzle if P_shock_exit < back pressure
            if (M_exit > 1 and P_exit * self.P2_over_P1_shock(M_exit) <= self.P_b): 
                # location of shock in the nozzle
                x_shock = self.shock_position_in_nozzle_solver(P_0)
                self.x_shock = x_shock
                A_shock = self.Area_vs_position(x_shock)
                # Mach number at postion x
                M_x1 = self.mach_number_solver(A_shock/self.A_throat, supersonic=True)
                # New stagnation pressure after the shock
                P_02 = P_0 * self.P02_over_P01_shock(M_x1)
                # exit mach number based on new stagnation pressure
                M_exit = np.sqrt(2/(self.gamma-1) * ((self.P_b/P_02) ** (-(self.gamma-1)/self.gamma) - 1))
                # Pressure at exit of nozzle based on loss of stagnation pressure over shock then isentopic contraction
                P_exit = P_02 * self.P_over_P0(M_exit)
                # No loss in stagnation temperature over a shock
                T_exit = T_0 * self.T_over_T0(M_exit)
                speed_of_sound_exit = np.sqrt(self.gamma * self.R * T_exit)
            
            else: x_shock = 0
            self.x_shock = x_shock

            self.P_exit = P_exit
            self.T_exit = T_exit
            self.M_exit = M_exit
            self.v_eff = M_exit * speed_of_sound_exit
        
        else:
            # print("Warning: Low chamber pressure!")
            self.choked = False
            self.flowrate = 0
            self.P_exit = P_0
            self.T_exit = T_0
            self.M_exit = 0
            self.M_throat = 0
            self.v_eff = 0
            self.x_shock = 0

    def is_choked(self) -> bool:
        if(self.choked == None): raise InputError("Define the flowrate first!")
        return self.choked

    def get_thrust_Isp(self) -> tuple[float, float]:
        '''
            Calculate the thrust and specific impulse (Isp) of the engine
        '''
        if self.flowrate == 0: return 0,0
        thrust = self.flowrate * self.v_eff + (self.P_exit - self.P_b)*self.A_exit
        Isp = thrust / self.flowrate / 9.81
        return thrust, Isp
    
    def get_thrust_coeff_characteristic_velocity(self) -> tuple[float, float]:
        '''
            Calculate the thrust coefficient, c_F and the characteristic velocity, c* of the engine
        '''
        if self.flowrate == 0: return 0,0
        thrust = self.flowrate * self.v_eff + (self.P_exit - self.P_b)*self.A_exit
        c_F = thrust / (self.P_0 * self.A_throat)
        c_star = (self.P_0 * self.A_throat) / self.flowrate
        return c_F, c_star
    
    def T_over_T0(self, M):
        return ((1 + (self.gamma-1)/2 * M**2))**(-1)
        
    def P_over_P0(self, M):
        return ((1 + (self.gamma-1)/2 * M**2))**(-(self.gamma)/(self.gamma-1))
        
    def T_over_Tstar(self, M):
        return ((1 + (self.gamma-1)/2 * M**2) / (1 + (self.gamma-1)/2))**(-1)
    
    def P_over_Pstar(self, M):
        return ((1 + (self.gamma-1)/2 * M**2) / (1 + (self.gamma-1)/2))**(-(self.gamma)/(self.gamma-1))
        
    def A_over_Astar(self, M):
        return ((((self.gamma+1)/2)**(-(self.gamma+1)/(2*(self.gamma-1)))) * ((1 + ((self.gamma-1)/2) * M**2)**((self.gamma+1)/(2*(self.gamma-1)))) / M)
        
    def M2_shock(self, M1):
        if M1 < 1: raise NozzleError("Shocks only exist supersonic -> subsonic!")
        return np.sqrt(((self.gamma-1)*M1**2 + 2) / (2*self.gamma*M1**2 - (self.gamma-1)))
        
    def T2_over_T1_shock(self, M1):
        if M1 < 1: raise NozzleError("Shocks only exist supersonic -> subsonic!")
        return ((((self.gamma-1)*M1**2 + 2) * (2*self.gamma*M1**2 - (self.gamma-1))) / ((self.gamma-1)*M1**2))
        
    def P2_over_P1_shock(self, M1):
        if M1 < 1: raise NozzleError("Shocks only exist supersonic -> subsonic!")
        return ((2*self.gamma*M1**2 - (self.gamma-1)) / (self.gamma+1))
    
    def P02_over_P01_shock(self, M1):
        if M1 < 1: raise NozzleError("Shocks only exist supersonic -> subsonic!")
        return ((((self.gamma+1)*M1**2) / ((self.gamma-1)*M1**2 + 2)) ** (self.gamma / (self.gamma-1))) * (((self.gamma+1) / (2*self.gamma*M1**2-(self.gamma-1))) ** (1 / (self.gamma-1)))
        
    def Area_vs_position(self,x):
        '''
        Linear expanding nozzle contour
        '''
        r_throat = np.sqrt(self.A_throat / np.pi)
        r_exit = np.sqrt(self.A_exit / np.pi)
        r_vs_x = r_throat + (r_exit - r_throat) * x / self.length
        return np.pi * r_vs_x ** 2

    def shock_position_in_nozzle_solver(self, P_0):
        '''
        Determines the location of a shock in the nozzle
        '''
        def exit_pressure_after_shock_at_x(x):
            # Area at position x
            A_x = self.Area_vs_position(x)
            # Mach number at postion x
            M_x1 = self.mach_number_solver(A_x/self.A_throat, supersonic=True)
            # New stagnation pressure after the shock
            P_02 = P_0 * self.P02_over_P01_shock(M_x1)
            # Mach number after shock
            M_x2 = self.M2_shock(M_x1)
            # Area required to choke new subsonic flow
            A_star = A_x / self.A_over_Astar(M_x2)
            # Exit mach number from area ratio using this new reference choking area
            M_exit = self.mach_number_solver(self.A_exit/A_star,supersonic=False)
            # Pressure at exit of nozzle based on loss of stagnation pressure over shock then isentopic contraction
            return P_02 * self.P_over_P0(M_exit)

        def func(x):
            return exit_pressure_after_shock_at_x(x) - self.P_b

        # golden search domain
        a = 0
        b = self.length
        guess = 0.5 * self.length

        return golden_search(func,a,b,tol=1e-8) #fsolve(func, guess)[0]
