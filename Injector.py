from algorithm import *

class InjectorError(Exception):
    pass

class Injector:

    def __init__(self, N_orifice:int, A_orifice:float):
        '''
        Oxidizer Injector class. Returns a flowrate parameter or returns a flowrate vs time function
        '''
        self.N_orifice = N_orifice
        self.A_orifice = A_orifice
        self.flowrate = None
        self.Cd = 0.66 # Coefficient of orifice dischage.

    def throttle(self,t):
        pass
    
    def calc_flowrate_spi(self, P_injector_1:float, P_injector_2:float) -> float:
        '''
        Calculate flowrate of injector assuming single phase incompressible fluid

        Parameters:
            P_injector_1: Pressure Upstream of Injector
            P_injector_2: Pressure Downstream of Injector
        '''
        # Upstream Injector Properties
        x_injector_1 = 0 # Assume: Quality = staturated liquid (x = 0)
        rho_injector_1 = CP.PropsSI("D", "P", P_injector_1, "Q", x_injector_1, NOX) # Density = rho_sat at P
        
        if (P_injector_1 < P_injector_2):
            m_dot_SPI = 0
        else:
            m_dot_SPI = self.Cd * self.N_orifice * self.A_orifice * np.sqrt(2 * rho_injector_1 * (P_injector_1 - P_injector_2))

        self.flowrate = m_dot_SPI
        return m_dot_SPI
    
    def calc_flowrate_hem(self, P_injector_1:float, P_injector_2:float) -> float:
        '''
        Calculate flowrate of injector assuming homogenous equilibrium phase change for compressible fluid

        Parameters:
            P_injector_1: Pressure Upstream of Injector
            P_injector_2: Pressure Downstream of Injector
        '''
        # Upstream Injector Properties
        x_injector_1 = 0 # Assume: Quality = staturated liquid (x = 0)
        s_injector_1 = CP.PropsSI("S", "P", P_injector_1, "Q", x_injector_1, NOX) # Entropy = s_sat at P
        h_injector_1 = CP.PropsSI("H", "P", P_injector_1, "Q", x_injector_1, NOX)   # Enthalpy = h_sat at P

        # Downstream Injector Properties
        s_injector_2 = s_injector_1 # Assume Isentropic flow through injector
        rho_injector_2 = CP.PropsSI("D", "P", P_injector_2, "S", s_injector_2, NOX) # Density = rho at P and s
        h_injector_2 = CP.PropsSI("H", "P", P_injector_2, "S", s_injector_2, NOX)   # Enthalpy = h at P and s

        if (P_injector_1 < P_injector_2) or (h_injector_1 < h_injector_2): 
            m_dot_HEM = 0
        else:
            m_dot_HEM = self.Cd * self.N_orifice * self.A_orifice * rho_injector_2 * np.sqrt(2 * (h_injector_1 - h_injector_2))

        self.flowrate = m_dot_HEM
        return m_dot_HEM

    def calc_flowrate_dyer(self, P_injector_1:float, P_injector_2:float) -> float:
        '''
        Calculate flowrate of injector assuming homogenous non-equilibrium phase change for compressible fluid,
        proposed by Dyer et al.

        Parameters:
            P_injector_1: Pressure Upstream of Injector
            P_injector_2: Pressure Downstream of Injector
            P_super: supercharged pressure - how much greater is the pressure than P_sat
        '''
        m_dot_SPI = self.calc_flowrate_spi(P_injector_1, P_injector_2)
        m_dot_HEM = self.calc_flowrate_hem(P_injector_1, P_injector_2)
        # Non equilibrium factor
        kappa = 1 # Assume saturated liquid => kappa = 1
        # set the flowrate based on Dyer model
        m_dot_dyer = (kappa / (kappa + 1)) * m_dot_SPI + (1 / (kappa + 1)) * m_dot_HEM
        self.flowrate = m_dot_dyer
        return m_dot_dyer

    def flowrate_step_func(self, t:float, t0:float) -> float:
        '''
        Flowrate step down function
        '''
        if t < t0:
            return self.flowrate
        else:
            return 0
    
    def sinusoidal_flowrate(self, t:float) -> float:
        '''
        Flowrate sinusoidal function, for testing throttle
        '''
        amplitude = self.flowrate
        frequency = 4 # Hz
        offset = 2 * amplitude

        flowrate = -amplitude * np.cos(2*np.pi * frequency * t) + offset
        if flowrate < 0: raise InjectorError("negative flowrate should not be specified. Bad flowrate definition!")
        return amplitude * np.sin(2*np.pi * frequency * t) + offset
    