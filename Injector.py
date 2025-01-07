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
        Flowrate sinusoidal function
        '''
        amplitude = self.flowrate
        frequency = 4 # Hz
        offset = 2 * amplitude

        flowrate = -amplitude * np.cos(2*np.pi * frequency * t) + offset
        if flowrate < 0: raise InjectorError("negative flowrate should not be specified. Bad flowrate definition!")
        return amplitude * np.sin(2*np.pi * frequency * t) + offset
    
