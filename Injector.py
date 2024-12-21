from algorithm import *

class Injector:

    def __init__(self, flowrate):
        '''
        Oxidizer Injector class. Returns a flowrate parameter or returns a flowrate vs time function
        '''
        self.flowrate = flowrate

    def set_flowrate(self, flowrate:float):
        self.flowrate = flowrate
    
    def flowrate_step_func(self, t, t0) -> float:
        if t < t0:
            return self.flowrate
        else:
            return 0
    
