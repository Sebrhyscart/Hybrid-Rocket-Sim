from algorithm import *

class Injector: # oxidizer injector

    def __init__(self, flowrate):
        self.flowrate = flowrate

    def set_flowrate(self, flowrate:float):
        self.flowrate = flowrate
    
    def flowrate_step_func(self, t, t0) -> float:
        if t < t0:
            return self.flowrate
        else:
            return 0
    
