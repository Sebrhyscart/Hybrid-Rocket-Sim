from algorithm import *

class Injector: # oxidizer injector

    def __init__(self, flowrate):
        self.flowrate = flowrate

    def set_flowrate(self, flowrate:float):
        self.flowrate = flowrate

    def get_flowrate(self) -> float:
        return self.flowrate
    
