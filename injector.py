from algorithm import *

class Injector: # oxidizer injector

    def __init__(self):
        self.flowrate = None

    def set_flowrate(self, flowrate:float):
        self.flowrate = flowrate

    def get_flowrate(self) -> float:
        if (self.flowrate == None): raise ImportError("Injector Flow rate has not been defined yet!")
        return self.flowrate
    
