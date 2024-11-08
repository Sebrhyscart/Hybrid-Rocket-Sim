from algorithm import *

class Injector: # oxidizer injector

    def __init__(self):
        self.flowrate = None

    def set_flowrate(self, flowrate):
        self.flowrate = flowrate

    def get_flowrate(self):
        if (self.flowrate == None): raise ImportError("Injector Flow rate has not been defined yet!")
        return self.flowrate
    
