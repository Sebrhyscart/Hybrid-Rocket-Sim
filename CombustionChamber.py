from algorithm import *

class CombustionChamber:

    def __init__(self, r_wall: float, r_port: float, l_fuel: float, l_pre: float, l_post: float):
        # GEOMETRY DEFINITION
        self.r_wall = r_wall
        self.r_port = r_port
        self.l_fuel = l_fuel
        self.l_pre = l_pre
        self.l_post = l_post
        self.a_regress = 1 # regression rate coefficient
        self.n_regress = 1 # regression rate coefficient

    def get_r_port(self):
        return self.r_port
    
    def regression_rate(self, G_ox: float):
        return self.a_regress * (G_ox ** self.n_regress)
    