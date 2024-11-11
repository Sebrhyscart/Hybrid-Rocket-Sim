from algorithm import *

class CombustionChamber:

    def __init__(self, r_wall: float, r_port: float, l_fuel: float, l_pre: float, l_post: float):
        # GEOMETRY DEFINITION
        self.r_wall = r_wall
        self.r_port = r_port
        self.l_fuel = l_fuel
        self.l_pre = l_pre
        self.l_post = l_post
        self.a_regress = 1.688e-4 # regression rate coefficient (BASE SI UNITS)
        self.n_regress = 0.600 # regression rate coefficient (BASE SI UNITS)

    def set_r_port(self, r_port:float):
        self.r_port = r_port

    def get_r_port(self) -> float:
        return self.r_port
    
    def regression_rate(self, G_ox: float) -> float:
        return self.a_regress * (G_ox ** self.n_regress)
    
    def get_volume(self) -> float:
        return (np.pi * self.r_port**2 * self.l_fuel) + (np.pi * self.r_wall**2 * (self.l_pre + self.l_post))
    