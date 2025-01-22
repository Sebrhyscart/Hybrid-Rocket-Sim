from algorithm import *

class CombustionChamber:

    def __init__(self, r_wall: float, r_port: float, l_fuel: float, l_pre: float, l_post: float):
        # GEOMETRY DEFINITION
        self.r_wall = r_wall
        self.r_port = r_port
        self.l_fuel = l_fuel
        self.l_pre = l_pre
        self.l_post = l_post

    def set_r_port(self, r_port:float):
        self.r_port = r_port

    def get_r_port(self) -> float:
        return self.r_port
    
    def paraffin_regression_rate(self, G_ox: float) -> float:
        '''
        Returns the regression rate based on an experimental coorliation: r_dot = a * G_ox ^ n

        if given G_ox in g/cm^2/s and a in mm/s then to use G_ox kg/m^2/s and a in m/s, put in a (m/s) = 10^n /1000 * a (mm/s)

        Parameters:
            G_ox: Oxidizer mass flux
        '''
        a_regress = 1.688e-4 # regression rate coefficient (BASE SI UNITS)
        n_regress = 0.600 # regression rate coefficient (BASE SI UNITS)
        return a_regress * (G_ox ** n_regress)
    
    def HTPB_regression_rate(self, G_ox: float) -> float:
        '''
        Returns the regression rate based on an experimental coorliation: r_dot = a * G_ox ^ n

        Parameters:
            G_ox: Oxidizer mass flux
        '''
        a_regress = 9.03e-5 # regression rate coefficient (BASE SI UNITS)
        n_regress = 0.527 # regression rate coefficient (BASE SI UNITS)
        return a_regress * (G_ox ** n_regress)
    
    def HTPB2_regression_rate(self, G_ox: float) -> float:
        '''
        Returns the regression rate based on an experimental coorliation: r_dot = a * G_ox ^ n

        Parameters:
            G_ox: Oxidizer mass flux
        '''
        a_regress = 4.148e-5 # regression rate coefficient (BASE SI UNITS)
        n_regress = 0.670 # regression rate coefficient (BASE SI UNITS)
        return a_regress * (G_ox ** n_regress)
    
    def HTPB20percentAl_regression_rate(self, G_ox: float) -> float:
        '''
        Returns the regression rate based on an experimental coorliation: r_dot = a * G_ox ^ n

        Parameters:
            G_ox: Oxidizer mass flux
        '''
        a_regress = 1.959e-5 # regression rate coefficient (BASE SI UNITS)
        n_regress = 0.956 # regression rate coefficient (BASE SI UNITS)
        return a_regress * (G_ox ** n_regress)
    
    def get_volume(self) -> float:
        return (np.pi * self.r_port**2 * self.l_fuel) + (np.pi * self.r_wall**2 * (self.l_pre + self.l_post))
    