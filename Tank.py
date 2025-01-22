class Tank:

    def __init__(self, m_oxidizer:float):
        self.m_oxidizer = m_oxidizer
        self.pressure = None

    def linear_pressure_func(self, t):
        P0 = 3.447e6
        P1 = 2.413e6
        dt = 6
        pressure = ((P1 - P0)/dt) * t + P0
        self.pressure = pressure
        return pressure
