class Tank:

    def __init__(self, m_oxidizer:float):
        self.m_oxidizer = m_oxidizer
        self.pressure = None

    def linear_pressure_func(self, t):
        P0 = 550 * 6894.76
        P1 = 350 * 6894.76
        dt = 8
        pressure = ((P1 - P0)/dt) * t + P0
        self.pressure = pressure
        return pressure
