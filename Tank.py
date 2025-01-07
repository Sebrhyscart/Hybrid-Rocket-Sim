import CoolProp.CoolProp as CP
NOX = "NitrousOxide"   

class Tank:

    def __init__(self, volume, m_oxidizer:float):
        if (m_oxidizer / volume > CP.PropsSI("D", "T", 300, "Q", 0, NOX)): raise AttributeError("The Oxidizer tank is over-full! Decrease the mass of oxidizer!")
        self.volume = volume
        self.m_oxidizer = m_oxidizer
        self.x_tank = CP.PropsSI("Q", "T", 300, "D", m_oxidizer / volume, NOX)

    def set_m_oxidizer(self, m_oxidizer:float):
        self.m_oxidizer = m_oxidizer


