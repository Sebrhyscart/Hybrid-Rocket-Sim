
class Tank:

    def __init__(self,m_oxidizer:float):
        self.m_oxidizer = m_oxidizer

    def set_m_oxidizer(self, m_oxidizer:float):
        self.m_oxidizer = m_oxidizer

    def get_m_oxidizer(self) -> float:
        return self.m_oxidizer