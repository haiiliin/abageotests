from .AbaqusVariable import AbaqusVariable


class AbaqusMaterialElastic(AbaqusVariable):
    modulus: float = None
    poisson_ratio: float = None

    def __init__(self, modulus: float = None, poisson_ratio: float = None):
        self.setParameters(modulus, poisson_ratio)

    def setModulus(self, modulus: float):
        self.modulus = modulus

    def setPoissonRatio(self, poisson_ratio: float):
        self.poisson_ratio = poisson_ratio

    def setParameters(self, modulus: float, poisson_ratio: float):
        self.modulus, self.poisson_ratio = modulus, poisson_ratio


class AbaqusMaterial(AbaqusVariable):
    density: float = None
    elastic: AbaqusMaterialElastic

    def __init__(self, modulus: float = None, poisson_ratio: float = None, density: float = None):
        self.elastic = AbaqusMaterialElastic()
        self.setParameters(modulus, poisson_ratio, density)

    def setParameters(self, modulus: float, poisson_ratio: float, density: float):
        self.setDensityParameters(density)
        self.setElasticParameters(modulus, poisson_ratio)

    def setDensityParameters(self, density: float):
        self.density = density

    def setElasticParameters(self, modulus: float, poisson_ratio: float):
        self.elastic.setParameters(modulus, poisson_ratio)