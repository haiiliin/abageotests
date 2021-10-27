from .AbaqusVariable import AbaqusVariable


class AbaqusMaterialMohrCoulomb(AbaqusVariable):
    friction_angle: float = None
    dilation_angle: float = None
    cohesion_yield_stress: float = None
    abs_plastic_strain: float = None

    def __init__(self, friction_angle: float = None, dilation_angle: float = None,
                 cohesion_yield_stress: float = None, abs_plastic_strain: float = None):
        self.setParameters(friction_angle, dilation_angle, cohesion_yield_stress, abs_plastic_strain)

    def setParameters(self, friction_angle: float = None, dilation_angle: float = None,
                      cohesion_yield_stress: float = None, abs_plastic_strain: float = None):
        self.friction_angle, self.dilation_angle = friction_angle, dilation_angle
        self.cohesion_yield_stress, self.abs_plastic_strain = cohesion_yield_stress, abs_plastic_strain


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
    mohr_coulomb: AbaqusMaterialMohrCoulomb

    def __init__(self, modulus: float = None, poisson_ratio: float = None, density: float = None,
                 friction_angle: float = None, dilation_angle: float = None,
                 cohesion_yield_stress: float = None, abs_plastic_strain: float = None):
        self.elastic = AbaqusMaterialElastic()
        self.mohr_coulomb = AbaqusMaterialMohrCoulomb()
        self.setParameters(modulus, poisson_ratio, density, friction_angle, dilation_angle, cohesion_yield_stress,
                           abs_plastic_strain)

    def setParameters(self, modulus: float = None, poisson_ratio: float = None, density: float = None,
                      friction_angle: float = None, dilation_angle: float = None, cohesion_yield_stress: float = None,
                      abs_plastic_strain: float = None):
        self.setDensityParameters(density)
        self.setElasticParameters(modulus, poisson_ratio)

    def setDensityParameters(self, density: float):
        self.density = density

    def setElasticParameters(self, modulus: float, poisson_ratio: float):
        self.elastic.setParameters(modulus, poisson_ratio)

    def setMohrCoulombParameters(self, friction_angle: float, dilation_angle: float,  cohesion_yield_stress: float,
                                 abs_plastic_strain: float):
        self.mohr_coulomb.setParameters(friction_angle, dilation_angle, cohesion_yield_stress, abs_plastic_strain)
