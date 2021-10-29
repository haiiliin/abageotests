import typing
from enum import IntEnum

from .AbaqusVariable import AbaqusVariable
from .abaqusConstants import *


class AbaqusNormalContactType(IntEnum):
    Hard = 0
    Exponential = 1
    Linear = 2
    Tabular = 3
    ScaleFactor = 4


class AbaqusNormalContact(AbaqusVariable):
    contactStiffness: typing.Union[SymbolicConstant, float] = DEFAULT
    pressureOverclosure: SymbolicConstant = HARD
    allowSeparation: bool = ON
    maxStiffness: float = None
    table: tuple = ()
    constraintEnforcementMethod: SymbolicConstant = DEFAULT
    overclosureFactor: float = 0.0
    overclosureMeasure: float = 0.0
    contactStiffnessScaleFactor: float = 1.0
    initialStiffnessScaleFactor: float = 1.0
    clearanceAtZeroContactPressure: float = 0.0
    stiffnessBehavior: SymbolicConstant = LINEAR
    stiffnessRatio: float = 0.01
    upperQuadraticFactor: float = 0.03
    lowerQuadraticRatio: float = 0.3333

    def __init__(self, *args, **kwargs):
        self.setParameters(*args, **kwargs)

    def setParameters(self, *args, **kwargs):
        pass


class AbaqusHardNormalContact(AbaqusNormalContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
    """
    constraintEnforcementMethod: SymbolicConstant = DEFAULT  # AUGMENTED-LAGRANGE, PENALTY (STANDARD), DIRECT (STANDARD)
    allowSeparation: Boolean = ON

    pressureOverclosure = HARD

    def __init__(self, constraintEnforcementMethod: SymbolicConstant = DEFAULT, allowSeparation: Boolean = ON,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(constraintEnforcementMethod, allowSeparation)

    def setParameters(self, constraintEnforcementMethod: SymbolicConstant = DEFAULT, allowSeparation: Boolean = ON,
                      *args, **kwargs):
        self.constraintEnforcementMethod, self.allowSeparation = constraintEnforcementMethod, allowSeparation


class AbaqusExponentialNormalContact(AbaqusNormalContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=EXPONENTIAL, table=((1.0, 0.0), (0.0, 2.0)),
        maxStiffness=None, constraintEnforcementMethod=DEFAULT)
    constraintEnforcementMethod = DEFAULT ONLY
    table = ((pressure, 0.0), (0.0, clearance))
    """
    maxStiffness: float = None  # EXPLICIT ONLY
    # Fit exponential curve through this points
    pressure: float = None  # point 1: (pressure, 0)
    clearance: float = None  # point 2: (0, clearance)

    pressureOverclosure = EXPONENTIAL

    def __init__(self, maxStiffness: float = None, pressure: float = None, clearance: float = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(maxStiffness, pressure, clearance)

    def setParameters(self, maxStiffness: float, pressure: float, clearance: float, *args, **kwargs):
        self.maxStiffness, self.pressure, self.clearance = maxStiffness, pressure, clearance
        self.table = ((pressure, 0.0), (0.0, clearance),)


class AbaqusLinearNormalContact(AbaqusNormalContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=LINEAR, contactStiffness=2.0, constraintEnforcementMethod=DEFAULT)
    """
    contactStiffness: typing.Union[SymbolicConstant, float] = DEFAULT

    pressureOverclosure = LINEAR

    def __init__(self, contactStiffness: typing.Union[SymbolicConstant, float] = DEFAULT, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(contactStiffness)

    def setParameters(self, contactStiffness: typing.Union[SymbolicConstant, float] = DEFAULT, *args, **kwargs):
        self.contactStiffness = contactStiffness


class AbaqusTabularNormalContact(AbaqusNormalContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=TABULAR, table=((0.0, 1.0), (1.0, 2.0)),
        constraintEnforcementMethod=DEFAULT)
    constraintEnforcementMethod = DEFAULT ONLY
    table = ((pressure[0], overclosure[0]), (pressure[1], overclosure[1]), ...)
    """
    # Provide data in order descending of overclosure, a negative overclosure is a positive clearance
    pressure: typing.Iterable = None
    overclosure: typing.Iterable = None

    pressureOverclosure = TABULAR

    def __init__(self, pressure: typing.Iterable = None, overclosure: typing.Iterable = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(pressure, overclosure)

    def setParameters(self, pressure: typing.Iterable, overclosure: typing.Iterable, *args, **kwargs):
        self.pressure, self.overclosure = pressure, overclosure
        self.table = tuple((list(pressure)[i], list(overclosure)[i],) for i in range(len(list(pressure))))


class AbaqusScaleFactorNormalContact(AbaqusNormalContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
        pressureOverclosure=SCALE_FACTOR, contactStiffnessScaleFactor=4.0, initialStiffnessScaleFactor=1.0,
        overclosureFactor=0.0, overclosureMeasure=1.0, constraintEnforcementMethod=DEFAULT)
    constraintEnforcementMethod = DEFAULT ONLY
    """
    overclosureFactor: float = 0.0  # overclosureFactor or overclosureMeasure, one of this two must be set to 0.0
    overclosureMeasure: float = 0.0
    contactStiffnessScaleFactor: float = 1.0
    initialStiffnessScaleFactor: float = 1.0

    pressureOverclosure = SCALE_FACTOR

    def __init__(self, contactStiffnessScaleFactor: float = 1.0, initialStiffnessScaleFactor: float = 1.0,
                 overclosure_method: str = 'factor', overclosure: float = 0.0, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(contactStiffnessScaleFactor, initialStiffnessScaleFactor, overclosure_method,
                           overclosure)

    def setParameters(self, contactStiffnessScaleFactor: float, initialStiffnessScaleFactor: float,
                      overclosure_method: str = 'factor',  # or 'measure'
                      overclosure: float = 0.0, *args, **kwargs):
        self.contactStiffnessScaleFactor = contactStiffnessScaleFactor
        self.initialStiffnessScaleFactor = initialStiffnessScaleFactor
        if overclosure_method == 'factor':
            self.overclosureFactor = overclosure
        else:
            self.overclosureMeasure = overclosure


class AbaqusTangentialContactType(IntEnum):
    Frictionless = 0
    Penalty = 1
    ExponentialDecay = 2
    Rough = 3
    LagrangeMultiplier = 4
    UserDefined = 5


class AbaqusTangentialContact(AbaqusVariable):
    formulation: SymbolicConstant = FRICTIONLESS
    directionality: SymbolicConstant = ISOTROPIC
    slipRateDependency: Boolean = OFF
    pressureDependency: Boolean = OFF
    temperatureDependency: Boolean = OFF
    dependencies: int = 0
    exponentialDecayDefinition: SymbolicConstant = COEFFICIENTS
    table: tuple = ()
    shearStressLimit: float = None
    maximumElasticSlip: SymbolicConstant = FRACTION
    fraction: float = 0.0
    absoluteDistance: float = 0.0
    elasticSlipStiffness: float = None
    nStateDependentVars: int = 0
    useProperties: Boolean = OFF

    def __init__(self, *args, **kwargs):
        pass

    def setParameters(self, *args, **kwargs):
        pass


class AbaqusFrictionlessTangentialContact(AbaqusTangentialContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(formulation=FRICTIONLESS)
    """
    formulation = FRICTIONLESS

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(*args, **kwargs)

    def setParameters(self, *args, **kwargs):
        pass


class AbaqusPenaltyTangentialContact(AbaqusTangentialContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=ON,
        pressureDependency=ON, temperatureDependency=ON, dependencies=0,
        table=((0.1, 0.2, 0.3, 0.4), ), shearStressLimit=None,
        maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)

    table = (friction_coefficient [, friction_coefficient_2], slip_rate, contact_pressure, temp, field [, 2, 3 ...])
    maximumElasticSlip=FRACTION, fraction=0.005
    maximumElasticSlip=ABSOLUTE_DISTANCE, absoluteDistance=0.01
    """
    # FRICTION
    directionality: SymbolicConstant = ISOTROPIC
    frictionCoefficients: typing.Union[float, typing.Iterable] = ()  # 1 value for ISOTROPIC, 2 values ANISOTROPIC
    slipRate: float = None
    contactPressure: float = None
    temp: float = None
    dependencies: int = 0
    field: typing.Iterable = ()  # dependencies values
    # SHEAR STRESS
    shearStressLimit: float = None  # None for no limit
    # ELASTIC SLIP
    elasticSlipStiffness: float = None  # None for infinite
    maximumElasticSlip: SymbolicConstant = FRACTION  # OR absoluteDistance STANDARD ONLY
    fraction: float = 0.0
    absoluteDistance: float = 0.0

    formulation = PENALTY

    @typing.overload
    def __init__(self, directionality: SymbolicConstant = ISOTROPIC,
                 frictionCoefficients: typing.Union[float, typing.Iterable] = (), slipRate: float = None,
                 contactPressure: float = None, temp: float = None, dependencies: int = 0,
                 field: typing.Iterable = None, shearStressLimit: float = None,
                 elasticSlipStiffness: float = None,
                 maximumElasticSlip: SymbolicConstant = FRACTION, friction: float = 0.005):
        pass

    @typing.overload
    def __init__(self, directionality: SymbolicConstant = ISOTROPIC,
                 frictionCoefficients: typing.Union[float, typing.Iterable] = (), slipRate: float = None,
                 contactPressure: float = None, temp: float = None, dependencies: int = 0,
                 field: typing.Iterable = None, shearStressLimit: float = None,
                 elasticSlipStiffness: float = None,
                 maximumElasticSlip: SymbolicConstant = absoluteDistance, absoluteDistance: float = 0.0):
        pass

    def __init__(self, directionality: SymbolicConstant = ISOTROPIC,
                 frictionCoefficients: typing.Union[float, typing.Iterable] = (), slipRate: float = None,
                 contactPressure: float = None, temp: float = None, dependencies: int = 0,
                 field: typing.Iterable = None, shearStressLimit: float = None,
                 elasticSlipStiffness: float = None,
                 maximumElasticSlip: SymbolicConstant = FRACTION, maximumElasticSlip_value: float = 0.005,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(directionality, frictionCoefficients, slipRate, contactPressure, temp, dependencies,
                           field, shearStressLimit, elasticSlipStiffness, maximumElasticSlip,
                           maximumElasticSlip_value)

    @typing.overload
    def setParameters(self, directionality: SymbolicConstant = ISOTROPIC,
                      frictionCoefficients: typing.Union[float, typing.Iterable] = (), slipRate: float = None,
                      contactPressure: float = None, temp: float = None, dependencies: int = 0,
                      field: typing.Iterable = None, shearStressLimit: float = None,
                      elasticSlipStiffness: float = None,
                      maximumElasticSlip: SymbolicConstant = FRACTION, friction: float = 0.005):
        pass

    @typing.overload
    def setParameters(self, directionality: SymbolicConstant = ISOTROPIC,
                      frictionCoefficients: typing.Union[float, typing.Iterable] = (), slipRate: float = None,
                      contactPressure: float = None, temp: float = None, dependencies: int = 0,
                      field: typing.Iterable = None, shearStressLimit: float = None,
                      elasticSlipStiffness: float = None,
                      maximumElasticSlip: SymbolicConstant = absoluteDistance, absoluteDistance: float = 0.0):
        pass

    def setParameters(self, directionality: SymbolicConstant = ISOTROPIC,
                      frictionCoefficients: typing.Union[float, typing.Iterable] = (), slipRate: float = None,
                      contactPressure: float = None, temp: float = None, dependencies: int = 0,
                      field: typing.Iterable = None, shearStressLimit: float = None,
                      elasticSlipStiffness: float = None,
                      maximumElasticSlip: SymbolicConstant = FRACTION, maximumElasticSlip_value: float = 0.005,
                      *args, **kwargs):
        self.directionality, self.frictionCoefficients = directionality, tuple(frictionCoefficients)
        self.slipRate, self.contactPressure, self.temp = slipRate, contactPressure, temp
        self.dependencies, self.frictionCoefficients = dependencies, tuple(field)
        self.shearStressLimit = shearStressLimit
        self.elasticSlipStiffness = elasticSlipStiffness

        if slipRate is not None:
            self.slipRateDependency = ON
        if contactPressure is not None:
            self.pressureDependency = ON
        if temp is not None:
            self.temperatureDependency = ON

        self.maximumElasticSlip = maximumElasticSlip
        if maximumElasticSlip == FRACTION:
            self.fraction = maximumElasticSlip_value
        else:
            self.absoluteDistance = maximumElasticSlip_value

        self.table = (*self.frictionCoefficients, slipRate, contactPressure, temp, *self.field)


class AbaqusExponentialDecayTangentialContact(AbaqusTangentialContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=EXPONENTIAL_DECAY, exponentialDecayDefinition=COEFFICIENTS,
        table=((1.0, 1.0, 1.0), ), maximumElasticSlip=FRACTION, fraction=0.005,
        elasticSlipStiffness=None)
    table = ((static_coefficient, kinetic_coefficient, decay_coefficient), )
    """
    exponentialDecayDefinition: SymbolicConstant = COEFFICIENTS  # OR 'TEST_DATA'
    staticCoefficient: float = None
    kineticCoefficient: float = None
    decayCoefficient: float = None
    frictionCoefficients: list = None
    slipRate: float = None
    # SHEAR STRESS
    elasticSlipStiffness: float = None  # None for infinite
    # ELASTIC SLIP
    maximumElasticSlip: SymbolicConstant = FRACTION  # OR absoluteDistance
    fraction: float = 0.005
    absoluteDistance: float = 0.0

    formulation = EXPONENTIAL_DECAY

    @typing.overload
    def __init__(self, exponentialDecayDefinition: SymbolicConstant = COEFFICIENTS,
                 staticCoefficient: float = None, kineticCoefficient: float = None,
                 decayCoefficient: float = None, frictionCoefficients: list = None,
                 slipRate: float = None, elasticSlipStiffness: float = None,
                 maximumElasticSlip: SymbolicConstant = FRACTION, friction: float = 0.005):
        pass

    @typing.overload
    def __init__(self, exponentialDecayDefinition: SymbolicConstant = COEFFICIENTS,
                 staticCoefficient: float = None, kineticCoefficient: float = None,
                 decayCoefficient: float = None, frictionCoefficients: list = None,
                 slipRate: float = None, elasticSlipStiffness: float = None,
                 maximumElasticSlip: SymbolicConstant = absoluteDistance, absoluteDistance: float = 0.0):
        pass

    def __init__(self, exponentialDecayDefinition: SymbolicConstant = COEFFICIENTS,
                 staticCoefficient: float = None, kineticCoefficient: float = None,
                 decayCoefficient: float = None, frictionCoefficients: list = None,
                 slipRate: float = None, elasticSlipStiffness: float = None,
                 maximumElasticSlip: SymbolicConstant = FRACTION, maximumElasticSlip_value: float = 0.005,
                 *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(exponentialDecayDefinition, staticCoefficient, kineticCoefficient, decayCoefficient,
                           frictionCoefficients, slipRate, elasticSlipStiffness, maximumElasticSlip,
                           maximumElasticSlip_value)

    @typing.overload
    def setParameters(self, exponentialDecayDefinition: SymbolicConstant = COEFFICIENTS,
                      staticCoefficient: float = None, kineticCoefficient: float = None,
                      decayCoefficient: float = None, frictionCoefficients: list = None,
                      slipRate: float = None, elasticSlipStiffness: float = None,
                      maximumElasticSlip: SymbolicConstant = FRACTION, friction: float = 0.005):
        pass

    @typing.overload
    def setParameters(self, exponentialDecayDefinition: SymbolicConstant = COEFFICIENTS,
                      staticCoefficient: float = None, kineticCoefficient: float = None,
                      decayCoefficient: float = None, frictionCoefficients: list = None,
                      slipRate: float = None, elasticSlipStiffness: float = None,
                      maximumElasticSlip: SymbolicConstant = absoluteDistance, absoluteDistance: float = 0.0):
        pass

    def setParameters(self, exponentialDecayDefinition: SymbolicConstant = COEFFICIENTS,
                      staticCoefficient: float = None, kineticCoefficient: float = None,
                      decayCoefficient: float = None, frictionCoefficients: list = None,
                      slipRate: float = None, elasticSlipStiffness: float = None,
                      maximumElasticSlip: SymbolicConstant = FRACTION, maximumElasticSlip_value: float = 0.005,
                      *args, **kwargs):
        self.exponentialDecayDefinition = exponentialDecayDefinition

        self.staticCoefficient = staticCoefficient
        self.kineticCoefficient = kineticCoefficient
        self.decayCoefficient = decayCoefficient

        self.frictionCoefficients, self.slipRate = frictionCoefficients, slipRate

        self.maximumElasticSlip = maximumElasticSlip
        if maximumElasticSlip == FRACTION:
            self.fraction = maximumElasticSlip_value
        else:
            self.absoluteDistance = maximumElasticSlip_value

        if exponentialDecayDefinition == COEFFICIENTS:
            self.table = ((staticCoefficient, kineticCoefficient, decayCoefficient),)
        else:
            self.table = (frictionCoefficients[0], frictionCoefficients[1], slipRate, frictionCoefficients[2])


class AbaqusRoughTangentialContact(AbaqusTangentialContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(formulation=ROUGH)
    """
    formulation = ROUGH

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def setParameters(self, *args, **kwargs):
        pass


class AbaqusLagrangeMultiplierTangentialContact(AbaqusTangentialContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=LAGRANGE, directionality=ISOTROPIC, slipRateDependency=ON,
        pressureDependency=ON, temperatureDependency=ON, dependencies=0, table=((
        0.1, 0.2, 0.3, 0.4), ), shearStressLimit=None)

    table = ((friction_coefficient [, 2], slip_rate, contact_pressure, temp, field [, 2, 3 ...]), )
    maximumElasticSlip=FRACTION, fraction=0.005
    maximumElasticSlip=ABSOLUTE_DISTANCE, absoluteDistance=0.01
    """
    # FRICTION
    directionality: SymbolicConstant = ISOTROPIC
    frictionCoefficients: typing.Union[float, typing.Iterable] = ()  # 1 value for ISOTROPIC, 2 values ANISOTROPIC
    slipRate: float = None
    contactPressure: float = None
    temp: float = None
    dependencies: int = 0
    field: typing.Iterable = ()  # dependencies values
    # SHEAR STRESS
    shearStressLimit: float = None  # None for no limit

    formulation = LAGRANGE

    def __init__(self, directionality: SymbolicConstant = ISOTROPIC,
                 frictionCoefficients: typing.Union[float, typing.Iterable] = (), slipRate: float = None,
                 contactPressure: float = None, temp: float = None, dependencies: int = 0,
                 field: typing.Iterable = (), shearStressLimit: float = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(directionality, frictionCoefficients, slipRate, contactPressure, temp, dependencies,
                           field, shearStressLimit)

    def setParameters(self, directionality: SymbolicConstant = ISOTROPIC,
                      frictionCoefficients: typing.Union[float, typing.Iterable] = (), slipRate: float = None,
                      contactPressure: float = None, temp: float = None, dependencies: int = 0,
                      field: typing.Iterable = (), shearStressLimit: float = None, *args, **kwargs):
        self.directionality, self.frictionCoefficients = directionality, frictionCoefficients
        self.slipRate, self.contactPressure, self.temp = slipRate, contactPressure, temp
        self.dependencies, self.frictionCoefficients = dependencies, field
        self.shearStressLimit = shearStressLimit

        self.table = (*self.frictionCoefficients, slipRate, contactPressure, temp, *self.field)


class AbaqusUserDefinedTangentialContact(AbaqusTangentialContact):
    """
    mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
        formulation=USER_DEFINED, nStateDependentVars=N_STATE_DEPENDENT_VARS,
        useProperties=ON, table=FRICTION_PARAMETERS)
    """
    nStateDependentVars: int = 0
    useProperties: Boolean = ON

    cohesion: float = 0.0
    G0: float = 250.0
    R: float = 5.19
    e_cref: float = 0.89
    Lambda: float = 0.147
    xi: float = 0.424
    phi: float = 31.2
    dd: float = 7.57
    n_p: float = 2.06
    n_d: float = 0.46
    e0: float = 0.5771
    thickness: float = 5 * 0.23 * 1e-3

    properties: typing.Iterable = None

    formulation = USER_DEFINED

    def __init__(self, nStateDependentVars: int = 0, cohesion: float = 0.0, G0: float = 250.0, R: float = 5.19,
                 e_cref: float = 0.89, Lambda: float = 0.147, xi: float = 0.424, phi: float = 31.2, dd: float = 7.57,
                 n_p: float = 2.06, n_d: float = 0.46, e0: float = 0.5771, thickness: float = 5 * 0.23 * 1e-3, *args,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.setParameters(nStateDependentVars, cohesion, G0, R, e_cref, Lambda, xi, phi, dd, n_p, n_d, e0, thickness)

    @typing.overload
    def setParameters(self, nStateDependentVars: int = 0,
                      properties: typing.Iterable = (0.0, 250.0, 5.19, 0.89, 0.147, 0.424,
                                                     31.2, 7.57, 2.06, 0.46, 0.5771, 0.00115), *args, **kwargs):
        pass

    @typing.overload
    def setParameters(self, nStateDependentVars: int = 0, cohesion: float = 0.0, G0: float = 250.0, R: float = 5.19,
                      e_cref: float = 0.89, Lambda: float = 0.147, xi: float = 0.424, phi: float = 31.2,
                      dd: float = 7.57, n_p: float = 2.06, n_d: float = 0.46, e0: float = 0.5771,
                      thickness: float = 5 * 0.23 * 1e-3, *args, **kwargs):
        pass

    def setParameters(self, nStateDependentVars: int = 0, *args, **kwargs):
        if len(args) == 1:
            properties = tuple(arg for arg in args[0]) + tuple(value for value in kwargs.values())
        else:
            properties = tuple(arg for arg in args) + tuple(value for value in kwargs.values())
        cohesion, G0, R, e_cref, Lambda, xi, phi, dd, n_p, n_d, e0, thickness = properties
        self.nStateDependentVars = nStateDependentVars
        self.cohesion, self.G0, self.R, self.e_cref, self.Lambda = cohesion, G0, R, e_cref, Lambda
        self.xi, self.phi, self.dd, self.n_p, self.n_d, self.e0, self.thickness = xi, phi, dd, n_p, n_d, e0, thickness
        self.properties = [cohesion, G0, R, e_cref, Lambda, xi, phi, dd, n_p, n_d, e0, thickness]
        self.table = ((cohesion,), (G0,), (R,), (e_cref,), (Lambda,), (xi,), (phi,), (dd,), (n_p,),
                      (n_d,), (e0,), (thickness,),)


class AbaqusContact(AbaqusVariable):
    normal_contact: typing.Union[
        AbaqusNormalContact,
        AbaqusHardNormalContact,
        AbaqusExponentialNormalContact,
        AbaqusLinearNormalContact,
        AbaqusTabularNormalContact,
        AbaqusScaleFactorNormalContact
    ] = AbaqusHardNormalContact()
    tangential_contact: typing.Union[
        AbaqusTangentialContact,
        AbaqusFrictionlessTangentialContact,
        AbaqusPenaltyTangentialContact,
        AbaqusExponentialDecayTangentialContact,
        AbaqusRoughTangentialContact,
        AbaqusLagrangeMultiplierTangentialContact,
        AbaqusUserDefinedTangentialContact
    ] = AbaqusUserDefinedTangentialContact()
    normal_contact_type = AbaqusNormalContactType.Hard
    tangential_contact_type = AbaqusTangentialContactType.UserDefined

    def __init__(self, normal_contact_type=AbaqusNormalContactType.Hard,
                 tangential_contact_type=AbaqusTangentialContactType.UserDefined):
        self.setParameters(normal_contact_type, tangential_contact_type)

    def setParameters(self, normal_contact_type=AbaqusNormalContactType.Hard,
                      tangential_contact_type=AbaqusTangentialContactType.UserDefined):
        self._create_normal_contact(normal_contact_type)
        self._create_tangential_contact(tangential_contact_type)

    def NormalContact(self, contact_type: AbaqusNormalContactType, *args, **kwargs):
        if self.normal_contact_type == contact_type:
            self.normal_contact.setParameters(*args, **kwargs)
        else:
            self._create_normal_contact(contact_type, *args, **kwargs)

    def TangentialContact(self, contact_type: AbaqusTangentialContactType, *args, **kwargs):
        if self.tangential_contact_type == contact_type:
            self.tangential_contact.setParameters(*args, **kwargs)
        else:
            self._create_tangential_contact(contact_type, *args, **kwargs)

    def NormalContactType(self, contact_type: AbaqusNormalContactType):
        pass

    def _create_normal_contact(self, contact_type: AbaqusNormalContactType, *args, **kwargs):
        if self.normal_contact_type == contact_type:
            return

        self.normal_contact_type = contact_type
        if contact_type == AbaqusNormalContactType.Hard:
            self.normal_contact = AbaqusHardNormalContact(*args, **kwargs)
        elif contact_type == AbaqusNormalContactType.Exponential:
            self.normal_contact = AbaqusExponentialNormalContact(*args, **kwargs)
        elif contact_type == AbaqusNormalContactType.Linear:
            self.normal_contact = AbaqusLinearNormalContact(*args, **kwargs)
        elif contact_type == AbaqusNormalContactType.Tabular:
            self.normal_contact = AbaqusTabularNormalContact(*args, **kwargs)
        elif contact_type == AbaqusNormalContactType.ScaleFactor:
            self.normal_contact = AbaqusScaleFactorNormalContact(*args, **kwargs)

    def _create_tangential_contact(self, contact_type: AbaqusTangentialContactType, *args, **kwargs):
        if self.tangential_contact_type == contact_type:
            return

        self.tangential_contact_type = contact_type
        if contact_type == AbaqusTangentialContactType.Frictionless:
            self.tangential_contact = AbaqusFrictionlessTangentialContact(*args, **kwargs)
        elif contact_type == AbaqusTangentialContactType.Penalty:
            self.tangential_contact = AbaqusPenaltyTangentialContact(*args, **kwargs)
        elif contact_type == AbaqusTangentialContactType.ExponentialDecay:
            self.tangential_contact = AbaqusExponentialDecayTangentialContact(*args, **kwargs)
        elif contact_type == AbaqusTangentialContactType.Rough:
            self.tangential_contact = AbaqusRoughTangentialContact(*args, **kwargs)
        elif contact_type == AbaqusTangentialContactType.LagrangeMultiplier:
            self.tangential_contact = AbaqusLagrangeMultiplierTangentialContact(*args, **kwargs)
        elif contact_type == AbaqusUserDefinedTangentialContact:
            self.tangential_contact = AbaqusUserDefinedTangentialContact(*args, **kwargs)
