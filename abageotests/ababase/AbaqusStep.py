import typing
from enum import IntEnum

from .AbaqusVariable import AbaqusVariable


class AbaqusCalculationMethod(IntEnum):
    Standard = 0
    Explicit = 1


class AbaqusBaseStepIncrementationType(IntEnum):
    Automatic = 0
    Fixed = 1


class AbaqusBaseIncrementation(AbaqusVariable):
    incrementation_type = AbaqusBaseStepIncrementationType.Automatic

    def incrementationType(self):
        if self.incrementation_type == AbaqusBaseStepIncrementationType.Automatic:
            return 'AUTOMATIC'
        else:
            return 'FIXED'


class AbaqusStep(AbaqusVariable):
    calculation_method = AbaqusCalculationMethod.Standard
    time_period: float = 1
    description: str = ''

    def __init__(self, calculation_method=AbaqusCalculationMethod.Standard, time_period: float = 1,
                 description: str = ''):
        self.calculation_method = calculation_method
        self.time_period, self.description = time_period, description

    def setParameters(self, calculation_method: AbaqusCalculationMethod, time_period: float, description: str):
        self.calculation_method = calculation_method
        self.time_period, self.description = time_period, description


class AbaqusGeneralStaticStepIncrementation(AbaqusBaseIncrementation):
    # incrementation_type = AbaqusBaseStepIncrementationType.Automatic
    maximal_increments: int = 100
    initial_increment_size: float = 1.0
    minimal_increment_size: float = 1e-5
    maximal_increment_size: float = 1.0
    fixed_increment_size: float = None

    def __init__(self, incrementation_type=AbaqusBaseStepIncrementationType.Automatic, maximal_increments: int = 100,
                 initial_increment_size: float = 1.0, minimal_increment_size: float = 1e-5,
                 maximal_increment_size: float = 1.0, fixed_increment_size: float = None):
        if incrementation_type == AbaqusBaseStepIncrementationType.Automatic:
            self.setParameters(incrementation_type, maximal_increments, initial_increment_size, minimal_increment_size,
                               maximal_increment_size)
        elif incrementation_type == AbaqusBaseStepIncrementationType.Fixed:
            self.setParameters(incrementation_type, maximal_increments, fixed_increment_size)

    @typing.overload
    def setParameters(self, incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                      maximal_increments: int = 100, initial_increment_size: float = 1.0,
                      minimal_increment_size: float = 1e-5, maximal_increment_size: float = 1.0):
        pass

    @typing.overload
    def setParameters(self, incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                      maximal_increments: int = 100, fixed_increment_size: float = 0.1):
        pass

    def setParameters(self, incrementation_type=AbaqusBaseStepIncrementationType.Automatic, *args, **kwargs):
        if incrementation_type == AbaqusBaseStepIncrementationType.Fixed:
            keys = ['maximal_increments', 'fixed_increment_size']
            for arg, key in zip(args, keys[:len(args)]):
                kwargs[key] = arg
            maximal_increments = kwargs['maximal_increments'] if 'maximal_increments' in kwargs.keys() else 100
            fixed_increment_size = kwargs['fixed_increment_size'] if 'fixed_increment_size' in kwargs.keys() else 0.1
            self.incrementation_type = incrementation_type
            self.maximal_increments, self.fixed_increment_size = maximal_increments, fixed_increment_size
        else:
            keys = ['maximal_increments', 'initial_increment_size', 'minimal_increment_size', 'maximal_increment_size']
            for arg, key in zip(args, keys[:len(args)]):
                kwargs[key] = arg
            maximal_increments = kwargs['maximal_increments'] if 'maximal_increments' in kwargs.keys() else 100
            initial_increment_size = kwargs['initial_increment_size'] if 'initial_increment_size' in kwargs.keys() else 1.0
            minimal_increment_size = kwargs['minimal_increment_size'] if 'minimal_increment_size' in kwargs.keys() else 1e-5
            maximal_increment_size = kwargs['maximal_increment_size'] if 'maximal_increment_size' in kwargs.keys() else 1.0
            self.incrementation_type = incrementation_type
            self.maximal_increments, self.initial_increment_size = maximal_increments, initial_increment_size
            self.minimal_increment_size, self.maximal_increment_size = minimal_increment_size, maximal_increment_size


class AbaqusGeneralStaticStep(AbaqusStep):
    incrementation: AbaqusGeneralStaticStepIncrementation

    def __init__(self, time_period: float = 1.0, description: str = '',
                 incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                 maximal_increments: int = 100, initial_increment_size: float = 1.0,
                 minimal_increment_size: float = 1e-5, maximal_increment_size: float = 1.0,
                 fixed_increment_size: float = None):
        super(AbaqusGeneralStaticStep, self).__init__(AbaqusCalculationMethod.Standard, time_period, description)
        self.incrementation = AbaqusGeneralStaticStepIncrementation()
        if incrementation_type == AbaqusBaseStepIncrementationType.Automatic:
            self.incrementation.setParameters(incrementation_type, maximal_increments, initial_increment_size,
                                              minimal_increment_size, maximal_increment_size)
        elif incrementation_type == AbaqusBaseStepIncrementationType.Fixed:
            self.incrementation.setParameters(incrementation_type, maximal_increments, fixed_increment_size)

    @typing.overload
    def setParameters(self, time_period: float, description: str,
                      incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                      maximal_increments: int = 100, initial_increment_size: float = 0.1,
                      minimal_increment_size: float = 1e-5, maximal_increment_size: float = 1.0):
        pass

    @typing.overload
    def setParameters(self, time_period: float, description: str,
                      incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                      maximal_increments: int = 100, fixed_increment_size: float = 0.1):
        pass

    def setParameters(self, time_period: float, description: str,
                      incrementation_type=AbaqusBaseStepIncrementationType.Automatic, *args, **kwargs):
        super(AbaqusGeneralStaticStep, self).setParameters(AbaqusCalculationMethod.Standard, time_period, description)
        self.incrementation.setParameters(incrementation_type, *args, **kwargs)


class AbaqusExplicitDynamicStepIncrementation(AbaqusBaseIncrementation):
    # incrementation_type = AbaqusBaseStepIncrementationType.Automatic
    maximal_time_increment: typing.Union[float, str] = 'UNLIMITED'
    increment_size: typing.Union[float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'
    scaling_factor: float = 1.0

    def __init__(self, scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                 maximal_time_increment: typing.Union[float, str] = 'UNLIMITED',
                 increment_size: typing.Union[float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        if incrementation_type == AbaqusBaseStepIncrementationType.Automatic:
            self.setParameters(scaling_factor, incrementation_type, maximal_time_increment)
        elif incrementation_type == AbaqusBaseStepIncrementationType.Fixed:
            self.setParameters(scaling_factor, incrementation_type, increment_size)

    @typing.overload
    def setParameters(self, scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                      maximal_time_increment: typing.Union[float, str] = 'UNLIMITED'):
        pass

    @typing.overload
    def setParameters(self, scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                      increment_size: typing.Union[float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        pass

    def setParameters(self, scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                      *args, **kwargs):
        self.scaling_factor = scaling_factor
        if incrementation_type == AbaqusBaseStepIncrementationType.Automatic:
            self.incrementation_type = incrementation_type
            self.maximal_time_increment = args[0] if len(args) > 0 else \
                kwargs['maximal_time_increment'] if 'maximal_time_increment' in kwargs.keys() else 'UNLIMITED'
        elif incrementation_type == AbaqusBaseStepIncrementationType.Fixed:
            self.incrementation_type = incrementation_type
            self.increment_size = args[0] if len(args) > 0 else \
                kwargs['increment_size'] if 'increment_size' in kwargs.keys() else 'ELEMENT_BY_ELEMENT_INCREMENTATION'


class AbaqusExplicitDynamicStep(AbaqusStep):
    linear_bulk_viscosity: float = 0.06
    quad_bulk_viscosity: float = 1.2

    incrementation: AbaqusExplicitDynamicStepIncrementation

    def __init__(self, time_period: float, description: str,
                 linear_bulk_viscosity: float, quad_bulk_viscosity: float, scaling_factor: float = 1.0,
                 incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                 maximal_time_increment: typing.Union[float, str] = 'UNLIMITED',
                 increment_size: typing.Union[float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        super(AbaqusExplicitDynamicStep, self).__init__(AbaqusCalculationMethod.Explicit, time_period, description)
        self.incrementation = AbaqusExplicitDynamicStepIncrementation()
        if incrementation_type == AbaqusBaseStepIncrementationType.Automatic:
            self.setParameters(time_period, description, linear_bulk_viscosity,
                               quad_bulk_viscosity, scaling_factor, incrementation_type, maximal_time_increment)
        elif incrementation_type == AbaqusBaseStepIncrementationType.Fixed:
            self.setParameters(time_period, description, linear_bulk_viscosity,
                               quad_bulk_viscosity, scaling_factor, incrementation_type, increment_size)

    @typing.overload
    def setParameters(self, time_period: float, description: str,
                      linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2, scaling_factor: float = 1.0,
                      incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                      maximal_time_increment: typing.Union[float, str] = 'UNLIMITED'):
        pass

    @typing.overload
    def setParameters(self, time_period: float, description: str,
                      linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2, scaling_factor: float = 1.0,
                      incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                      increment_size: typing.Union[float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        pass

    def setParameters(self, time_period: float, description: str,
                      linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                      scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Automatic, 
                      *args, **kwargs):
        super(AbaqusExplicitDynamicStep, self).setParameters(AbaqusCalculationMethod.Explicit, time_period, description)
        self.incrementation.setParameters(scaling_factor, incrementation_type, *args, **kwargs)
        self.linear_bulk_viscosity, self.quad_bulk_viscosity = linear_bulk_viscosity, quad_bulk_viscosity
