import typing

from .AbaqusVariable import AbaqusVariable


class AbaqusLoad(AbaqusVariable):
    pass


class AbaqusPressure(AbaqusLoad):
    amplitude: typing.Union[float, typing.Iterable] = None

    def __init__(self, amplitude: typing.Union[float, typing.Iterable] = None):
        self.amplitude = amplitude

    def setParameters(self, amplitude: typing.Union[float, typing.Iterable] = None):
        self.amplitude = amplitude


class AbaqusPredefinedStress(AbaqusLoad):
    predefined_stress: typing.Iterable = None

    @typing.overload
    def __init__(self, predefined_stress: typing.Iterable = None):
        pass

    @typing.overload
    def __init__(self, sig11: float, sig22: float, sig33: float, sig12: float, sig13: float, sig23: float):
        pass

    def __init__(self, *args):
        self.setParameters(*args)

    @typing.overload
    def setParameters(self, predefined_stress: typing.Iterable):
        pass

    @typing.overload
    def setParameters(self, sig11: float, sig22: float, sig33: float, sig12: float, sig13: float, sig23: float):
        pass

    def setParameters(self, *args):
        if len(args) == 6:
            self.predefined_stress = [args[0], args[1], args[2], args[3], args[4], args[5]]
        else:
            self.predefined_stress = args[0]


class AbaqusDisplacement(AbaqusLoad):
    displacement: typing.Iterable = None

    @typing.overload
    def __init__(self, displacement: typing.Iterable = None):
        pass

    @typing.overload
    def __init__(self, u1: float = 0, u2: float = 0, u3: float = 0, ur1: float = 0, ur2: float = 0, ur3: float = 0):
        pass

    def __init__(self, *args):
        self.setParameters(*args)

    @typing.overload
    def setParameters(self, displacement: typing.Iterable):
        pass

    @typing.overload
    def setParameters(self, u1: float = 0, u2: float = 0, u3: float = 0, ur1: float = 0, ur2: float = 0, ur3: float = 0):
        pass

    def setParameters(self, *args):
        if len(args) == 6:
            self.displacement = [args[0], args[1], args[2], args[3], args[4], args[5]]
        else:
            self.displacement = args[0]


class AbaqusVelocity(AbaqusLoad):
    velocity: typing.Iterable = None

    @typing.overload
    def __init__(self, velocity: typing.Iterable = None):
        pass

    @typing.overload
    def __init__(self, v1: float = 0, v2: float = 0, v3: float = 0, vr1: float = 0, vr2: float = 0, vr3: float = 0):
        pass

    def __init__(self, *args):
        self.setParameters(*args)

    @typing.overload
    def setParameters(self, velocity: typing.Iterable):
        pass

    @typing.overload
    def setParameters(self, v1: float = 0, v2: float = 0, v3: float = 0, vr1: float = 0, vr2: float = 0, vr3: float = 0):
        pass

    def setParameters(self, *args):
        if len(args) == 6:
            self.velocity = [args[0], args[1], args[2], args[3], args[4], args[5]]
        else:
            self.velocity = args[0]


