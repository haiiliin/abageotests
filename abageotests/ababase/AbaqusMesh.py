import typing
from enum import IntEnum

from .AbaqusVariable import AbaqusVariable


class AbaqusMeshMethod(IntEnum):
    ByNumber = 0
    BySize = 1


class AbaqusMeshDenseMethod(IntEnum):
    General = 0
    Dense = 1


class AbaqusMeshBase(AbaqusVariable):
    method = AbaqusMeshMethod.ByNumber
    seeds_number: int = None
    seeds_size: float = None

    @typing.overload
    def __init__(self, method=AbaqusMeshMethod.ByNumber, seeds_number: int = None):
        pass

    @typing.overload
    def __init__(self, method=AbaqusMeshMethod.BySize, seeds_size: float = None):
        pass

    def __init__(self, method=AbaqusMeshMethod.ByNumber, seeds: typing.Union[int, float] = None):
        self.setParameters(method, seeds)

    def methodString(self):
        if self.method == AbaqusMeshMethod.ByNumber:
            return 'NUMBER'
        else:
            return 'SIZE'

    @typing.overload
    def setParameters(self, method=AbaqusMeshMethod.ByNumber, seeds_number: int = None):
        pass

    @typing.overload
    def setParameters(self, method=AbaqusMeshMethod.BySize, seeds_size: float = None):
        pass

    def setParameters(self, method=AbaqusMeshMethod.ByNumber, seeds: typing.Union[int, float] = None):
        self.method = method
        if method == AbaqusMeshMethod.BySize:
            self.seeds_size = seeds
        else:
            self.seeds_number = seeds


class AbaqusMesh(AbaqusVariable):

    dense = AbaqusMeshBase()
    general = AbaqusMeshBase()

    @typing.overload
    def setGeneralParameters(self, method=AbaqusMeshMethod.ByNumber, seeds_number: int = None):
        pass

    @typing.overload
    def setGeneralParameters(self, method=AbaqusMeshMethod.BySize, seeds_size: float = None):
        pass

    def setGeneralParameters(self, method=AbaqusMeshMethod.ByNumber, seeds: typing.Union[int, float] = None):
        self.general.method = method
        if method == AbaqusMeshMethod.BySize:
            self.general.seeds_size = seeds
        else:
            self.general.seeds_number = seeds

    @typing.overload
    def setDenseParameters(self, method=AbaqusMeshMethod.ByNumber, seeds_number: int = None):
        pass

    @typing.overload
    def setDenseParameters(self, method=AbaqusMeshMethod.BySize, seeds_size: float = None):
        pass

    def setDenseParameters(self, method=AbaqusMeshMethod.ByNumber, seeds: typing.Union[int, float] = None):
        self.dense.method = method
        if method == AbaqusMeshMethod.BySize:
            self.dense.seeds_size = seeds
        else:
            self.dense.seeds_number = seeds
