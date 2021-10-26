import typing
from enum import IntEnum

from .AbaqusVariable import AbaqusVariable


class AbaqusMeshMethod(IntEnum):
    ByNumber = 0
    BySize = 1


class AbaqusMesh(AbaqusVariable):
    method = AbaqusMeshMethod.ByNumber
    seeds_number: int = 1
    seeds_size: float = 0.01

    @typing.overload
    def __init__(self, method=AbaqusMeshMethod.ByNumber, seeds_number: int = 1):
        pass

    @typing.overload
    def __init__(self, method=AbaqusMeshMethod.BySize, seeds_size: float = 0.01):
        pass

    def __init__(self, method=AbaqusMeshMethod.ByNumber, seeds: typing.Union[int, float] = 1):
        self.method = method
        if method == AbaqusMeshMethod.BySize:
            self.setParameters(method, seeds)
        else:
            self.setParameters(method, seeds)

    def methodString(self):
        if self.method == AbaqusMeshMethod.ByNumber:
            return 'NUMBER'
        else:
            return 'SIZE'

    @typing.overload
    def setParameters(self, method=AbaqusMeshMethod.ByNumber, seeds_number: int = 1):
        pass

    @typing.overload
    def setParameters(self, method=AbaqusMeshMethod.BySize, seeds_size: float = 0.01):
        pass

    def setParameters(self, method=AbaqusMeshMethod.ByNumber, seeds: typing.Union[int, float] = 1):
        self.method = method
        if method == AbaqusMeshMethod.BySize:
            self.seeds_size = seeds
        else:
            self.seeds_number = seeds
