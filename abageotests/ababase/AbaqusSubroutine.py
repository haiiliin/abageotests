import typing
from enum import IntEnum

from .AbaqusVariable import AbaqusVariable


class AbaqusSubroutineType(IntEnum):
    UMAT = 0
    VUMAT = 1
    FRIC = 2
    VFRIC = 3


class AbaqusSubroutine(AbaqusVariable):
    subroutine_type = AbaqusSubroutineType.FRIC
    subroutine_filepath: str = ''

    n_state_dependent_variables: int = 0
    properties: typing.Iterable = None

    def __init__(self, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = '',
                 n_state_dependent_variables: int = 0, properties: typing.Iterable = 0):
        self.setParameters(subroutine_type, subroutine_filepath, n_state_dependent_variables, properties)

    def setParameters(self, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = '',
                      n_state_dependent_variables: int = 0, properties: typing.Iterable = None):
        if properties is not None:
            properties = tuple((P, ) for P in properties)
        self.subroutine_type, self.subroutine_filepath = subroutine_type, subroutine_filepath
        self.n_state_dependent_variables, self.properties = n_state_dependent_variables, properties

    def subroutineString(self):
        if self.subroutine_type == AbaqusSubroutineType.UMAT:
            return 'UMAT'
        elif self.subroutine_type == AbaqusSubroutineType.VUMAT:
            return 'VUMAT'
        elif self.subroutine_type == AbaqusSubroutineType.FRIC:
            return 'FRIC'
        elif self.subroutine_type == AbaqusSubroutineType.VFRIC:
            return 'VFRIC'
