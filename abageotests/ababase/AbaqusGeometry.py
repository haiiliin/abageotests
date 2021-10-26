import typing

from .AbaqusVariable import AbaqusVariable


class AbaqusGeometry(AbaqusVariable):
    values: tuple[typing.Any] = None

    def __init__(self, *values: typing.Any):
        self.setParameters(*values)

    def setParameters(self, *values: typing.Any):
        self.values = values
