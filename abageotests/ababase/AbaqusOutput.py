from .AbaqusVariable import AbaqusVariable


class AbaqusOutput(AbaqusVariable):
    variables: list[str] = []

    def __init__(self, variables: list[str]):
        self.setParameters(variables)

    def setParameters(self, variables: list[str]):
        self.variables = variables

    def addOutputVariable(self, variable: str):
        self.variables.append(variable)

    def addOutputVariables(self, variables: list[str]):
        self.variables.extend(variables)


class AbaqusFieldOutput(AbaqusOutput):
    frequency: int = 1

    def __init__(self, variables: list[str], frequency: int = 1):
        super(AbaqusFieldOutput, self).__init__(variables)
        self.frequency = frequency

    def setParameters(self, variables: list[str], frequency: int = 1):
        super(AbaqusFieldOutput, self).setParameters(variables)
        self.frequency = frequency


class AbaqusHistoryOutput(AbaqusOutput):
    num_intervals: int = 20

    def __init__(self, variables: list[str], num_intervals: int = 1):
        super(AbaqusHistoryOutput, self).__init__(variables)
        self.num_intervals = num_intervals

    def setParameters(self, variables: list[str], num_intervals: int = 1):
        super(AbaqusHistoryOutput, self).setParameters(variables)
        self.num_intervals = num_intervals
