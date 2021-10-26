import os
import typing

from .AbaqusGeometry import AbaqusGeometry
from .AbaqusLoad import AbaqusPressure, AbaqusPredefinedStress, AbaqusDisplacement, AbaqusVelocity
from .AbaqusMaterial import AbaqusMaterial
from .AbaqusMesh import AbaqusMesh, AbaqusMeshMethod
from .AbaqusStep import (AbaqusBaseStepIncrementationType, AbaqusGeneralStaticStep, AbaqusExplicitDynamicStep,
                         AbaqusCalculationMethod)
from .AbaqusSubroutine import AbaqusSubroutine, AbaqusSubroutineType
from.AbaqusOutput import AbaqusFieldOutput, AbaqusHistoryOutput


class AbaqusModelBase:
    calculation_method = AbaqusCalculationMethod.Standard
    geometries: dict[str, AbaqusGeometry] = {}
    materials: dict[str, AbaqusMaterial] = {}
    steps: dict[str, typing.Union[AbaqusGeneralStaticStep, AbaqusExplicitDynamicStep]] = {}
    loads: dict[str, typing.Union[AbaqusPressure, AbaqusPredefinedStress, AbaqusDisplacement, AbaqusVelocity]] = {}
    subroutines: dict[str, AbaqusSubroutine] = {}
    meshes: dict[str, AbaqusMesh] = {}
    outputs: dict[str, typing.Union[AbaqusFieldOutput, AbaqusHistoryOutput]] = {}

    old_work_directory = os.getcwd()
    work_directory: str = os.getcwd()
    model_name = 'dst'
    output_name = 'output'

    def __init__(self):
        pass

    def WorkDirectory(self, work_directory: str):
        if work_directory == '':
            return

        if not os.path.exists(work_directory):
            os.mkdir(work_directory)
        os.chdir(work_directory)
        self.work_directory = os.getcwd()

    def ModelName(self, model_name: str):
        self.model_name = model_name

    def OutputName(self, output_name: str):
        self.output_name = output_name

    def setCalculationMethod(self, calculation_method=AbaqusCalculationMethod.Standard):
        self.calculation_method = calculation_method

    def calculationMethodString(self):
        return 'STANDARD' if self.calculation_method == AbaqusCalculationMethod.Standard else 'EXPLICIT'

    def Geometry(self, name: str, *values: typing.Any):
        if name in self.geometries.keys:
            self.geometries[name].setParameters(*values)
        else:
            self.geometries[name] = AbaqusGeometry(values)

    def Material(self, name: str, modulus: float = None, poisson_ratio: float = None, density: float = None):
        if name in self.materials.keys():
            self.materials[name].setParameters(modulus, poisson_ratio, density)
        else:
            self.materials[name] = AbaqusMaterial(modulus, poisson_ratio, density)

    def Pressure(self, name: str, amplitude: typing.Union[float, typing.Iterable]):
        if name in self.loads.keys() and type(self.loads[name]) == AbaqusPressure:
            self.loads[name].setParameters(amplitude)
        else:
            self.loads[name] = AbaqusPressure(amplitude)

    @typing.overload
    def PredefinedStress(self, name: str, predefined_stress: typing.Iterable):
        pass

    @typing.overload
    def PredefinedStress(self, name: str, sig11: float, sig22: float, sig33: float, sig12: float, sig13: float,
                         sig23: float):
        pass

    def PredefinedStress(self, name: str, *args):
        if name in self.loads.keys() and type(self.loads[name]) == AbaqusPredefinedStress:
            self.loads[name].setParameters(*args)
        else:
            self.loads[name] = AbaqusPredefinedStress(*args)

    @typing.overload
    def Displacement(self, name: str, displacement: typing.Iterable = None):
        pass

    @typing.overload
    def Displacement(self, name: str, u1: float, u2: float, u3: float, ur1: float, ur2: float, ur3: float):
        pass

    def Displacement(self, name: str, *args):
        if name in self.loads.keys() and type(self.loads[name]) == AbaqusDisplacement:
            self.loads[name].setParameters(*args)
        else:
            self.loads[name] = AbaqusDisplacement(*args)

    @typing.overload
    def Velocity(self, name: str, displacement: typing.Iterable):
        pass

    @typing.overload
    def Velocity(self, name: str, v1: float, v2: float, v3: float, vr1: float, vr2: float, vr3: float):
        pass

    def Velocity(self, name: str, *args):
        if name in self.loads.keys() and type(self.loads[name]) == AbaqusVelocity:
            self.loads[name].setParameters(*args)
        else:
            self.loads[name] = AbaqusVelocity(*args)

    @typing.overload
    def Mesh(self, name: str, method: AbaqusMeshMethod.ByNumber, seeds_number: int):
        pass

    @typing.overload
    def Mesh(self, name: str, method: AbaqusMeshMethod.BySize, seeds_size: float):
        pass

    def Mesh(self, name: str, method: AbaqusMeshMethod, seeds: typing.Union[int, float]):
        if name in self.meshes.keys():
            self.meshes[name].setParameters(method, seeds)
        else:
            self.meshes[name] = AbaqusMesh(method, seeds)

    @typing.overload
    def GeneralStaticStep(self, name: str, time_period: float = 1.0, description: str = '',
                          incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                          maximal_increments: int = 100, initial_increment_size: float = 1.0,
                          minimal_increment_size: float = 1e-5, maximal_increment_size: float = 1.0):
        pass

    @typing.overload
    def GeneralStaticStep(self, name: str, time_period: float = 1.0, description: str = '',
                          incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                          maximal_increments: int = 100, fixed_increment_size: float = 0.1):
        pass

    def GeneralStaticStep(self, name: str, time_period: float = 1.0, description: str = '',
                          incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                          *args, **kwargs):
        if name in self.steps.keys() and type(self.steps[name]) == AbaqusGeneralStaticStep:
            self.steps[name].setParameters(time_period, description, incrementation_type, *args, **kwargs)
        else:
            self.steps[name] = AbaqusGeneralStaticStep(time_period, description, incrementation_type, *args, **kwargs)

    @typing.overload
    def ExplicitDynamicStep(self, name: str, time_period: float = 1.0, description: str = '',
                            linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                            scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                            maximal_time_increment: typing.Union[float, str] = 'UNLIMITED'):
        pass

    @typing.overload
    def ExplicitDynamicStep(self, name: str, time_period: float = 1.0, description: str = '',
                            linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                            scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                            increment_size: typing.Union[float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        pass

    def ExplicitDynamicStep(self, name: str, time_period: float = 1.0, description: str = '',
                            linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                            scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                            *args, **kwargs):
        if name in self.steps.keys() and type(self.steps[name]) == AbaqusExplicitDynamicStep:
            self.steps[name].setParameters(time_period, description, linear_bulk_viscosity, quad_bulk_viscosity,
                                           scaling_factor, incrementation_type, *args, **kwargs)
        else:
            self.steps[name] = AbaqusExplicitDynamicStep(time_period, description, linear_bulk_viscosity,
                                                         quad_bulk_viscosity, scaling_factor, incrementation_type,
                                                         *args, **kwargs)

    def Subroutine(self, name: str, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = '',
                   n_state_dependent_variables: int = 0, properties: typing.Iterable = None):
        if name in self.subroutines.keys():
            self.subroutines[name].setParameters(subroutine_type, subroutine_filepath, n_state_dependent_variables,
                                                 properties)
        else:
            self.subroutines[name] = AbaqusSubroutine(subroutine_type, subroutine_filepath, n_state_dependent_variables,
                                                      properties)

    def FieldOutput(self, name: str, variables: list[str], frequency: int = 1):
        if name in self.outputs.keys() and type(self.outputs[name]) == AbaqusFieldOutput:
            self.outputs[name].setParameters(variables, frequency)
        else:
            self.outputs[name] = AbaqusFieldOutput(variables, frequency)

    def DefaultFieldOutput(self, variables: list[str], frequency: int = 1):
        self.FieldOutput('F-Output-1', variables, frequency)

    def HistoryOutput(self, name: str, variables: list[str], num_intervals: int = 10):
        if name in self.outputs.keys() and type(self.outputs[name]) == AbaqusHistoryOutput:
            self.outputs[name].setParameters(variables, num_intervals)
        else:
            self.outputs[name] = AbaqusHistoryOutput(variables, num_intervals)

    def DefaultHistoryOutput(self, variables: list[str], num_intervals: int = 10):
        self.HistoryOutput('H-Output-1', variables, num_intervals)

    def generateAbaqusPythonScript(self):
        pass

    def submit(self):
        pass
    
    def resetWorkDirectory(self):
        """Reset the work wirectory to your own outer directory, since the work directory has been set to 
        AbaqusModelBase.work_directory or AbaqusDirectShearTest.work_directory that is set by the method 
        WorkDirectory or specified by the work_driectory parameter.
        """
        os.chdir(self.old_work_directory)
