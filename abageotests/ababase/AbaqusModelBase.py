import typing

from .AbaqusContact import AbaqusContact, AbaqusNormalContactType, AbaqusTangentialContactType
from .AbaqusGeometry import AbaqusGeometry
from .AbaqusLoad import AbaqusPressure, AbaqusPredefinedStress, AbaqusDisplacement, AbaqusVelocity
from .AbaqusMaterial import AbaqusMaterial
from .AbaqusMesh import AbaqusMesh, AbaqusMeshMethod
from .AbaqusStep import (AbaqusBaseStepIncrementationType, AbaqusGeneralStaticStep, AbaqusExplicitDynamicStep,
                         AbaqusCalculationMethod)
from .AbaqusSubroutine import AbaqusSubroutine, AbaqusSubroutineType
from .AbaqusOutput import AbaqusFieldOutput, AbaqusHistoryOutput
from .abaqusConstants import *


class AbaqusModelBase:
    calculation_method = AbaqusCalculationMethod.Standard
    geometries: dict[str, AbaqusGeometry] = {}
    materials: dict[str, AbaqusMaterial] = {}
    steps: dict[str, typing.Union[AbaqusGeneralStaticStep, AbaqusExplicitDynamicStep]] = {}
    loads: dict[str, typing.Union[AbaqusPressure, AbaqusPredefinedStress, AbaqusDisplacement, AbaqusVelocity]] = {}
    contacts: dict[str, AbaqusContact] = {}
    subroutines: dict[str, AbaqusSubroutine] = {}
    meshes: dict[str, AbaqusMesh] = {}
    outputs: dict[str, typing.Union[AbaqusFieldOutput, AbaqusHistoryOutput]] = {}

    _old_work_directory = os.getcwd()
    work_directory: str = os.getcwd()
    model_name = 'dst'
    output_name = 'output'

    def __init__(self):
        pass

    def WorkDirectory(self, work_directory: str):
        if work_directory == '':
            return

        if not os.path.exists(work_directory):
            os.makedirs(work_directory)
        os.chdir(work_directory)
        self.work_directory = os.getcwd()

    def ModelName(self, model_name: str):
        self.model_name = model_name

    def OutputName(self, output_name: str):
        self.output_name = output_name

    def CalculationMethod(self, calculation_method=AbaqusCalculationMethod.Standard):
        self.calculation_method = calculation_method

    def _calculationMethodString(self):
        return 'STANDARD' if self.calculation_method == AbaqusCalculationMethod.Standard else 'EXPLICIT'

    def _Geometry(self, name: str = 'Geometry', *values: typing.Any):
        if name in self.geometries.keys:
            self.geometries[name].setParameters(*values)
        else:
            self.geometries[name] = AbaqusGeometry(values)

    def _Material(self, name: str, modulus: float = None, poisson_ratio: float = None, density: float = None,
                  friction_angle: float = None, dilation_angle: float = None, cohesion_yield_stress: float = None,
                  abs_plastic_strain: float = None):
        if name in self.materials.keys():
            self.materials[name].setParameters(modulus, poisson_ratio, density, friction_angle, dilation_angle,
                                               cohesion_yield_stress, abs_plastic_strain)
        else:
            self.materials[name] = AbaqusMaterial(modulus, poisson_ratio, density, friction_angle, dilation_angle,
                                                  cohesion_yield_stress, abs_plastic_strain)

    def _Pressure(self, name: str, amplitude: typing.Union[float, typing.Iterable]):
        if name in self.loads.keys() and type(self.loads[name]) == AbaqusPressure:
            self.loads[name].setParameters(amplitude)
        else:
            self.loads[name] = AbaqusPressure(amplitude)

    @typing.overload
    def _PredefinedStress(self, name: str, predefined_stress: typing.Iterable):
        pass

    @typing.overload
    def _PredefinedStress(self, name: str, sig11: float, sig22: float, sig33: float, sig12: float, sig13: float,
                          sig23: float):
        pass

    def _PredefinedStress(self, name: str, *args):
        if name in self.loads.keys() and type(self.loads[name]) == AbaqusPredefinedStress:
            self.loads[name].setParameters(*args)
        else:
            self.loads[name] = AbaqusPredefinedStress(*args)

    @typing.overload
    def _Displacement(self, name: str, displacement: typing.Iterable = None):
        pass

    @typing.overload
    def _Displacement(self, name: str, u1: float, u2: float, u3: float, ur1: float, ur2: float, ur3: float):
        pass

    def _Displacement(self, name: str, *args):
        if name in self.loads.keys() and type(self.loads[name]) == AbaqusDisplacement:
            self.loads[name].setParameters(*args)
        else:
            self.loads[name] = AbaqusDisplacement(*args)

    @typing.overload
    def _Velocity(self, name: str, displacement: typing.Iterable):
        pass

    @typing.overload
    def _Velocity(self, name: str, v1: float, v2: float, v3: float, vr1: float, vr2: float, vr3: float):
        pass

    def _Velocity(self, name: str, *args):
        if name in self.loads.keys() and type(self.loads[name]) == AbaqusVelocity:
            self.loads[name].setParameters(*args)
        else:
            self.loads[name] = AbaqusVelocity(*args)

    @typing.overload
    def _GeneralMesh(self, name: str, method: AbaqusMeshMethod.ByNumber, seeds_number: int):
        pass

    @typing.overload
    def _GeneralMesh(self, name: str, method: AbaqusMeshMethod.BySize, seeds_size: float):
        pass

    def _GeneralMesh(self, name: str, method: AbaqusMeshMethod, seeds: typing.Union[int, float]):
        if name in self.meshes.keys():
            self.meshes[name].general.setParameters(method, seeds)
        else:
            self.meshes[name] = AbaqusMesh()
            self.meshes[name].general.setParameters(method, seeds)

    @typing.overload
    def _DenseMesh(self, name: str, method: AbaqusMeshMethod.ByNumber, seeds_number: int):
        pass

    @typing.overload
    def _DenseMesh(self, name: str, method: AbaqusMeshMethod.BySize, seeds_size: float):
        pass

    def _DenseMesh(self, name: str, method: AbaqusMeshMethod, seeds: typing.Union[int, float]):
        if name in self.meshes.keys():
            self.meshes[name].dense.setParameters(method, seeds)
        else:
            self.meshes[name] = AbaqusMesh()
            self.meshes[name].dense.setParameters(method, seeds)

    @typing.overload
    def _GeneralStaticStep(self, name: str, time_period: float = 1.0, description: str = '',
                           incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                           maximal_increments: int = 100, initial_increment_size: float = 1.0,
                           minimal_increment_size: float = 1e-5, maximal_increment_size: float = 1.0):
        pass

    @typing.overload
    def _GeneralStaticStep(self, name: str, time_period: float = 1.0, description: str = '',
                           incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                           maximal_increments: int = 100, fixed_increment_size: float = 0.1):
        pass

    def _GeneralStaticStep(self, name: str, time_period: float = 1.0, description: str = '',
                           incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                           *args, **kwargs):
        if name in self.steps.keys() and type(self.steps[name]) == AbaqusGeneralStaticStep:
            self.steps[name].setParameters(time_period, description, incrementation_type, *args, **kwargs)
        else:
            self.steps[name] = AbaqusGeneralStaticStep(time_period, description, incrementation_type, *args, **kwargs)

    @typing.overload
    def _ExplicitDynamicStep(self, name: str, time_period: float = 1.0, description: str = '',
                             linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                             scaling_factor: float = 1.0,
                             incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                             maximal_time_increment: typing.Union[float, str] = 'UNLIMITED'):
        pass

    @typing.overload
    def _ExplicitDynamicStep(self, name: str, time_period: float = 1.0, description: str = '',
                             linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                             scaling_factor: float = 1.0, incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                             increment_size: typing.Union[float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        pass

    def _ExplicitDynamicStep(self, name: str, time_period: float = 1.0, description: str = '',
                             linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                             scaling_factor: float = 1.0,
                             incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                             *args, **kwargs):
        if name in self.steps.keys() and type(self.steps[name]) == AbaqusExplicitDynamicStep:
            self.steps[name].setParameters(time_period, description, linear_bulk_viscosity, quad_bulk_viscosity,
                                           scaling_factor, incrementation_type, *args, **kwargs)
        else:
            self.steps[name] = AbaqusExplicitDynamicStep(time_period, description, linear_bulk_viscosity,
                                                         quad_bulk_viscosity, scaling_factor, incrementation_type,
                                                         *args, **kwargs)

    def _Contact(self, name: str, normal_contact_type=AbaqusNormalContactType.Hard,
                 tangential_contact_type=AbaqusTangentialContactType.UserDefined):
        if name in self.contacts.keys():
            self.contacts[name].setParameters(normal_contact_type, tangential_contact_type)
        else:
            self.contacts[name] = AbaqusContact(normal_contact_type, tangential_contact_type)

    def _NormalContact(self, name: str, contact_type: AbaqusNormalContactType, *args, **kwargs):
        if name in self.contacts.keys():
            self.contacts[name].NormalContact(contact_type, *args, **kwargs)
        else:
            self.contacts[name] = AbaqusContact()
            self.contacts[name].NormalContact(contact_type, *args, **kwargs)

    @typing.overload
    def NormalContact(self, contact_type=AbaqusNormalContactType.Hard, constraintEnforcementMethod=DEFAULT,
                      allowSeparation: Boolean = ON):
        """Hard normal contact

        Parameters
        ----------
        contact_type: AbaqusNormalContactType
            contact type, AbaqusNormalContactType.Hard
        constraintEnforcementMethod: str
            constraint enforcement method, possible values are
            'DEFAULT', 'AUGMENTED-LAGRANGE', 'PENALTY' (ONLY FOR STANDARD), 'DIRECT' (ONLY FOR STANDARD)
        allowSeparation: Boolean
            allow separation after contact

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def NormalContact(self, contact_type=AbaqusNormalContactType.Exponential, maxStiffness: float = None,
                      pressure: float = None, clearance: float = None):
        """Exponential normal contact,
        Fit exponential curve through 2 points

        Parameters
        ----------
        contact_type: AbaqusNormalContactType
            contact type, AbaqusNormalContactType.Exponential
        maxStiffness: float
            maximal stiffness, None for infinite
        pressure: float
            the first point is (pressure, 0)
        clearance: float
            the second point is (0, clearance)

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def NormalContact(self, contact_type=AbaqusNormalContactType.Linear,
                      contactStiffness: typing.Union[SymbolicConstant, float] = DEFAULT):
        """Linear normal contact

        Parameters
        ----------
        contact_type: AbaqusNormalContactType
            contact type, AbaqusNormalContactType.Linear
        contactStiffness: typing.Union[SymbolicConstant, float]
            contact stiffness

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def NormalContact(self, contact_type=AbaqusNormalContactType.Tabular, pressure: typing.Iterable = None,
                      overclosure: typing.Iterable = None):
        """Tabular normal contact,
        Provide data in order descending of overclosure, a negative overclosure is a positive clearance

        Parameters
        ----------
        contact_type: AbaqusNormalContactType
            contact type, AbaqusNormalContactType.Tabular
        pressure: typing.Iterable
            pressure data points
        overclosure: typing.Iterable
            overclosure data points

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def NormalContact(self, contact_type=AbaqusNormalContactType.ScaleFactor,
                      contactStiffnessScaleFactor: float = 1.0, initialStiffnessScaleFactor: float = 1.0,
                      overclosure_method: str = 'factor',  # or 'measure'
                      overclosure: float = 0.0):
        """Scale factor normal contact, General contact (EXPLICIT ONLY)

        Parameters
        ----------
        contact_type: AbaqusNormalContactType
            contact type, AbaqusNormalContactType.ScaleFactor
        contactStiffnessScaleFactor: float
            contact stiffness scale factor
        initialStiffnessScaleFactor: float
            initial stiffness scale factor
        overclosure_method: str
            overclosure method, 'factor' or 'measure'
        overclosure: float
            overclosure value, overclosure factor or overclosure measure

        Returns
        -------
        None
        """
        pass

    def NormalContact(self, contact_type=AbaqusNormalContactType.Hard, *args, **kwargs):
        self._NormalContact('Contact', contact_type, *args, **kwargs)

    def _TangentialContact(self, name: str, contact_type: AbaqusTangentialContactType, *args, **kwargs):
        if name in self.contacts.keys():
            self.contacts[name].TangentialContact(contact_type, *args, **kwargs)
        else:
            self.contacts[name] = AbaqusContact()
            self.contacts[name].TangentialContact(contact_type, *args, **kwargs)

    @typing.overload
    def TangentialContact(self, contact_type=AbaqusTangentialContactType.Frictionless):
        """Frictionless tangential contact

        Parameters
        ----------
        contact_type: AbaqusTangentialContactType
            contact type, AbaqusTangentialContactType.Frictionless

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def TangentialContact(self, contact_type=AbaqusTangentialContactType.Penalty, directionality=ISOTROPIC,
                          frictionCoefficient: typing.Union[float, typing.Iterable] = None, slipRate: float = None,
                          contactPressure: float = None, temp: float = None, dependencies: int = 0,
                          field: typing.Iterable = None, shearStressLimit: float = None,
                          elasticSlipStiffness: float = None,
                          maximumElasticSlip=FRACTION, maximumElasticSlipValue: float = 0.005):
        """Penalty tangential contact

        Parameters
        ----------
        contact_type: AbaqusTangentialContactType
            contact type, AbaqusTangentialContactType.Penalty
        directionality: SymbolicConstant
            directionality, ISOTROPIC or ANISOTROPIC
        frictionCoefficient: float
            friction coefficient
        slipRate: float
            use slip-rate-dependent data, None fot don't use this data
        contactPressure: float
            use contact-pressure-dependent data, None fot don't use this data
        temp: float
            use temperature-dependent data None fot don't use this data
        dependencies: int
            number of field variables
        field: typing.Iterable
            field variables
        shearStressLimit: float
            shear stress limit, None for infinite
        elasticSlipStiffness: float
            elastic slip stiffness, None for infinite
        maximumElasticSlip: SymbolicConstant
            specify maximum elastic slip, STANDARD ONLY, FRACTION or ABSOLUTE_DISTANCE
        maximumElasticSlipValue: float
            fraction or absolute distance

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def TangentialContact(self, contact_type=AbaqusTangentialContactType.ExponentialDecay,
                          exponential_decay_definition=COEFFICIENTS,
                          static_coefficient: float = None, kinetic_coefficient: float = None,
                          decay_coefficient: float = None, friction_coefficient: typing.Iterable = None,
                          slip_rate: float = None, elasticSlipStiffness: float = None,
                          maximumElasticSlip=FRACTION, maximumElasticSlipValue: float = 0.005):
        """Static-Kinetic Exponential Decay tangential contact

        Parameters
        ----------
        contact_type: AbaqusTangentialContactType
            contact type, AbaqusTangentialContactType.ExponentialDecay
        exponential_decay_definition: SymbolicConstant
            exponential decay definition method, 'COEFFICIENTS' or TEST_DATA
        static_coefficient: float
            static coefficient, when exponential_decay_definition = 'COEFFICIENTS'
        kinetic_coefficient: float
            kinetic coefficient, when exponential_decay_definition = 'COEFFICIENTS'
        decay_coefficient: float
            decay coefficient, when exponential_decay_definition = 'COEFFICIENTS'
        friction_coefficient: typing.Iterable
            friction coefficient data, when exponential_decay_definition = 'TEST_DATA'
        slip_rate: float
            slip rate data, (0, slip_rate, ∞), when exponential_decay_definition = 'TEST_DATA'
        elasticSlipStiffness: float
            elastic slip stiffness, None for infinite
        maximumElasticSlip: SymbolicConstant
            specify maximum elastic slip, STANDARD ONLY, 'FRACTION' or 'ABSOLUTE_DISTANCE'
        maximumElasticSlipValue: float
            fraction or absolute distance

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def TangentialContact(self, contact_type=AbaqusTangentialContactType.Rough):
        """Rough tangential contact

        Parameters
        ----------
        contact_type: AbaqusTangentialContactType
            contact type, AbaqusTangentialContactType.Rough

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def TangentialContact(self, contact_type=AbaqusTangentialContactType.UserDefined):
        """Rough tangential contact

        Parameters
        ----------
        contact_type: AbaqusTangentialContactType
            contact type, AbaqusTangentialContactType.UserDefined

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def TangentialContact(self, contact_type=AbaqusTangentialContactType.LagrangeMultiplier,
                          directionality=ISOTROPIC, friction_coefficient: typing.Union[float, typing.Iterable] = None,
                          slip_rate: float = None, contact_pressure: float = None, temp: float = None,
                          dependencies: int = 0, field: typing.Iterable = None, shearStressLimit: float = None):
        """Lagrange Multiplier tangential contact

        Parameters
        ----------
        contact_type: AbaqusTangentialContactType
            contact type, AbaqusTangentialContactType.LagrangeMultiplier
        directionality: str
            directionality, 'ISOTROPIC' or 'ANISOTROPIC'
        friction_coefficient: float
            friction coefficient
        slip_rate: float
            use slip-rate-dependent data, None fot don't use this data
        contact_pressure: float
            use contact-pressure-dependent data, None fot don't use this data
        temp: float
            use temperature-dependent data None fot don't use this data
        dependencies: int
            number of field variables
        field: typing.Iterable
            field variables
        shearStressLimit: float
            shear stress limit, None for infinite

        Returns
        -------

        """
        pass

    @typing.overload
    def TangentialContact(self, contact_type=AbaqusTangentialContactType.UserDefined,
                          nStateDependentVars: int = 0, cohesion: float = 0.0, G0: float = 250.0, R: float = 5.19,
                          e_cref: float = 0.89, Lambda: float = 0.147, xi: float = 0.424, phi: float = 31.2,
                          dd: float = 7.57, n_p: float = 2.06, n_d: float = 0.46, e0: float = 0.5771,
                          thickness: float = 5 * 0.23 * 1e-3):
        """Exponential FRIC/VFRIC tangential contact subroutine

        Parameters
        ----------
        contact_type: AbaqusTangentialContactType
            contact type, AbaqusTangentialContactType.UserDefined
        nStateDependentVars: int
            number of dependent variables
        cohesion: float
            cohesion of the soil
        G0: float
            initial shear modulus
        R: float
            a material constant
        e_cref: float
            reference void ratio
        Lambda: float
            parameter controlling the shape of CSL line
        xi: float
            parameter controlling the nonlinearity of CSL line
        phi: float
            friction angle
        dd: float
            a dilatancy parameter
        n_p: float
            a parameter in the critical state concept
        n_d: float
            a parameter in the critical state concept
        e0: float
            initial void ratio
        thickness: float
            thickness of the soil-structure interface

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def TangentialContact(self, contact_type=AbaqusTangentialContactType.UserDefined,
                          nStateDependentVars: int = 0,
                          properties: typing.Iterable = (0.0, 250.0, 5.19, 0.89, 0.147, 0.424,
                                                         31.2, 7.57, 2.06, 0.46, 0.5771, 0.00115)):
        """Exponential FRIC/VFRIC tangential contact subroutine

        Parameters
        ----------
        contact_type: AbaqusTangentialContactType
            contact type, AbaqusTangentialContactType.UserDefined
        nStateDependentVars: int
            number of dependent variables
        properties: typing.Iterable
            properties of the FRIC/VFRIC subroutine

        Returns
        -------
        None
        """
        pass

    def TangentialContact(self, contact_type=AbaqusTangentialContactType.Frictionless, *args, **kwargs):
        self._TangentialContact('Contact', contact_type, *args, **kwargs)

    def _Subroutine(self, name: str, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = '',
                    n_state_dependent_variables: int = 0, properties: typing.Iterable = None):
        if name in self.subroutines.keys():
            self.subroutines[name].setParameters(subroutine_type, subroutine_filepath, n_state_dependent_variables,
                                                 properties)
        else:
            self.subroutines[name] = AbaqusSubroutine(subroutine_type, subroutine_filepath, n_state_dependent_variables,
                                                      properties)

    @typing.overload
    def FrictionSubroutine(self, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = '',
                           n_state_dependent_variables: int = 0, friction_properties: typing.Iterable = None):
        """Set properties of the FRIC/VFRIC subroutine

        Parameters
        ----------
        subroutine_type: AbaqusSubroutineType
            subroutine type, either AbaqusSubroutineType.FRIC or AbaqusSubroutineType.VFRIC
        subroutine_filepath: str
            filepath of the subroutine
        n_state_dependent_variables: int
            number of dependent variables
        friction_properties: typing.Iterable
            friction properties in the FRIC/VFRIC subroutine
        Returns
        -------
        None
        """
        pass

    def FrictionSubroutine(self, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = '',
                           n_state_dependent_variables: int = 0, friction_properties: typing.Iterable = None):
        """Set properties of the FRIC/VFRIC subroutine

        Parameters
        ----------
        subroutine_type: AbaqusSubroutineType
            subroutine type, either AbaqusSubroutineType.FRIC or AbaqusSubroutineType.VFRIC
        subroutine_filepath: str
            filepath of the subroutine
        n_state_dependent_variables: int
            number of dependent variables
        friction_properties
            friction properties in the FRIC/VFRIC subroutine, if you use non-key-value way to call this function, please
            provide in correct order, if you use the key-value way to call this function, please provide all the
            properties, if you call this function in a mixed way, please provide the properties in correct order in
            non-key-value and it must be the first several properties, please provide the rest properties in key-value
            way. Whatever in what way you use, all parameters must be specified, otherwise an error would occur.
        Returns
        -------
        None
        """
        self._Subroutine('Friction', subroutine_type, subroutine_filepath, n_state_dependent_variables,
                         friction_properties)

    def _FieldOutput(self, name: str, variables: list[str], frequency: int = 1):
        if name in self.outputs.keys() and type(self.outputs[name]) == AbaqusFieldOutput:
            self.outputs[name].setParameters(variables, frequency)
        else:
            self.outputs[name] = AbaqusFieldOutput(variables, frequency)

    def FieldOutput(self, variables: list[str], frequency: int = 1):
        self._FieldOutput('F-Output-1', variables, frequency)

    def _HistoryOutput(self, name: str, variables: list[str], num_intervals: int = 10):
        if name in self.outputs.keys() and type(self.outputs[name]) == AbaqusHistoryOutput:
            self.outputs[name].setParameters(variables, num_intervals)
        else:
            self.outputs[name] = AbaqusHistoryOutput(variables, num_intervals)

    def HistoryOutput(self, variables: list[str], num_intervals: int = 10):
        self._HistoryOutput('H-Output-1', variables, num_intervals)

    def generateAbaqusPythonScript(self):
        pass

    def submit(self):
        pass

    def resetWorkDirectory(self):
        """Reset the work directory to your own outer directory, since the work directory has been set to
        AbaqusModelBase.work_directory or AbaqusDirectShearTest.work_directory that is set by the method 
        WorkDirectory or specified by the work_directory parameter.
        """
        os.chdir(self._old_work_directory)

    @staticmethod
    def _setParameters(text: str, **kwargs):
        for key, value in kwargs.items():
            pattern = '\n{}\s?=\s?.+\n'.format(key)
            if type(value) != str:
                text = re.sub(pattern, '\n{} = {}\n'.format(key, value), text)
            else:
                text = re.sub(pattern, "\n{} = '{}'\n".format(key, value), text)
        return text

    @staticmethod
    def _setParameters_string_not_sensitive(text: str, **kwargs):
        for key, value in kwargs.items():
            pattern = '\n{}\s?=\s?.+\n'.format(key)
            text = re.sub(pattern, '\n{} = {}\n'.format(key, value), text)
        return text
