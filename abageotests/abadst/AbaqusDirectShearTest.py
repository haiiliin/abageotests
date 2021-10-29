import os
import re
import shutil
import typing

import numpy as np
from matplotlib import pyplot as plt

from ..ababase.AbaqusGeometry import AbaqusGeometry
from ..ababase.AbaqusMesh import AbaqusMeshMethod
from ..ababase.AbaqusStep import AbaqusCalculationMethod, AbaqusBaseStepIncrementationType
from ..ababase.AbaqusSubroutine import AbaqusSubroutineType
from ..abageoabstract.AbaGeoTestBase import AbaGeoTestBase
from ..ababase.abaqusConstants import *


class AbaqusDirectShearTestGeometry(AbaqusGeometry):
    length: float = None
    width: float = None
    height: float = None

    def __init__(self, length: float = None, width: float = None, height: float = None):
        super().__init__(length, width, height)
        self.setParameters(length, width, height)

    def setParameters(self, length: float, width: float, height: float):
        super().setParameters(length, width, height)
        self.length, self.width, self.height = length, width, height


class AbaqusDirectShearTest(AbaGeoTestBase):
    """
    This is a python class to manage the information of the direct shear test, and to generate python script to build
    the abaqus model. You can use this class to

    - generate python script to build abaqus direct shear test model
    - submit the abaqus job using the abaqus command
    - generate python script to obtain the output data of the direct shear test model
    - plot some figures using the output data
    """
    geometries: dict[str, AbaqusDirectShearTestGeometry] = {}

    def __init__(self, calculation_method=AbaqusCalculationMethod.Standard,
                 soil_geometry: tuple[float, float, float] = (0.025, 0.025, 0.025),
                 solid_geometry: tuple[float, float, float] = (0.05, 0.05, 0.05),
                 soil_material: tuple[float, float, typing.Optional[float]] = (4e4, 0.0, None),
                 solid_material: tuple[float, float, typing.Optional[float]] = (2e8, 0.0, None),
                 soil_mesh_method=AbaqusMeshMethod.ByNumber, soil_seeds: typing.Union[int, float] = 1,
                 solid_mesh_method=AbaqusMeshMethod.ByNumber, solid_seeds: typing.Union[int, float] = 1,
                 vertical_pressure: float = 310,
                 predefined_stress: tuple[float, float, float, float, float, float] = (-310, -310, -310, 0, 0, 0),
                 displacement: tuple[float, float, float, float, float, float] = (0.006, 0, 0, 0, 0, 0),
                 velocity: tuple[float, float, float, float, float, float] = (1e-5, 0, 0, 0, 0, 0),
                 initial_stress_time_period: float = 1.0, shear_time_period: float = 10.0,
                 initial_stress_description: str = '', shear_description: str = '', subroutine_filepath: str = '',
                 n_state_dependent_variables: int = 0, friction_properties: tuple = None):

        super().__init__()
        self.CalculationMethod(calculation_method)
        self.SoilGeometry(*soil_geometry)
        self.SolidGeometry(*solid_geometry)
        self.SoilMaterial(*soil_material)
        self.SolidMaterial(*solid_material)
        self.SoilMesh(soil_mesh_method, soil_seeds)
        self.SolidMesh(solid_mesh_method, solid_seeds)
        self.Pressure(vertical_pressure)
        self.PredefinedStress(*predefined_stress)
        self.HistoryOutput(['U1', 'RF1'], 10)

        if calculation_method == AbaqusCalculationMethod.Standard:
            self.Displacement(*displacement)
            self.InitialStressGeneralStaticStep(time_period=initial_stress_time_period,
                                                description=initial_stress_description)
            self._GeneralStaticStep('Initial-Stress', time_period=initial_stress_time_period,
                                    description=initial_stress_description)
            self.ShearGeneralStaticStep(time_period=shear_time_period, description=shear_description)
            self.FrictionSubroutine(AbaqusSubroutineType.FRIC, subroutine_filepath=subroutine_filepath,
                                    n_state_dependent_variables=n_state_dependent_variables,
                                    friction_properties=friction_properties)
            self.FieldOutput(['S', 'E', 'LE', 'U', 'RF', 'RT', 'RM', 'P',
                              'CSTRESS', 'CDISP', 'CFORCE', 'CNAREA', 'CSTATUS'])
        else:
            self.Velocity(*velocity)
            self.InitialStressExplicitDynamicStep(time_period=initial_stress_time_period,
                                                  description=initial_stress_description)
            self.ShearExplicitDynamicStep(time_period=shear_time_period, description=shear_description)
            self.FrictionSubroutine(AbaqusSubroutineType.VFRIC, subroutine_filepath=subroutine_filepath,
                                    n_state_dependent_variables=n_state_dependent_variables,
                                    friction_properties=friction_properties)
            self.FieldOutput(['S', 'E', 'LE', 'U', 'V', 'A', 'RF', 'P',
                              'CSTRESS', 'CFORCE', 'FSLIPR', 'FSLIP', 'PPRESS'])

    def _Geometry(self, name: str = 'Soil', length: float = 1.0, width: float = 1.0, height: float = 1.0):
        """Set geometry of the soil and solid

        Parameters
        ----------
        name: str
            name of the geometry
        length: float
            length of the soil/solid
        width: float
            width of the soil/solid
        height: float
            height of the soil/solid
        Returns
        -------
        None
        """
        if name in self.geometries.keys():
            self.geometries[name].setParameters(length, width, height)
        else:
            self.geometries[name] = AbaqusDirectShearTestGeometry(length, width, height)

    def SoilGeometry(self, length: float, width: float, height: float):
        """Set geometry of the soil, identical to Geometry('Soil', length, width, height)

        Parameters
        ----------
        length: float
            length of the soil
        width: float
            width of the soil
        height: float
            height of the soil

        Returns
        -------
        None
        """
        self._Geometry('Soil', length, width, height)

    def SolidGeometry(self, length: float, width: float, height: float):
        """Set geometry of the solid, identical to Geometry('Solid', length, width, height)

        Parameters
        ----------
        length: float
            length of the solid
        width: float
            width of the solid
        height: float
            height of the solid

        Returns
        -------
        None
        """
        self._Geometry('Solid', length, width, height)

    def SoilMaterial(self, modulus: float, poisson_ratio: float, density: float = None):
        """Set soil properties, identical to Material('Soil', modulus, poisson_ratio, density)

        Parameters
        ----------
        modulus: float
            elastic modulus of the soil
        poisson_ratio: float
            poisson's ratio of the soil, should be in the range of (0, 0.5)
        density: float
            density of the soil
        Returns
        -------
        None
        """
        self._Material('Soil', modulus, poisson_ratio, density)

    def SolidMaterial(self, modulus: float, poisson_ratio: float, density: float = None):
        """Set solid properties, identical to Material('Solid', modulus, poisson_ratio, density)

        Parameters
        ----------
        modulus: float
            elastic modulus of the solid
        poisson_ratio: float
            poisson's ratio of the solid, should be in the range of (0, 0.5)
        density: float
            density of the solid
        Returns
        -------
        None
        """
        self._Material('Solid', modulus, poisson_ratio, density)

    @typing.overload
    def SoilMesh(self, method=AbaqusMeshMethod.ByNumber, seeds_number: int = 1):
        """Set mesh parameter of the soil by number of seeds, identical to Mesh('Soil', method, seeds_number)

        Parameters
        ----------
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.ByNumber
        seeds_number: int
            number of seeds in each edge
        Returns
        -------
        None
        """
        pass

    @typing.overload
    def SoilMesh(self, method=AbaqusMeshMethod.BySize, seeds_size: float = 0.01):
        """Set mesh parameter of the soil by size of seeds, identical to Mesh('Soil', method, seeds_size)

        Parameters
        ----------
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.BySize
        seeds_size: float
            spacing size of seeds in each edge
        Returns
        -------
        None
        """
        pass

    def SoilMesh(self, method=AbaqusMeshMethod.ByNumber, seeds: typing.Union[int, float] = 1):
        """Set mesh parameter of the soil by number/size of seeds, identical to Mesh('Soil', method, seeds)

        Parameters
        ----------
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.BySize
        seeds: float or int
            number of seeds or spacing size of seeds in each edge
        Returns
        -------
        None
        """
        self._GeneralMesh('Soil', method, seeds)

    @typing.overload
    def SolidMesh(self, method=AbaqusMeshMethod.ByNumber, seeds_number: int = 1):
        """Set mesh parameter of the solid by number of seeds, identical to Mesh('Solid', method, seeds_number)

        Parameters
        ----------
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.ByNumber
        seeds_number: int
            number of seeds in each edge
        Returns
        -------
        None
        """
        pass

    @typing.overload
    def SolidMesh(self, method=AbaqusMeshMethod.BySize, seeds_size: float = 0.01):
        """Set mesh parameter of the solid by size of seeds, identical to Mesh('Solid', method, seeds_size)

        Parameters
        ----------
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.BySize
        seeds_size: float
            spacing size of seeds in each edge
        Returns
        -------
        None
        """
        pass

    def SolidMesh(self, method=AbaqusMeshMethod.ByNumber, seeds: typing.Union[int, float] = 1):
        """Set mesh parameter of the solid by number/size of seeds, identical to Mesh('Solid', method, seeds)

        Parameters
        ----------
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.BySize
        seeds: float or int
            number of seeds or spacing size of seeds in each edge
        Returns
        -------
        None
        """
        self._GeneralMesh('Solid', method, seeds)

    @typing.overload
    def InitialStressGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                                       incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                       maximal_increments: int = 100, initial_increment_size: float = 1.0,
                                       minimal_increment_size: float = 1e-5, maximal_increment_size: float = 1.0):
        pass

    @typing.overload
    def InitialStressGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                                       incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                                       maximal_increments: int = 100, fixed_increment_size: float = 0.1):
        pass

    def InitialStressGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                                       incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                       *args, **kwargs):
        self._GeneralStaticStep('Initial-Stress', time_period, description, incrementation_type, *args, **kwargs)

    @typing.overload
    def ShearGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                               incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                               maximal_increments: int = 100, initial_increment_size: float = 1.0,
                               minimal_increment_size: float = 1e-5, maximal_increment_size: float = 1.0):
        pass

    @typing.overload
    def ShearGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                               incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                               maximal_increments: int = 100, fixed_increment_size: float = 0.1):
        pass

    def ShearGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                               incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                               *args, **kwargs):
        self._GeneralStaticStep('Shear', time_period, description, incrementation_type, *args, **kwargs)

    @typing.overload
    def InitialStressExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                         linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                         scaling_factor: float = 1.0,
                                         incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                         maximal_time_increment: typing.Union[float, str] = 'UNLIMITED'):
        pass

    @typing.overload
    def InitialStressExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                         linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                         scaling_factor: float = 1.0,
                                         incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                                         increment_size: typing.Union[
                                             float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        pass

    def InitialStressExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                         linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                         scaling_factor: float = 1.0,
                                         incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                         *args, **kwargs):
        self._ExplicitDynamicStep('Initial-Stress', time_period, description, linear_bulk_viscosity,
                                  quad_bulk_viscosity,
                                  scaling_factor, incrementation_type, *args, **kwargs)

    @typing.overload
    def ShearExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                 linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                 scaling_factor: float = 1.0,
                                 incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                 maximal_time_increment: typing.Union[float, str] = 'UNLIMITED'):
        pass

    @typing.overload
    def ShearExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                 linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                 scaling_factor: float = 1.0,
                                 incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                                 increment_size: typing.Union[
                                     float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        pass

    def ShearExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                 linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                 scaling_factor: float = 1.0,
                                 incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                 *args, **kwargs):
        self._ExplicitDynamicStep('Shear', time_period, description, linear_bulk_viscosity, quad_bulk_viscosity,
                                  scaling_factor, incrementation_type, *args, **kwargs)

    def generateAbaqusPythonScript(self):
        """Generate abaqus-python script file to the work directory

        Returns
        -------
        None
        """
        # TODO ABAQUS MAIN MODEL SCRIPT
        output_template_name = 'dst.py'
        template_path = os.path.join(os.path.split(os.path.abspath(__file__))[0], output_template_name)
        with open(template_path, 'r+') as file:
            template_text = file.read()

        subroutine_filename = os.path.split(self.subroutines['Friction'].subroutine_filepath)[-1]
        template_text = self._setParameters(
            template_text,
            CALCULATION_METHOD=self._calculationMethodString(),
            SOIL_LENGTH=self.geometries['Soil'].length,  # SOIL GEOMETRY
            SOIL_WIDTH=self.geometries['Soil'].width,
            SOIL_HEIGHT=self.geometries['Soil'].height,
            SOLID_LENGTH=self.geometries['Solid'].length,  # SOLID GEOMETRY
            SOLID_WIDTH=self.geometries['Solid'].width,
            SOLID_HEIGHT=self.geometries['Solid'].height,
            SOIL_ELASTIC_MODULUS=self.materials['Soil'].elastic.modulus,  # SOIL MATERIAL PROPERTIES
            SOIL_POISSON_RATIO=self.materials['Soil'].elastic.poisson_ratio,
            SOIL_DENSITY=self.materials['Soil'].density,
            SOLID_ELASTIC_MODULUS=self.materials['Solid'].elastic.modulus,  # SOLID MATERIAL PROPERTIES
            SOLID_POISSON_RATIO=self.materials['Solid'].elastic.poisson_ratio,
            SOLID_DENSITY=self.materials['Solid'].density,
            SOIL_MESH_METHOD=self.meshes['Soil'].general.methodString(),  # MESH PARAMETERS
            SOIL_SEEDS_NUMBER=self.meshes['Soil'].general.seeds_number,
            SOIL_SEEDS_SIZE=self.meshes['Soil'].general.seeds_size,
            SOLID_MESH_METHOD=self.meshes['Solid'].general.methodString(),
            SOLID_SEEDS_NUMBER=self.meshes['Solid'].general.seeds_number,
            SOLID_SEEDS_SIZE=self.meshes['Solid'].general.seeds_size,
            VERTICAL_PRESSURE=self.loads['Vertical-Pressure'].amplitude,  # VERTICAL PRESSURE
            PREDEFINED_STRESS=self.loads['Predefined-Stress'].predefined_stress,  # PREDEFINED STRESS
            SUBROUTINE=self.subroutines['Friction'].subroutineString(),  # FRIC/VFRIC SUBROUTINE
            SUBROUTINE_PATH=subroutine_filename,
            # N_STATE_DEPENDENT_VARS=self.subroutines['Friction'].n_state_dependent_variables,
            # FRICTION_PARAMETERS=self.subroutines['Friction'].properties,
            FIELD_OUTPUT_VARIABLES=self.outputs['F-Output-1'].variables,
            HISTORY_OUTPUT_VARIABLES=self.outputs['H-Output-1'].variables,
            EXPLICIT_NUM_INTERVALS=self.outputs['H-Output-1'].num_intervals,
        )
        template_text = self._setParameters_string_not_sensitive(
            template_text,
            # NORMAL CONTACT
            normal_contactStiffness=self.contacts['Contact'].normal_contact.contactStiffness,
            normal_pressureOverclosure=self.contacts['Contact'].normal_contact.pressureOverclosure,
            normal_allowSeparation=self.contacts['Contact'].normal_contact.allowSeparation,
            normal_maxStiffness=self.contacts['Contact'].normal_contact.maxStiffness,
            normal_table=self.contacts['Contact'].normal_contact.table,
            normal_constraintEnforcementMethod=self.contacts['Contact'].normal_contact.constraintEnforcementMethod,
            normal_overclosureFactor=self.contacts['Contact'].normal_contact.overclosureFactor,
            normal_overclosureMeasure=self.contacts['Contact'].normal_contact.overclosureMeasure,
            normal_contactStiffnessScaleFactor=self.contacts['Contact'].normal_contact.contactStiffnessScaleFactor,
            normal_initialStiffnessScaleFactor=self.contacts['Contact'].normal_contact.initialStiffnessScaleFactor,
            normal_clearanceAtZeroContactPressure=self.contacts[
                'Contact'].normal_contact.clearanceAtZeroContactPressure,
            normal_stiffnessBehavior=self.contacts['Contact'].normal_contact.stiffnessBehavior,
            normal_stiffnessRatio=self.contacts['Contact'].normal_contact.stiffnessRatio,
            normal_upperQuadraticFactor=self.contacts['Contact'].normal_contact.upperQuadraticFactor,
            normal_lowerQuadraticRatio=self.contacts['Contact'].normal_contact.lowerQuadraticRatio,
            # TANGENTIAL CONTACT
            tangential_formulation=self.contacts['Contact'].tangential_contact.formulation,
            tangential_directionality=self.contacts['Contact'].tangential_contact.directionality,
            tangential_slipRateDependency=self.contacts['Contact'].tangential_contact.slipRateDependency,
            tangential_pressureDependency=self.contacts['Contact'].tangential_contact.pressureDependency,
            tangential_temperatureDependency=self.contacts['Contact'].tangential_contact.temperatureDependency,
            tangential_dependencies=self.contacts['Contact'].tangential_contact.dependencies,
            tangential_exponentialDecayDefinition=self.contacts[
                'Contact'].tangential_contact.exponentialDecayDefinition,
            tangential_table=self.contacts['Contact'].tangential_contact.table,
            tangential_shearStressLimit=self.contacts['Contact'].tangential_contact.shearStressLimit,
            tangential_maximumElasticSlip=self.contacts['Contact'].tangential_contact.maximumElasticSlip,
            tangential_fraction=self.contacts['Contact'].tangential_contact.fraction,
            tangential_absoluteDistance=self.contacts['Contact'].tangential_contact.absoluteDistance,
            tangential_elasticSlipStiffness=self.contacts['Contact'].tangential_contact.elasticSlipStiffness,
            tangential_nStateDependentVars=self.contacts['Contact'].tangential_contact.nStateDependentVars,
            tangential_useProperties=self.contacts['Contact'].tangential_contact.useProperties,
        )
        if self.calculation_method == AbaqusCalculationMethod.Standard:
            template_text = self._setParameters(
                template_text, RP_DISPLACEMENT=self.loads['Displacement'].displacement,
                INITIAL_STRESS_TIME_PERIOD=self.steps['Initial-Stress'].time_period,
                INITIAL_STRESS_DESCRIPTION=self.steps['Initial-Stress'].description,
                INITIAL_STRESS_INCREMENTATION_TYPE=self.steps['Initial-Stress'].incrementation.incrementationType(),
                INITIAL_STRESS_MAXIMAL_INCREMENTS=self.steps['Initial-Stress'].incrementation.maximal_increments,
                INITIAL_STRESS_INITIAL_INCREMENT_SIZE=self.steps[
                    'Initial-Stress'].incrementation.initial_increment_size,
                INITIAL_STRESS_MINIMAL_INCREMENT_SIZE=self.steps[
                    'Initial-Stress'].incrementation.minimal_increment_size,
                INITIAL_STRESS_MAXIMAL_INCREMENT_SIZE=self.steps[
                    'Initial-Stress'].incrementation.maximal_increment_size,
                INITIAL_STRESS_FIXED_INCREMENT_SIZE=self.steps['Initial-Stress'].incrementation.fixed_increment_size,
                SHEAR_TIME_PERIOD=self.steps['Shear'].time_period,
                SHEAR_DESCRIPTION=self.steps['Shear'].description,
                SHEAR_INCREMENTATION_TYPE=self.steps['Shear'].incrementation.incrementationType(),
                SHEAR_MAXIMAL_INCREMENTS=self.steps['Shear'].incrementation.maximal_increments,
                SHEAR_INITIAL_INCREMENT_SIZE=self.steps['Shear'].incrementation.initial_increment_size,
                SHEAR_MINIMAL_INCREMENT_SIZE=self.steps['Shear'].incrementation.minimal_increment_size,
                SHEAR_MAXIMAL_INCREMENT_SIZE=self.steps['Shear'].incrementation.maximal_increment_size,
                SHEAR_FIXED_INCREMENT_SIZE=self.steps['Shear'].incrementation.fixed_increment_size,
            )
        else:
            template_text = self._setParameters(
                template_text, RP_VELOCITY=self.loads['Velocity'].velocity,
                INITIAL_STRESS_TIME_PERIOD=self.steps['Initial-Stress'].time_period,
                INITIAL_STRESS_DESCRIPTION=self.steps['Initial-Stress'].description,
                INITIAL_STRESS_MAXIMAL_TIME_INCREMENT=self.steps[
                    'Initial-Stress'].incrementation.maximal_time_increment,
                INITIAL_STRESS_INCREMENT_SIZE=self.steps['Initial-Stress'].incrementation.increment_size,
                INITIAL_STRESS_SCALING_FACTOR=self.steps['Initial-Stress'].incrementation.scaling_factor,
                INITIAL_STRESS_LINEAR_BULK_VISCOSITY=self.steps['Initial-Stress'].linear_bulk_viscosity,
                INITIAL_STRESS_QUAD_BULK_VISCOSITY=self.steps['Initial-Stress'].quad_bulk_viscosity,
                SHEAR_TIME_PERIOD=self.steps['Shear'].time_period,
                SHEAR_DESCRIPTION=self.steps['Shear'].description,
                SHEAR_MAXIMAL_TIME_INCREMENT=self.steps['Shear'].incrementation.maximal_time_increment,
                SHEAR_INCREMENT_SIZE=self.steps['Shear'].incrementation.increment_size,
                SHEAR_SCALING_FACTOR=self.steps['Shear'].incrementation.scaling_factor,
                SHEAR_LINEAR_BULK_VISCOSITY=self.steps['Shear'].linear_bulk_viscosity,
                SHEAR_QUAD_BULK_VISCOSITY=self.steps['Shear'].quad_bulk_viscosity
            )
        template_text = template_text.replace('new directory', self.work_directory)
        pattern = '\n{}.+\n'.format('WORK_DIRECTORY\s?=\s?')
        template_text = re.sub(pattern, '\n{}\n'.format("WORK_DIRECTORY = r'new directory'"), template_text)
        template_text = template_text.replace('new directory', self.work_directory)
        template_text = template_text.replace('dst.cae', self.model_name + '.cae')

        with open(os.path.join(self.work_directory, self.model_name + '.py'), 'w+') as file:
            file.write(template_text)

        # TODO ABAQUS OUTPUT SCRIPT FILE
        if self.calculation_method == AbaqusCalculationMethod.Standard:
            output_template_name = 'dst_output_fric.py'
        else:
            output_template_name = 'dst_output_vfric.py'
        output_template_path = os.path.join(os.path.split(os.path.abspath(__file__))[0], output_template_name)
        with open(output_template_path, 'r+') as file:
            output_template_text = file.read()

        output_template_text = output_template_text.replace('new directory', self.work_directory)
        pattern = '\n{}.+\n'.format('WORK_DIRECTORY\s?=\s?')
        output_template_text = re.sub(pattern, '\n{}\n'.format("WORK_DIRECTORY = r'new directory'"),
                                      output_template_text)
        output_template_text = output_template_text.replace('new directory', self.work_directory)
        output_template_text = output_template_text.replace('dst.cae', self.model_name + '.cae')

        output_template_text = self._setParameters(
            output_template_text,
            SUBROUTINE=self.subroutines['Friction'].subroutineString(),
        )

        with open(os.path.join(self.work_directory, self.output_name + '.py'), 'w+') as file:
            file.write(output_template_text)

        try:
            shutil.copy(os.path.join(self._old_work_directory, self.subroutines['Friction'].subroutine_filepath),
                        self.work_directory)

        except FileNotFoundError:
            pass
        except shutil.SameFileError:
            pass

    def plot(self):
        """Plot some figures when everything has been done properly, more specifically, the output data have been saved
        to the DATA folder in the work directory

        Returns
        -------

        """
        os.chdir(self.work_directory)
        if self.contacts['Contact'].tangential_contact_type != AbaqusTangentialContactType.UserDefined:
            return

        friction_parameters = tuple(self.contacts['Contact'].tangential_contact.properties)
        (cohesion, G0, R, e_ref0, Lambda, xi, phi_u, dd, n_p, n_d, e0, thickness) = friction_parameters

        if self.calculation_method == AbaqusCalculationMethod.Standard:
            # TODO LOAD FRIC DATA
            GAMMA_NORM = np.loadtxt('DATA/GAMMA_NORM.csv', skiprows=1, delimiter=',')
            CSHEAR1 = np.loadtxt('DATA/CSHEAR1.csv', skiprows=1, delimiter=',')
            CSLIP1 = np.loadtxt('DATA/CSLIP1.csv', skiprows=1, delimiter=',')
            U1_RP = np.loadtxt('DATA/U1_RP.csv', skiprows=1, delimiter=',')
            EPSN_FRIC = np.loadtxt('DATA/EPSN_FRIC.csv', skiprows=1, delimiter=',')
            U3_TOP = np.loadtxt('DATA/U3_TOP.csv', skiprows=1, delimiter=',')
            RF1_RP = np.loadtxt('DATA/RF1_RP.csv', skiprows=1, delimiter=',')

            FRIC_TANGENTIAL_DISPLACEMENT = GAMMA_NORM[:, 1] * thickness * 1e3
            FRIC_SHEAR_STRESS = CSHEAR1[:, 1]
            FRIC_U1_RP = abs(U1_RP[:, 1]) * 1e3
            FRIC_CSLIP1 = abs(CSLIP1[:, 1]) * 1e3
            FRIC_RF1_RP = RF1_RP[:, 1] / (self.geometries['Soil'].length * self.geometries['Soil'].width)

            FRIC_EPSN = EPSN_FRIC[:, 1] * thickness * 1e3 * (-1)
            FRIC_UTOP = (U3_TOP[:, 1] - U3_TOP[0, 1]) * 1e3

            WORK_DIRECTORY = os.getcwd()
            if not os.path.exists(os.path.join(WORK_DIRECTORY, 'Figures')):
                os.mkdir(os.path.join(WORK_DIRECTORY, 'Figures'))

            plt.figure(figsize=(8, 5))
            plt.plot(FRIC_TANGENTIAL_DISPLACEMENT, FRIC_SHEAR_STRESS, 'c-', label='FRIC (GAMMA_N - CSHEAR1)')
            plt.plot(FRIC_U1_RP, FRIC_SHEAR_STRESS, 'r:', label='FRIC (RPU1 - CSHEAR1)')
            plt.plot(FRIC_U1_RP, FRIC_RF1_RP, 'k--', label='FRIC (RPU1 - RPRF1)')
            plt.plot(FRIC_CSLIP1, FRIC_SHEAR_STRESS, 'm-', label='FRIC (CSLIP1 - CSHEAR1)')
            plt.xlabel('Tangential displacement (mm)')
            plt.ylabel('Shear stress (kPa)')
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig('Figures/TANGENTIAL-DISPLACEMENT-SHEAR-STRESS.pdf')

            plt.figure(figsize=(8, 5))
            plt.plot(FRIC_TANGENTIAL_DISPLACEMENT, FRIC_EPSN, 'c-', label='FRIC (GAMMA_N - EPSN)')
            plt.plot(FRIC_U1_RP, FRIC_EPSN, 'r:', label='FRIC (RPU1 - EPSN)')
            plt.plot(FRIC_U1_RP, FRIC_UTOP, 'm-', label='FRIC (RPU1 - UTOP)')
            plt.plot(FRIC_CSLIP1, FRIC_EPSN, 'k--', label='FRIC (CSLIP1 - EPSN)')
            plt.xlabel('Tangential displacement (mm)')
            plt.ylabel('Normal displacement (mm)')
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig('Figures/TANGENTIAL-DISPLACEMENT-NORMAL-DISPLACEMENT.pdf')
        else:
            # TODO LOAD VFRIC DATA
            # GAMMA_NORM = np.loadtxt('DATA/GAMMA_NORM.csv', skiprows=1, delimiter=',')
            CSHEAR1 = np.loadtxt('DATA/CSHEAR1.csv', skiprows=1, delimiter=',')
            # CSLIP1 = np.loadtxt('DATA/CSLIP1.csv', skiprows=1, delimiter=',')
            U1_RP = np.loadtxt('DATA/U1_RP.csv', skiprows=1, delimiter=',')
            # EPSN_FRIC = np.loadtxt('DATA/EPSN_FRIC.csv', skiprows=1, delimiter=',')
            # U3_TOP = np.loadtxt('DATA/U3_TOP.csv', skiprows=1, delimiter=',')
            RF1_RP = np.loadtxt('DATA/RF1_RP.csv', skiprows=1, delimiter=',')
            FSLIPEQ = np.loadtxt('DATA/FSLIPEQ.csv', skiprows=1, delimiter=',')

            # VFRIC_TANGENTIAL_DISPLACEMENT = GAMMA_NORM[:, 1] * thickness * 1e3
            VFRIC_SHEAR_STRESS = CSHEAR1[:, 1]
            VFRIC_U1_RP = abs(U1_RP[:, 1]) * 1e3
            # VFRIC_CSLIP1 = abs(CSLIP1[:, 1]) * 1e3
            VFRIC_RF1_RP = RF1_RP[:, 1] / (self.geometries['Soil'].length * self.geometries['Soil'].width)
            VFRIC_FSLIPEQ = FSLIPEQ[:, 1] * 1000

            # VFRIC_EPSN = EPSN_FRIC[:, 1] * thickness * 1e3 * (-1)
            # VFRIC_UTOP = (U3_TOP[:, 1] - U3_TOP[0, 1]) * 1e3

            WORK_DIRECTORY = os.getcwd()
            if not os.path.exists(os.path.join(WORK_DIRECTORY, 'Figures')):
                os.mkdir(os.path.join(WORK_DIRECTORY, 'Figures'))

            plt.figure(figsize=(8, 5))
            # plt.plot(FRIC_TANGENTIAL_DISPLACEMENT, VFRIC_SHEAR_STRESS, 'c-', label='VFRIC (GAMMA_N - CSHEAR1)')
            plt.plot(VFRIC_U1_RP, VFRIC_SHEAR_STRESS, 'r:', label='VFRIC (RPU1 - CSHEAR1)')
            plt.plot(VFRIC_U1_RP, VFRIC_RF1_RP, 'k--', label='VFRIC (RPU1 - RPRF1)')
            plt.plot(VFRIC_FSLIPEQ, VFRIC_SHEAR_STRESS, 'm-', label='VFRIC (FSLIPQ - CSHEAR1)')
            # plt.plot(FRIC_CSLIP1, VFRIC_SHEAR_STRESS, 'm-', label='VFRIC (CSLIP1 - CSHEAR1)')
            plt.xlabel('Tangential displacement (mm)')
            plt.ylabel('Shear stress (kPa)')
            plt.legend()
            plt.grid()
            plt.tight_layout()
            plt.savefig('Figures/TANGENTIAL-DISPLACEMENT-SHEAR-STRESS.pdf')

            plt.figure(figsize=(8, 5))
            # plt.plot(VFRIC_TANGENTIAL_DISPLACEMENT, VFRIC_EPSN, 'c-', label='FRIC (GAMMA_N - EPSN)')
            # plt.plot(VFRIC_U1_RP, VFRIC_EPSN, 'r:', label='VFRIC (RPU1 - EPSN)')
            # plt.plot(VFRIC_FSLIPQ, VFRIC_EPSN, 'r:', label='VFRIC (RPU1 - EPSN)')
            # plt.plot(VFRIC_U1_RP, VFRIC_UTOP, 'm-', label='VFRIC (RPU1 - UTOP)')
            # plt.plot(VFRIC_CSLIP1, VFRIC_EPSN, 'k--', label='VFRIC (CSLIP1 - EPSN)')
            plt.xlabel('Tangential displacement (mm)')
            plt.ylabel('Normal displacement (mm)')
            plt.legend()
            plt.grid()
            plt.tight_layout()
            # plt.savefig('Figures/TANGENTIAL-DISPLACEMENT-NORMAL-DISPLACEMENT.pdf')

            # plt.show()
