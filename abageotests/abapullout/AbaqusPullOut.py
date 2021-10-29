import os
import re
import shutil
import typing

from abageotests import AbaqusSubroutineType

from ..ababase.AbaqusGeometry import AbaqusGeometry
from ..ababase.AbaqusMesh import AbaqusMeshMethod, AbaqusMeshDenseMethod
from ..ababase.AbaqusStep import AbaqusBaseStepIncrementationType, AbaqusCalculationMethod
from ..abageoabstract.AbaGeoTestBase import AbaGeoTestBase


class AbaqusPullOutNailGeometry(AbaqusGeometry):
    radius: float = None
    length: float = None

    def __init__(self, radius: float = None, length: float = None):
        super().__init__(radius, length)

    def setParameters(self, radius: float, length: float):
        super().setParameters(radius, length)
        self.radius, self.length = radius, length


class AbaqusPullOutSoilGeometry(AbaqusGeometry):
    length: float = None
    width: float = None
    height: float = None

    def __init__(self, length: float = None, width: float = None, height: float = None):
        super().__init__(length, width, height)
        self.setParameters(length, width, height)

    def setParameters(self, length: float, width: float, height: float):
        super().setParameters(length, width, height)
        self.length, self.width, self.height = length, width, height


class AbaqusPulloutNailOffsetGeometry(AbaqusGeometry):
    offset: float = None

    def __init__(self, offset: float = None):
        super().__init__(offset)

    def setParameters(self, offset: float):
        super().setParameters(offset)
        self.offset = offset


class AbaqusPullOut(AbaGeoTestBase):
    geometries: dict[str, typing.Union[AbaqusPullOutNailGeometry,
                                       AbaqusPullOutSoilGeometry,
                                       AbaqusPulloutNailOffsetGeometry]] = {}

    def __init__(self, calculation_method=AbaqusCalculationMethod.Standard):
        super().__init__()
        self.CalculationMethod(calculation_method)
        self.FrictionSubroutine(AbaqusSubroutineType.FRIC, subroutine_filepath='')

    @typing.overload
    def _Geometry(self, name: str = 'Nail', radius: float = 0.05, length: float = 1.2):
        """Set geometry of the nail

        Parameters
        ----------
        name: str
            name of the nail
        radius: float
            radius of the nail
        length: float
            length of the nail

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def _Geometry(self, name: str = 'Soil', length: float = 1.0, width: float = 0.3, height: float = 0.8):
        """Set geometry of the soil

        Parameters
        ----------
        name: str
            name of the soil
        length: float
            length of the soil, along with the nail
        width: float
            width of the soil
        height: float
            height of the soil, a vertical pressure will be applied on the top of the soil

        Returns
        -------
        None
        """
        pass

    @typing.overload
    def _Geometry(self, name: str = 'Offset', offset: float = 0.4):
        """Set geometry of the offset of the nail in the height dimension of the soil, usually half of the height of
        the soil

        Parameters
        ----------
        name: str
            name of the offset
        offset: float
            offset of the nail
        Returns
        -------
        None
        """
        pass

    def _Geometry(self, name: str = 'Nail', *values):
        """Set geometry of the nail, soil, and the offset of the nail

        Parameters
        ----------
        name: str
            name of the geometry
        values: float [, float [, float]]
            parameters of the geometry, radius, length for nail, length, width, height for soil, offset for offset of
            the nail
        Returns
        -------
        None
        """
        if name in self.geometries.keys():
            self.geometries[name].setParameters(*values)
        else:
            if name == 'Nail' or len(values) == 2:
                self.geometries[name] = AbaqusPullOutNailGeometry(*values)
            elif name == 'Soil' or len(values) == 3:
                self.geometries[name] = AbaqusPullOutSoilGeometry(*values)
            elif name == 'Offset' or len(values) == 1:
                self.geometries[name] = AbaqusPulloutNailOffsetGeometry(*values)

    def NailGeometry(self, radius: float, length: float):
        """Set geometry of the nail

        Parameters
        ----------
        radius: float
            radius of the nail
        length: float
            length of the nail

        Returns
        -------
        None
        """
        self._Geometry('Nail', radius, length)
        if 'Soil' in self.geometries.keys():
            self.geometries['Nail'].setParameters(radius, length)
        else:
            self.geometries['Nail'] = AbaqusPullOutNailGeometry(radius, length)

    def SoilGeometry(self, length: float, width: float, height: float):
        """Set geometry of the soil

        Parameters
        ----------
        length: float
            length of the soil, along with the nail
        width: float
            width of the soil
        height: float
            height of the soil, a vertical pressure will be applied on the top of the soil

        Returns
        -------
        None
        """
        self._Geometry('Soil', length, width, height)

    def NailOffsetGeometry(self, offset: float):
        """Set geometry of the offset of the nail in the height dimension of the soil, usually half of the height of
        the soil

        Parameters
        ----------
        offset: float
            offset of the nail
        Returns
        -------
        None
        """
        self._Geometry('Offset', offset)

    def NailMaterial(self, modulus: float, poisson_ratio: float, density: float = None):
        """Set nail properties, identical to Material('Nail', modulus, poisson_ratio, density)

        Parameters
        ----------
        modulus: float
            elastic modulus of the nail
        poisson_ratio: float
            poisson's ratio of the nail, should be in the range of (0, 0.5)
        density: float
            density of the nail
        Returns
        -------
        None
        """
        self._Material('Nail', modulus, poisson_ratio, density)

    def SoilMaterial(self, modulus: float, poisson_ratio: float, density: float = None,
                     friction_angle: float = None, dilation_angle: float = None, cohesion_yield_stress: float = None,
                     abs_plastic_strain: float = None):
        """Set soil properties, identical to Material('Soil, modulus, poisson_ratio, density, friction_angle,
        friction_angle, dilation_angle, cohesion_yield_stress, abs_plastic_strain)

        Parameters
        ----------
        modulus: float
            elastic modulus of the soil
        poisson_ratio: float
            poisson's ratio of the soil, should be in the range of (0, 0.5)
        density: float
            density of the soil
        friction_angle: float
            friction angle of the soil
        dilation_angle: float
            dilation angle of the soil
        cohesion_yield_stress: float
            cohesion yield stress of the soil
        abs_plastic_strain: float
            absolute plastic strain of the soil

        Returns
        -------

        """
        self._Material('Soil', modulus, poisson_ratio, density, friction_angle, dilation_angle,
                       cohesion_yield_stress, abs_plastic_strain)
        pass

    @typing.overload
    def SoilMesh(self, mesh_type: AbaqusMeshDenseMethod, method=AbaqusMeshMethod.ByNumber, seeds_number: int = None):
        """Set mesh parameter of the soil by number of seeds, identical to Mesh('Soil', method, seeds_number)

        Parameters
        ----------
        mesh_type: AbaqusMeshDenseMethod
            type of parameters, AbaqusMeshDenseMethod.General or AbaqusMeshDenseMethod.Dense
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
    def SoilMesh(self, mesh_type: AbaqusMeshDenseMethod, method=AbaqusMeshMethod.BySize, seeds_size: float = None):
        """Set mesh parameter of the soil by size of seeds, identical to Mesh('Soil', method, seeds_size)

        Parameters
        ----------
        mesh_type: AbaqusMeshDenseMethod
            type of parameters, AbaqusMeshDenseMethod.General or AbaqusMeshDenseMethod.Dense
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.BySize
        seeds_size: float
            spacing size of seeds in each edge
        Returns
        -------
        None
        """
        pass

    def SoilMesh(self, mesh_type: AbaqusMeshDenseMethod, method=AbaqusMeshMethod.ByNumber,
                 seeds: typing.Union[int, float] = None):
        """Set mesh parameter of the soil by number/size of seeds, identical to Mesh('Soil', method, seeds)

        Parameters
        ----------
        mesh_type: AbaqusMeshDenseMethod
            type of parameters, AbaqusMeshDenseMethod.General or AbaqusMeshDenseMethod.Dense
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.BySize
        seeds: float or int
            number of seeds or spacing size of seeds in each edge
        Returns
        -------
        None
        """
        if mesh_type == AbaqusMeshDenseMethod.General:
            self._GeneralMesh('Soil', method, seeds)
        else:
            self._DenseMesh('Soil', method, seeds)

    @typing.overload
    def NailMesh(self, mesh_type: AbaqusMeshDenseMethod, method=AbaqusMeshMethod.ByNumber, seeds_number: int = None):
        """Set mesh parameter of the solid by number of seeds, identical to Mesh('Solid', method, seeds_number)

        Parameters
        ----------
        mesh_type: AbaqusMeshDenseMethod
            type of parameters, AbaqusMeshDenseMethod.General or AbaqusMeshDenseMethod.Dense
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
    def NailMesh(self, mesh_type: AbaqusMeshDenseMethod, method=AbaqusMeshMethod.BySize, seeds_size: float = None):
        """Set mesh parameter of the solid by size of seeds, identical to Mesh('Solid', method, seeds_size)

        Parameters
        ----------
        mesh_type: AbaqusMeshDenseMethod
            type of parameters, AbaqusMeshDenseMethod.General or AbaqusMeshDenseMethod.Dense
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.BySize
        seeds_size: float
            spacing size of seeds in each edge
        Returns
        -------
        None
        """
        pass

    def NailMesh(self, mesh_type: AbaqusMeshDenseMethod, method=AbaqusMeshMethod.ByNumber,
                 seeds: typing.Union[int, float] = None):
        """Set mesh parameter of the solid by number/size of seeds, identical to Mesh('Solid', method, seeds)

        Parameters
        ----------
        mesh_type: AbaqusMeshDenseMethod
            type of parameters, AbaqusMeshDenseMethod.General or AbaqusMeshDenseMethod.Dense
        method: AbaqusMeshMethod
            mesh method, should be AbaqusMeshMethod.BySize
        seeds: float or int
            number of seeds or spacing size of seeds in each edge
        Returns
        -------
        None
        """
        if mesh_type == AbaqusMeshDenseMethod.General:
            self._GeneralMesh('Nail', method, seeds)
        else:
            self._DenseMesh('Nail', method, seeds)

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
    def PulloutGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                                 incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                 maximal_increments: int = 100, initial_increment_size: float = 1.0,
                                 minimal_increment_size: float = 1e-5, maximal_increment_size: float = 1.0):
        pass

    @typing.overload
    def PulloutGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                                 incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                                 maximal_increments: int = 100, fixed_increment_size: float = 0.1):
        pass

    def PulloutGeneralStaticStep(self, time_period: float = 1.0, description: str = '',
                                 incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                 *args, **kwargs):
        self._GeneralStaticStep('Pullout', time_period, description, incrementation_type, *args, **kwargs)

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
                                  quad_bulk_viscosity, scaling_factor, incrementation_type, *args, **kwargs)

    @typing.overload
    def PulloutExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                   linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                   scaling_factor: float = 1.0,
                                   incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                   maximal_time_increment: typing.Union[float, str] = 'UNLIMITED'):
        pass

    @typing.overload
    def PulloutExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                   linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                   scaling_factor: float = 1.0,
                                   incrementation_type=AbaqusBaseStepIncrementationType.Fixed,
                                   increment_size: typing.Union[
                                       float, str] = 'ELEMENT_BY_ELEMENT_INCREMENTATION'):
        pass

    def PulloutExplicitDynamicStep(self, time_period: float = 1.0, description: str = '',
                                   linear_bulk_viscosity: float = 0.06, quad_bulk_viscosity: float = 1.2,
                                   scaling_factor: float = 1.0,
                                   incrementation_type=AbaqusBaseStepIncrementationType.Automatic,
                                   *args, **kwargs):
        self._ExplicitDynamicStep('Pullout', time_period, description, linear_bulk_viscosity, quad_bulk_viscosity,
                                  scaling_factor, incrementation_type, *args, **kwargs)

    def generateAbaqusPythonScript(self):
        """Generate abaqus-python script file to the work directory

        Returns
        -------
        None
        """
        # TODO ABAQUS MAIN MODEL SCRIPT
        output_template_name = 'pullout.py'
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
            SOIL_NAIL_OFFSET_HEIGHT=self.geometries['Offset'].offset,
            NAIL_RADIUS=self.geometries['Nail'].radius,  # NAIL GEOMETRY
            NAIL_LENGTH=self.geometries['Nail'].length,
            SOIL_ELASTIC_MODULUS=self.materials['Soil'].elastic.modulus,  # SOIL MATERIAL PROPERTIES
            SOIL_POISSON_RATIO=self.materials['Soil'].elastic.poisson_ratio,
            SOIL_DENSITY=self.materials['Soil'].density,
            SOIL_MOHR_COULOMB_FRICTION_ANGLE=self.materials['Soil'].mohr_coulomb.friction_angle,
            SOIL_MOHR_COULOMB_DILATION_ANGLE=self.materials['Soil'].mohr_coulomb.dilation_angle,
            SOIL_MOHR_COULOMB_COHESION_YIELD_STRESS=self.materials['Soil'].mohr_coulomb.cohesion_yield_stress,
            SOIL_MOHR_COULOMB_ABS_PLASTIC_STRAIN=self.materials['Soil'].mohr_coulomb.abs_plastic_strain,
            NAIL_ELASTIC_MODULUS=self.materials['Nail'].elastic.modulus,  # NAIL MATERIAL PROPERTIES
            NAIL_POISSON_RATIO=self.materials['Nail'].elastic.poisson_ratio,
            NAIL_DENSITY=self.materials['Nail'].density,
            SOIL_MESH_METHOD=self.meshes['Soil'].general.methodString(),  # MESH PARAMETERS
            SOIL_SEEDS_GENERAL_NUMBER=self.meshes['Soil'].general.seeds_number,
            SOIL_SEEDS_DENSE_NUMBER=self.meshes['Soil'].dense.seeds_number,
            SOIL_SEEDS_GENERAL_SIZE=self.meshes['Soil'].general.seeds_size,
            SOIL_SEEDS_DENSE_SIZE=self.meshes['Soil'].dense.seeds_size,
            NAIL_MESH_METHOD=self.meshes['Nail'].general.methodString(),
            NAIL_SEEDS_GENERAL_NUMBER=self.meshes['Nail'].general.seeds_number,
            NAIL_SEEDS_DENSE_NUMBER=self.meshes['Nail'].dense.seeds_number,
            NAIL_SEEDS_GENERAL_SIZE=self.meshes['Nail'].general.seeds_size,
            NAIL_SEEDS_DENSE_SIZE=self.meshes['Nail'].dense.seeds_size,
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
                INITIAL_STRESS_INITIAL_INCREMENT_SIZE=self.steps['Initial-Stress'].incrementation.initial_increment_size,
                INITIAL_STRESS_MINIMAL_INCREMENT_SIZE=self.steps['Initial-Stress'].incrementation.minimal_increment_size,
                INITIAL_STRESS_MAXIMAL_INCREMENT_SIZE=self.steps['Initial-Stress'].incrementation.maximal_increment_size,
                INITIAL_STRESS_FIXED_INCREMENT_SIZE=self.steps['Initial-Stress'].incrementation.fixed_increment_size,
                PULLOUT_TIME_PERIOD=self.steps['Pullout'].time_period,
                PULLOUT_DESCRIPTION=self.steps['Pullout'].description,
                PULLOUT_INCREMENTATION_TYPE=self.steps['Pullout'].incrementation.incrementationType(),
                PULLOUT_MAXIMAL_INCREMENTS=self.steps['Pullout'].incrementation.maximal_increments,
                PULLOUT_INITIAL_INCREMENT_SIZE=self.steps['Pullout'].incrementation.initial_increment_size,
                PULLOUT_MINIMAL_INCREMENT_SIZE=self.steps['Pullout'].incrementation.minimal_increment_size,
                PULLOUT_MAXIMAL_INCREMENT_SIZE=self.steps['Pullout'].incrementation.maximal_increment_size,
                PULLOUT_FIXED_INCREMENT_SIZE=self.steps['Pullout'].incrementation.fixed_increment_size,
            )
        else:
            template_text = self._setParameters(
                template_text, RP_VELOCITY=self.loads['Velocity'].velocity,
                INITIAL_STRESS_TIME_PERIOD=self.steps['Initial-Stress'].time_period,
                INITIAL_STRESS_DESCRIPTION=self.steps['Initial-Stress'].description,
                INITIAL_STRESS_MAXIMAL_TIME_INCREMENT=self.steps['Initial-Stress'].incrementation.maximal_time_increment,
                INITIAL_STRESS_INCREMENT_SIZE=self.steps['Initial-Stress'].incrementation.increment_size,
                INITIAL_STRESS_SCALING_FACTOR=self.steps['Initial-Stress'].incrementation.scaling_factor,
                INITIAL_STRESS_LINEAR_BULK_VISCOSITY=self.steps['Initial-Stress'].linear_bulk_viscosity,
                INITIAL_STRESS_QUAD_BULK_VISCOSITY=self.steps['Initial-Stress'].quad_bulk_viscosity,
                PULLOUT_TIME_PERIOD=self.steps['Pullout'].time_period,
                PULLOUT_DESCRIPTION=self.steps['Pullout'].description,
                PULLOUT_MAXIMAL_TIME_INCREMENT=self.steps['Pullout'].incrementation.maximal_time_increment,
                PULLOUT_INCREMENT_SIZE=self.steps['Pullout'].incrementation.increment_size,
                PULLOUT_SCALING_FACTOR=self.steps['Pullout'].incrementation.scaling_factor,
                PULLOUT_LINEAR_BULK_VISCOSITY=self.steps['Pullout'].linear_bulk_viscosity,
                PULLOUT_QUAD_BULK_VISCOSITY=self.steps['Pullout'].quad_bulk_viscosity
            )
        template_text = template_text.replace('new directory', self.work_directory)
        pattern = '\n{}.+\n'.format('WORK_DIRECTORY\s?=\s?')
        template_text = re.sub(pattern, '\n{}\n'.format("WORK_DIRECTORY = r'new directory'"), template_text)
        template_text = template_text.replace('new directory', self.work_directory)
        template_text = template_text.replace('dst.cae', self.model_name + '.cae')

        with open(os.path.join(self.work_directory, self.model_name + '.py'), 'w+') as file:
            file.write(template_text)

        """# TODO ABAQUS OUTPUT SCRIPT FILE
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
            file.write(output_template_text)"""

        try:
            shutil.copy(os.path.join(self._old_work_directory, self.subroutines['Friction'].subroutine_filepath),
                        self.work_directory)

        except FileNotFoundError:
            pass
        except shutil.SameFileError:
            pass
