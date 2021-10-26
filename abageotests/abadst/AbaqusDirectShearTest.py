import os
import re
import shutil
import typing

import numpy as np
from matplotlib import pyplot as plt

from ..ababase.AbaqusGeometry import AbaqusGeometry
from ..ababase.AbaqusMesh import AbaqusMeshMethod
from ..ababase.AbaqusModelBase import AbaqusModelBase
from ..ababase.AbaqusStep import AbaqusCalculationMethod
from ..ababase.AbaqusSubroutine import AbaqusSubroutineType


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


class AbaqusDirectShearTest(AbaqusModelBase):
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
        self.setCalculationMethod(calculation_method)
        self.SoilGeometry(*soil_geometry)
        self.SolidGeometry(*solid_geometry)
        self.SoilMaterial(*soil_material)
        self.SolidMaterial(*solid_material)
        self.SoilMesh(soil_mesh_method, soil_seeds)
        self.SolidMesh(solid_mesh_method, solid_seeds)
        self.VerticalPressure(vertical_pressure)
        self.PredefinedConfiningStress(*predefined_stress)
        self.DefaultHistoryOutput(['U1', 'RF1'], 10)

        if calculation_method == AbaqusCalculationMethod.Standard:
            self.SoilDisplacement(*displacement)
            self.GeneralStaticStep('Initial-Stress', time_period=initial_stress_time_period,
                                   description=initial_stress_description)
            self.GeneralStaticStep('Shear', time_period=shear_time_period, description=shear_description)
            self.FrictionSubroutine(AbaqusSubroutineType.FRIC, subroutine_filepath=subroutine_filepath,
                                    n_state_dependent_variables=n_state_dependent_variables,
                                    friction_properties=friction_properties)
            self.DefaultFieldOutput(['S', 'E', 'LE', 'U', 'RF', 'RT', 'RM', 'P',
                                     'CSTRESS', 'CDISP', 'CFORCE', 'CNAREA', 'CSTATUS'])
        else:
            self.SoilVelocity(*velocity)
            self.ExplicitDynamicStep('Initial-Stress', time_period=initial_stress_time_period,
                                     description=initial_stress_description)
            self.ExplicitDynamicStep('Shear', time_period=shear_time_period, description=shear_description)
            self.FrictionSubroutine(AbaqusSubroutineType.VFRIC, subroutine_filepath=subroutine_filepath,
                                    n_state_dependent_variables=n_state_dependent_variables,
                                    friction_properties=friction_properties)
            self.DefaultFieldOutput(['S', 'E', 'LE', 'U', 'V', 'A', 'RF', 'P',
                                     'CSTRESS', 'CFORCE', 'FSLIPR', 'FSLIP', 'PPRESS'])

    def Geometry(self, name: str, length: float = 1.0, width: float = 1.0, height: float = 1.0):
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
        self.Geometry('Soil', length, width, height)

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
        self.Geometry('Solid', length, width, height)

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
        self.Material('Soil', modulus, poisson_ratio, density)

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
        self.Material('Solid', modulus, poisson_ratio, density)

    def VerticalPressure(self, amplitude: typing.Union[float, typing.Iterable]):
        """Set vertical pressure, identical to Pressure('Vertical-Pressure', amplitude)

        Parameters
        ----------
        amplitude: float
            amplitude of the vertical pressure

        Returns
        -------
        None
        """
        self.Pressure('Vertical-Pressure', amplitude)

    @typing.overload
    def PredefinedConfiningStress(self, predefined_stress: typing.Iterable):
        """Set predefined stress, identical to PredefinedStress('Predefined-Confining-Stress', predefined_stress)

        Parameters
        ----------
        predefined_stress: typing.Iterable
            predefined stress, an array with the dimension of 6.
        Returns
        -------
        None
        """
        pass

    @typing.overload
    def PredefinedConfiningStress(self, sig11: float, sig22: float, sig33: float, sig12: float, sig13: float,
                                  sig23: float):
        """Set predefined stress, identical to
        PredefinedStress('Predefined-Confining-Stress', sig11, sig22, sig33, sig12, sig13, sig23)

        Parameters
        ----------
        sig11: float
            component sig11 of the predefined stress
        sig22: float
            component sig22 of the predefined stress
        sig33: float
            component sig33 of the predefined stress
        sig12: float
            component sig12 of the predefined stress
        sig13: float
            component sig13 of the predefined stress
        sig23: float
            component sig23 of the predefined stress

        Returns
        -------
        None
        """
        pass

    def PredefinedConfiningStress(self, *args):
        """Set predefined stress, identical to PredefinedStress('Predefined-Confining-Stress', *args)

        Parameters
        ----------
        args: typing.Iterable or float, float, float, float, float, float
            arguments of the predefined stress, either an array with the dimension of 6, or 6 components of the predefined
            stress

        Returns
        -------
        None
        """
        self.PredefinedStress('Predefined-Confining-Stress', *args)

    @typing.overload
    def SoilDisplacement(self, displacement: typing.Iterable):
        """Set displacement of the shearing, identical to Displacement('Displacement', displacement)

        Parameters
        ----------
        displacement: typing.Iterable
            shearing displacement, an array with the dimension of 6
        Returns
        -------
        None
        """
        pass

    @typing.overload
    def SoilDisplacement(self, u1: float = 0, u2: float = 0, u3: float = 0, ur1: float = 0, ur2: float = 0,
                         ur3: float = 0):
        """Set displacement of the shearing, identical to Displacement('Displacement', u1, u2, u3, ur1, ur2, ur3)

        Parameters
        ----------
        u1: float
            component u1 of the displacement
        u2: float
            component u2 of the displacement
        u3: float
            component u3 of the displacement
        ur1: float
            component ur1 of the displacement
        ur2: float
            component ur2 of the displacement
        ur3: float
            component ur3 of the displacement

        Returns
        -------
        None
        """

        pass

    def SoilDisplacement(self, *args):
        """Set displacement of the shearing, identical to Displacement('Displacement', *args)

        Parameters
        ----------
        args: typing.Iterable or float, float, float, float, float, float
            arguments of the displacement, either an array with the dimension of 6, or 6 components of the displacement

        Returns
        -------
        None
        """
        self.Displacement('Displacement', *args)

    @typing.overload
    def SoilVelocity(self, velocity: typing.Iterable):
        """Set velocity of the shearing, identical to Velocity('Velocity', velocity)

        Parameters
        ----------
        velocity: typing.Iterable
            shearing velocity, an array with the dimension of 6
        Returns
        -------
        None
        """
        pass

    @typing.overload
    def SoilVelocity(self, v1: float = 0, v2: float = 0, v3: float = 0, vr1: float = 0, vr2: float = 0,
                     vr3: float = 0):
        """Set velocity of the shearing, identical to Velocity('Velocity', v1, v2, v3, vr1, vr2, vr3)

        Parameters
        ----------
        v1: float
            component v1 of the velocity
        v2: float
            component v2 of the velocity
        v3: float
            component v3 of the velocity
        vr1: float
            component vr1 of the velocity
        vr2: float
            component vr2 of the velocity
        vr3: float
            component vr3 of the velocity

        Returns
        -------
        None
        """
        pass

    def SoilVelocity(self, *args):
        """Set velocity of the shearing, identical to Velocity('Velocity', *args)

        Parameters
        ----------
        args: typing.Iterable or float, float, float, float, float, float
            arguments of the displacement, either an array with the dimension of 6, or 6 components of the velocity

        Returns
        -------
        None
        """
        self.Velocity('Velocity', *args)

    @typing.overload
    def SoilMesh(self, method: AbaqusMeshMethod.ByNumber, seeds_number: int):
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
    def SoilMesh(self, method: AbaqusMeshMethod.BySize, seeds_size: float):
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

    def SoilMesh(self, method: AbaqusMeshMethod, seeds: typing.Union[int, float]):
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
        self.Mesh('Soil', method, seeds)

    @typing.overload
    def SolidMesh(self, method: AbaqusMeshMethod.ByNumber, seeds_number: int):
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
    def SolidMesh(self, method: AbaqusMeshMethod.BySize, seeds_size: float):
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

    def SolidMesh(self, method: AbaqusMeshMethod, seeds: typing.Union[int, float]):
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
        self.Mesh('Solid', method, seeds)

    @typing.overload
    def FrictionSubroutine(self, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = '',
                           n_state_dependent_variables: int = 0, cohesion: float = 0.0, G0: float = 250.0,
                           R: float = 5.19, e_cref: float = 0.89, Lambda: float = 0.147, xi: float = 0.424,
                           phi: float = 31.2, dd: float = 7.57, n_p: float = 2.06, n_d: float = 0.46,
                           e0: float = 0.5771, thickness: float = 5 * 0.23 * 1e-3):
        """Set properties of the exponential interface model (FRIC/VFRIC) subroutine,

        Common properties

        - Dense: [0.0,  250.0, 5.19, 0.89, 0.147, 0.424, 31.2, 7.57, 2.06,  0.46, 0.5771, 0.00115]
        - Loose: [0.0,  290.0, 2.69, 1.20, 0.750, 0.161, 29.5, 14.2, 0.10,  0.10, 0.7700, 0.00115]
        - C:     [0.0, 2000.0, 6.96, 0.90, 0.413, 1.000, 36.4, 0.10, 5.70, 14.01, 0.5771, 0.00500]

        Parameters
        ----------
        subroutine_type: AbaqusSubroutineType
            subroutine type, either AbaqusSubroutineType.FRIC or AbaqusSubroutineType.VFRIC
        subroutine_filepath: str
            filepath of the subroutine
        n_state_dependent_variables: int
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
                           n_state_dependent_variables: int = 0, *args, **kwargs):
        """Set properties of the FRIC/VFRIC subroutine

        Parameters
        ----------
        subroutine_type: AbaqusSubroutineType
            subroutine type, either AbaqusSubroutineType.FRIC or AbaqusSubroutineType.VFRIC
        subroutine_filepath: str
            filepath of the subroutine
        n_state_dependent_variables: int
            number of dependent variables
        args, kwargs
            friction properties in the FRIC/VFRIC subroutine, if you use non-key-value way to call this function, please
            provide in correct order, if you use the key-value way to call this function, please provide all the
            properties, if you call this function in a mixed way, please provide the properties in correct order in
            non-key-value and it must be the first several properties, please provide the rest properties in key-value
            way. Whatever in what way you use, all parameters must be specified, otherwise an error would occur.
        Returns
        -------
        None
        """
        if len(args) + len(kwargs) == 1:
            friction_properties = args[0] if len(args) == 1 else kwargs['friction_properties']
        else:
            keys = ['cohesion', 'G0', 'R', 'e_cref', 'Lambda', 'xi', 'phi', 'dd', 'n_p', 'n_d', 'e0', 'thickness']
            for arg, key in zip(args, keys[:len(args)]):
                kwargs[key] = arg
            friction_properties = [kwargs[key] for key in keys]
        self.Subroutine('Friction', subroutine_type, subroutine_filepath, n_state_dependent_variables,
                        friction_properties)

    @staticmethod
    def _setParameters(text: str, **kwargs):
        for key, value in kwargs.items():
            pattern = '\n{}\s?=\s?.+\n'.format(key)
            if type(value) != str:
                text = re.sub(pattern, '\n{} = {}\n'.format(key, value), text)
            else:
                text = re.sub(pattern, "\n{} = '{}'\n".format(key, value), text)
        return text

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

        template_text = self._setParameters(
            template_text,
            CALCULATION_METHOD=self.calculationMethodString(),
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
            SOIL_MESH_METHOD=self.meshes['Soil'].methodString(),  # MESH PARAMETERS
            SOIL_SEEDS_NUMBER=self.meshes['Soil'].seeds_number,
            SOIL_SEEDS_SIZE=self.meshes['Soil'].seeds_size,
            SOLID_MESH_METHOD=self.meshes['Solid'].methodString(),
            SOLID_SEEDS_NUMBER=self.meshes['Solid'].seeds_number,
            SOLID_SEEDS_SIZE=self.meshes['Solid'].seeds_size,
            VERTICAL_PRESSURE=self.loads['Vertical-Pressure'].amplitude,  # VERTICAL PRESSURE
            PREDEFINED_STRESS=self.loads['Predefined-Confining-Stress'].predefined_stress,  # PREDEFINED STRESS
            SUBROUTINE=self.subroutines['Friction'].subroutineString(),  # FRIC/VFRIC SUBROUTINE
            SUBROUTINE_PATH=self.subroutines['Friction'].subroutine_filepath,
            N_STATE_DEPENDENT_VARS=self.subroutines['Friction'].n_state_dependent_variables,
            FRICTION_PARAMETERS=self.subroutines['Friction'].properties,
            FIELD_OUTPUT_VARIABLES=self.outputs['F-Output-1'].variables,
            HISTORY_OUTPUT_VARIABLES=self.outputs['H-Output-1'].variables, 
            EXPLICIT_NUM_INTERVALS=self.outputs['H-Output-1'].num_intervals,  # FIXME
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
                INITIAL_STRESS_MAXIMAL_TIME_INCREMENT=self.steps['Initial-Stress'].incrementation.maximal_time_increment,
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
        output_template_text = re.sub(pattern, '\n{}\n'.format("WORK_DIRECTORY = r'new directory'"), output_template_text)
        output_template_text = output_template_text.replace('new directory', self.work_directory)
        output_template_text = output_template_text.replace('dst.cae', self.model_name + '.cae')

        output_template_text = self._setParameters(
            output_template_text,
            SUBROUTINE=self.subroutines['Friction'].subroutineString(),
        )

        with open(os.path.join(self.work_directory, self.output_name + '.py'), 'w+') as file:
            file.write(output_template_text)

        try:
            shutil.copy(os.path.join(self.old_work_directory, self.subroutines['Friction'].subroutine_filepath),
                        self.work_directory)
        except FileNotFoundError:
            pass
        except shutil.SameFileError:
            pass

    def submit(self):
        """Submit the abaqus job, 3 steps will be done in this function

        - use abaqus command to run the model script file the build the abaqus model
        - use abaqus command to submit the job
        - use abaqus command to ru the output script file to obtain the output data and save it to the work directory

        Returns
        -------

        """
        os.system('abaqus cae noGUI={}.py'.format(self.model_name))
        os.system('abaqus job=Job-1 user={} cpus=1 int double'.format(self.subroutines['Friction'].subroutine_filepath))
        os.system('abaqus cae noGUI={}.py'.format(self.output_name))

        for filename in os.listdir():
            if os.path.splitext(filename)[-1] in ['.com', '.dat', '.log', '.msg', '.prt', '.rpy', '.sim', '.sta', '.1']:
                os.remove(filename)

    def plot(self):
        """Plot some figures when everything has been done properly, more specifically, the output data have been saved
        to the DATA folder in the work directory

        Returns
        -------

        """
        os.chdir(self.work_directory)

        friction_parameters = (parameter[0] for parameter in self.subroutines['Friction'].properties)
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
