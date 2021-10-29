import typing

from ..ababase.AbaqusContact import AbaqusTangentialContactType
from ..ababase.AbaqusModelBase import AbaqusModelBase
from ..ababase.AbaqusSubroutine import AbaqusSubroutineType
from ..ababase.abaqusConstants import *


class AbaGeoTestBase(AbaqusModelBase):

    def Pressure(self, amplitude: typing.Union[float, typing.Iterable]):
        """Set vertical pressure, identical to Pressure('Vertical-Pressure', amplitude)

        Parameters
        ----------
        amplitude: float
            amplitude of the vertical pressure

        Returns
        -------
        None
        """
        self._Pressure('Vertical-Pressure', amplitude)

    @typing.overload
    def PredefinedStress(self, predefined_stress: typing.Iterable):
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
    def PredefinedStress(self, sig11: float, sig22: float, sig33: float, sig12: float, sig13: float,
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

    def PredefinedStress(self, *args):
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
        self._PredefinedStress('Predefined-Stress', *args)

    @typing.overload
    def Displacement(self, displacement: typing.Iterable):
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
    def Displacement(self, u1: float = 0, u2: float = 0, u3: float = 0, ur1: float = 0, ur2: float = 0,
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

    def Displacement(self, *args):
        """Set displacement of the shearing, identical to Displacement('Displacement', *args)

        Parameters
        ----------
        args: typing.Iterable or float, float, float, float, float, float
            arguments of the displacement, either an array with the dimension of 6, or 6 components of the displacement

        Returns
        -------
        None
        """
        self._Displacement('Displacement', *args)

    @typing.overload
    def Velocity(self, velocity: typing.Iterable):
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
    def Velocity(self, v1: float = 0, v2: float = 0, v3: float = 0, vr1: float = 0, vr2: float = 0,
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

    def Velocity(self, *args):
        """Set velocity of the shearing, identical to Velocity('Velocity', *args)

        Parameters
        ----------
        args: typing.Iterable or float, float, float, float, float, float
            arguments of the displacement, either an array with the dimension of 6, or 6 components of the velocity

        Returns
        -------
        None
        """
        self._Velocity('Velocity', *args)

    @typing.overload
    def ExponentialFrictionSubroutine(self, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = ''):
        pass

    @typing.overload
    def ExponentialFrictionSubroutine(self, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = ''):
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

        Returns
        -------
        None
        """
        self.FrictionSubroutine(subroutine_type, subroutine_filepath)

    def ExponentialFrictionSubroutine(self, subroutine_type=AbaqusSubroutineType.FRIC, subroutine_filepath: str = ''):
        """Set properties of the FRIC/VFRIC subroutine
        Parameters
        ----------
        subroutine_type: AbaqusSubroutineType
            subroutine type, either AbaqusSubroutineType.FRIC or AbaqusSubroutineType.VFRIC
        subroutine_filepath: str
            filepath of the subroutine
        Returns
        -------
        None
        """
        self.FrictionSubroutine(subroutine_type, subroutine_filepath)

    def generateAbaqusCaeModel(self):
        """Generate abaqus cae model

        Returns
        -------
        None
        """
        os.system('abaqus cae noGUI={}.py'.format(self.model_name))

    def submit(self):
        """Submit the abaqus job, 3 steps will be done in this function

        - use abaqus command to run the model script file the build the abaqus model
        - use abaqus command to submit the job
        - use abaqus command to ru the output script file to obtain the output data and save it to the work directory

        Returns
        -------
        None
        """
        subroutine_filename = os.path.split(self.subroutines['Friction'].subroutine_filepath)[-1]
        if subroutine_filename == '':  # TODO ADD OUTPUT SCRIPT FOR BUILD-IN SUBROUTINE
            os.system('abaqus job=Job-1 cpus=1 int double')
        else:
            os.system('abaqus job=Job-1 user={} cpus=1 int double'.format(subroutine_filename))

    def extractOutputData(self):
        """Extract output data

        Returns
        -------
        None
        """
        os.chdir(self.work_directory)
        if self.contacts['Contact'].tangential_contact_type != AbaqusTangentialContactType.UserDefined:
            return
        os.system('abaqus cae noGUI={}.py'.format(self.output_name))

        for filename in os.listdir():
            if os.path.splitext(filename)[-1] in ['.com', '.dat', '.log', '.msg', '.prt', '.rpy', '.sim', '.sta', '.1',
                                                  '.abq', '.mdl', '.pac', '.res', '.sel', '.stt']:
                os.remove(filename)
