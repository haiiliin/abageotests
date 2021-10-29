# abageotests
A package to generate and execute python scripts of abaqus models of geotechnical tests.

## Dependencies

In order to use this package, you need to install the `numpy` and `matplotlib` package:

```bash
pip install numpy
pip install matplotlib
```

And furthermore,  install abaqus properly, and link the FORTRAN compiler to abaqus, and add the abaqus commands directory to the system path environment.

And then install this package:

```bash
pip install abageotests
```

You can also find this project in [abageotests · PyPI](https://pypi.org/project/abageotests/).

## Examples

## An example of direct shear test

```python
from abageotests import *

dst = AbaqusDirectShearTest(AbaqusCalculationMethod.Standard)
dst.InitialStressGeneralStaticStep(
    time_period=1.0, initial_increment_size=0.1,
    maximal_increment_size=1.0, maximal_increments=10000)
dst.ShearGeneralStaticStep(
    time_period=1.0, initial_increment_size=0.01,
    maximal_increment_size=0.01, maximal_increments=1000000)

dst.Displacement(0.006, 0, 0, 0, 0, 0)
dst.NormalContact(AbaqusNormalContactType.Hard, constraintEnforcementMethod=DEFAULT)
dst.TangentialContact(AbaqusTangentialContactType.Frictionless)

dst.WorkDirectory('dst-example')

dst.SoilGeometry(0.06, 0.06, 0.02)
dst.SolidGeometry(0.1, 0.1, 0.02)
dst.SoilMesh(AbaqusMeshMethod.BySize, 0.01)
dst.SolidMesh(AbaqusMeshMethod.BySize, 0.01)

dst.SoilMaterial(4e4, 0.0, 2e8)
dst.SolidMaterial(2e8, 0.0, 2.0)

dst.FieldOutput(['S', 'E', 'LE', 'U', 'RF', 'RT', 'RM', 'P', 
                 'CSTRESS', 'CDISP', 'CFORCE', 'CNAREA', 'CSTATUS'])
dst.HistoryOutput(['U1', 'RF1'], 50)

dst.Pressure(310)
dst.PredefinedStress(-310, -310, -310, 0, 0, 0)

dst.ModelName('dst')
dst.OutputName('output')
dst.generateAbaqusPythonScript()
dst.generateAbaqusCaeModel()
dst.submit()
dst.extractOutputData()
dst.plot()
dst.resetWorkDirectory()
```

## An example of pullout test

```python
from abageotests import *


pullout = AbaqusPullOut(AbaqusCalculationMethod.Standard)
pullout.NailGeometry(0.05, 1.2)
pullout.SoilGeometry(1.0, 0.3, 0.8)
pullout.NailOffsetGeometry(0.4)

pullout.SoilMaterial(4e4, 0.3, None, 30.0, 3.0, 5.0, 0.0)
pullout.NailMaterial(3.2e7, 0.2)

pullout.Displacement(0.0, 0.0, 8e-4, 0.0, 0.0, 0.0)
pullout.Pressure(310.0)
pullout.PredefinedStress(-310.0, -310.0, -310.0, 0.0, 0.0, 0.0)

pullout.NormalContact(AbaqusNormalContactType.Hard, constraintEnforcementMethod=DEFAULT)
pullout.TangentialContact(AbaqusTangentialContactType.Frictionless)

pullout.SoilMesh(AbaqusMeshDenseMethod.General, AbaqusMeshMethod.BySize, 0.025)
pullout.SoilMesh(AbaqusMeshDenseMethod.Dense, AbaqusMeshMethod.BySize, 0.01)
pullout.NailMesh(AbaqusMeshDenseMethod.General, AbaqusMeshMethod.BySize, 0.025)
pullout.NailMesh(AbaqusMeshDenseMethod.Dense, AbaqusMeshMethod.BySize, 0.01)

pullout.InitialStressGeneralStaticStep(
    time_period=1.0, initial_increment_size=0.01, maximal_increment_size=0.1,
    minimal_increment_size=0.001, maximal_increments=1000)
pullout.PulloutGeneralStaticStep(
    time_period=1.0, initial_increment_size=0.01, maximal_increment_size=0.1,
    minimal_increment_size=0.001, maximal_increments=1000)

pullout.FieldOutput(['S', 'E', 'LE', 'U', 'V', 'A', 'RF', 'P', 
                     'CSTRESS', 'CFORCE', 'FSLIPR', 'FSLIP', 'PPRESS'])
pullout.HistoryOutput(['U1', 'RF1'], 50)

pullout.WorkDirectory('pullout-example')
pullout.ModelName('pullout')
pullout.OutputName('output')
pullout.generateAbaqusPythonScript()
pullout.generateAbaqusCaeModel()
pullout.submit()
pullout.extractOutputData()
pullout.resetWorkDirectory()
```

