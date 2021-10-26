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

This is a minimal example:

```python
from abageotests import *

dst = AbaqusDirectShearTest(AbaqusCalculationMethod.Standard)

# geometries of the model
dst.SoilGeometry(0.025, 0.025, 0.025)
dst.SolidGeometry(0.05, 0.05, 0.05)
# material properties
dst.SoilMaterial(4e4, 0.0, 2e8)
dst.SolidMaterial(2e8, 0.0, 2.0)
# steps
dst.GeneralStaticStep('Initial-Stress', time_period=1.0, initial_increment_size=0.1, 
                      maximal_increment_size=1.0, maximal_increments=10000)
dst.GeneralStaticStep('Shear', time_period=1.0, initial_increment_size=0.01, 
                      maximal_increment_size=0.01, maximal_increments=1000000)
# mesh parameters
dst.SoilMesh(AbaqusMeshMethod.ByNumber, 1)
dst.SolidMesh(AbaqusMeshMethod.ByNumber, 1)
# output variables
dst.DefaultFieldOutput(['S', 'E', 'LE', 'U', 'RF', 'RT', 'RM', 'P', 
                        'CSTRESS', 'CDISP', 'CFORCE', 'CNAREA', 'CSTATUS'])
dst.DefaultHistoryOutput(['U1', 'RF1'], 50)
# loads
dst.SoilDisplacement(0.006, 0, 0, 0, 0, 0)
dst.VerticalPressure(310)
dst.PredefinedConfiningStress(-310, -310, -310, 0, 0, 0)
# friction subroutines
dst.FrictionSubroutine(AbaqusSubroutineType.FRIC, '', 28, [])
# other parameters
dst.WorkDirectory('dst-sample')
dst.ModelName('dst')
dst.OutputName('output')
dst.generateAbaqusPythonScript()
dst.submit()
dst.plot()
dst.resetWorkDirectory()
```

