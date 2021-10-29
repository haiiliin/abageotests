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
