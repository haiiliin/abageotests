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
