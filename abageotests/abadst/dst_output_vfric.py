import os

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

import numpy as np

# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2021 replay file
# Internal Version: 2020_03_06-22.50.37 167380
# Run by Hailin on Wed Oct 13 16:22:48 2021
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
    
SUBROUTINE = 'FRIC'

executeOnCaeStartup()

WORK_DIRECTORY = r'C:\Users\user\OneDrive\Documents\GitHub\CH-Model\CH_model_Python\Model\SoilStructureInterface\Exponential\Abaqus\DirectShear\FRIC'
mdbPath = os.path.join(WORK_DIRECTORY, 'dst.cae')
odbPath = os.path.join(WORK_DIRECTORY, 'Job-1.odb')
os.chdir(WORK_DIRECTORY)

# OPEN OUTPUT DATABASE FILE
odb = session.openOdb(name=odbPath)
session.viewports['Viewport: 1'].setValues(displayedObject=odb)

# TODO CREATE XYDATA
# TODO CONTACT SHEAR1
session.xyDataListFromField(odb=odb, outputPosition=ELEMENT_NODAL, variable=((
    'CSHEAR1', ELEMENT_NODAL), ), elementSets=("INSTANCE-SOIL.SET-SOIL-Z0", ))

# TODO SDV8: gamma_norm
"""session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=((
    'SDV8     ASSEMBLY_INSTANCE-SOIL_SURFACE-SOIL-Z0/ASSEMBLY_INSTANCE-SOLID_SURFACE-SOLID-Z1', 
    NODAL), ), nodeSets=("INSTANCE-SOIL.SET-SOIL-Z0", ))"""

# TODO U1_RP
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
    NODAL, ((COMPONENT, 'U1'), )), ), nodeSets=("INSTANCE-SOLID.SET-RP", ))

# TODO RF1_RP
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
    NODAL, ((COMPONENT, 'RF1'), )), ), nodeSets=("INSTANCE-SOLID.SET-RP", ))

# TODO CSLIP1
"""session.xyDataListFromField(odb=odb, outputPosition=ELEMENT_NODAL, variable=((
    'CSLIP1', ELEMENT_NODAL), ), elementSets=("INSTANCE-SOIL.SET-SOIL-Z0", ))"""

# TODO EPSN_FRIC
"""session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=((
    'SDV19    ASSEMBLY_INSTANCE-SOIL_SURFACE-SOIL-Z0/ASSEMBLY_INSTANCE-SOLID_SURFACE-SOLID-Z1', 
    NODAL), ), nodeSets=("INSTANCE-SOIL.SET-SOIL-Z0", ))"""

# TODO U3_TOP
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
    NODAL, ((COMPONENT, 'U3'), )), ), nodeSets=("INSTANCE-SOIL.SET-SOIL-Z1", ))

# TODO FSLIPQ
session.xyDataListFromField(odb=odb, outputPosition=ELEMENT_NODAL, variable=((
    'FSLIPEQ  ASSEMBLY_INSTANCE-SOIL_SURFACE-SOIL-Z0/ASSEMBLY_INSTANCE-SOLID_SURFACE-SOLID-Z1', 
    ELEMENT_NODAL), ), elementSets=("INSTANCE-SOIL.SET-SOIL-Y0", ))

CSHEAR1, GAMMA_NORM, U1_RP, RF1_RP, CSLIP1, EPSN_FRIC, U3, FSLIPEQ = [], [], [], [], [], [], [], []
for key, xyDataObject in session.xyDataObjects.items():
    if key.startswith('CSHEAR1'):
        CSHEAR1.append(xyDataObject)
    elif key.startswith('SDV8'):
        GAMMA_NORM.append(xyDataObject)
    elif key.startswith('U:U1'):
        U1_RP.append(xyDataObject)
    elif key.startswith('RF:RF1'):
        RF1_RP.append(xyDataObject)
    elif key.startswith('CSLIP1'):
        CSLIP1.append(xyDataObject)
    elif key.startswith('SDV19'):
        EPSN_FRIC.append(xyDataObject)
    elif key.startswith('U:U3'):
        U3.append(xyDataObject)
    elif key.startswith('FSLIPEQ'):
        FSLIPEQ.append(xyDataObject)


def SAVE_DATA_TO_CSV(DATA_LIST, name, ABS = False):
    DATA = np.array(DATA_LIST)
    if ABS:
        DATA_AVG = abs(np.mean(DATA, axis=0))
    else:
        DATA_AVG = np.mean(DATA, axis=0)
    if not os.path.exists(os.path.join(WORK_DIRECTORY, 'DATA')):
        os.mkdir(os.path.join(WORK_DIRECTORY, 'DATA'))
    np.savetxt(os.path.join(WORK_DIRECTORY, ('DATA\\' + name + '.csv')), DATA_AVG, delimiter=',', header='TIME, ' + name)
    
SAVE_DATA_TO_CSV(CSHEAR1, 'CSHEAR1', ABS=True)
# SAVE_DATA_TO_CSV(GAMMA_NORM, 'GAMMA_NORM')
SAVE_DATA_TO_CSV(U1_RP, 'U1_RP')
SAVE_DATA_TO_CSV(RF1_RP, 'RF1_RP')
# SAVE_DATA_TO_CSV(CSLIP1, 'CSLIP1')
# SAVE_DATA_TO_CSV(EPSN_FRIC, 'EPSN_FRIC')
SAVE_DATA_TO_CSV(U3, 'U3_TOP')
SAVE_DATA_TO_CSV(FSLIPEQ, 'FSLIPEQ')
