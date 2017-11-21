"""
    AMPLE_gui.py: CCP4 GUI Project
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

from CCP4TaskWidget import CTaskWidget
# David's handy debug console
from phil2etree import debug_console

from PyQt4 import QtCore
from multiprocessing import cpu_count

class SIMBAD_gui(CTaskWidget):
    """
    Draw the SIMBAD gui
    """

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'SIMBAD' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'molecular_replacement' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='SIMBAD Molecular Replacement Pipeline'
    TASKTITLE='Sequence Free Molecular Replacement - SIMBAD'
    DESCRIPTION = '''This task is for running Molecular Replacement without a sequence'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild']
    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def drawContents(self):
#         self.createLine(['subtitle','Use SHELXE'])
#         x = self.container.inputData.AMPLE_USE_SHELXE.qualifiers()['guiLabel']
#         self.createLine( ['label', x, 'widget', 'AMPLE_USE_SHELXE'])

        self.openFolder(folderFunction='inputData',followFrom=False)
        
        self.createLine(['subtitle','Input reflections'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'tip','Input reflections','widget','SIMBAD_F_SIGF' ] )
        self.closeSubFrame()

        self.createLine(['subtitle','Search level:', 'widget', 'SIMBAD_SEARCH_LEVEL'])
        self.createLine(['subtitle','Organism:', 'widget', 'SIMBAD_ORGANISM'])
        self.drawOptions()
    
    def drawOptions(self):
        folder = self.openFolder(folderFunction='inputData',title='Advanced Options')
