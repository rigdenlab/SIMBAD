"""
    AMPLE.py: CCP4 GUI Project
    
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

from lxml import etree
import os
import shutil

# CCP4 imports
from CCP4PluginScript import CPluginScript

AMPLE_LOG_NODE = 'LogText'
LOGFILE_NAME = 'log.txt'
#LOGFILE_NAME = os.path.join('AMPLE_0','AMPLE.log')

class SIMBAD(CPluginScript):
    TASKNAME = 'SIMBAD'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'jens.thomas@liv.ac.uk'
    ERROR_CODES = { 1 : {'description' : 'Something not very good has happened.' },
                    }
    WHATNEXT = ['prosmart_refmac','buccaneer_build_refine','coot_rebuild']
#     PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
#                        ['log_mtzjoin.txt', 0]
#                        ]
    TASKCOMMAND="simbad-foo"

# Andre's stuff for a clean shutdown - a file caled INTERUPT will be created.
#     INTERRUPTABLE = True
#     INTERRUPTLABEL = 'Stop and keep current best solution'
    
    def __init__(self, *args, **kws):
        super(SIMBAD, self).__init__(*args, **kws)

    def processInputFiles(self):
        #Preprocess reflections to generate an "HKLIN" file
        '''
        #makeHklin0 takes as arguments a list of sublists
        #Each sublist comprises 1) A reflection data object identifier (one of those specified in the inputData container 
        #                           the task in the corresponding .def.xml
        #                       2) The requested data representation type to be placed into the file that is generated
        #
        #makeHklin0 returns a tuple comprising:
        #                       1) the file path of the file that has been created
        #                       2) a list of strings, each of which contains a comma-separated list of column labels output from
        #                       the input data objects
        #                       3) A CCP4 Error object       
        ''' 
        import CCP4XtalData
        import CCP4ErrorHandling
        # No idea why we need the 'SIMBAD_F_SIGF' bit...
        self.hklin, self.columns, error = self.makeHklin0([
                                                           ['SIMBAD_F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]
        ])
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        if self.hklin is None: return CPluginScript.FAILED
        
        self.F, self.SIGF = self.columns.split(',')
        
        #Preprocess coordinates to extract a subset
        '''
        # The method "getSelectedAtomsPdbFile" applied to a coordinate data object
        # selects those atoms declared in the objects "selectionString" property and writes them into
        # a pruned down file, the name of which is provided in the argument
        self.selectedCoordinatesPath = os.path.join(self.getWorkDirectory(), "selected_xyzin.pdb")
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.selectedCoordinatesPath)
        '''
        
        return self.SUCCEEDED

    def makeCommandAndScript(self):
        params = self.container.inputData
        if params.SIMBAD_SEARCH_LEVEL == 'Lattice':
            self.TASKCOMMAND = 'simbad-lattice'
        else: assert False
            
        # General flags
        #self.appendCommandLine(['-nproc', str(params.AMPLE_NPROC)])
        self.appendCommandLine(['-ccp4i2_xml', self.makeFileName('PROGRAMXML')])
        self.appendCommandLine(['-F', self.F])
        self.appendCommandLine(['-SIGF', self.SIGF])
        self.appendCommandLine([self.hklin])
        return self.SUCCEEDED

    def handleLogChanged(self, filename):
        with open(os.path.join(self.getWorkDirectory(),'foo.txt'),'a') as w:
            w.write('flushXML: {0}\n'.format(self.makeFileName('PROGRAMXML')))
        for ampleTxtNode in self.xmlroot.xpath(AMPLE_LOG_NODE):
            self.xmlroot.remove(ampleTxtNode)
        element = etree.SubElement(self.xmlroot,AMPLE_LOG_NODE)
        with open (filename,'r') as logFile:
            element.text = etree.CDATA(logFile.read())
        self.flushXML()

    def flushXML(self):
        tmpFilename = self.makeFileName('PROGRAMXML')+'_tmp'
        with open(tmpFilename,'w') as xmlFile:
            xmlFile.write(etree.tostring(self.xmlroot,pretty_print=True))
        if os.path.exists(self.makeFileName('PROGRAMXML')):
            os.remove(self.makeFileName('PROGRAMXML'))
        os.rename(tmpFilename, self.makeFileName('PROGRAMXML'))

    def processOutputFiles(self):
        '''
        Associate the tasks output coordinate file with the output coordinate object XYZOUT:
        '''
        
        #debug_console()
        
        # Split an MTZ file into minimtz data objects
        '''
        outputFilesToMake = ['FPHIOUT','DIFFPHIOUT']
        columnsToTake = ['FWT,PHWT','DELFWT,PHDELWT']
        infile = os.path.join(self.workDirectory,'final.mtz')
        error = self.splitHklout(outputFilesToMake, columnsToTake, infile=infile)
        import CCP4ErrorHandling
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        '''
        #Create (dummy) PROGRAMXML, which basically contains only the log text of the job
        #without this, a report will not be generated
#         logfilePath = os.path.join(self.getWorkDirectory(),LOGFILE_NAME)
#         with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
#             xmlStructure = etree.Element(AMPLE_ROOT_NODE)
#             logText = etree.SubElement(xmlStructure,AMPLE_LOG_NODE)
#             with open(logfilePath,"r") as logFile:
#                 logText.text = etree.CDATA(logFile.read())
#             programXMLFile.write(etree.tostring(xmlStructure))

        # programXML file is generated by pyrvapi so we only handle the results specific for I2 here.
        
#         if os.path.isfile(logPath):
#             with open(logPath, 'r') as logFile:
#                 element = etree.SubElement(self.xmlroot,AMPLE_LOG_NODE)
#                 element.text = etree.CDATA(logFile.read())
#         self.flushXML()

        return self.SUCCEEDED
