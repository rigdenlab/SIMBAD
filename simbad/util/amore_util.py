"""
Code to run the amore rotation search

@author: hlasimpk
"""

import copy_reg
import os
import iotbx.pdb
from iotbx.pdb import mining
import logging
import multiprocessing
import types

import simbad_util


_logger = logging.getLogger(__name__)

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

class amore(object):
    def __init__(self, optd=None):
        self.amore_exe = os.path.join(os.environ["CCP4"], "bin", 'amoreCCB2.exe')
        self.optd = optd
        self.work_dir = os.path.relpath(optd.d['work_dir'])
        return

    def amore_start(self):
        """A function to prepare files for AMORE rotation function"""

        _logger.info("Preparing files for AMORE rotation function")

        cmd = ['amore',
               'hklin', self.optd.d['mtz'],
               'hklpck0', os.path.join(self.optd.d['work_dir'], 'spmipch.hkl')]

        command_line = os.linesep.join(map(str, cmd))

        key = """TITLE   ** spmi  packing h k l F for crystal**
SORTFUN RESOL 100.  2.5
LABI FP={0}  SIGFP={1}""".format(self.optd.d['F'],
                                 self.optd.d['SIGF'])

        logfile = os.path.join(self.optd.d['work_dir'], 'SORTFUN.log')
        simbad_util.run_job(command_line, logfile, key)
        return

    def amore_run(self, models_dir):
        job_queue = multiprocessing.Queue()

        def run(job_queue):
            """processes element of job queue if queue not empty"""
            TIME_OUT_IN_SECONDS = 60

            while not job_queue.empty():
                model = job_queue.get(timeout=TIME_OUT_IN_SECONDS)
                self._amore_run(model)

            return

        for e in os.walk(models_dir):
            for model in e[2]:
                relpath = os.path.relpath(models_dir)
                job_queue.put(os.path.join(relpath, model))

        processes = []
        for i in range(self.optd.d['nproc']):
            process = multiprocessing.Process(target=run, args=(job_queue,))
            process.start()
            processes.append(process)

    def _amore_run(self, model):

        if self.optd.d['mode'] == 'CONTAM_ROT':
            self.name = os.path.basename(model)[0:6]
        elif self.optd.d['mode'] == 'FULL_ROT':
            self.name = os.path.basename(model)[0:7]

        _logger.info("Running AMORE rotation function on {0}".format(self.name))

        if self.optd.d['mode'] == 'CONTAM_ROT':
            # Set up variables for the run
            x, y, z, intrad = self.calculate_intr_box(model)

            # Run tabfun
            self.amore_tabfun(model, x, y, z)

            # Run rotfun 1
            self.amore_rotfun(intrad)

        return

    def amore_tabfun(self, model, x, y, z):
        """Function to perform AMORE table function,
        note: this is not needed if spherical harmonics pre-calculated"""

        cmd = [self.amore_exe,
               'xyzin1', model,
               'xyzout1', os.path.join(self.work_dir, 'output', '{0}_rot.pdb'.format(self.name)),
               'table1', os.path.join(self.work_dir, 'output', '{0}_sfs.tab'.format(self.name))]
        command_line = os.linesep.join(map(str, cmd))

        key = """TITLE: Produce table for MODEL FRAGMENT
TABFUN
CRYSTAL {0} {1} {2} 90 90 120 ORTH 1
MODEL 1 BTARGET 23.5
SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0""".format(x, y, z)

        logfile = os.path.join(self.work_dir, '{0}_tabfun.log'.format(self.name))

        simbad_util.run_job(command_line, logfile, key)

        # Clean up
        os.remove(logfile)

        return

    def amore_rotfun(self, intrad):
        """Function to perform first AMORE rotation function,
        note: this is not needed if spherical harmonics pre-calculated"""

        cmd = [self.amore_exe,
               'table1', os.path.join(self.work_dir, 'output', '{0}_sfs.tab'.format(self.name)),
               'HKLPCK1', os.path.join(self.work_dir, 'output', '{0}.hkl'.format(self.name)),
               'hklpck0', os.path.join(self.work_dir, 'spmipch.hkl'),
               'clmn1', os.path.join(self.work_dir, 'output', '{0}.clmn'.format(self.name)),
               'clmn0', os.path.join(self.work_dir, 'output', '{0}_spmipch.clmn'.format(self.name)),
               'MAPOUT', os.path.join(self.work_dir, 'output', 'amore_cross.map')]
        command_line = os.linesep.join(map(str, cmd))

        key = """ROTFUN
TITLE: Generate HKLPCK1 from MODEL FRAGMENT 1
GENE 1   RESO 100.0 {0}  CELL_MODEL 80 75 65
CLMN CRYSTAL ORTH  1 RESO  20.0  {0}  SPHERE   {1}
CLMN MODEL 1     RESO  20.0  {0} SPHERE   {1}
ROTA  CROSS  MODEL 1  PKLIM {2}  NPIC {3} STEP {4}""".format(self.optd.d['SHRES'],
                                                             intrad,
                                                             self.optd.d['PKLIM'],
                                                             self.optd.d['NPIC'],
                                                             self.optd.d['ROTASTEP'])

        logfile = os.path.join(self.work_dir, 'clogs', '{0}.log'.format(self.name))

        simbad_util.run_job(command_line, logfile, key)

        return

    def calculate_intr_box(self, model):
        """Function to calculate the integration radius or minimal box for an input PDB"""
        pdb_input = iotbx.pdb.pdb_input(file_name=model)
        hierarchy = pdb_input.construct_hierarchy()

        # Get resolution
        x = pdb_input.extract_remark_iii_records(2)
        resolution = mining.extract_best_resolution(x)

        # Set a default resolution if mining fails
        if resolution == None:
            resolution = 2.0

        # Get a list of all xyz coordinates
        x, y, z = [], [], []

        for residue_group in hierarchy.models()[0].chains()[0].residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms():
                    x.append(atom.xyz[0])
                    y.append(atom.xyz[1])
                    z.append(atom.xyz[2])

        # Get the smallest box containing the model
        xdiff = max(x) - min(x)
        ydiff = max(y) - min(y)
        zdiff = max(z) - min(z)

        # Get integration radius (note, for spherical structure)
        intrad = min(xdiff, ydiff, zdiff) * 0.75

        # Add together for each coordinate
        x = xdiff + intrad + resolution
        y = ydiff + intrad + resolution
        z = zdiff + intrad + resolution

        return x, y, z, intrad

    def matthews_coef(self, model):
        """Function to run matthews coefficient to decide if the model can fit in the unit cell"""

        # Get the molecular weight of the input model
        molecular_weight = self.rwcontents(model)

        cmd = ["matthews_coef"]
        key = """CELL {0}
        symm {1}
        molweight {2}
        auto""".format(self.optd.d['cell_paramaters'],
                       self.optd.d['space_group'],
                       molecular_weight)
        name = os.basename(model)[0:3]
        logfile = os.path.join(self.optd.d['work_dir'], 'matt_coef_{0}.log'.format(name))
        simbad_util.run_job(cmd, logfile, key)

        # Determine if the model can fit in the unit cell
        with open(logfile, 'r') as f:
            for line in f:
                if line.startswith('  1'):
                    solvent_content = float(line.split()[2])
                    if solvent_content >= 30:
                        result = True
                    else:
                        result = False

        # Clean up
        os.remove(logfile)

        return result

    def rwcontents(self, model):
        """Function to run rwcontents to get the molecular weight of a model"""

        cmd = ['rwcontents',
               'xyzin', model]

        name = os.basename(model)[0:3]
        logfile = 'rwcontents_{0}.log'.format(name)
        simbad_util.run_job(cmd, logfile)

        # Exctract molecular weight from log file
        molecular_weight = None
        with open(logfile, 'r') as f:
            for line in f:
                if line.startswith(" Molecular Weight of protein"):
                    molecular_weight = float(line.split()[-1])
        if not molecular_weight:
            msg = "Cannot find Molecular weight in logfile {0]".format(logfile)
            _logger.debug(msg)
            raise RuntimeError(msg)

        # Clean up
        os.remove(logfile)

        return molecular_weight
