"""Module to run the AMORE rotation search"""

__author__ = "Adam Simpkin & Felix Simkovic"
__date__ = "21 May 2017"
__version__ = "0.2"

import base64
import logging
import glob
import numpy
import os
import shutil
import zlib

from pyjob import Job, cexec
from pyjob.misc import make_script, tmp_file

from simbad.parsers import rotsearch_parser
from simbad.util import mtz_util
from simbad.util import molecular_weight

import cctbx.crystal
import iotbx.pdb
import iotbx.pdb.mining
import mmtbx.scaling.matthews

logger = logging.getLogger(__name__)

EXPORT = "SET" if os.name == "nt" else "export"


class _AmoreRotationScore(object):
    """An amore rotation scoring class"""

    __slots__ = ("pdb_code", "ALPHA", "BETA", "GAMMA", "CC_F", "RF_F", "CC_I", "CC_P", "Icp",
                 "CC_F_Z_score", "CC_P_Z_score", "Number_of_rotation_searches_producing_peak")

    def __init__(self, pdb_code, ALPHA, BETA, GAMMA, CC_F, RF_F, CC_I, CC_P, Icp,
                 CC_F_Z_score, CC_P_Z_score, Number_of_rotation_searches_producing_peak):
        self.pdb_code = pdb_code
        self.ALPHA = ALPHA
        self.BETA = BETA
        self.GAMMA = GAMMA
        self.CC_F = CC_F
        self.RF_F = RF_F
        self.CC_I = CC_I
        self.CC_P = CC_P
        self.Icp = Icp
        self.CC_F_Z_score = CC_F_Z_score
        self.CC_P_Z_score = CC_P_Z_score
        self.Number_of_rotation_searches_producing_peak = Number_of_rotation_searches_producing_peak

    def __repr__(self):
        string = "{name}(pdb_code={pdb_code} ALPHA={ALPHA} BETA={BETA} GAMMA={GAMMA} CC_F=CC_F RF_F={RF_F} " \
                 "CC_I={CC_I} CC_P={CC_P} Icp={Icp} CC_F_Z_score={CC_F_Z_score} CC_P_Z_score={CC_P_Z_score} " \
                 "Number_of_rotation_searches_producing_peak={Number_of_rotation_searches_producing_peak})"
        return string.format(name=self.__class__.__name__, **{k: getattr(self, k) for k in self.__slots__})

    def _as_dict(self):
        """Convert the :obj:`_AmoreRotationScore <simbad.rotsearch.amore_search._AmoreRotationScore>`
        object to a dictionary"""
        return {k: getattr(self, k) for k in self.__slots__}


class AmoreRotationSearch(object):
    """A class to perform the amore rotation search

    Attributes
    ----------
    amore_exe : str
        The path to the amore executable
    mtz : str
        The path to the input MTZ
    work_dir : str
        The path to the working directory
    max_to_keep : int
        The maximum number of results to keep [default: 20]

    Examples
    --------
    >>> from simbad.rotsearch.amore_search import AmoreRotationSearch
    >>> rotation_search = AmoreRotationSearch('<amore_exe>', '<mtz>', '<work_dir>', '<max_to_keep>')
    >>> rotation_search.sortfun()
    >>> rotation_search.run_pdb(
    ...     '<models_dir>', '<output_dir>', '<nproc>', '<shres>', '<pklim>', '<npic>', '<rotastep>',
    ...     '<min_solvent_content>', '<submit_qtype>', '<submit_queue>', '<monitor>', '<chunk_size>'
    ... )
    >>> rotation_search.summarize()
    >>> search_results = rotation_search.search_results


    If any results are found, an object is returned containing the pdb_code, and the various associated scores
    from amore.

    """

    def __init__(self, amore_exe, mtz, work_dir, max_to_keep=20):
        """Initialise a new amore rotation search class

        Parameters
        ----------
        amore_exe : str
            The path to the amore executable
        mtz : str
            The path to the input MTZ
        work_dir : str
            The path to the working directory
        max_to_keep : int
            The maximum number of results to keep [default: 20]
        """
        self._amore_exe = None
        self._max_to_keep = 0
        self._mtz = None
        self._search_results = None
        self._work_dir = None

        self.amore_exe = amore_exe
        self.max_to_keep = max_to_keep
        self.mtz = mtz
        self.work_dir = work_dir

    @property
    def amore_exe(self):
        """The amore exectutable"""
        return self._amore_exe

    @amore_exe.setter
    def amore_exe(self, amore_exe):
        """Define the amore executable"""
        self._amore_exe = amore_exe

    @property
    def max_to_keep(self):
        """The maximum number of results to keep"""
        return self._max_to_keep

    @max_to_keep.setter
    def max_to_keep(self, max_to_keep):
        """Define the maximum number of results to keep"""
        self._max_to_keep = max_to_keep

    @property
    def mtz(self):
        """The input MTZ file"""
        return self._mtz

    @mtz.setter
    def mtz(self, mtz):
        """Define the input MTZ file"""
        self._mtz = mtz

    @property
    def search_results(self):
        """The results from the amore rotation search"""
        return sorted(self._search_results, key=lambda x: float(x.CC_F_Z_score), reverse=True)[:self._max_to_keep]

    @property
    def work_dir(self):
        """The path to the working directory"""
        return self._work_dir

    @work_dir.setter
    def work_dir(self, work_dir):
        """Define the working directory"""
        self._work_dir = work_dir

    @staticmethod
    def cleanup(logfile):
        """Simple function to clean up log files after a run"""
        os.remove(logfile)

    @staticmethod
    def calculate_integration_box(model):
        """Function to calculate the integration radius or minimal box for an input PDB

        Parameters
        ----------
        model : str
            Path to input model

        Returns
        -------
        float
            The X coordinate
        float
            The Y coordinate
        float
            The Z coordinate
        float
            The integration radius for spherical structure

        """
        pdb_input = iotbx.pdb.pdb_input(file_name=model)
        hierarchy = pdb_input.construct_hierarchy()

        # Get resolution
        resolution = iotbx.pdb.mining.extract_best_resolution(
            pdb_input.extract_remark_iii_records(2)
        )

        # Set a default resolution if mining fails
        if resolution is None:
            resolution = 2.0

        # Get a list of all xyz coordinates
        chain = hierarchy.models()[0].chains()[0]
        xyz = numpy.zeros((chain.atoms_size(), 3))
        for i, atom in enumerate(chain.atoms()):
            xyz[i] = atom.xyz

        # Get the smallest box containing the model
        #   numpy.ptp() ==> "Range of values (maximum - minimum) along an axis"
        diffs = numpy.asarray([
            numpy.ptp(xyz[:, 0]),
            numpy.ptp(xyz[:, 1]),
            numpy.ptp(xyz[:, 2])
        ])
        # Get integration radius (note, for spherical structure)
        intrad = diffs.min() * 0.75

        # Add together for each coordinate
        x, y, z = diffs + intrad + resolution

        return x.item(), y.item(), z.item(), intrad.item()

    @staticmethod
    def rotfun(amore_exe, table1, hklpck1, clmn1,  shres, intrad, hklpck0=None,
               clmn0=None, mapout=None, pklim=None, npic=None, rotastep=None):
        """Function to perform amore rotation function,

        Parameters
        ----------
        amore_exe : str
           The path to the AMORE executable
        table1 : str
           The path to the TABLE1 file
        hklpck1 : str
           The path to the HKLPCK1 file
        hklpck0 : str
           The path to the HKLPCK0 file
        clmn1 : str
           The path to the CLMN1 file
        clmn0 : str
           The path to the CLMN0 file
        mapout : str
           The path to the MAPOUT file
        shres : int, float
           Spherical harmonic resolution
        intrad : int, float
           The tolerance
        pklim : int, float
           Peak limit, output all peaks above :obj:`float`
        npic : int, float
           Number of peaks to output from the translation function map for each orientation
        rotastep : int, float
           Size of rotation step

        Returns
        -------
        cmd : list
            rotation function command
        stdin : str
            rotation function standard input

        """
        if table1 and hklpck1 and clmn1 and not hklpck0 and not clmn0 and not mapout:
            cmd = [amore_exe, 'table1', table1,
                   'HKLPCK1', hklpck1, 'clmn1', clmn1]
            stdin = os.linesep.join([
                "ROTFUN",
                "VERB",
                "TITLE : Generate HKLPCK1 from MODEL FRAGMENT   1",
                "GENE 1   RESO 100.0 {0}  CELL_MODEL 80 75 65",
                "CLMN MODEL 1     RESO  20.0  {0} SPHERE   {1}"
            ])
            stdin = stdin.format(shres, intrad)

        elif table1 and hklpck1 and clmn1 and hklpck0 and clmn0 and mapout:
            if os.path.isfile(clmn1) and os.path.isfile(hklpck1):
                cmd = [amore_exe, 'table1', table1, 'HKLPCK1', hklpck1, 'hklpck0', hklpck0,
                       'clmn1', clmn1, 'clmn0', clmn0, 'MAPOUT', mapout]
                stdin = os.linesep.join([
                    "ROTFUN",
                    "TITLE : Generate HKLPCK1 from MODEL FRAGMENT   1",
                    "CLMN CRYSTAL ORTH  1 RESO  20.0  {0}  SPHERE   {1}",
                    "ROTA  CROSS  MODEL 1  PKLIM {2}  NPIC {3} STEP {4}",
                ])
                stdin = stdin.format(shres, intrad, pklim, npic, rotastep)
            else:
                cmd = [amore_exe, 'table1', table1, 'HKLPCK1', hklpck1, 'hklpck0', hklpck0,
                       'clmn1', clmn1, 'clmn0', clmn0, 'MAPOUT', mapout]
                stdin = os.linesep.join([
                    "ROTFUN",
                    "TITLE: Generate HKLPCK1 from MODEL FRAGMENT 1",
                    "GENE 1   RESO 100.0 {0}  CELL_MODEL 80 75 65",
                    "CLMN CRYSTAL ORTH  1 RESO  20.0  {0}  SPHERE   {1}",
                    "CLMN MODEL 1     RESO  20.0  {0} SPHERE   {1}",
                    "ROTA  CROSS  MODEL 1  PKLIM {2}  NPIC {3} STEP {4}"
                ])
                stdin = stdin.format(shres, intrad, pklim, npic, rotastep)
        else:
            msg = "Incorrect combination of input arguments"
            logger.critical(msg)
            raise RuntimeError(msg)

        return cmd, stdin

    @staticmethod
    def solvent_content(pdbin, cell_parameters, space_group):
        """Get the solvent content for an input pdb
        
        Parameters
        ----------
        pdbin : str
            Path to input PDB file
        cell_parameters : str
            The parameters describing the unit cell of the crystal
        space_group : str
            The space group of the crystal

        Returns
        -------
        float
            The solvent content
        """
        crystal_symmetry = cctbx.crystal.symmetry(
            unit_cell=cell_parameters, space_group_symbol=space_group)
        dens_calc = mmtbx.scaling.matthews.density_calculator(crystal_symmetry)
        return dens_calc.solvent_fraction(molecular_weight(pdbin), 0.74) * 100

    @staticmethod
    def submit_chunk(chunk_scripts, output_dir, nproc, job_name, submit_qtype, submit_queue, monitor):
        """Submit jobs in small chunks to avoid using too much disk space
        
        Parameters
        ----------
        chunk_scripts : list
            List of scripts for each chunk
        nproc : int, optional
            The number of processors to run the job on
        job_name : str
            The name of the job to submit
        submit_qtype : str
            The cluster submission queue type - currently support SGE and LSF
        submit_queue : str
            The queue to submit to on the cluster

        """
        j = Job(submit_qtype)
        j.submit(chunk_scripts, directory=output_dir, name=job_name,
                 nproc=nproc, queue=submit_queue, permit_nonzero=True)
        interval = int(numpy.log(len(chunk_scripts)) / 3)
        interval_in_seconds = interval if interval >= 5 else 5
        j.wait(interval=interval_in_seconds, monitor=monitor)

    @staticmethod
    def tabfun(amore_exe, xyzin1, xyzout1, table1, x=200, y=200, z=200, a=90, b=90, c=120):
        """Function to perform amore table function

        Parameters
        ----------
        amore_exe : str
           The path to the AMORE executable
        xyzin1 : str
           Path to input model
        xyzout1 : str
           Path to the output model
        table1 : str
           Path to the output table file
        x : int, float, optional
           x value for minimal box [default: 200]
        y : int, float, optional
           y value for minimal box [default: 200]
        z : int, float, optional
           z value for minimal box [default: 200]
        a : int, float, optional
           alpha value for minimal box [default: 90]
        b : int, float, optional
           beta value for minimal box [default: 90]
        c : int, float, optional
           gamma value for minimal box [default: 120]

        Returns
        -------
        file
            Output PDB for use in rotfun
        file
            Output table file for use in rotfun

        """
        cmd = [amore_exe, 'xyzin1', xyzin1,
               'xyzout1', xyzout1, 'table1', table1]
        stdin = os.linesep.join([
            "TITLE: Produce table for MODEL FRAGMENT",
            "TABFUN",
            "CRYSTAL {0} {1} {2} {3} {4} {5} ORTH 1",
            "MODEL 1 BTARGET 23.5",
            "SAMPLE 1 RESO 2.5 SHANN 2.5 SCALE 4.0",
        ])
        stdin = stdin.format(x, y, z, a, b, c)
        return cmd, stdin

    def run_pdb(self, models_dir, output_model_dir, nproc=2, shres=3.0, pklim=0.5, npic=50, rotastep=1.0,
                min_solvent_content=20, submit_qtype=None, submit_queue=None, monitor=None, chunk_size=5000):
        """Run amore rotation function on a directory of models

        Parameters
        ----------
        models_dir : str
            The directory containing the models to run the rotation search on
        output_model_dir : str
            Path to the directory to move top ranking models from the rotation search
        nproc : int, optional
            The number of processors to run the job on
        shres : int, float, optional
            Spherical harmonic resolution [default 3.0]
        pklim : int, float, optional
            Peak limit, output all peaks above <float> [default: 0.5]
        npic : int, optional
            Number of peaks to output from the translation function map for each orientation [default: 50]
        rotastep : int, float, optional
            Size of rotation step [default : 1.0]
        min_solvent_content : int, float, optional
            The minimum solvent content present in the unit cell with the input model [default: 30]
        submit_qtype : str
            The cluster submission queue type - currently support SGE and LSF
        submit_queue : str
            The queue to submit to on the cluster
        monitor
        chunk_size : int, optional
            The number of jobs to submit at the same time

        Returns
        -------
        file
            log file for each model in the models_dir

        """
        simbad_dat_files = [
            os.path.join(root, filename) for root, _, files in os.walk(models_dir)
            for filename in files if filename.endswith('.dat')
        ]

        # Creating temporary output directory
        output_dir = os.path.join(self.work_dir, 'output')
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        # Get the space group and cell parameters for the input mtz
        space_group, _, cell_parameters = mtz_util.crystal_data(self.mtz)
        cell_parameters = " ".join(map(str, cell_parameters))

        # Save some data for populating results later on
        rotation_data = []
        total_chunk_cycles = len(
            simbad_dat_files) // chunk_size + (len(simbad_dat_files) % 5 > 0)
        for cycle, i in enumerate(range(0, len(simbad_dat_files), chunk_size)):
            logger.info("Working on chunk %d out of %d",
                        cycle + 1, total_chunk_cycles)
            chunk_dat_files = simbad_dat_files[i:i + chunk_size]

            # Generate the relevant scripts and data
            amore_files, to_delete = [], []
            for dat_model in chunk_dat_files:
                root = dat_model.replace('.dat', '')
                name = os.path.basename(root)

                # Convert .dat to .pdb
                input_model = tmp_file(
                    directory=output_dir, prefix="", stem=name, suffix='.pdb')
                with open(dat_model, 'rb') as f_in, open(input_model, 'w') as f_out:
                    f_out.write(zlib.decompress(base64.b64decode(f_in.read())))

                # Compute the solvent content and decide if we trial this structure
                try:
                    solvent_content = self.solvent_content(
                        input_model, cell_parameters, space_group)
                except:
                    logger.critical(
                        "Error calculating solvent content for %s", name)
                    continue
                if solvent_content < min_solvent_content:
                    msg = "Skipping %s: solvent content is predicted to be less than %.2f"
                    logger.debug(msg, name, min_solvent_content)
                    continue

                logger.debug(
                    "Generating script to perform AMORE rotation function on %s", name)

                # Set up variables for the __TAB__ run
                x, y, z, intrad = AmoreRotationSearch.calculate_integration_box(
                    input_model)
                output_model = os.path.join(
                    self.work_dir, 'output', '{0}.pdb'.format(name))
                table1 = os.path.join(
                    self.work_dir, 'output', '{0}_sfs.tab'.format(name))

                # Get the command and stdin
                tab_cmd, tab_key = AmoreRotationSearch.tabfun(
                    self.amore_exe, input_model, output_model, table1, x, y, z)

                # Set up variables for the __ROT__ run
                hklpck1 = os.path.join(
                    self.work_dir, 'output', '{0}.hkl'.format(name))
                hklpck0 = os.path.join(self.work_dir, 'spmipch.hkl')
                clmn1 = os.path.join(
                    self.work_dir, 'output', '{0}.clmn'.format(name))
                clmn0 = os.path.join(
                    self.work_dir, 'output', '{0}_spmipch.clmn'.format(name))
                mapout = os.path.join(
                    self.work_dir, 'output', '{0}_amore_cross.map'.format(name))

                # Get the command and stdin
                rot_cmd, rot_key = AmoreRotationSearch.rotfun(self.amore_exe, table1, hklpck1, clmn1, shres, intrad,
                                                              hklpck0, clmn0, mapout, pklim, npic, rotastep)

                # Set up script, log and stdin for the amore table function
                prefix, stem = "tabfun_", name
                tab_stdin = tmp_file(directory=output_dir,
                                     prefix=prefix, stem=stem, suffix=".stdin")
                with open(tab_stdin, 'w') as f_out:
                    f_out.write(tab_key)

                # Set up script, log and stdin for the amore rotation function
                prefix, stem = "rotfun_", name
                rot_stdin = tmp_file(directory=output_dir,
                                     prefix=prefix, stem=stem, suffix=".stdin")
                with open(rot_stdin, 'w') as f_out:
                    f_out.write(rot_key)

                # Generate script
                amore_script = make_script(
                    [[EXPORT, "CCP4_SCR=" + output_dir], tab_cmd + ["<", tab_stdin],
                     os.linesep, rot_cmd + ["<", rot_stdin], os.linesep,
                     ["rm", clmn0, clmn1, hklpck1, table1, mapout]],
                    directory=output_dir, prefix=prefix, stem=stem
                )
                amore_log = amore_script.rsplit(".", 1)[0] + '.log'

                # Save a copy of the files we need to run
                amore_files += [(amore_script, tab_stdin,
                                 rot_stdin, amore_log)]

                # Save the data
                rotation_data += [(input_model, amore_log)]

            results = []
            if len(amore_files) > 0:
                # Run the AMORE tab/rot function on chunk
                logger.info("Running AMORE tab/rot functions")
                amore_scripts, _, _, _ = zip(*amore_files)
                self.submit_chunk(amore_scripts, output_dir, nproc,
                                  'simbad_amore', submit_qtype, submit_queue, monitor)

                # Remove some files to clear disk space
                amore_tmp_files = glob.glob(os.path.join(
                    output_dir, "{}_".format(os.path.basename(self.amore_exe))
                ))
                map(os.remove, amore_tmp_files)

                # Populate the results
                for input_model, rot_log in rotation_data:
                    pdb_code = os.path.basename(rot_log).replace(
                        "rotfun_", "").replace(".log", "")
                    RP = rotsearch_parser.RotsearchParser(rot_log)
                    score = _AmoreRotationScore(pdb_code, RP.alpha, RP.beta, RP.gamma, RP.cc_f, RP.rf_f, RP.cc_i,
                                                RP.cc_p, RP.icp, RP.cc_f_z_score, RP.cc_p_z_score, RP.num_of_rot)
                    if RP.cc_f_z_score is not None:
                        results += [score]
                        # Need to move input models to specific directory
                        if os.path.isfile(input_model):
                            shutil.move(input_model, output_model_dir)

            else:
                msg = "No structures to be trialled"
                logger.critical(msg)

        # Save the results
        self._search_results = results

        # Remove the large temporary tmp directory
        shutil.rmtree(output_dir)

        return

    def sortfun(self):
        """A function to prepare files for amore rotation function

        Parameters
        ----------
        self.mtz : str
            mtz file input to :obj:`AmoreRotationSearch`
        self.work_dir : str
            working directory input to :obj:`AmoreRotationSearch`

        Returns
        -------
        file
            spmipch.hkl file needed for rotfun and tabfun

        """
        logger.info("Preparing files for AMORE rotation function")
        f, sigf, _, _, _, _, _ = mtz_util.get_labels(self.mtz)
        cmd = [self.amore_exe, 'hklin', self.mtz, 'hklpck0',
               os.path.join(self.work_dir, 'spmipch.hkl')
               ]
        stdin = os.linesep.join([
            "TITLE   ** spmi  packing h k l F for crystal**",
            "SORTFUN RESOL 100.  2.5",
            "LABI FP={0}  SIGFP={1}",
        ])
        stdin = stdin.format(f, sigf)
        cexec(cmd, stdin=stdin)

    def summarize(self, csv_file):
        """Summarize the search results

        Parameters
        ----------
        csv_file : str
           The path for a backup CSV file

        Raises
        ------
            No results found
        """
        from simbad.util import summarize_result
        columns = [
            "ALPHA", "BETA", "GAMMA", "CC_F", "RF_F", "CC_I", "CC_P", "Icp",
            "CC_F_Z_score", "CC_P_Z_score", "Number_of_rotation_searches_producing_peak"
        ]
        summarize_result(self.search_results,
                         csv_file=csv_file, columns=columns)
