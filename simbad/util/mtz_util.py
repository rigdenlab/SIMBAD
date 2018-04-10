"""Module for MTZ file I/O and manipulation"""

__author__ = "Adam Simpkin & Jens Thomas"
__date__ = "17 May 2017"
__version__ = "0.2"

from cctbx import crystal
from cctbx import miller
from cctbx import sgtbx
from cctbx.xray import observation_types
from iotbx import reflection_file_reader
from iotbx.reflection_file_utils import looks_like_r_free_flags_info

import logging

logger = logging.getLogger(__name__)


class ExperimentalData(object):
    """Class to create a temporary mtz containing all the columns needed for SIMBAD from input reflection file

    Attributes
    ----------
    input_reflection_file : str
        Path to the input reflection file in ccp4 mtz format
    output_mtz_file : str
        Path to the output mtz file

    Example 1
    ---------
    >>> from simbad.util import mtz_util
    >>> ED = mtz_util.ExperimentalData("<input_reflection_file>")
    >>> ED.process_miller_arrays()
    >>> ED.output_mtz("<output_mtz_file>")

    Example 2
    ---------
    >>> from simbad.util import mtz_util
    >>> ED = mtz_util.ExperimentalData("<input_reflection_file>")
    >>> ED.change_space_group('<new space group>')
    >>> ED.output_mtz("<output_mtz_file>")
    """

    def __init__(self, input_reflection_file):

        self.amplitude_array = None
        self.anomalous_amplitude_array = None
        self.reconstructed_amplitude_array = None
        self.intensity_array = None
        self.anomalous_intensity_array = None
        self.free_array = None
        self.mtz_dataset = None

        reflection_file = reflection_file_reader.any_reflection_file(file_name=input_reflection_file)
        if not reflection_file.file_type() == "ccp4_mtz":
            msg = "File is not of type ccp4_mtz: {0}".format(input_reflection_file)
            logging.critical(msg)
            raise RuntimeError(msg)

        self.all_miller_arrays = reflection_file.as_miller_arrays()
        self.get_array_types()

    def add_array_to_mtz_dataset(self, miller_array, column_root_label=None):
        """Function to add cctbx miller array obj to cctbx mtz dataset obj

        Parameters
        ----------
        miller_array : cctbx :obj:
            Input cctbx obj containing a miller array
        column_root_label : str
            The root for the label of the output column

        Returns
        -------
        self.mtz_dataset : cctbx :obj:
            cctbx mtz obj containing the input miller array
        """
        if self.mtz_dataset:
            if column_root_label:
                self.mtz_dataset.add_miller_array(miller_array, column_root_label=column_root_label)
            else:
                column_root_label = miller_array.info().labels[0]
                self.mtz_dataset.add_miller_array(miller_array, column_root_label=column_root_label)
        else:
            if column_root_label:
                self.mtz_dataset = miller_array.as_mtz_dataset(column_root_label=column_root_label)
            else:
                column_root_label = miller_array.info().labels[0]
                self.mtz_dataset = miller_array.as_mtz_dataset(column_root_label=column_root_label)
        return

    @staticmethod
    def check_anomalous_arrays(miller_array):
        """Function to check intensity/amplitude arrays from cctbx to ensure that they are not anomalous arrays

        Parameters
        ----------
        miller_array : cctbx :obj:
            A cctbx :obj: containing a miller array of either intensities or amplitudes

        Returns
        -------
        bool
            True/False
        """

        if miller_array.anomalous_flag():
            return True
        elif miller_array.info().type_hints_from_file == "anomalous_difference":
            return True
        elif len(miller_array.info().labels) > 2 and any("+" in label for label in miller_array.info().labels):
            return True
        return False

    def change_space_group(self, new_space_group):
        """Change space group of input mtz

        Parameters
        ----------
        new_space_group : str
            The new space group
        """

        for miller_array in self.all_miller_arrays:
            array_info = miller_array.info()
            if not looks_like_r_free_flags_info(miller_array.info()):
                new_space_group_info = sgtbx.space_group_info(symbol=new_space_group)
                new_crystal_symmetry = crystal.symmetry(unit_cell=miller_array.unit_cell(),
                                                        space_group_info=new_space_group_info,
                                                        assert_is_compatible_unit_cell=False)
                miller_array = miller_array.customized_copy(crystal_symmetry=new_crystal_symmetry,
                                                            info=array_info)

            self.add_array_to_mtz_dataset(miller_array)

    def create_amplitude_array(self, intensity_array):
        """Function to create a cctbx amplitude array from an cctbx intensity array

        Parameters
        ----------
        intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of intensities

        Returns
        -------
        self.amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of amplitudes
        """
        array_info = miller.array_info()
        array_info.labels = ['F', 'SIGF']
        self.amplitude_array = intensity_array.customized_copy(
            observation_type=observation_types.amplitude(),
            info=array_info
        )
        return

    def create_anomalous_amplitude_array(self, anomalous_intensity_array):
        """Function to create a cctbx anomalous amplitude array from a cctbx anomalous intensity array

        Parameters
        ----------
        anomalous_intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous intensities

        Returns
        -------
        self.anomalous_amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous amplitudes
        """
        array_info = miller.array_info()
        array_info.labels = ['F(+)', 'F(-)', 'SIGF(+)', 'SIGF(-)']
        self.anomalous_amplitude_array = anomalous_intensity_array.customized_copy(
            observation_type=observation_types.amplitude(),
            info=array_info
        )
        return

    def create_anomalous_intensity_array(self, anomalous_amplitude_array):
        """Function to create a cctbx anomalous intensity array from a cctbx anomalous amplitude array

        Parameters
        ----------
        anomalous_amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous amplitudes

        Returns
        -------
        self.anomalous_intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous intensities
        """
        array_info = miller.array_info()
        array_info.labels = ['I(+)', 'I(-)', 'SIGI(+)', 'SIGI(-)']
        self.anomalous_intensity_array = anomalous_amplitude_array.customized_copy(
            observation_type=observation_types.intensity(),
            info=array_info
        )
        return

    def create_merged_intensity_array(self, anomalous_intensity_array):
        """Function to create a cctbs intensity array from a cctbx anomalous intensity array

        Parameters
        ----------
        anomalous_intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous intensities

        Returns
        -------
        self.intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of intensities
        """
        array_info = miller.array_info()
        array_info.labels = ['I', 'SIGI']
        merged_intensity_array = anomalous_intensity_array.copy().as_non_anomalous_array().merge_equivalents()
        self.intensity_array = merged_intensity_array.array().customized_copy(
            observation_type=observation_types.intensity(),
            info=array_info
        )
        return

    def create_reconstructed_amplitude_array(self, anomalous_amplitude_array):
        """Function to create a cctbx reconstructed amplitude array from a cctbx anomalous amplitude array

        Parameters
        ----------
        anomalous_amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous amplitudes

        Returns
        -------
        self.reconstructed_amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of reconstructed amplitudes
        """
        array_info = miller.array_info()
        array_info.labels = ['F', 'SIGF', 'DANO', 'SIGDANO', 'ISYM']
        self.reconstructed_amplitude_array = anomalous_amplitude_array.customized_copy(
            observation_type=observation_types.reconstructed_amplitude(),
            info=array_info
        )
        return

    def get_array_types(self):
        """Function to assign array types contained within cctbx obj, in cases where there are multiple instances of
        a type of array, only the first will be considered. Due to limitations in cctbx

        Returns
        -------
        self.free_array : cctbx :obj:
            A cctbx :obj: containing a miller array containing the Free R data
        self.amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of amplitudes
        self.anomalous_amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous amplitudes
        self.reconstructed_amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of reconstructed amplitudes
        self.intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of intensities
        self.anomalous_intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous intensities
        """
        for miller_array in self.all_miller_arrays:
            if miller_array.observation_type() is None:
                if looks_like_r_free_flags_info(miller_array.info()):
                    if not self.free_array:
                        self.free_array = miller_array

            if miller_array.is_xray_amplitude_array() and not self.check_anomalous_arrays(miller_array):
                if not self.amplitude_array:
                    self.amplitude_array = miller_array
            elif miller_array.info().type_hints_from_file == 'amplitude' and \
                    not self.check_anomalous_arrays(miller_array):
                if not self.amplitude_array:
                    self.amplitude_array = miller_array
            elif miller_array.is_xray_reconstructed_amplitude_array():
                if not self.reconstructed_amplitude_array:
                    self.reconstructed_amplitude_array = miller_array
            elif miller_array.is_xray_amplitude_array() and self.check_anomalous_arrays(miller_array):
                if not self.anomalous_amplitude_array:
                    self.anomalous_amplitude_array = miller_array
            elif miller_array.is_xray_intensity_array() and not self.check_anomalous_arrays(miller_array):
                if not self.intensity_array:
                    self.intensity_array = miller_array
            elif miller_array.info().type_hints_from_file == 'intensity' and \
                    not self.check_anomalous_arrays(miller_array):
                if not self.intensity_array:
                    self.intensity_array = miller_array
            elif miller_array.is_xray_intensity_array() and self.check_anomalous_arrays(miller_array):
                if not self.anomalous_intensity_array:
                    self.anomalous_intensity_array = miller_array
        return

    def output_mtz(self, output_mtz_file):
        """Function to output an mtz file from processed miller arrays

        Parameters
        ----------
        output_mtz_file : str
            Path to output mtz file

        Returns
        -------
        file
            mtz file containing all the columns needed to run SIMBAD
        """
        mtz_object = self.mtz_dataset.mtz_object()
        mtz_object.write(file_name=output_mtz_file)
        return

    def process_miller_arrays(self):
        """Function to process the miller arrays needed for SIMBAD

        Parameters
        ----------
        self.free_array : cctbx :obj:
            A cctbx :obj: containing a miller array containing the Free R data
        self.amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of amplitudes
        self.anomalous_amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous amplitudes
        self.reconstructed_amplitude_array : cctbx :obj:
            A cctbx :obj: containing a miller array of reconstructed amplitudes
        self.intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of intensities
        self.anomalous_intensity_array : cctbx :obj:
            A cctbx :obj: containing a miller array of anomalous intensities

        Returns
        -------
        self.mtz_dataset : cctbx :obj:
            cctbx mtz obj containing all the miller arrays needed to run SIMBAD
        """

        # Add amplitudes
        if self.reconstructed_amplitude_array:
            self.add_array_to_mtz_dataset(self.reconstructed_amplitude_array)
        elif self.anomalous_amplitude_array:
            self.create_reconstructed_amplitude_array(self.anomalous_amplitude_array)
            self.add_array_to_mtz_dataset(self.reconstructed_amplitude_array)
        elif self.amplitude_array:
            self.add_array_to_mtz_dataset(self.amplitude_array)
        elif self.intensity_array:
            self.create_amplitude_array(self.intensity_array)
            self.add_array_to_mtz_dataset(self.amplitude_array)
        elif self.anomalous_intensity_array:
            self.create_anomalous_amplitude_array(self.anomalous_intensity_array)
            self.create_reconstructed_amplitude_array(self.anomalous_amplitude_array)
            self.add_array_to_mtz_dataset(self.reconstructed_amplitude_array)

        # Add intensities
        if self.intensity_array:
            self.add_array_to_mtz_dataset(self.intensity_array)
        elif self.anomalous_intensity_array:
            self.create_merged_intensity_array(self.anomalous_intensity_array)
            self.add_array_to_mtz_dataset(self.intensity_array)

        # if not self.intensity_array or not self.amplitude_array:
        #     msg = "No amplitudes or intensities found in input reflection file"
        #     logging.critical(msg)
        #     raise RuntimeError(msg)

        # Add free flag
        if self.free_array:
            try:
                self.add_array_to_mtz_dataset(self.free_array)
            except RuntimeError:
                pass
        else:
            self.free_array = self.amplitude_array.generate_r_free_flags(format='ccp4')
            self.add_array_to_mtz_dataset(self.free_array, "FreeR_flag")
        return


def crystal_data(mtz_file):
    """Set crystallographic parameters from mtz file

    Parameters
    ----------
    mtz_file : str
       The path to the mtz file

    Returns
    -------
    space_group : str
       The space group
    resolution : str
       The resolution
    cell_parameters : tuple
       The cell parameters

    """

    reflection_file = reflection_file_reader.any_reflection_file(file_name=mtz_file)
    content = reflection_file.file_content()
    space_group = content.space_group_name().replace(" ", "")
    resolution = content.max_min_resolution()[1]
    cell_parameters = content.crystals()[0].unit_cell_parameters()

    return space_group, resolution, cell_parameters


def get_labels(mtz_file):
    """Function to get the column labels for input mtz file

    Parameters
    ----------
    mtz_file : str
       The path to the mtz file

    Returns
    -------
    f : str
        f column label
    fp : str
        fp column label
    i : str
        i column label
    sigi : str
        sigi column label
    dano : str
        dano column label
    sigdano : str
        sigdano column label
    free : str
        free column label
    """

    reflection_file = reflection_file_reader.any_reflection_file(file_name=mtz_file)
    if not reflection_file.file_type() == "ccp4_mtz":
        msg="File is not of type ccp4_mtz: {0}".format(mtz_file)
        logging.critical(msg)
        raise RuntimeError(msg)

    content = reflection_file.file_content()
    ctypes = content.column_types()
    clabels = content.column_labels()
    ftype = 'F'
    jtype = 'J'
    dtype = 'D'

    if ftype not in ctypes:
        msg = "Cannot find any structure amplitudes in: {0}".format(mtz_file)
        raise RuntimeError(msg)
    f = clabels[ctypes.index(ftype)]

    # FP derived from F
    fp = 'SIG' + f
    fp_alt = 'PHI' + f
    if fp in clabels:
        pass
    elif fp_alt in clabels:
        fp = fp_alt
        pass
    else:
        msg = "Cannot find label {0} or {1} in file: {2}".format(fp, fp_alt, mtz_file)
        logging.warning(msg)

    i, sigi = None, None
    if jtype in ctypes:
        i = clabels[ctypes.index(jtype)]

        # SIGI derired from I
        sigi = 'SIG' + i
        if sigi not in clabels:
            msg = "Cannot find label {0} in file: {1}".format(sigi, mtz_file)
            raise RuntimeError(msg)

    try:
        if dtype not in ctypes:
            msg = "Cannot find any structure amplitudes in: {0}".format(mtz_file)
            raise RuntimeError(msg)
        dano = clabels[ctypes.index(dtype)]

        # SIGDANO derived from DANO
        sigdano = 'SIG' + dano
        if sigdano not in clabels:
            msg = "Cannot find label {0} in file: {1}".format(sigdano, mtz_file)
            raise RuntimeError(msg)
    except RuntimeError:
        dano, sigdano = None, None

    free = None
    for label in clabels:
        for word in ["free", "test", "cross", "status", "flag"]:
            if label.lower().find(word) >= 0:
                if free:
                    logger.warning("FOUND >1 R FREE label in file!")
                free = label
                break

    return f, fp, i, sigi, dano, sigdano, free

