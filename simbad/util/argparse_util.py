"""
Code to add arguments to SIMBAD

@author: hlasimpk

Largely borrowed from AMPLE

"""

from simbad.util import version

def add_core_options(parser):
    """"Function to add any argument required by all runtypes"""

    parser.add_argument('-config_file', help="user configuration file")

    parser.add_argument('-mtz', metavar='MTZ in', type=str,
                        help='Path to the MTZ file with the reflection data.')

    parser.add_argument('-sf_cif', metavar='CIF in', type=str,
                        help='Path to the CIF file with the reflection data.')

    parser.add_argument('-nproc', type=int,
                        help="Number of processors [1]. For local, serial runs the jobs will be split across nproc processors. " + \
                             "For cluster submission, this should be the number of processors on a node.")

    parser.add_argument('-work_dir', type=str,
                        help='Path to the directory where SIMBAD will run (will be created if it doesn\'t exist)')
    return

def add_cluster_submit_options(parser):
    """Add the options for submission to a cluster queuing system"""

    submit_group = parser.add_argument_group('Cluster queue submission options')

    submit_group.add_argument('-submit_array', metavar='True/False', type=str,
                              help='Submit SGE jobs as array jobs')

    submit_group.add_argument('-submit_cluster', metavar='True/False', type=str,
                              help='Submit jobs to a cluster - need to set -submit_qtype flag to specify the batch queue system.')

    submit_group.add_argument('-submit_max_array', type=int,
                              help='The maximum number of jobs to run concurrently with SGE array job submission')

    submit_group.add_argument('-submit_num_array_jobs', type=int,
                              help='The number of jobs to run concurrently with SGE array job submission')

    submit_group.add_argument('-submit_pe_lsf',
                              help='Cluster submission: string to set number of processors for LSF queueing system')

    submit_group.add_argument('-submit_pe_sge',
                              help='Cluster submission: string to set number of processors for SGE queueing system')

    submit_group.add_argument('-submit_queue', type=str,
                              help='The queue to submit to on the cluster.')

    submit_group.add_argument('-submit_qtype', type=str,
                              help='cluster submission queue type - currently support SGE and LSF')
    return

def add_general_options(parser):

    # Add core options

    add_core_options(parser)


    # Exectutable options

    parser.add_argument('-amore', metavar='amore exe', type=str,
                        help='Path to amore executable')

    # Database options

    parser.add_argument('-pdb_db', metavar='PDB database', type=str,
                        help='Path to local installation of the PDB')

    parser.add_argument('-cont_db', metavar='Contaminant database', type=str,
                        help='Path to local installation of the contaminant database')

    parser.add_argument('-morda_db', metavar='MoRDa database', type=str,
                        help='Path to local installation of the MoRDa database')

    # MTZ options

    parser.add_argument('-F', type=str,
                        help='Flag for F column in the MTZ')

    parser.add_argument('-SIGF', type=str,
                        help='Flag for SIGF column in the MTZ')

    parser.add_argument('-FREE', type=str,
                        help='Flag for FREE column in the MTZ')

    parser.add_argument('-DANO', type=str,
                        help='Flag for the DANO column in the MTZ')

    parser.add_argument('-SIGDANO', type=str,
                        help='Flag for the SIGDANO column in the MTZ')

    # Other options

    parser.add_argument('-name', metavar='job_name',
                        help='4-letter identifier for job [simb]')

    parser.add_argument('--version', action='version', version='%(prog)s {0}'.format(version.__version__))

    parser.add_argument('-webserver_uri',
                        help='URI of the webserver directory - also indicates we are running as a webserver')

    return

def add_anomalous_options(parser):

    anomalous_group = parser.add_argument_group('Anomalous signal options')

    anomalous_group.add_argument('-hatom_num', type=int,
                                 help="Input the number of heavy atoms to search for")

    anomalous_group.add_argument('-hatom_type', type=str,
                                 help="Input the type of heavy atoms to search for")

    return