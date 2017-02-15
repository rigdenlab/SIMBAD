# SIMBAD
Sequence Independent Molecular Replacement Based on Available Database

### Installation

    git clone https://github.com/rigdenlab/SIMBAD.git
    cd SIMBAD
    ccp4-python setup.py install --prefix $CCP4 --install-scripts $CCP4/bin

This will install SIMBAD into your CCP4 installation. On top of installing the source code, an executable script ``simbad`` should be automatically installed.

### Usage

You can run SIMBAD by executing the following command:

    simbad -mtz <MTZ FILE>

### Development

If you would like to install SIMBAD in development mode, execute the following:

    git clone https://github.com/rigdenlab/SIMBAD.git
    cd SIMBAD
    ccp4-python setup.py develop --prefix $CCP4 --install-scripts $CCP4/bin

This will link this version into the root of your CCP4 installation. Therefore, you do not need to re-install it after every change.

### Contributors

- Adam Simpkin
- Jens Thomas
- Ronan Keegan
- Felix Simkovic
- Daniel Rigden
