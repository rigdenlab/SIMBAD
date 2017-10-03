cmake_minimum_required(VERSION 2.8)
project(simbad NONE)

install (
    DIRECTORY
        simbad
    DESTINATION
        lib/py2
)

install (
    DIRECTORY
        static
    DESTINATION    
        share/simbad    
)    

install (
    PROGRAMS
        bin/simbad bin/simbad-contaminant bin/simbad-database bin/simbad-full bin/simbad-lattice bin/simbad-morda
    DESTINATION     
        bin   
)    
