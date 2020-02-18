cmake_minimum_required(VERSION 2.8)
project(simbad NONE)

install(
    DIRECTORY
        simbad
    DESTINATION
        lib/py2
)

install(
    DIRECTORY
        static
    DESTINATION
        share/simbad
)

if (WIN32)
  install(
      PROGRAMS
          bin/simbad.bat bin/simbad-contaminant.bat bin/simbad-database.bat bin/simbad-full.bat bin/simbad-lattice.bat bin/simbad-morda.bat
      DESTINATION
          bin
  )
else()
  install(
      PROGRAMS
          bin/simbad bin/simbad-contaminant bin/simbad-database bin/simbad-full bin/simbad-lattice bin/simbad-morda
      DESTINATION
          bin
)
endif()

install(DIRECTORY ccp4i DESTINATION share)
