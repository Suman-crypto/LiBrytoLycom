#!/bin/bash

# SETTINGS
#-------------------------------

para=1

#-------------------------------

a0=Ntest06_paleo

if [ -e "libry.x" ]
then
  rm libry.x
fi

if [ -e "libry*.mod" ]
then
  rm libry*.mod
fi

if [[ $para == 1 ]]
then

  # DKRZ
  #source /sw/rhel6-x64/etc/profile.mistral 	 	# make module command available
  #module load intel/18.0.2 intelmpi/2018.1.163
  #module add intel/18.0.2 intelmpi/2018.1.163
  
  #mpiifort -fp-model strict -free -real-size 64 -O0 \
  #         $(/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-intel14/bin/nf-config --fflags --flibs) \
  #         -g -traceback -ftrapuv -debug all -gen-interfaces \
  #         libry_par.f90 libry_common.f90 libry_global.f90 libry_nolocal.f90 libry_noopt.f90 libry.f90 libry_main.f90 -o libry.x
  #         #libry_par.f90 libry_common.f90 libry_global.f90 libry_local.f90 libry_opt.f90 libry.f90 libry_main.f90 -o libry.x

  ########

  # Hummel
  source /sw/modules/init.sh

  module switch env env/2015Q2-intel16-impi

  mpiifort -fp-model strict -free -real-size 64 -O0 \
           $(/sw/env/intel-16.0.1_impi-5.1.2/pkgsrc/2015Q2/bin/nf-config --fflags --flibs) \
           -g -traceback -ftrapuv -debug all -gen-interfaces \
           libry_par.f90 libry_common.f90 libry_global.f90 libry_nolocal.f90 libry_noopt.f90 libry.f90 libry_land.f90 libry_interface.f90 libry_main.f90 -o libry.x

           #libry_par.f90 libry_common.f90 libry_global.f90 libry_local.f90 libry_opt.f90 libry.f90 libry_main.f90 -o libry.x

  # -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -diag-enable=sc \
else

  gfortran -ffree-line-length-256 \
           libry_par.f90 libry_common.f90 libry_noglobal.f90 libry_local.f90 libry_noopt.f90 libry.f90 libry_land.f90 libry_interface.f90 libry_main.f90 -o libry.x

  # -fcheck=all 
fi

rm libry*.mod

if [[ $para == 1 ]]
then
  
  if [ -e "binaries/binary_${a0}" ]
  then
    if [[ $1 == "c" ]]
    then
      rm -r binaries/binary_${a0}
    else
      echo "binary directory exists!"
      exit
    fi
  fi

  mkdir binaries/binary_${a0}

  mv libry.x binaries/binary_${a0}/libry_${a0}.x
  cp libry*f90 binaries/binary_${a0}
  cp compile_libry.sh binaries/binary_${a0}

else

  if [ -e "binariesL/binary_${a0}" ]
  then
    if [[ $1 == "c" ]]
    then
      rm -r binariesL/binary_${a0}
    else
      echo "binary directory exists!"
      exit
    fi
  fi

  mkdir binariesL/binary_${a0}

  mv libry.x binariesL/binary_${a0}/libry_${a0}.x
  cp libry*f90 binariesL/binary_${a0}
  cp compile_libry.sh binariesL/binary_${a0}

fi

exit
