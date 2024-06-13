#!/bin/bash

################################
# Settings
################################
#
# name of binary
#
vname="Ntest06_paleo" #4b" #..
expname="Ntest_t06_PALEO" #4b_long" #-n8c36
#
model="libry_${vname}.x"
#
# initial number of species
#
nspecies=200
#
# simulation length
#
simlength=9 #600  # years
startout=8 #581   
endout=9 #600
#
# output interval
#
intvout=24 #168      # hours
#
# run control ( WATCH:year0=1958; ERA5:year0=1979 )
#
firstyear=1
firstCyear=1958
tscurrent=1
accumts=1
tposit=1
runlength=250           # restart interval in years
lrestart0=".false."     # needs to be false at 1st call!
#
setupR="libry_restart0"
setupRX="libry_restartX"
#
# switch (0/1) for submission to queuing system
#
lsf=1
#
# switches
#
lBSC=".false."
lNOHONO=".false."
#
# PATHS
#
if [[ $lsf == 1 ]]
then
  basedir=/home/bay2464/libry

  bindir=${basedir}/binaries

  rundir=/work/bay2464/simulations/${expname}

  ##############################################################inputdir=/work/bay2464/group_data/input_global
  inputdir=/work/bay2464/group_data/input_global/Paleo_new
else
  basedir="/home/philipp/work/projects/LIBRY-2"

  bindir=${basedir}/binariesL

  rundir=${basedir}/simulations/${expname}
  #rundir=/home/philipp/work/supervision/elevatedCO2_bryophytes/simulations/${expname}

 # inputdir=${basedir}/input_local/Hamburg_BOTshadow
 # inputdir=${basedir}/input_local/Hamburg_IPMopen
  inputdir=${basedir}/input_local/Hamburg_ERA5
 # inputdir=/home/philipp/Desktop/testrun_Hamburg_4diana/input_ERA5_Hamburg
 # inputdir=${basedir}/input_local/CostaRica_WATCH
 # inputdir=${basedir}/input_local/Almeria_Station
 # inputdir=${basedir}/input_local/Sardinia2
 # inputdir=${basedir}/input_local/Gradient_Italy
 # inputdir=/home/philipp/work/supervision/elevatedCO2_bryophytes/input_CR_Maaike
 # inputdir=/home/philipp/work/supervision/elevatedCO2_bryophytes/input_CR_global
 # inputdir=/home/philipp/work/cooperations/epiphyte_saturation/data_Soltis
fi
#
# lsf settings
#
ptn="std"  	# [ compute / compute2 ]
nn=16           # nodes ???
tpn=16          # max, tasks p. node [ (24 / 36) ] ???
#
runID=${expname} #${rundir##*/}
#
################################

# make run directory

if [ -e ${rundir} ]
then

  if [ -e ${rundir}/${setupR} ]
  then

    i=0
    while read line
    do
      a0=${line#*=}
      a[${i}]=${a0%,}
      i=$[ $i + 1 ]
    done < ${rundir}/${setupR}

    firstyear=${a[0]}
    firstCyear=${a[1]}
    tscurrent=${a[2]}
    accumts=${a[3]}
    tposit=${a[4]}
    lrestartX=${a[5]}

    if [[ $lrestartX == "T" ]]
    then
      lrestart0=".true."
    else
      lrestart0=".false"
    fi

    rm ${rundir}/${setupR}
    rm ${rundir}/libry_namelist
    rm ${rundir}/libry_batch
  else
    echo "run directory already exists, but no restart set!"
    exit
  fi
else

  mkdir ${rundir}

  cp ${basedir}/run_libry.sh ${rundir}
fi


# copy/link input files to run directory

if [[ $lsf == 1 ]]
then
  if [[ $lrestart0 == ".false." ]]
  then

  # ln -s ${inputdir}/srad.nc ${rundir}/srad.nc
  # ln -s ${inputdir}/lrad.nc ${rundir}/lrad.nc
  # ln -s ${inputdir}/tair.nc ${rundir}/tair.nc
  # ln -s ${inputdir}/rhum.nc ${rundir}/rhum.nc
  # ln -s ${inputdir}/rain.nc ${rundir}/rain.nc
  # ln -s ${inputdir}/snow.nc ${rundir}/snow.nc
  # ln -s ${inputdir}/wind.nc ${rundir}/wind.nc
  # ln -s ${inputdir}/libry_specpar ${rundir}/libry_specpar
  # ln -s ${inputdir}/landsea.nc ${rundir}/landsea.nc
  # ln -s ${inputdir}/biome.nc ${rundir}/biome.nc
  # ln -s ${inputdir}/LAI.nc ${rundir}/LAI.nc
  # ln -s ${inputdir}/SAI.nc ${rundir}/SAI.nc
  # ln -s ${inputdir}/SSA.nc ${rundir}/SSA.nc

    ln -s ${inputdir}/srad_new.nc ${rundir}/srad.nc
    ln -s ${inputdir}/lrad_new.nc ${rundir}/lrad.nc
    ln -s ${inputdir}/tair_new.nc ${rundir}/tair.nc
    ln -s ${inputdir}/rhum_new.nc ${rundir}/rhum.nc
    ln -s ${inputdir}/rain_new.nc ${rundir}/rain.nc
    ln -s ${inputdir}/snow_new.nc ${rundir}/snow.nc
    ln -s ${inputdir}/wind_new.nc ${rundir}/wind.nc
    ln -s /work/bay2464/group_data/input_global/libry_specpar ${rundir}/libry_specpar
    ln -s ${inputdir}/landsea_new.nc ${rundir}/landsea.nc
    ln -s ${inputdir}/biome_new.nc ${rundir}/biome.nc
    ln -s ${inputdir}/LAI_new.nc ${rundir}/LAI.nc
    ln -s ${inputdir}/SAI_new.nc ${rundir}/SAI.nc
    ln -s ${inputdir}/landsea_new.nc ${rundir}/SSA.nc

    if [[ \${lNOHONO} == .true. ]]; then ln -s ${inputdir}/NO_HONO_fSatBins ${rundir}/NO_HONO_fSatBins ; fi

    cp ${bindir}/binary_${vname}/${model} ${rundir}
  fi
else

  cd ${inputdir}

  cp srad lrad tair rhum rain snow wind libry_specpar LAI SAI SSA biome ${rundir}

  if [[ $lNOHONO == .true. ]]
  then
    cp NO_HONO_fSatBins ${rundir}
  fi

  cp ${bindir}/binary_${vname}/${model} ${rundir}
fi


# start simulation

cd ${rundir}

cat > libry_namelist << EOF
&librypar
year0=${firstyear},
cyear0=${firstCyear},
tsindata0=${tscurrent}
accts0=${accumts}
tpos0=${tposit}
lastyear=${simlength},
runperiod=${runlength},
tsl=3600,
yearout1=${startout},
yearoutX=${endout},
outint=${intvout},
nSites=1,
p_nspec=${nspecies},
fracTH_s_init=0.1,
lrestart=${lrestart0},
NOHONO=${lNOHONO},
BSCtypes=${lBSC},
llevels=3,
ltraits=0
&end
EOF
#outdiurnal=.true.,

if [[ $lsf == 1 ]]      # parallel run
then

cat > libry_batch << EOF
#!/bin/bash
#SBATCH --job-name=run_${runID}         # Specify job name
#SBATCH --partition=${ptn}              # Specify partition name (compute/compute2)
#SBATCH --nodes=${nn}                   # Specify number of nodes
#SBATCH --ntasks-per-node=${tpn}        # Specify number of tasks per node (max. 24/36)
#SBATCH --mail-type=FAIL                # Notify user by email in case of job failure
#SBATCH --mail-user=philipp.porada@uni-hamburg.de	# email address
#SBATCH --export=NONE					# ?
#SBATCH --output=job_${runID}.o%j       # File name for standard output
#SBATCH --error=job_${runID}.e%j        # File name for standard error output

source /sw/batch/init.sh

# Environment settings

module switch env env/2015Q2-intel16-impi

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

# launch program

mpirun ./${model}

echo "waiting for restart setup - start year was:" $firstyear

while [ ! -e $setupR ]
do
  sleep 0.1
done

echo "restart setup found"

if [ ! -e $setupRX ]
then
  bash ${rundir}/run_libry.sh
fi

EOF

  # submit script to cluster

  sbatch --partition=${ptn} libry_batch

else

  ./${model}

fi

exit

# rm srad.nc
# rm lrad.nc
# rm tair.nc
# rm rhum.nc
# rm rain.nc
# rm snow.nc
# rm wind.nc
# rm libry_specpar
# rm landsea.nc
# rm biome.nc
# rm LAI.nc
# rm SAI.nc
# rm SSA.nc
# if [[ \${lNOHONO} == .true. ]]; then rm NO_HONO_fSatBins; fi

