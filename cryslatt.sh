#!/bin/bash

# Perform lattice parameter optimization for HF and DFT, and for ACM functionals in a post-SCF fashion
# The script uses external programs such as acmxc, fit.py
# The enviroment variables CRY23_UTILS, PYTHONBIN must be defined
# The variable STORE must be modified
cryslattdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#set -x

# initial minimization search
do_min="yes"

#---------------rerun--------------------------------------
# if SCF calculation is already done, with rerun="no" the output is read
# with rerun="yes" calculation is rerun using GUESSP
rerun="no"

readonly="no"

#-------------no do_min--------------
npnt_nomin=13
step_nomin=0.0136
#------------with do_min-----------
npnt=12
step=0.0034

#--shells for cluster/bsse
ashell=3

#--------INPUT FILES DIRECTORY---------------------
STORE="$cryslattdir/REFDIR"
#--------------------------------------------------

systrange="Na"
crystalpath=$CRY23_UTILS

# if the MPI parallel version of CRYSTAL is installed, the variable MPI_ROOT
# should be defined in order for the following condition to work
# parallel
if [[ -z "${MPI_ROOT}" ]]; then
  crystalexe="$crystalpath/runcry23OMP"
  crystalprop=""
  atomseq="no"
else
  crystalexe="$crystalpath/runPcry23"
  crystalprop="$crystalpath/runprop23"
  atomseq="no"
fi
# serial
crystalseq="$crystalpath/runcry23"

ncpu=10

# the environment variable PYTHONBIN, containing the path
# of the python executable, must be present
fitexe="$PYTHONBIN $cryslattdir/fit.py"
eos="Birch-Murnaghan"
acmexe="$cryslattdir/acmxc.py -p crystal  "

wfunc="hpc"
#------------scf calculation--------------
func="HF"
#------------total energy------------------
# usually total-energy is scf, efunc=func, but it may not be scf
efunc=$func
finalconv=9
#
acmformula="genisi"
forcedmetal="no"
#----------------------
copyprev="yes"

#-----options and arguments---------------
while [ "$1" != "" ]; do
  if [ "$1" == "-f" ]; then
    shift
    func=$1
    efunc=$1
  elif [ "$1" == "-readonly" ]; then
    readonly="yes"
    rerun="no"
  elif [ "$1" == "-rerun" ]; then
    rerun="yes"
  elif [ "$1" == "-r" ]; then
    shift
    systrange=$1
  elif [ "$1" == "-n" ]; then
    shift
    ncpu=$1
  elif [ "$1" == "-conv" ]; then
    shift
    finalconv=$1
  elif [ "$1" == "-copyprev" ]; then
    copyprev="yes"
  elif [ "$1" == "-forcedmetal" ]; then
    forcedmetal="yes"
  elif [ "$1" == "-nocopyprev" ]; then
    copyprev="no"
  elif [ "$1" == "-do_min" ]; then
    do_min="yes"
  elif [ "$1" == "-nodo_min" ]; then
    do_min="no"
  elif [ "$1" == "-center" ]; then
    shift
    center=$1
  elif [ "$1" == "-npnt" ]; then
    shift
    npnt=$1
  elif [ "$1" == "-npnt_nomin" ]; then
    shift
    npnt_nomin=$1
  elif [ "$1" == "-ashell" ]; then
    shift
    ashell=$1
  elif [ "$1" == "-wfunc" ] || [ "$1" == "-w" ]; then
    shift
    wfunc=$1
  elif [ "$1" == "-acm" ]; then
    shift
    func="HF"
    efunc="acm"
    acmformula=$1
  elif [ "$1"  == "-h" ] || [ "$1" == "-help" ]; then
    echo
    echo "        ====== cryslatt.sh ======"
    echo "      driver for the cryslatt.sh script"
    echo
    echo "Possible options:"
    echo
    echo "-f THEORY: THEORY indicated the level of theory. Possible values: pbe, HF (default: HF)"
    echo
    echo "-readonly: make the script read results of existing calculations only"
    echo
    echo "-rerun: run again all calculations"
    echo
    echo "-r DIR: name of the system directory to perform calculations on"
    echo
    echo "-n NUM: NUM is the number of processors to be used (default: 10)"
    echo
    echo "-conv NUM: SCF tolerance on total energy (10^(-NUM)$) (default: 9)"
    echo
    echo "-copyprev: use the density matrix of previous calculations"
    echo
    echo "-nocopyprev: do not use the density matrix of previous calculations"
    echo
    echo "-forcedmetal: assume the system is metallic (the MP2 correlation energy is not computed)"
    echo
    echo "-do_min: use the initial minimization procedure"
    echo
    echo "-nodo_min: do not use the initial minimization procedure"
    echo
    echo "-center LAT: LAT is the starting lattice constant when no initial minimization procedure is used. It must be used along with -nodo_min"
    echo
    echo "-npnt NUM: NUM is the number of lattice constants to perform calculations on when the initial minimization procedure is used"
    echo
    echo "-npnt_nomin NUM: NUM is the number of lattice constants to perform calculations on when the initial minimization procedure is not used"
    echo
    echo "-ashell NUM: NUM is the number of shells of ghost atoms to be used"
    echo
    echo "-w METHOD, -wfunc METHOD: METHOD is the correlation energy approximation to be used. Possible values: pc, hpc, lda"
    echo
    echo "-acm FUNC: FUNC is the ACM functional to be used"
    echo
    echo "-h, -help: display help menu"
    echo
    exit 0
  fi
  shift
done



#
#-----------------------------------------------------------------#
#                          FUNCTIONS
#-----------------------------------------------------------------#

function make_scaling_factors_nomin {
#   prepare steps to be added
#   to the starting lattice constant
#   uses: $npnt_nomin
#   output: $factrange
  factrange=""
  ch=$(echo "$npnt_nomin % 2" | bc)
  if [ "$ch" == "1" ]; then
    npnt_nomin=$(echo "$npnt_nomin + 1" | bc)
  fi
# more points for large cell
  nnn=$(echo "$npnt_nomin / 2 - 2" | bc)
  nnp=$(echo "$npnt_nomin / 2 + 2" | bc)
  echo "nnn nnp" $nnn $nnp &>> $LOGFILE
  for i in $(seq $nnn); do
    val=$(echo "-$step_nomin*($nnn - ($i-1))" | bc -l)
    factrange=$(echo "$factrange $val")
  done
  factrange=$(echo "$factrange 0.0")
  for i in $(seq $nnp); do
    val=$(echo "$step_nomin*$i" | bc -l)
    factrange=$(echo "$factrange $val")
  done
}


function make_scaling_factors {
# create a list of scaling factors
# these will be added to reflat to create the lattice constants
# to be scanned
# uses: $npnt $step
  factrange=""
  ch=$(echo "$npnt % 2" | bc)
  if [ "$ch" == "1" ]; then
    npnt=$(echo "$npnt + 1" | bc)
  fi
  nn=$(echo "$npnt / 2" | bc)
  for i in $(seq $nn); do
    val=$(echo "-$step*($nn - ($i-1))" | bc -l)
    factrange=$(echo "$factrange $val")
  done
  factrange=$(echo "$factrange 0.0")
  for i in $(seq $nn); do
    val=$(echo "$step*$i" | bc -l)
    factrange=$(echo "$factrange $val")
  done
}


function prepare_crystal_input() {
#   set the lattice constant value
#   uses: $copyprev $lat $func $conv $ashell
#   output: $lastene $f9file $wasbadconv
# INPUT:
  name=$1
  prevdir=$2
  sed s/LLLL/$lat/g $name.d12 > tmp
  mv tmp $name.d12
#   select the right functional
  if [ "$func" == "HF" ]; then
    sed '/DFT/d' $name.d12 | sed '/FFFF/d' > tmp
    sed '/XXLGRID/d' tmp | sed '/#END/d' > $name.d12
    rm tmp &> /dev/null
  else
    sed s/FFFF/$func/g $name.d12 | sed s/#END/END/g > tmp
    mv tmp $name.d12
  fi

#   set convergence threshold
  sed s/CCCC/$conv/g $name.d12 > tmp
  mv tmp $name.d12

#   set number of shells of ghost atoms
  sed s/AAAA/$ashell/g $name.d12 > tmp
  mv tmp $name.d12

#   density matrix from a previous run, if present,
#   is used as SCF guess; if not present, the program
#   tries to locate the density matrix of a previously
#   simulated directory; if this also fails,
#   a calculation from scratch is run
  f9file="yes"
  if [ ! -f $name.f9 ]; then
    if [ "$copyprev" == "yes" ] && [ "$prevdir" != "" ]; then
      echo " using $prevdir as starting guess" | tee -a $LOGFILE
      cp $prevdir/$name.f9 .
#                diss active with copyprev
      sed '/ANDERSON/d' $name.d12  > tmp
      mv tmp $name.d12
    else
#   calculation from scratch: diis active
      f9file="no"
      sed '/GUESSP/d' $name.d12 | sed '/NODIIS/d' |  sed '/ANDERSON/d'  > tmp
      mv tmp $name.d12
    fi
  else
    lastene=$(grep "DETOT" scf.out | tail -n 1 | awk '{printf "%d",sqrt($6*$6)*1000000}')
    wasbadconv="no"
#    rerun a badly converged calculation
    if [ "$lastene" != "" ]; then
      if [ "$lastene" -gt "10" ]; then
        echo "bad convergence" $lastene "*10^-6"    | tee -a $LOGFILE
        sed '/ANDERSON/d' $name.d12  > tmp
        wasbadconv="yes"
        mv tmp $name.d12
       fi
    fi

# rerun: GUESSP and NODISS
  fi
}


function run_guessp() {
# run crystal program using GUESSP
# uses: $f9file $mycry $rrcpu $name(from run_program) $crystalprop
  isrun="yes"
  if [ "$f9file" == "yes" ]; then
    echo "$mycry $rrcpu $name $name ===>" >> $LOGFILE
    $mycry $rrcpu $name $name &>> $LOGFILE
#        if f9file is present then GUESSP is set
#        thus if you don't use name name, the file will not be COPIED in the parallel directory
  else
    echo "$mycry $rrcpu $name ===> " >> $LOGFILE
    $mycry $rrcpu $name &>> $LOGFILE
  fi
  if [ "$crystalprop" != "" ]; then
    $crystalprop $name $name &>> $LOGFILE
  fi
}

function run_program() {
# run crystal using its script or acmxc
# uses: $mycry $efunc $energy(from get_ene_cyc_dee) $forcedmetal
# uses: $rerun $toldee(from get_ene_cyc_dee) $conv $acmexe $acmformula $wfunc $cc $nn(from get_ene_cyc_dee)

  name=$1
  rcpu=$2
  exiterr=$3
  echo "run_program:",$1,$2,$3 >> $LOGFILE
  if [ "$rcpu" == "1" ]; then
    mycry=$crystalseq
    rrcpu=""
  else
    mycry=$crystalexe
    rrcpu=$rcpu
  fi
  isrun="no"

#--------------------acm----------------

  if [ "$efunc" == "acm" ]; then
    echo &>> $LOGFILE
    echo " Calculation in" $PWD | tee -a  $LOGFILE
    get_ene_cyc_dee scf.out
    wasnew="no"
    if [ "$energy" == "NA" ]; then
      echo " New calculation!"  |  tee -a  $LOGFILE
      wasnew="yes"
    fi

    ssmm=""
    if [ -e "METAL" ] || [ "$forcedmetal" == "yes" ]; then
      ssmm="--metal"
    fi

    ssrr=""
    if [ "$rerun" == "yes" ]; then # the -rerun flag works only with an old version of acmxc; now acmxc does not support rerun
      ssrr="--rerun"
    fi

    if  [ "$toldee" != "NA" ] && [ "$toldee" -lt "$conv" ]; then
      echo " TOLDEE differs " $toldee $conv  | tee -a  $LOGFILE
      ssrr="--rerun"
    fi

#    echo " $acmexe -f $acmformula -i $name $ssmm $ssrr -w $wfunc -n $rrcpu" | tee -a $LOGFILE
    echo " $acmexe -f $acmformula -i $name $ssmm -w $wfunc -n $rrcpu" | tee -a $LOGFILE

# ATOM contains both calculation with --metal and without
# - yesmp2.out is the exact computed MP2 value
# - mp2.out contains the current one (computed or -inf)
    if [ -e mp2.out ]; then
      kmp2=$(grep inf mp2.out | wc -l)

      if [ "$kmp2" == "1" ]; then # mp2=-inf
        if [ "$ssmm" != "--metal" ]; then
# mp2 calculation is needed, but file contains mp2=-inf
          rm mp2.out
          cp yesmp2.out mp2.out
# do mp2 again
        fi
      else #mp2 calc done
        if [ "$ssmm" == "--metal" ]; then
# caso metal , file mp2=-values
          if [ ! -e yesmp2.out ]; then
# save for future run
            mv mp2.out yesmp2.out
# redo just -inf
            rm mp2.out
          fi
        else
# caso mp2, file mp2=-valus
        if [ ! -e yesmp2.out ]; then
          cp mp2.out yesmp2.out
        fi
      fi
    fi
 # mp2.out exists
    fi
    if [ "$name" == "inputatom" ]; then
      realout=acmxc.$namemethoda.out
    else
      realout=acmxc.$namemethod.out
    fi

    echo 'VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV' >> $LOGFILE

    if [ "$readonly" == "no" ]; then
      if [ "$rrcpu" == "" ]; then
        $acmexe -f $acmformula -i $name           $ssmm $ssrr -w $wfunc | tee $realout >> $LOGFILE
      else
        $acmexe -f $acmformula -i $name -n $rrcpu $ssmm $ssrr -w $wfunc | tee $realout >> $LOGFILE
      fi
      acmstatus="$?"

      echo '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' >> $LOGFILE

    else
      echo " Readonly!" | tee -a $LOGFILE
    fi

    if [ ${acmstatus} != 0 ]; then
      echo "ACM failed!"    | tee -a $LOGFILE
      if [ "$wasbadconv" == "yes" ]; then
        echo "Remove" $PWD "and restart ! "
      fi
      echo " see file $PWD/$realout"    | tee -a $LOGFILE
      if [ "exiterr" == "yes" ];  then
        exit 1
      fi
    fi
  else
#----------------------dft or hf-----------
    echo &>> $LOGFILE
    echo " Calculation in" $PWD | tee -a  $LOGFILE
    if [ -e $name.out ]; then
      get_ene_cyc_dee $name.out
      if [ "$energy" != "NONE" ]; then
        echo " already done and converged in " $cycles "cycles" | tee -a $LOGFILE
        if [ "$rerun" == "yes" ] || [ "$toldee" -lt "$conv" ]; then
          if  [ "$toldee" -lt "$conv" ]; then
            echo " TOLDEE differs " $toldee $conv                  | tee -a  $LOGFILE
          fi
          echo " rerun:" | tee -a  $LOGFILE
          run_guessp
        fi
	cp $name.out scf.out
      else
        echo " calculation error; toomany=$cc energy=$nn => rerun" | tee -a $LOGFILE
        if [ "$readonly" == "no" ]; then
          run_guessp
        else
          echo " Readonly!" | tee -a $LOGFILE
        fi
	cp $name.out scf.out
      fi
    else  # name.out not present
      echo " New calculation!"
      run_guessp
      cp $name.out scf.out
    fi
  fi
}

function get_ene_cyc_dee() {
# extract data from CRYSTAL SCF output file
# output: $nn $cc

  filename=$1
  energy="NA"
  cycles="NA"
  toldee="NA"
  if [ -e $filename ]; then
    nn=$(grep "TOTAL ENERGY(" $filename | wc -l)
    cc=$(grep "SCF ENDED - TOO MANY CYCLES" $filename | wc -l)

    if [ "$nn" -gt "0" ] && [ "$cc" -eq "0" ]; then
#energy present and toomany not present
      energy=$(grep "TOTAL ENERGY(" $filename | tail -n 1 | awk '{print $4}')
      cycles=$(grep "TOTAL ENERGY(" $filename | tail -n 1 | awk -F\) '{print $3}' | awk '{print $2}')
      toldee=$(grep "**** TOLDEE" scf.out | awk '{print $12}')
    elif [ "$cc" -eq "0" ] && [ "$cc" -eq "0" ]; then
# energy not present, there must be an error
      grep "ERROR \*" $filename | tee -a $LOGFILE
      energy="NONE"
    else
# done but not converged
      energy="NONE"
    fi
  fi
}


function extract_results {
# extract results from acmxc output file or CRYSTAL output file
# uses: $energy $efunc $realout $isrun

# filename is fixed here
  if [ ! -e scf.out ];  then
    echo "$PWD/scf.out does not exist"
    echo "switch off readonly or use -rerun"
    exit 1
  fi
  get_ene_cyc_dee scf.out
  if [ "$energy" != "NONE" ]; then
    if [ "$efunc" == "acm" ]; then #uses acmxc output
# reads total energy and exchange energy from acmxc output
      volume=$(grep "PRIMITIVE CELL - CENTRING CODE" scf.out | tail -n 1 | awk '{print $8}')
      nn=$(grep "Total energy:" $realout | wc -l)
      if [ "$nn" -gt "0" ]; then
        energy=$(grep "Total energy:" $realout | awk '{print $3}') # acmxc was changed: now the output contains "Total energy", not "ACM Total energy", so $3 is right, $4 is wrong
      else
        energy="NONE"
      fi
      cycles=$(grep "TOTAL ENERGY(" scf.out | tail -n 1 | awk -F\) '{print $3}' | awk '{print $2}')
      nn=$(grep "HF-exchange" $realout | wc -l)
      if [ "$nn" -gt "0" ]; then
        exx=$(grep "HF-exchange" $realout | awk '{print $3}')
      else
        exx="NONE"
      fi
    else #use crystal output
      volume=$(grep "PRIMITIVE CELL - CENTRING CODE" scf.out | tail -n 1 | awk '{print $8}')
      get_ene_cyc_dee scf.out
      exx=$(grep " FOCK EXCHANGE" scf.out | tail -n 1 | awk '{print $4}')
    fi

    if [ "$isrun" == "yes" ]; then
      echo " Calculation done in " $cycles "cycles"
    fi

  fi
}


function make_fit () {
# use the python script fit.py to fit the curve to EOS
# uses; $fitexe

  datatofit=$1
# ./store.$namemethod.dat
  fitname=$2
# fitname=fit.$namemethod.out
  printcoe=$3

  echo >> $LOGFILE
  echo "makefit()" $1 ">" $2 | tee -a $LOGFILE

#   read reference values
  reflat=$(head -n 1 reference | awk '{printf "%12.8f",$1}')
  refbmod=$(head -n 2 reference | tail -n 1 | awk '{printf "%12.8f",$1}')
  refcoe=$(head -n 3 reference | tail -n 1 | awk '{printf "%12.8f",$1}')
  latyp=$(tail -n 1 reference | awk '{print $1}')

  if [ -s $datatofit ]; then
    len=$(wc -l $datatofit | awk '{print $1}')
  else
    len=0
  fi

  if [ "$len" -ge "3" ]; then
    echo "FIT" $datatofit "REF:" $reflat $refbmod $latyp | tee -a $LOGFILE
    echo "$datatofit :" >> $LOGFILE
    cat $datatofit >> $LOGFILE
    echo " $fitexe $datatofit $reflat $refbmod $latyp &> $fitname " >> $LOGFILE
    $fitexe $datatofit $reflat $refbmod $latyp &> $fitname

    grep ERROR $fitname | tee -a $LOGFILE
    errcount=$(grep ERROR $fitname | wc -l)
    echo $errcount >> $LOGFILE
    if [ "$errcount" -gt "0" ]; then
      echo "WARNING: Errors  in the FIT !!! " | tee -a $LOGFILE
#            exit
    fi
    grep Warning $fitname | tee -a $LOGFILE
    warncount=$(grep Warning $fitname | wc -l)

    if [ "$warncount" -gt "0" ]; then
      echo "WARNING: warning  in the FIT !!! " | tee -a $LOGFILE
# hai commentato la riga seguente
#            exit
    fi
#    mv wfit0.dat wfit0.$namemethod.dat
#    mv fit0.dat fit0.$namemethod.dat
#   read volume and compute lattice constant
#   head  -n 1 means the not weighted results
#
# $fitname file:
#Method                E_0 (Ha)    V_0 (A^3)    B_0 (GPa)     B' (Ha/A^6)   Error
#--------------------------------------------------------------------------------
#Murnaghan                -0.0208925      25.769        8.73         3.01      4.73e-12
    if  [ "$printcoe" == "yes" ]; then
#
# transform in eV
#
      coe=$(grep "$eos" $fitname | head -n 1 | awk '{printf "%15.8f", -$2*27.2113961}')
    else
      coe=$(grep "$eos" $fitname | head -n 1 | awk '{printf "%15.8f", $2}')
    fi

    vol=$(grep "$eos" $fitname | head -n 1 | awk '{print $3}')
    echo "vol coe" $vol $coe | tee -a $LOGFILE
    if [ "$latyp" == "A1" ] || [ "$latyp" == "A4" ] || [ "$latyp" == "B1" ] || [ "$latyp" == "B3" ]; then
# ccp: Ag, Al, Cu, Ca
      tmp1=$(echo "$vol/2" | bc -l)
      tmp2=$(echo "$tmp1" | awk '{print $1^0.33333}')
      lat=$(echo "2*$tmp2" | bc -l | awk '{printf "%8.3f",$1}')
    elif [ "$latyp" == "A2" ] ; then
# bcc: Li, Na
      tmp1=$(echo "$vol/4" | bc -l)
      tmp2=$(echo "$tmp1" | awk '{print $1^0.33333}')
      lat=$(echo "2*$tmp2" | bc -l | awk '{printf "%8.3f",$1}')
    fi
#   read bulk modulus
    bulk=$(grep "$eos" $fitname | head -n 1 | awk '{print $4}')
# save results of the fit
    echo "$lat"   > results_$fitname
    echo "$bulk" >> results_$fitname
    echo "$coe"  >> results_$fitname
    echo
    lattdiff=$(echo "$lat - $reflat" | bc -l | awk '{printf "%8.3f",$1}')
    bmdiff=$(echo "$bulk - $refbmod" | bc -l | awk '{printf "%8.2f",$1}')
    if [ "$printcoe" == "yes" ]; then
      coediff=$(echo "$coe - $refcoe" | bc -l | awk '{printf "%14.8f",$1}')
    fi

    echo "Lattice constant: $lat   reference: $reflat   error: $lattdiff" | tee -a $LOGFILE lattopt_$fitname
    echo "Bulk modulus: $bulk  reference: $refbmod  error: $bmdiff"       | tee -a $LOGFILE lattopt_$fitname
    if  [ "$printcoe" == "yes" ]; then
      echo "Cohesive energy: $coe   reference: $refcoe   error: $coediff"   | tee -a $LOGFILE lattopt_$fitname
    else
      echo "Total energy: $coe"                                            | tee -a $LOGFILE lattopt_$fitname
    fi
    echo
  else
    echo "only" $len "point: cannot fit !"   | tee -a $LOGFILE
    echo | tee -a $LOGFILE
  fi
}


function run_at_lat (){
# run calculation at specific lattice
# uses: $latold $ncpu $energy $volume $cycles $exx
#INPUT
  lattic=$1
  exiterr=$2
# OUTPUT:  energydiff energybulk energyatom
#                     exxbulk    exxatom
#                     volume
  lat_or=$lat
  lat=$lattic

  pdir=""
#latold is "" for the first simulation or if the previous simulation fails
  if [ "$latold" != "" ]; then
    pdir=$PWD/LAT_$latold
  fi

  echo "run_at_lat:" $lattic,$lattic,$exiterr,$latold,$pdir >> $LOGFILE

# bulk calculation
  mkdir LAT_$lattic &> /dev/null
  cd LAT_$lattic
  cp -r $STORE/${syst}.SAVE/* .
  prepare_crystal_input  input $pdir
  run_program            input $ncpu $exiterr
  extract_results
  energybulk=$energy
  volumebulk=$volume
  cyclesbulk=$cycles
  exxbulk=$exx
  echo "energy vol cycles exx:" $energybulk   $volumebulk    $cyclesbulk $exxbulk >> $LOGFILE

# atom calculation
  cd ATOM
  prepare_crystal_input inputatom ""
  run_program           inputatom $ncpu $exiterr
  extract_results

  energyatom=$energy
  exxatom=$exx
  cyclesatom=$cycles
  echo "energy cycles exx:" $energyatom    $cyclesatom $exxatom >> $LOGFILE
  cd ..
  if [ "$energybulk" != "NONE" ] && [ "$energyatom" != "NONE" ]; then
    energydiff=$(echo $energybulk $energyatom | awk '{printf "%12.8f", $1-$2}')
  else
    energydiff="NONE"
  fi
  lat=$lat_or
  cd ..
}

function run_lat_range() {
# run the program on set of lattice constants
# uses: $latrange $namemethod $namemethoda $efunc $energydiff
# uses: $energybulk $energyatom $volumebulk $exxbulk $exxatom
# uses: $lattic
#   outputfile : $1
#   exit on error : $2
  outputfilestore=$1
  exiterr=$2
  echo "run_lat_range:" $1,$2 >> $LOGFILE
  rm tmplatene.$namemethoda
  if [ "$efunc" == "acm" ]; then
    rm $1.acmdat
  else
    rm $1.xcdat
  fi
  for lat in $latrange; do
    echo >> $LOGFILE
    echo  ">>>>> lat = $lat A" | tee -a $LOGFILE

    run_at_lat $lat $exiterr

    if [ "$efunc" == "acm" ]; then
      wscf=$(grep "SCF energy" LAT_$lat/acmxc.$namemethod.out | awk '{print $3}')
      awscf=$(grep "SCF energy" LAT_$lat/ATOM/acmxc.$namemethoda.out | awk '{print $3}')
      ww=$(grep W_inf         LAT_$lat/acmxc.$namemethod.out | awk '{print $3}')
      w1=$(grep W1_inf        LAT_$lat/acmxc.$namemethod.out | awk '{print $3}')
      aww=$(grep W_inf        LAT_$lat/ATOM/acmxc.$namemethoda.out | awk '{print $3}')
      aw1=$(grep W1_inf       LAT_$lat/ATOM/acmxc.$namemethoda.out | awk '{print $3}')
      zmp2=$(grep "MP2 corr"  LAT_$lat/ATOM/acmxc.$namemethoda.out | awk '{print $4}')
      cc=$(grep "Correlation energy" LAT_$lat/acmxc.$namemethod.out | awk '{print $3}') # acmxc was changed: now the output contains "Correlation energy", not "ACM Correlation energy", so $3 is right, $4 is wrong
      acc=$(grep "Correlation energy" LAT_$lat/ATOM/acmxc.$namemethoda.out | awk '{print $3}') # acmxc was changed: now the output contains "Correlation energy", not "ACM Correlation energy", so $3 is right, $4 is wrong
      echo $lat $wscf $ww $w1 $cc  "   "  $awscf $aww $aw1 $zmp2 $acc >> $1.acmdat
    elif [ "$efunc" != "HF" ]; then # only for pbe
      cc=$(grep "CORRELATION ENERGY" LAT_$lat/input.out | awk '{print $6}')
      acc=$(grep "CORRELATION ENERGY"  LAT_$lat/ATOM/inputatom.out | awk '{print $6}')
      echo $lat $cc $acc >> $1.xcdat
    fi
    echo $lat $energydiff >> tmplatene.$namemethoda
    echo " vol= $volumebulk A^3   E= $energydiff au  cycles= $cyclesbulk"   | tee -a $LOGFILE
    if [ "$energydiff" != "NONE" ]; then
      echo " $volumebulk  $energydiff $energybulk $energyatom $exxbulk $exxatom $lat" >> $outputfilestore
      echo " $volumebulk  $energydiff $energybulk $energyatom $exxbulk $exxatom $lat" > LAT_$lattic/$outputfilestore
    else
      if [ "$exiterr" == "yes" ]; then
        echo "Error in scf procedure"
        exit 1
      fi
    fi

    if  [ "$energybulk" != "NONE" ]; then
      echo " $volumebulk  $energybulk $lat"                                           >> $outputfilestore.tot
      echo " $volumebulk  $energybulk $lat"                                           >> LAT_$lattic/$outputfilestore.tot
    else
      if [ "$exiterr" == "yes" ]; then
        echo "Error in scf procedure" | tee -a $LOGFILE
        exit 1
      fi
    fi

    echo " lat = $lat  vol= $volumebulk A^3   E= $energydiff au" >> $LOGFILE

#non converged, do not use in prevdir
    if  [ "$energydiff" != "NONE" ]; then
      latold=$lat
    else
      latold=""
    fi

  done
}

function find_parabolic_minimum() {
# use a parabolic fit to estimate the starting lattice constant
# from the initial minimization procedure
#INPUT : $1
  inputfile=$1
#OUTPUT: $minvol

  sort -n -k 2 $inputfile | uniq | head -n 3 > tmpdata.dat
  x1=$(head -n 1 tmpdata.dat | awk '{print $1}')
  y1=$(head -n 1 tmpdata.dat | awk '{print $2}')
  x2=$(head -n 2 tmpdata.dat | tail -n 1 | awk '{print $1}')
  y2=$(head -n 2 tmpdata.dat | tail -n 1 | awk '{print $2}')
  x3=$(head -n 3 tmpdata.dat | tail -n 1 | awk '{print $1}')
  y3=$(head -n 3 tmpdata.dat | tail -n 1 | awk '{print $2}')
  echo "parabola" >> $LOGFILE
  echo "point 1: $x1 $y1" >> $LOGFILE
  echo "point 2: $x2 $y2" >> $LOGFILE
  echo "point 3: $x3 $y3" >> $LOGFILE
#   tmp1=denom
  tmp1=$(echo "($x1-($x2))*($x1-($x3))*($x2-($x3))" | bc -l)
#   tmp2=a*denom
  tmp2=$(echo "$x3*($y2-($y1)) + $x2*($y1-($y3)) + $x1*($y3-($y2))" | bc -l)
#   tmp3=b*denom
  tmp3=$(echo "($y2-($y3))*($x1)*($x1) + ($y3-($y1))*($x2)*($x2) + ($y1-($y2))*($x3)*($x3)" | bc -l)
#   c is not required
  a=$(echo "$tmp2/($tmp1)" | bc -l)
  b=$(echo "$tmp3/($tmp1)" | bc -l)
  minvol=$(echo "-($b / (2 * $a))" | bc -l | awk '{printf "%10.6f",$1}')

  rm tmpdata.dat
}


function check_parabolic_convergence {
# check if the parabolic fit reaches convergence
# uses: $minvol $oldvol

  echo "check_parabolic_convergence" $minvol $oldvol >> $LOGFILE
#   ch=1 OK, ch=0 not converged
  tol=0.7
  ch=$(echo "sqrt(($minvol - $oldvol)^2) < $tol " | bc -l)
}

function lat2vol () {
# convert lattice constant to volume
# uses: $latyp
# INPUT : $1
  lattic=$1
# OUTPUT : $newvol

  if [ "$latyp" == "A1" ] || [ "$latyp" == "A4" ] || [ "$latyp" == "B1" ] || [ "$latyp" == "B3" ]; then
    newvol=$(echo "2*(($lattic/2)^3)" | bc -l | awk '{printf "%10.6f",$1}')
  elif [ "$latyp" == "A2" ] ; then
    newvol=$(echo "4*(($lattic/2)^3)" | bc -l | awk '{printf "%10.6f",$1}')
  fi

}

function vol2lat () {
# convert volume to lattice constant
# uses: $latyp
#INPUT : $1
  volum=$1
#OUTPUT: $newlat

  if [ "$latyp" == "A1" ] || [ "$latyp" == "A4" ] || [ "$latyp" == "B1" ] || [ "$latyp" == "B3" ]; then
    tmp1=$(echo "$volum/2" | bc -l)
    tmp2=$(echo "$tmp1" | awk '{print $1^0.33333}')
#        newlat=$(echo "2*$tmp2" | bc -l | awk '{printf "%8.4f",$1}')
    newlat=$(echo "2*$tmp2" | bc -l | awk '{print $1}')
  elif [ "$latyp" == "A2" ] ; then
    tmp1=$(echo "$volum/4" | bc -l)
    tmp2=$(echo "$tmp1" | awk '{print $1^0.33333}')
#        newlat=$(echo "2*$tmp2" | bc -l | awk '{printf "%8.4f",$1}')
    newlat=$(echo "2*$tmp2" | bc -l | awk '{print $1}')
  fi
}


function make_latrange () {
# obtain set of lattice constants
# uses: $do_min $center $factrange
#INPUT : $1
  startl=$1
#OUTPUT : $latrange
  if [ "$do_min" == "no" ] && [ -z "$center" ]; then
    make_scaling_factors_nomin
  else
    make_scaling_factors
  fi
  echo "FACTRANGE:" $factrange &>> $LOGFILE
  latrange=""
  for fact in $factrange; do
    ll=$(echo "$startl*(1 + $fact)" | bc -l | awk '{printf "%7.2f",$1}')
    latrange=$(echo "$latrange $ll")
  done
  echo "LATRANGE:" $latrange &>> $LOGFILE
}

function shrink_latrange() {
  ls -d LAT*/ | cut -f1 -d'/' | awk -F_ '{print $2}' > tmplist
  tmplen=$(wc -l tmplist | awk '{print $1}')
  declare -a myArray
  myArray=($(cat tmplist))
# myArray is a list of simulation already done
  finalrange=""
  for la in $latrange; do
    ffound="no"
    length=${#myArray[@]}
    for (( j=0; j<length; j++ )); do
      str=${myArray[$j]}
      if [ "$la" != "$str" ]; then
        dk=$(echo $str $la | awk '{printf "%d\n",sqrt(($1-$2)*($1-$2))*100}')
        if [ "$dk" -eq "0" ]; then
          finalrange=$(echo "$finalrange $str")
          ffound="yes"

          break
        fi
      fi
    done
    if [ "$ffound" == "no" ]; then
      echo $la $ffound   &>> $LOGFILE
      finalrange=$(echo "$finalrange $la")
    else
      echo $la $ffound $str  &>> $LOGFILE
    fi

  done
  echo "NEW LATTICES:" $finalrange
  echo $finalrange | tr ' ' '\n' | uniq > frfile.tmp
  finalrange=$(cat frfile.tmp | tr ' ' '\n' | sort | tr '\n' ' ')
  echo "FINAL RANGE:" $finalrange &>> $LOGFILE
}

#-----------------------------------------------------------------#
#-----------------------------------------------------------------#
#                            MAIN CODE                            #
#-----------------------------------------------------------------#
#-----------------------------------------------------------------#

echo "-----[ [ [ CRYSTAL LATTICE OPTIMIZATION ] ] ] ----"
echo "        by F. Della Sala and E. Fabiano            "
echo

# read options and arguments
if [ "$efunc" == "acm" ]; then
  namemethod=$(echo $acmformula $wfunc   | awk '{ printf "ACM-%s-%s",$1,$2}')
  namemethoda=$(echo $acmformula $wfunc   | awk '{ printf "ACM-%s-%s",$1,$2}')
  if [ "$forcedmetal" == "yes" ]; then
    namemethoda=$(echo $acmformula $wfunc   | awk '{ printf "ACM-%s-metal-%s",$1,$2}')
  fi
elif [ "$efunc" == "$func" ]; then
  namemethod=$(echo $func    | awk '{ printf "%s",$1,$2}')
  namemethoda=$namemethod
else
  namemethod=$(echo $efunc $func   | awk '{ printf "%s[%s]",$1,$2}')
  namemethoda=$namemethod
fi

echo "NAME METHOD:" $namemethod $namemethoda
echo "copy previous:" $copyprev
if [ "$efunc" == "acm" ]; then
 echo "forced metal:" $forcedmetal
fi
echo "readonly" $readonly
echo "ncpu:" $ncpu
echo "final conv:" $finalconv
echo "shell in atom:" $ashell

# start loop on systems
for syst in $systrange; do
  namedir=$syst"_"$func
  echo "==== $namedir ===="

# check if reference directories are present
  mkdir $namedir &> /dev/null
  if [ ! -d $STORE/${syst}.SAVE/ ]; then
    echo "ERROR! directory" $STORE/${syst}.SAVE/ "is not present"
    exit 1
  fi
  if [ ! -d $STORE/${syst}.SAVE/ATOM ]; then
    echo "ERROR! directory" $STORE/${syst}.SAVE/ATOM "is not present"
    exit 1
  fi

# the reference directories are copied
  cp -r $STORE/${syst}.SAVE/* $namedir
  cd $namedir

  echo "--------------------------CRYSLATT------------------------------" >> cryslatt.$namemethoda.log
  LOGFILE=$(realpath cryslatt.$namemethoda.log)
  echo $LOGFILE
  date                                         >> $LOGFILE
  echo "cryslatt $1 $2 $3 $4 $5 $6 "           >> $LOGFILE
  echo $namemethoda                            >> $LOGFILE

#   read reference lattice constant and lattice type
  if [ -s  $STORE/${syst}.SAVE/reference ]; then
    reflat=$(head -n 1 $STORE/${syst}.SAVE/reference | awk '{printf "%12.8f",$1}')
    latyp=$(tail -n 1 $STORE/${syst}.SAVE/reference | awk '{print $1}')
    echo "reflat, latytp:  $reflat ,  $latyp" | tee -a $LOGFILE
  else
    echo "Error! reference data in " $STORE/${syst}.SAVE/reference " not found"
    exit 1
  fi

#---------------optional search of the parabolic minimum---------------------------------
  if [ "$do_min" == "yes" ]; then
    conv=7 #calculations for the parabolic minimum are run with a lower tolerance
    echo  | tee -a $LOGFILE
    echo "Minimization Procedure: yes " | tee -a $LOGFILE

#--------------- initial run----------------------------

    larguess="no"

    if [ "$efunc" == "HF" ]; then
      larguess="yes"
      echo "HF force LargeGuess=true" | tee -a $LOGFILE
    fi

# if an error is found in a restart file, use a larger interval
    if [ -e  tmplatene.restart.$namemethoda ]; then
      echo "Restart file found:" tmplatene.restart.$namemethoda
      cat tmplatene.restart.$namemethoda
      nerr=$(grep NONE tmplatene.restart.$namemethoda | wc -l )
      if [ "$nerr" != "0" ]; then
        larguess="yes"
      fi
    fi

    echo "LargeLattice Guess: $larguess" | tee -a $LOGFILE

# values of the guessing procedure are different for HF and ACM
#----------------first 3 run are fixed => all digits----
#       start with
#         DFT        -5%   0  +10%
#         HF         +5% +10% +15%

    if [ "$larguess" == "no" ]; then
      latrange=$(echo "$reflat*(1.0-0.05)" | bc -l )
      val=$(echo "$reflat*(1.0)" | bc -l )
      latrange=$(echo "$latrange $val")
    else
      val=$(echo "$reflat*(1.0+0.05)" | bc -l )
      latrange=$(echo "$latrange $val")
    fi

    val=$(echo "$reflat*(1.0+0.1)" | bc -l )
    latrange=$(echo "$latrange $val")

    if  [ "$larguess" == "yes" ]; then
      val=$(echo "$reflat*(1.0+0.15)" | bc -l )
      latrange=$(echo "$latrange $val")
    fi

    filestore=storemin.$namemethoda.dat
    rm -rf $filestore &> /dev/null

    echo "Initial Guess using: " $latrange | tee -a  $LOGFILE
#                            do not exit on error
#                            |
    run_lat_range $filestore no

#----try to correct initial range------------
    nerr=$(grep NONE tmplatene.$namemethoda | wc -l )
    echo $nerr >> $LOGFILE
    if [ "$nerr" != "0" ]; then
      cat tmplatene.$namemethoda
      cp tmplatene.$namemethoda tmplatene.restart.$namemethoda
      exit 1
    else
# there are no errors, thus restart is not needed (?)
# it is needed for large guess
      cat tmplatene.restart.$namemethoda tmplatene.$namemethoda | sort | uniq   > tmpk
      mv tmpk tmplatene.restart.$namemethoda
    fi
#------------ minimization procedure-------------
#------------all points must converge--------
    latold=""  # first point from scratch
    lat2vol $reflat
    refvol=$newvol
    echo "refvol" $refvol >> $LOGFILE
#       output $newvol

    oldvol=$newvol
    ch=0

    while [ "$ch" == "0" ]; do

      find_parabolic_minimum $filestore
                #               output is $minvol

      check_parabolic_convergence
      if [ "$ch" == "1" ]; then break; fi

      oldvol=$minvol

      vol2lat $minvol
                    #           output is $newlat

#  new lat may depend on convergence issue and methods !

      latrange=$(echo $newlat | awk '{printf "%7.2f",$1}')
      shrink_latrange
      latrange=$finalrange

# if a problem occurs and lattice constant equals 0, exit
      if [ $(echo "$latrange" | xargs) == "0.00" ]; then
        echo "lattice=zero"
        exit 1
      fi

      if [ "$latold" == "$latrange" ]; then break; fi

      run_lat_range $filestore yes

      latold=$latrange
    done

# set the lattice constant corresponding to minvol as starting one
    vol2lat $minvol
    startlat=$newlat

# when a starting lattice constant is given as an argument
    elif [ "$do_min" == "no" ] && [ ! -z "$center" ]; then
      startlat="$center"
      echo $startlat
    else
# if the minimization is skipped the starting lattice constant is the reference one
    startlat=$reflat
    if [ "$efunc" == "HF" ]; then
      startlat=$(echo $startlat | awk '{printf "%7.2f",$1*1.05}')
    fi
  fi

#----------final run for the fit--------------
  latold=""
  conv=$finalconv
  echo  | tee -a $LOGFILE
  echo "FINAL RUN: with $startlat" | tee -a $LOGFILE

# values are stored in a file
  filestore=store.$namemethoda.dat
  rm -rf $filestore &> /dev/null
  rm -rf $filestore.* &> /dev/null


  make_latrange $startlat
  echo "GUESS RANGE:" $latrange
  shrink_latrange
# some calculations may have been done with toldee=7
  latrange=$finalrange
  echo "FINAL RANGE:" $latrange
#                          do not exit on error
  run_lat_range $filestore no


#---------   perform the fit and print the results -------------
  echo

  c1=$(echo $reflat) # reference lattice constant
#c2=`head -n 1 $filestore | awk '{print $2}'`
  c3=$(head -n 1 $filestore | awk '{print $3}') # bulk ACM total energy of the first point
  c4=$(head -n 1 $filestore | awk '{print $4}') # atom ACM total energy of the first point
  c5=$(head -n 1 $filestore | awk '{print $5}') # bulk exchange energy of the first point
  c6=$(head -n 1 $filestore | awk '{print $6}') # atom exchange energy of the first point
# for cmin you need the fit beforehand
  cmin=0
  echo $c1 $revcol $cmin
  awk -v cc1=$c1 -v cc2=$cmin -v cc3=$c3 -v cc4=$c4 -v cc5=$c5 -v cc6=$c6 '{print $7-cc1,$2-cc2,$3-cc3,$4-cc4,$5-cc5,$6-cc6}' $filestore > $filestore.0

  if [ "$efunc" == "acm" ]; then
#---------store zero acm------------
    c1=$(echo $reflat)
#bulk
    c2=$(head -n 1 $filestore.acmdat | awk '{print $2}') # bulk SCF energy
    c3=$(head -n 1 $filestore.acmdat | awk '{print $3}') # bulk W
    c4=$(head -n 1 $filestore.acmdat | awk '{print $4}') # bulk W1
    c5=$(head -n 1 $filestore.acmdat | awk '{print $5}') # bulk ACM correlation energy
# atom
    c6=$(head -n 1 $filestore.acmdat | awk '{print $6}') # atom SCF energy
    c7=$(head -n 1 $filestore.acmdat | awk '{print $7}') # atom W
    c8=$(head -n 1 $filestore.acmdat | awk '{print $8}') # atom W1
    c9=$(head -n 1 $filestore.acmdat | awk '{print $9}') # atom MP2 correlation energy
    c10=$(head -n 1 $filestore.acmdat | awk '{print $10}') # atom ACM correlation energy

#       bulk                                                   atom
    awk -v cc1=$c1 -v cc2=$c2 -v cc3=$c3 -v cc4=$c4 -v cc5=$c5 -v cc6=$c6 -v cc7=$c7 -v cc8=$c8 -v cc9=$c9 -v cc10=$c10 '{print $1-cc1,$2-cc2,$3-cc3,$4-cc4,$5-cc5,"    ",$6-cc6,$7-cc7,$8-cc8,$9-cc9,$10-cc10}' $filestore.acmdat > $filestore.acmdat.0
#-----------------------------------------HF     HF+corr
    awk '{printf "%10.5f %16.6f %16.6f \n       ",$1, $2-$6, $2-$6+$5-$10}' $filestore.acmdat.0 > $filestore.corrcomp.0
    awk '{printf "%10.5f %16.6f %16.6f %16.6f \n",$1, $5,   $10, $5-$10}' $filestore.acmdat > $filestore.abscorr
#                                     latt corrbulk corratom diff
    awk '{printf "%10.5f %16.6f %16.6f %16.6f \n",$1, $5,  $10,  $5-$10}' $filestore.acmdat.0 > $filestore.abscorr.0
  fi
#------------------------------------------
  make_fit $filestore      fitdiff.$namemethoda.out yes

#--
  cat fitdiff.$namemethoda.out &>> $LOGFILE

    # this is total energy, using namemethod not namemethoda
    #                       but we used for simplicity
  make_fit $filestore.tot  fittote.$namemethoda.out no
#                                                    donotprintcoe

  c1=$(echo $reflat)
  cmin=$(tail -n 1 results_$fitname | awk '{print $1}')
  echo $c1 $cmin
  awk -v cc1=$c1 -v cc2=$cmin  '{print $3-cc1,$2-cc2}' $filestore.tot > $filestore.tot.0
  echo " DONE!"
  echo
  cat LAT*/store*.$namemethoda.dat | sort -n > storeall.$namemethoda.dat
done
