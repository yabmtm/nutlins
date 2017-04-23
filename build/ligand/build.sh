#!/bin/bash

# This script will take a mol2 file for each of the 7 nutlins
# and generate GROMACS structure and topology files using the
# partial charges calculated via AM1BCC semi-empirical methods.

for NAME in 1MQ 20Q 20U I18 I31 NUT YIN; do
    cd $NAME
    rm -rf A* dlm* em.mdp leap.log md.mdp  *crd *edited* *off *GMX* *prmtop sqm* build.tleap *nut; cd ..; done
ls *; sleep 2

for NAME in 1MQ 20Q 20U I18 I31 NUT YIN; do
    cd $NAME
    python ~/research/repos/voelzlab/custom_topologies/scripts/Allsteps_VAV.py $NAME.mol2 $NAME CappingAtoms.dat B G 0 && sleep 5 # runs allsteps

## Build tleap script.
    echo source leaprc.gaff >> build.tleap
    echo gaff = loadamberparams gaff.dat >> build.tleap
    echo loadoff $NAME.off >> build.tleap
    echo $NAME = sequence { $NAME } >> build.tleap
    echo saveAmberParm $NAME $NAME.prmtop $NAME.crd >> build.tleap
    echo quit >> build.tleap

## Run the tleap script
    tleap -f build.tleap

## Run Acpype
    python ~/research/scripts/acpype.py -p $NAME.prmtop -x $NAME.crd -b $NAME

    cd ..

done

# organize this output into simulation-ready .top/.gro/.itp's
for i in 1MQ 20U 20Q I18 I31 NUT YIN; do
    mkdir $i/mdm2_$i
    cd ${i}/mdm2_${i}
    cp ../*GMX* .
    mv ${i}_GMX.top ${i}.itp
    mv ${i}_GMX.gro ${i}.gro
    cp ../../../protein/$i/conf.gro mdm2.gro
    ls; sleep 5

# writes an overall GROMACS topol.top that will include ligand/protein .itp's
    printf "; Include forcefield topology\n#include \"amber99sb-ildn.ff/forcefield.itp\"" > topol.top
    cat n18_GMX.top | grep -A 50 atomtypes | sed -e '/^$/,$d' >> topol.top # assuming <=50 atomtypes
    printf "\n\n; Include structural topologies\n#include \"mdm2.itp\"\n#include \"$i.itp\"" >> topol.top
    printf "\n\n; Include solvation topology\n#include \"amber99sb-ildn.ff/tip3p.itp\"" >> topol.top
    printf "\n\n#ifdef POSRES_WATER\n[ position_restraints ]" >> topol.top
    printf ";  i funct       fcx        fcy        fcz" >> topol.top
    printf "\n   1    1       1000       1000       1000\n#endif" >> topol.top
    printf "\n\n[ system ]\n$i_mdm2" >> topol.top
    printf "\n\n; Include ion topology\n#include \"amber99sb-ildn.ff/ions.itp\"" >> topol.top
    printf "\n\n[ molecules ]\nProtein_Chain_A    1\n$i    21" >> topol.top

# building mdm2.itp
    TOP_LEN=$(cat ~/research/nutlin/protein/$i/topol.top | wc -l)
    BOF=$(cat ~/research/nutlin/protein/$i/topol.top | grep -n "moleculetype" | sed "s/:.*$//")
    EOF=$(cat ~/research/nutlin/protein/20Q/topol.top | grep -n "Include Position restraint" | sed "s/:.*$//")
    head -n $((TOP_LEN-BOF)) ~/research/nutlin/protein/$i/topol.top | tail -n $((TOP_LEN-EOF)) ~/research/nutlin/protein/$i/topol.top > mdm2.itp

# putting tons of ligands into the structure (probably don't want to do this)
#    echo "Adding ligands for $i:"
#    gmx editconf -f mdm2.gro -o box.gro -d 1.0 -bt cubic
#    gmx insert-molecules -f box.gro -ci -try 100 $i.gro -nmol 21 -o ${i}_mdm2.gro | grep Added

    cd ../..; done

