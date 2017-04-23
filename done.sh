
# extract out each ligand into a separate file in a separate directory
for i in {111..130}; do count=$(($i-110)); cat npt.gro | grep "${i}1MQ" >> RUN$count/1MQ; done
for i in {106..125}; do count=$(($i-105)); cat npt.gro | grep "${i}20Q" >> RUN$count/20Q; done
for i in {106..125}; do count=$(($i-105)); cat npt.gro | grep "${i}20U" >> RUN$count/20U; done
for i in {106..125}; do count=$(($i-105)); cat npt.gro | grep "${i}I18" >> RUN$count/I18; done
for i in {106..125}; do count=$(($i-105)); cat npt.gro | grep "${i}I31" >> RUN$count/I31; done
for i in {106..125}; do count=$(($i-105)); cat npt.gro | grep "${i}NUT" >> RUN$count/NUT; done
for i in {112..131}; do count=$(($i-111)); cat npt.gro | grep "${i}YIN" >> RUN$count/YIN; done

# delete all ligands and solvent from previously equilibrated system and copy the remaining system over into each run, concatenating the single ligand onto it
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; sed -i "/$i/d" npt.gro; sed -i "/SOL/d" npt.gro; for j in {1..20}; do cp npt.gro RUN${j}/conf.gro; cat RUN${j}/$i >> RUN${j}/conf.gro; done; cd ..; done

# delete all old ions from the system as well, being careful not to delete the CL that naturally exist on the ligand
for i in 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do sed -i "/NA/d" $j/conf.gro; sed -i "/[0-9]CL/d" $j/conf.gro; done; cd ..; done

# move the size coordinates of the box down from the middle of the file to the end
for i in RUN*; do a=$(cat $i/conf.gro | grep "^\ *[0-9]\."); sed -i "/$a/d" $i/conf.gro; echo $a >> $i/conf.gro; done

# correct the number of atoms at the top of each conf.gro to reflect the absence of old water/ions/ligands: note 1MQ and YIN had to be fixed because they didn't have titles in their .gro files
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do lines=$(cat $j/conf.gro | wc -l); atoms=$(($lines-3)); sed -i "2s/.*/$atoms/" $j/conf.gro; done; cd ..; done

# manually edit a master topology file to include 1 ligand, no solvent or ions, and then copy it and other necessary files into each run directory
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do yes | cp * $j; done; cd ..; done

# solvate all simulations
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; gmx solvate -cp conf.gro -cs spc216.gro -o solv.gro -p topol.top; cd ..; done; cd ..; done
# grompp an ions.tpr
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr; cd ..; done; cd ..; done
# add ions to 100mmol and neutralize
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; echo SOL | gmx genion -s ions.tpr -o solv.gro -p topol.top -pname NA -nname CL -neutral -conc 0.1; cd ..; done; cd ..; done
# grompp for minimization
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; gmx grompp -f minim.mdp -c solv.gro -p topol.top -o em.tpr; cd ..; done; cd ..; done
# submit minimization
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; qsub gmx_minim.sh; cd ..; done; cd ..; done
# grompp for npt equilibration
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; gmx grompp -f npt.mdp -c solv.gro -p topol.top -o npt.tpr; cd ..; done; cd ..; done
# submit equilibration
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; qsub gmx_npt.sh; cd ..; done; cd ..; done

# gather all the output structures in a separate production directory
mkdir production; cd production
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do mkdir production/$i; cd $i; for j in RUN*; do mkdir ../production/$i/$j; cd $j; cp npt.gro ../../production/$i/$j/conf.gro

# populate with topologies
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cp *.itp ../production/$i/$j; cp *top ../production/$i/$j

# populate with index files, adding groups for harmonic restraints (these will automatically be groups 24 and 25)
cd 1MQ; for i in RUN*; do cd $i; echo -e "a1564\na678\nq" | gmx make_ndx -f npt.gro -o index.ndx; cd ..; done cd ..
cd 20Q; for i in RUN*; do cd $i; echo -e "a1420\na567\nq" | gmx make_ndx -f npt.gro -o index.ndx; cd ..; done cd ..
cd 20U; for i in RUN*; do cd $i; echo -e "a1403\na550\nq" | gmx make_ndx -f npt.gro -o index.ndx; cd ..; done cd ..
cd I18; for i in RUN*; do cd $i; echo -e "a1421\na567\nq" | gmx make_ndx -f npt.gro -o index.ndx; cd ..; done cd ..
cd I31; for i in RUN*; do cd $i; echo -e "a1422\na567\nq" | gmx make_ndx -f npt.gro -o index.ndx; cd ..; done cd ..
cd NUT; for i in RUN*; do cd $i; echo -e "a1446\na567\nq" | gmx make_ndx -f npt.gro -o index.ndx; cd ..; done cd ..
cd YIN; for i in RUN*; do cd $i; echo -e "a1612\na696\nq" | gmx make_ndx -f npt.gro -o index.ndx; cd ..; done cd ..

# change .mdp pull group names to reflect the names instead of group numbers
for i in 1MQ; do for j in {1..20}; do sed -i "s/pull_group1_name        = 24/pull_group1_name        = a_1564/" $i/RUN$j/production.mdp; done; done
for i in 1MQ; do for j in {1..20}; do sed -i "s/pull_group2_name        = 25/pull_group2_name        = a_678/" $i/RUN$j/production.mdp; done; done
for i in 20Q; do for j in {1..20}; do sed -i "s/pull_group1_name        = 24/pull_group1_name        = a_1420/" $i/RUN$j/production.mdp; done; done
for i in 20Q; do for j in {1..20}; do sed -i "s/pull_group2_name        = 25/pull_group2_name        = a_567/" $i/RUN$j/production.mdp; done; done
for i in 20U; do for j in {1..20}; do sed -i "s/pull_group1_name        = 24/pull_group1_name        = a_1403/" $i/RUN$j/production.mdp; done; done
for i in 20U; do for j in {1..20}; do sed -i "s/pull_group2_name        = 25/pull_group2_name        = a_550/" $i/RUN$j/production.mdp; done; done
for i in I18; do for j in {1..20}; do sed -i "s/pull_group1_name        = 24/pull_group1_name        = a_1421/" $i/RUN$j/production.mdp; done; done
for i in I18; do for j in {1..20}; do sed -i "s/pull_group2_name        = 25/pull_group2_name        = a_567/" $i/RUN$j/production.mdp; done; done
for i in I31; do for j in {1..20}; do sed -i "s/pull_group1_name        = 24/pull_group1_name        = a_1422/" $i/RUN$j/production.mdp; done; done
for i in I31; do for j in {1..20}; do sed -i "s/pull_group2_name        = 25/pull_group2_name        = a_567/" $i/RUN$j/production.mdp; done; done
for i in NUT; do for j in {1..20}; do sed -i "s/pull_group1_name        = 24/pull_group1_name        = a_1446/" $i/RUN$j/production.mdp; done; done
for i in NUT; do for j in {1..20}; do sed -i "s/pull_group2_name        = 25/pull_group2_name        = a_567/" $i/RUN$j/production.mdp; done; done
for i in YIN; do for j in {1..20}; do sed -i "s/pull_group1_name        = 24/pull_group1_name        = a_1612/" $i/RUN$j/production.mdp; done; done
for i in YIN; do for j in {1..20}; do sed -i "s/pull_group2_name        = 25/pull_group2_name        = a_696/" $i/RUN$j/production.mdp; done; done

# add FEP code to production mdp and change pulling geometry for large distances
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do
    for j in {1..20}; do
        echo "init-lambda-state        = 1" >> $i/RUN$j/production.mdp
        echo "fep-lambdas              = 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00" >> $i/RUN$j/production.mdp
        echo "calc-lambda-neighbors = -1" >> $i/RUN$j/production.mdp
        echo "separate_dhdl_file = yes" >> $i/RUN$j/production.mdp
        sed -i "s/pull-geometry    = distance/pull-geometry    = direction-periodic/" $i/production.mdp; done; done


# grompp a production run
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; gmx grompp -f production.mdp -c conf.gro -p topol.top -n index.ndx -o prod.tpr; cd ..; done; cd ..; done

# submit production
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do cd $i; for j in RUN*; do cd $j; qsub prod.sh; cd ..; done; cd ..; done

###############################

# let's get some distances out
for i in 1MQ 20Q 20U I18 I31 NUT YIN; do
    j=$(cat $i/production.mdp | grep a_[0-9[0-9] | sed "s/^.*= a_//" | head -n 1)
    x1=$(cat $i/npt.gro | " $j " | awk '{print $4}'
    y1=$(cat $i/npt.gro | " $j " | awk '{print $5}'
    z1=$(cat $i/npt.gro | " $j " | awk '{print $6}'
    j=$(cat $i/production.mdp | grep a_[0-9[0-9] | sed "s/^.*= a_//" | tail -n 1)
    x2=$(cat $i/npt.gro | " $j " | awk '{print $4}'
    y2=$(cat $i/npt.gro | " $j " | awk '{print $5}'
    z2=$(cat $i/npt.gro | " $j " | awk '{print $6}'
    echo "$i $(echo "sqrt((${x2}-(${x1}))^2+(${y2}-(${y1}))^2+(${z2}-(${z1}))^2)" | bc -l)"
