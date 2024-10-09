#!/bin/bash

#export OMP_NUM_THREADS=64

start_time="$(date -u +%s)"

f1='_exited_particles.bin'
f2='_particle_restart.bin'


# --- spinup ---
# Using James' spinup runs
#dspn='ecoslim.spinup'   # directory name
#mkdir $dspn
#cd $dspn
#cp ../slimin.txt.spinup ./slimin.txt
#cp ../eastRiver_neon_spinUp_55-57_particle_restart.bin ./spinup$f2

#echo 'Working on spinup...'
#cp ../ER_dem.pfb .
#EcoSLIM.exe
#cd ..


# --- wy2015 ---
dwy15='ecoslim.wy2015'
mkdir $dwy15
cd $dwy15
cp ../slimin.txt.wy2015 ./slimin.txt

echo 'Working on wy2015...'
cp ../ER_dem.pfb .
#cp ../$dspn/spinup$f2 wy2015$f2
cp ../eastRiver_neon_spinUp_55-57_particle_restart.bin ./wy2015$f2
EcoSLIM.exe
cd ..


# --- wy2016 ---
dwy16='ecoslim.wy2016'
mkdir $dwy16
cd $dwy16
cp ../slimin.txt.wy2016 ./slimin.txt

echo 'Working on wy2016...'
cp ../ER_dem.pfb .
cp ../$dwy15/wy2015$f2 wy2016$f2
EcoSLIM.exe
cd ..


# --- wy2017 ---
dwy17='ecoslim.wy2017'
mkdir $dwy17
cd $dwy17
cp ../slimin.txt.wy2017 ./slimin.txt

echo 'Working on wy2017...'
cp ../ER_dem.pfb .
cp ../$dwy16/wy2016$f2 wy2017$f2
EcoSLIM.exe
cd ..


# --- wy2018 ---
dwy18='ecoslim.wy2018'
mkdir $dwy18
cd $dwy18
cp ../slimin.txt.wy2018 ./slimin.txt

echo 'Working on wy2018...'
cp ../ER_dem.pfb .
cp ../$dwy17/wy2017$f2 wy2018$f2
EcoSLIM.exe
cd ..


# --- wy2019 ---
dwy19='ecoslim.wy2019'
mkdir $dwy19
cd $dwy19
cp ../slimin.txt.wy2019 ./slimin.txt

echo 'Working on wy2019...'
cp ../ER_dem.pfb .
cp ../$dwy18/wy2018$f2 wy2019$f2
EcoSLIM.exe
cd ..


# --- wy2020 ---
dwy20='ecoslim.wy2020'
mkdir $dwy20
cd $dwy20
cp ../slimin.txt.wy2020 ./slimin.txt

echo 'Working on wy2020...'
cp ../ER_dem.pfb .
cp ../$dwy19/wy2019$f2 wy2020$f2
EcoSLIM.exe
cd ..


# --- wy2021 ---
dwy21='ecoslim.wy2021'
mkdir $dwy21
cd $dwy21
cp ../slimin.txt.wy2021 ./slimin.txt

echo 'Working on wy2021...'
cp ../ER_dem.pfb .
cp ../$dwy20/wy2020$f2 wy2021$f2
EcoSLIM.exe
cd ..


# --- wy2022 ---
dwy22='ecoslim.wy2022'
mkdir $dwy22
cd $dwy22
cp ../slimin.txt.wy2022 ./slimin.txt

echo 'Working on wy2022...'
cp ../ER_dem.pfb .
cp ../$dwy21/wy2021$f2 wy2022$f2
EcoSLIM.exe
cd ..

end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Runtime: $elapsed seconds"

