#!/usr/bin/env bash

# Typical Usage:
#
# README_restart_from_lammpstrj_dump_files.sh traj.lammpstrj system.data \
#  > system_most_recent.data
#
# Use this script to create a new LAMMPS "DATA" file containing the
# coordinates from the most recent snapshot in the trajectory ("dump") file
# (indicated by the 1st argument to this script).
# The rest of the information will be extracted from the "ORIG_DATA_FILE"
# (2nd argument, typically "system.data")
#
# If your simulation gets interrupted, you would typically use LAMMPS'
# "read_restart" command to read the most restart file from your simulation
# to continue running the simulation from where it left off.
# This is useful if you forgot to use the LAMMPS "restart" command to generate 
# these restart files, (or if LAMMPS is refusing to read them for some reason).
# Once you use this script to create a new "data" file, 
# you can then use LAMMPS' "read_data" command to read this file to continue
# your simulation where it left off.  (...approximately speaking. The atom types
# and the bond topology will be copied from the "ORIG_DATA_FILE" script.
# Only the atom coordinates will be updated.  Hopefully this is good enough.)
#
# WARNING: This script assumes your DATA files use "atom_style full" (which is
# the default in LAMMPS)  See:  http://lammps.sandia.gov/doc/atom_style.html

#TRAJ_FILE="traj.lammpstrj"
TRAJ_FILE=$1
#ORIG_DATA_FILE="system.data"
ORIG_DATA_FILE=$2


Natoms=`grep atoms "$ORIG_DATA_FILE" | awk '{print $1}'`
tail -n $((Natoms+9)) < "$TRAJ_FILE" \
   | dump2data.py -last "$ORIG_DATA_FILE" \
   | extract_lammps_data.py Atoms \
   | sort -g \
   | cut -d' ' -f5- \
   > TMP_CRDS

extract_lammps_data.py Header < "$ORIG_DATA_FILE"

echo "Atoms"
echo ""

extract_lammps_data.py Atoms < "$ORIG_DATA_FILE" | sort -g | cut -d' ' -f1,2,3,4 > TMP_ATOMTYPES

paste -d' ' TMP_ATOMTYPES TMP_CRDS
echo ""

extract_lammps_data.py -n Atoms -n Header < "$ORIG_DATA_FILE"

rm -f TMP_ATOMTYPES TMP_CRDS
