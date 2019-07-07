#!/usr/bin/python2.7

import os,sys
from fractions import Fraction
from numpy import *

### PARAMETERS ###
sigma0  = 1.00
epsilon = 1.00

sigma_hh = 0.95 * sigma0
sigma_ht = 0.95 * sigma0
sigma_tt = 1.00 * sigma0

r_min    = 0.000001
r_minsq = r_min**2
Nentries = 128
##################

### INPUTS ###
if len(sys.argv) == 2:
  w_cut_str = sys.argv[1]
  w_cut = float(w_cut_str)
else:
  w_cut = 1.6
  w_cut_str = str(w_cut)
# 1.6 seems to be 'good' for vesicles, bilayers 1.4
##############

def WCA_energy(b, r):
# Calculate WCA energy
  E_pot = 0
  val1 = math.pow((b / r), 12)
  val2 = -math.pow((b / r), 6)
  E_pot = 4 * epsilon * (val1 + val2 + 0.25)
  return E_pot

def WCA_forces(b, r):
# Calculate WCA forces
  Force = 0
  val1 = 24 * math.pow(b, 6) / math.pow(r, 7)
  val2 = -48 * math.pow(b, 12) / math.pow(r, 13)
  Force = -(val1 + val2)
  return Force

def Tail_energy(b, r, r_cut):
# Calculate extra Tail energy
  E_pot = 0
  if (r < r_cut):
    E_pot = -1 * epsilon
  else:
    val1 = math.cos((math.pi * (r - r_cut)) / (2 * w_cut))
    E_pot = -1 * epsilon * math.pow(val1, 2)
  return E_pot

def Tail_forces(b, r, r_cut):
  Force = 0
  if (r < r_cut):
    Force = 0;
  else:
    val1 = math.sin((math.pi * (r - r_cut)) / w_cut)
    Force = -math.pi * val1 / (2 * w_cut)
  return Force


##############
ofile = sys.stdout
#ofile = open('tabulated_potential.dat', 'w')

tot_potential_hh = zeros((Nentries, 4))
tot_potential_ht = zeros((Nentries, 4))
tot_potential_tt = zeros((Nentries, 4))

# Setup up formatting & distances in all arrays
for i in range(0, Nentries):
  tot_potential_hh[i,0] = i+1
  tot_potential_ht[i,0] = i+1
  tot_potential_tt[i,0] = i+1
# The r^2 values (not r values) of the table are evenly distributed.
# (Hopefully) this makes the table lookup faster in LAMMPS.
r_max = sigma_hh * math.pow(2,1.0/6)
r_maxsq = r_max**2
for i in range(0, Nentries):
  rsq = r_minsq + i*(r_maxsq-r_minsq)/(Nentries-1)
  r = sqrt(rsq)
  tot_potential_hh[i,1] = r
r_max = sigma_ht * math.pow(2,1.0/6)
r_maxsq = r_max**2
for i in range(0, Nentries):
  rsq = r_minsq + i*(r_maxsq-r_minsq)/(Nentries-1)
  r = sqrt(rsq)
  tot_potential_ht[i,1] = r
r_max = sigma_tt * math.pow(2,1.0/6) + w_cut
r_maxsq = r_max**2
for i in range(0, Nentries):
  rsq = r_minsq + i*(r_maxsq-r_minsq)/(Nentries-1)
  r = sqrt(rsq)
  tot_potential_tt[i,1] = r



ofile.write("# Tabulated potential for Cooke 3-bead lipid model, Wc = "+w_cut_str+"\n\n")

### Calcaulte first potential, H-H
ofile.write("HEAD_HEAD\n")
r_cut = 2**Fraction('1/6') * sigma_hh
r_max = tot_potential_hh[-1,1]
ofile.write("N "+str(Nentries)+" RSQ "+str(r_min)+" "+str(r_max)+" FPRIME 0.0 0.0\n\n")

for i in range(0, Nentries):
  r = tot_potential_hh[i,1]
  tot_potential_hh[i,2] = WCA_energy(sigma_hh, r)
  tot_potential_hh[i,3] = WCA_forces(sigma_hh, r)

for i in range(0, Nentries):
  ofile.write(str(i+1)+" "+
              str(tot_potential_hh[i,1])+" "+
              str(tot_potential_hh[i,2])+" "+
              str(tot_potential_hh[i,3])+"\n")
ofile.write("\n")



### Calcaulte second potential, H-T
ofile.write("HEAD_TAIL\n")
r_cut = 2**Fraction('1/6') * sigma_ht
r_max = tot_potential_ht[-1,1]
ofile.write("N "+str(Nentries)+" RSQ "+str(r_min)+" "+str(r_max)+" FPRIME 0.0 0.0\n\n")

for i in range(0, Nentries):
  r = tot_potential_ht[i,1]
  tot_potential_ht[i,2] = WCA_energy(sigma_ht, r)
  tot_potential_ht[i,3] = WCA_forces(sigma_ht, r)

for i in range(0, Nentries):
  ofile.write(str(i+1)+" "+
              str(tot_potential_ht[i,1])+" "+
              str(tot_potential_ht[i,2])+" "+
              str(tot_potential_ht[i,3])+"\n")
ofile.write("\n")



### Calcaulte third potential, T-T
# Also include extra tail-tail attraction term
ofile.write("TAIL_TAIL_Wc="+w_cut_str+"\n")
r_cut = 2**Fraction('1/6') * sigma_tt
r_max = tot_potential_tt[-1,1]
ofile.write("N "+str(Nentries)+" RSQ "+str(r_min)+" "+str(r_max)+" FPRIME 0.0 0.0\n\n")

for i in range(0, Nentries):
  r = tot_potential_tt[i,1]
  if r < r_cut:
    tot_potential_tt[i,2] = WCA_energy(sigma_tt, r)
    tot_potential_tt[i,3] = WCA_forces(sigma_tt, r)
  else:
    tot_potential_tt[i,2] = Tail_energy(sigma_tt, r, r_cut)
    tot_potential_tt[i,3] = Tail_forces(sigma_tt, r, r_cut)

for i in range(0, Nentries):
  ofile.write(str(i+1)+" "+
              str(tot_potential_tt[i,1])+" "+
              str(tot_potential_tt[i,2])+" "+
              str(tot_potential_tt[i,3])+"\n")
ofile.write("\n")


sys.exit()
