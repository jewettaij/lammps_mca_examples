/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Andrew Jewett (jewett.aij@gmail.com)
------------------------------------------------------------------------- */

#include <cmath>
#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include "fix_bond_modify.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
//#define NDEBUG //(disable assert())
using namespace std;

using namespace LAMMPS_NS;
using namespace FixConst;

#define DEBUG_FIX_BOND_MODIFY



#define DELTA 16
#ifndef INFINITY         //defined by ISO C99
#define INFINITY 1.0+37;
#define isfinite(x) ((x != INFINITY) && (x != -INFINITY))
#endif
static const int MSGLENGTH=2048;

/* ---------------------------------------------------------------------- */

FixBondModify::FixBondModify(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  partner(NULL), finalpartner(NULL), distsq(NULL), probability(NULL), broken(NULL), copy(NULL), random(NULL)
{
  char msg[MSGLENGTH];
  if (narg < 5) error->all(FLERR,"Illegal fix bond/modify command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix bond/modify command");
  nevery_delay = 0;

  conserve_ke_flag = 1;
  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;

  bond_rminsq = -1.0;
  bond_rmaxsq = -1.0;
  prob_transition = 1.0;
  seed = 1; //(default seed value)

  // The following variables are initialized with impossible values by default
  btypeold_lo = 0;
  btypeold_hi = -1;
  btypenew = -1;
  atype1old_lo = 0;
  atype1old_hi = -1;
  atype1new = -1;
  atype2old_lo = 0;
  atype2old_hi = -1;
  atype2new = -1;
  qold1_lo = 0.0;
  qold1_hi = -1.0;
  qold2_lo = 0.0;
  qold2_hi = -1.0;
  qnew1 = INFINITY;
  qnew2 = INFINITY;
  ignore_tags_flag = 1;

  check_bonded_atoms = true;
  change_bond_properties = false;
  break_bond = false;
  prioritize_long_bonds = false;
  #ifdef ALLOW_BOND_MODIFY_SINGLE_ATOMS
  // The following variables are initialized with impossible values by default
  atypeold_lo = 0;
  atypeold_hi = -1;
  atypenew = 0;
  qold_lo = 0.0;
  qold_hi = -1.0;
  qnew = INFINITY;
  check_individual_atoms = false;
  #endif
  //iterate_until_stationary = false;


  bool read_new_params = false;

  int iarg = 4;
  while (iarg < narg) {

    // Debugging:
    //fprintf(stderr, "iarg=%d, narg=%d, arg[iarg]=\"%s\"\n", iarg, narg, arg[iarg]);
    //if (strcmp(arg[iarg],"iterate") == 0) {
    //  if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
    //  iterate_until_stationary = (strcmp("yes", arg[iarg+1]) == 0);
    //  iarg += 2;
    //}

    if (strcmp(arg[iarg],"->") == 0) {
      read_new_params = true;
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"then") == 0) {
      read_new_params = true;
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"if") == 0) {
      read_new_params = false;
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"and") == 0) {
      iarg += 1;
    }
    else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
      nevery_delay = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"prob") == 0) {
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
      prob_transition = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (prob_transition < 0.0 || prob_transition > 1.0)
        error->all(FLERR,"Illegal fix bond/modify command");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"seed") == 0) {
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
      seed = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (seed <= 0) error->all(FLERR,"Illegal fix bond/modify command");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"distance") == 0) {
      if (iarg+2 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
      if (strcmp(arg[iarg+1],">=") == 0) {
        bond_rminsq = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        bond_rminsq *= bond_rminsq;
        prioritize_long_bonds = true;
      }
      else if (strcmp(arg[iarg+1],"<=") == 0) {
        bond_rmaxsq = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        bond_rmaxsq *= bond_rmaxsq;
        prioritize_long_bonds = false;
      }
      else error->all(FLERR,"Illegal fix bond/modify command. (\"distance\" should be followed by >= or <=)");
      iarg += 3;
    }
    //else if (strcmp(arg[iarg],"bondrmin") == 0) {
    //  check_bonded_atoms = true;
    //  if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
    //  bond_rminsq = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    //  if (bond_rminsq < 0.0) error->all(FLERR,"Illegal fix bond/modify command");
    //  bond_rminsq *= bond_rminsq;
    //  iarg += 2;
    //}
    //else if (strcmp(arg[iarg],"bondrmax") == 0) {
    //  check_bonded_atoms = true;
    //  if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
    //  bond_rmaxsq = utils::numeric(FLERR,arg[iarg+1],false,lmp);
    //  if (bond_rmaxsq < 0.0) error->all(FLERR,"Illegal fix bond/modify command");
    //  bond_rmaxsq *= bond_rmaxsq;
    //  iarg += 2;
    //}
    else if (strcmp(arg[iarg],"bond") == 0) {
      check_bonded_atoms = true;
      change_bond_properties = true;
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
      if (read_new_params) {
        if ((strcmp(arg[iarg+1],"break") == 0) ||
            (strcmp(arg[iarg+1],"BREAK") == 0)) {
          break_bond = true;
          prioritize_long_bonds = true;
          //btypenew = BONDBREAK;
        }
        else
          btypenew = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      }
      else
        utils::bounds(FLERR,arg[iarg+1], 1,
                      atom->nbondtypes, btypeold_lo, btypeold_hi, error);
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"atoms") == 0) {
      check_bonded_atoms = true;
      if (iarg+2 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
      if (read_new_params) {
        if ((strcmp(arg[iarg+1],"SAME")==0)||(strcmp(arg[iarg+1],"same")==0)||
            (strcmp(arg[iarg+1],"NULL")==0)||(strcmp(arg[iarg+1],"*")==0))
          atype1new = -1; //(disables)
        else
          atype1new = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
        if ((strcmp(arg[iarg+2],"SAME")==0)||(strcmp(arg[iarg+2],"same")==0)||
            (strcmp(arg[iarg+2],"NULL")==0)||(strcmp(arg[iarg+2],"*")==0))
          atype2new = -1; //(disables)
        else
          atype2new = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      }
      else {
        utils::bounds(FLERR,arg[iarg+1], 1,
                      atom->ntypes, atype1old_lo, atype1old_hi, error);
        utils::bounds(FLERR,arg[iarg+2], 1,
                      atom->ntypes, atype2old_lo, atype2old_hi, error);
      }
      iarg += 3;
    }
    else if (strcmp(arg[iarg],"charges") == 0) {
      check_bonded_atoms = true;
      if (read_new_params) {
        if (iarg+2 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
        qnew1 = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        qnew2 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      }
      else {
        if (iarg+4 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
        qold1_lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        qold1_hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        qold2_lo = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        qold2_hi = utils::numeric(FLERR,arg[iarg+4],false,lmp);
        iarg += 5;
      }
    }
    else if (strcmp(arg[iarg],"ordered") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/modify command");
      if (strcmp(arg[iarg+1],"no") == 0) ignore_tags_flag = 1;
      else if (strcmp(arg[iarg+1],"yes") == 0) ignore_tags_flag = 0;
      else error->all(FLERR,"Illegal fix bond/modify command");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"ke") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix bond/modify command");
      if (strcmp(arg[iarg+1],"no") == 0) conserve_ke_flag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) conserve_ke_flag = 1;
      else error->all(FLERR,"Illegal fix bond/modify command");
      iarg += 2;
    }


    #ifdef ALLOW_BOND_MODIFY_SINGLE_ATOMS
    else if (strcmp(arg[iarg],"atom") == 0) {
      check_individual_atoms = true;
      check_bonded_atoms = false;
      if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
      if (read_new_params)
        atypenew = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      else
        utils::bounds(FLERR,arg[iarg+1], 1,
                      atom->ntypes, atypeold_lo, atypeold_hi, error);
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"charge") == 0) {
      check_individual_atoms = true;
      check_bonded_atoms = false;
      if (read_new_params) {
        if (iarg+1 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
        qnew = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        iarg += 2;
      }
      else {
        if (iarg+2 >= narg) error->all(FLERR,"Illegal fix bond/modify command");
        qold_lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        qold_hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        iarg += 3;
      }
    }
    #endif //#ifdef ALLOW_BOND_MODIFY_SINGLE_ATOMS

    else {
      snprintf(msg, MSGLENGTH,
               "Illegal fix bond/modify command. Unrecognized keyword: \"%s\"\n",
               arg[iarg]);
      error->all(FLERR, msg);
    }
  }


  // error check

  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use fix bond/modify with non-molecular systems");

  #ifdef ALLOW_BOND_MODIFY_SINGLE_ATOMS
  if (check_individual_atoms && check_bonded_atoms)
      if (iarg+2 >= narg) error->all(FLERR,"Illegal fix bond/modify command.\n"
                                     "fix bond/modify must either be used to change the properties of\n"
                                     "individual atoms OR bonded atom pairs, but not both.\n"
                                     "This typically happens if you typed \"atom\" or \"charge\"\n"
                                     "when you meant to use \"atoms\" or \"charges\"\n");
  #endif

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + me);

  // set comm sizes needed by this fix
  // forward is big due to comm of broken bonds and 1-2 neighbors

  //comm_forward = MAX(2,2+atom->maxspecial);
  comm_forward = MAX(4,4+atom->maxspecial);
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;

  maxbreak = 0;

  // copy = special list for one atom
  // size = ms^2 + ms is sufficient
  // b/c in rebuild_special_one() neighs of all 1-2s are added,
  //   then a dedup(), then neighs of all 1-3s are added, then final dedup()
  // this means intermediate size cannot exceed ms^2 + ms

  int maxspecial = atom->maxspecial;
  copy = new tagint[maxspecial*maxspecial + maxspecial];

  // zero out stats
  breakcount = 0;
  breakcounttotal = 0;     //<--DO WE NEED THIS? -A 2819-7-08

} //FixBondModify::FixBondModify()

/* ---------------------------------------------------------------------- */

FixBondModify::~FixBondModify()
{
  delete random;

  // delete locally stored arrays

  memory->destroy(partner);
  memory->destroy(finalpartner);
  memory->destroy(distsq);
  memory->destroy(broken);
  delete [] copy;
}

/* ---------------------------------------------------------------------- */

int FixBondModify::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondModify::init()
{
  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // enable angle/dihedral/improper breaking if any defined

  if (atom->nangles) angleflag = 1;
  else angleflag = 0;
  if (atom->ndihedrals) dihedralflag = 1;
  else dihedralflag = 0;
  if (atom->nimpropers) improperflag = 1;
  else improperflag = 0;

  if (force->improper) {
    if (force->improper_match("class2") || force->improper_match("ring"))
      error->all(FLERR,"Cannot yet use fix bond/modify with this "
                 "improper style");
  }

  lastcheck = -1;

  // DEBUG
  //print_bb();
}


/* ---------------------------------------------------------------------- */

#ifdef DEBUG_FIX_BOND_MODIFY
// detect MPI-related issues
void FixBondModify::
CheckAtomConsistency()
{
  tagint *atom_counts = NULL; // sum of number of atoms "known"
                              //by each processor (nlocal+nghost)
  tagint *atom_displacements = NULL; // cumulative number of atoms known
  tagint natoms_known_all = 0; //sum of the number all atoms "known" 
                               //by each proc (nlocal+nghost).
                               //This will exceed the total number of 
                               //atoms in the simulation
                               //because some atoms are known by
                               //multiple procs. (ie. ghost atoms)
  tagint *tags_known_unsorted = NULL; //the atom-IDs for all the atoms known by all
                          //procs (including ghost atoms seen by multiple procs)
                          //This array has length "natoms_known_all"
  int *types_all = NULL;
  int *types_known_unsorted = NULL;
  double *q_all = NULL;
  double *q_known_unsorted = NULL;
  if (me == 0) {
    atom_counts = new int[nprocs];
    atom_displacements = new tagint[nprocs + 1];
  }
  tagint n_known = atom->nlocal + atom->nghost;
  MPI_Gather(&n_known,1,MPI_INT,atom_counts,1,MPI_INT,0,world);
  if (me == 0) {
    atom_displacements[0] = 0;
    for (int i = 1; i <= nprocs; i++)
      atom_displacements[i] = atom_displacements[i-1] + atom_counts[i-1];
    natoms_known_all = atom_displacements[nprocs];
    tags_known_unsorted = new tagint[natoms_known_all];
    types_all = new int[atom->natoms + 1];
    types_known_unsorted = new int[natoms_known_all];
    if (atom->q_flag) {
      q_all = new double[atom->natoms + 1];
      q_known_unsorted = new double[natoms_known_all];
    }
    for (tagint i=0; i < natoms_known_all; i++) {
      tags_known_unsorted[i] = -1;         //(an impossible value)
      types_known_unsorted[i] = -1;        //(an impossible value)
      if (atom->q_flag)
        q_known_unsorted[i] = -111111.111; //(an unlikely and suspicious value)
    }
    for (tagint i=0; i <= atom->natoms; i++) {
      types_all[i] = -2;                   //(an impossible value)
      if (atom->q_flag)
        q_all[i] = -222222.222;            //(an unlikely and suspicious value)
    }
  } //if (me == 0)

  MPI_Gatherv(atom->tag, n_known, MPI_LMP_TAGINT,
              tags_known_unsorted, atom_counts, atom_displacements, MPI_LMP_TAGINT,
              0, world);     
  MPI_Gatherv(atom->type, n_known, MPI_INT,
              types_known_unsorted, atom_counts, atom_displacements, MPI_INT,
              0, world);
  if (atom->q_flag) {
    assert(q_known_unsorted);
    MPI_Gatherv(atom->q, n_known, MPI_DOUBLE,
                q_known_unsorted, atom_counts, atom_displacements, MPI_DOUBLE,
                0, world);
  }
  if (me == 0) {
    //for (tagint i=0; i < natoms_known_all; i++) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      for (int i_localatom=0; i_localatom < atom_counts[iproc]; i_localatom++) {
        tagint i = atom_displacements[iproc] + i_localatom;  // atom ID for proc iproc
        tagint tag = tags_known_unsorted[i];                 // corresponding global atom ID

        assert((1 <= tag) && (tag <= atom->natoms));

        if (types_all[tag] >= 0)
        {
          //if already specified, check for consistency
          stringstream err_ss;
          if (types_all[tag] != types_known_unsorted[i]) {
            err_ss << "(MPI) Inconsistent atom types on different processors\n"
                   << "Atom-ID: " << tag << ",  Types: " << types_all[tag]
                   << ", " << types_known_unsorted[i] << "\n";
            error->all(FLERR, err_ss.str().c_str());
          }
          if (atom->q_flag && (q_all[tag] != q_known_unsorted[i])) {
            err_ss <<"(MPI) Inconsistent atom charges on different processors\n"
                   << "Atom-ID: " << tag << ",  Charges: " << q_all[tag]
                   << ", " << q_known_unsorted[i] << "\n";
            error->all(FLERR, err_ss.str().c_str());
          }
        }
        else {
          types_all[tag] = types_known_unsorted[i];
          if (atom->q_flag) {
            assert(q_all);
            q_all[tag] = q_known_unsorted[i];
          }
        }
      } // (loop over local atoms, i_localatom)
    } // (loop over processors, iproc)
  } // if (me == 0)

  // cleanup:
  if (me == 0) {
    delete [] tags_known_unsorted;
    delete [] types_all;
    delete [] types_known_unsorted;
    if (q_all)
      delete [] q_all;
    if (q_known_unsorted)
      delete [] q_known_unsorted;
    delete [] atom_counts;
    delete [] atom_displacements;
  }
} //CheckAtomConsistency()


/* ---------------------------------------------------------------------- */

void FixBondModify::
CheckBondConsistency()
{
  // loop over bond list
  // setup possible partner list of bonds to break

  tagint *atom_counts = NULL; // sum of number of atoms "known"
                              //by each processor (nlocal+nghost)
  tagint *atom_displacements = NULL; // cumulative number of atoms known
  tagint natoms_known_all = 0; //sum of the number all atoms "known" 
                               //by each proc (nlocal+nghost).
                               //This will exceed the total number of 
                               //atoms in the simulation
                               //because some atoms are known by
                               //multiple procs. (ie. ghost atoms)
  tagint *tags_known_unsorted = NULL; //the atom-IDs for all the atoms
                               //known by all procs (including ghost atoms)
                               //This array has length "natoms_known_all"
  vector<set<int> > procs_that_know_atom; //which procs know about this atom?
  vector<set<tagint> > atoms_known_by(nprocs); //set of atoms known by each proc
  vector<set<pair<tagint,tagint> > >
    bonded_pairs_known_by(nprocs); //pairs of bonded atoms known by each proc

  if (me == 0) {
    atom_counts = new int[nprocs];
    atom_displacements = new tagint[nprocs + 1];
  }

  tagint n_known = atom->nlocal + atom->nghost;
  MPI_Gather(&n_known,1,MPI_INT,atom_counts,1,MPI_INT,0,world);
  if (me == 0) {
    atom_displacements[0] = 0;
    for (int i = 1; i <= nprocs; i++)
      atom_displacements[i] = atom_displacements[i-1] + atom_counts[i-1];
    natoms_known_all = atom_displacements[nprocs];
    tags_known_unsorted = new tagint[natoms_known_all];
    procs_that_know_atom.resize(atom->natoms + 1);
  }

  MPI_Gatherv(atom->tag, n_known, MPI_LMP_TAGINT,
              tags_known_unsorted,
              atom_counts, atom_displacements, MPI_LMP_TAGINT,
              0, world);

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      for (int i_localatom=0;
           i_localatom < atom_counts[iproc];
           i_localatom++) {
        tagint i = atom_displacements[iproc] + i_localatom;
        tagint tag = tags_known_unsorted[i];
        assert(tag <= atom->natoms);
        procs_that_know_atom[tag].insert(iproc);
        // -- for debugging only --
        //set<tagint> atoms_known_by_this_proc = atoms_known_by[iproc];
        //set<tagint> atoms_known_by_this_proc;
        //atoms_known_by_this_proc.insert(tag);
        // ----
        atoms_known_by[iproc].insert(tag);
      }
    }
  }

  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  tagint nbonds_known_all = 0;
  tagint *bond_counts_x3 = NULL;
  tagint *bond_displacements = NULL;
  int *bond_atoms_known_unsorted = NULL;
  if (me == 0) {
    bond_counts_x3 = new tagint[nprocs];
    bond_displacements = new tagint[nprocs + 1];
  }

  MPI_Gather(&nbondlist,1,MPI_INT,bond_counts_x3,1,MPI_INT,0,world);
  if (me == 0) {
    for (int i = 0; i < nprocs; i++)
      bond_counts_x3[i] *= 3; //entry stores {bond_atom1, bond_atom2, bond_type}
    bond_displacements[0] = 0;
    for (int i = 1; i <= nprocs; i++)
      bond_displacements[i] = bond_displacements[i-1] + bond_counts_x3[i-1];
    nbonds_known_all = bond_displacements[nprocs];
    bond_atoms_known_unsorted = new int[3*nbonds_known_all];
  }

  int *bond_atoms_local = new int[3*nbondlist];
  for (int i = 0; i < nbondlist; i++) {
    int i1 = bondlist[i][0];
    int i2 = bondlist[i][1];
    int btype = bondlist[i][2];
    bond_atoms_local[3*i + 0] = i1;
    bond_atoms_local[3*i + 1] = i2;
    bond_atoms_local[3*i + 2] = btype;
  }

  MPI_Gatherv(bond_atoms_local, 3*nbondlist, MPI_INT,
              bond_atoms_known_unsorted,
              bond_counts_x3,bond_displacements,MPI_INT,
              0, world);

  if (me == 0) {
    //look up one or more bonds between pairs of atoms
    map<pair<tagint,tagint>, set<int> > apair2btypes;

    // loop over processors
    for (int iproc = 0; iproc < nprocs; iproc++) {

      map<pair<tagint,tagint>, set<int> > apair2btypes_local;

      // loop over the bonds for each processor
      tagint bond_counts_x1 = bond_counts_x3[iproc] / 3;
      assert(bond_counts_x3[iproc] % 3 == 0);

      for (int ibond = 0; ibond < bond_counts_x1; ibond++) {

        tagint ibond_offset = bond_displacements[iproc] + 3*ibond;
        int i1  = bond_atoms_known_unsorted[ibond_offset + 0];
        int i2  = bond_atoms_known_unsorted[ibond_offset + 1];
        int btype=bond_atoms_known_unsorted[ibond_offset + 2];

        tagint iatom_offset = atom_displacements[iproc];
        tagint tag1 = tags_known_unsorted[iatom_offset + i1]; //global atom id
        tagint tag2 = tags_known_unsorted[iatom_offset + i2]; //global atom id
        if (tag2 < tag1) {
          // (ordering to avoid double-counting: insure that tag1 <= tag2)
          tagint tmp = tag1;
          tag1 = tag2;
          tag2 = tmp;
        }
        pair<tagint, tagint> pa(tag1, tag2);
        apair2btypes_local[pa].insert(btype);
        bonded_pairs_known_by[iproc].insert(pa);
      }

      // CHECK 1:
      // Loop over pairs of bond atoms known by processor "iproc".
      // Check to make sure that the type (and number) of bonds
      // connecting for this pair of atoms agrees with information
      // we have learned from other processors we came accross earlier.
      for (map<pair<tagint,tagint>, set<int> >::const_iterator
             p = apair2btypes_local.begin();
           p != apair2btypes_local.end();
           p++)
      {
        pair<tagint,tagint> patoms = p->first; //(key = a pair of atom-IDs)
        set<int> btypes = p->second;  //(value = bond(s) (types) for that pair)
        if (apair2btypes.find(patoms) != apair2btypes.end()) {
          if (apair2btypes[patoms] != btypes) {
            stringstream err_ss;
            err_ss << " (MPI) Inconsistent bond types between atoms: "
                   << patoms.first << "," << patoms.second << "\n";
            error->all(FLERR, err_ss.str().c_str());
          }
          else
            apair2btypes[patoms] = btypes;
        }
      }

    } //for (int iproc = 0; iproc < nprocs; iproc++)


    // CHECKS 2 & 3:
    // Now loop over all of the pairs of bonded atoms in the system.
    // Check whether each processor knows about the existence of bonds
    // between pairs of atoms which are known by that processor.
    for (map< pair<tagint,tagint>, set<int> >::const_iterator
           pbondpair = apair2btypes.begin();
         pbondpair != apair2btypes.end();
         pbondpair++)
    {
      pair<tagint, tagint> p = pbondpair->first;
      tagint i = p.first;
      tagint j = p.second;
      set<int> iprocs = procs_that_know_atom[i];
      set<int> jprocs = procs_that_know_atom[j];
      for (set<tagint>::const_iterator p_iproc = iprocs.begin();
           p_iproc != iprocs.end();
           p_iproc++)
      {
        int iproc = *p_iproc;
        stringstream err_ss;
        // CHECK 2:
        // Complain if any processor which is aware of atom i,
        // fails to know about the existence of the other atom (j),
        if (atoms_known_by[iproc].find(j) == atoms_known_by[iproc].end()) {
          err_ss << " (MPI) Inconsistent bonded neighbor-list for atoms:\n"
                 << "      " << i << ", " << j
               //<< "  (on procs: " << iproc << ", " << *(jprocs.begin()) << ")"
                 <<"\n";
          error->all(FLERR, err_ss.str().c_str());
        }
        // CHECK 3:
        // ...OR fails to know there is a bond between them.
        if (bonded_pairs_known_by[iproc].find(p) ==
            bonded_pairs_known_by[iproc].end()) {
          err_ss << "(MPI) Inconsistent bond between atoms:\n"
                 << "      " << i << ", " << j
            //<< "  (on procs: " << iproc << ", " << *(jprocs.begin()) << ")"
                 <<"\n";
          error->all(FLERR, err_ss.str().c_str());
        }
      } //loop over the processers which know about atom i

      for (set<tagint>::const_iterator p_jproc = jprocs.begin();
           p_jproc != jprocs.end();
           p_jproc++)
      {
        int jproc = *p_jproc;
        stringstream err_ss;
        // CHECK 2:
        // Complain if any processor which is aware of atom j,
        // fails to know about the existence of the other atom (i),
        if (atoms_known_by[jproc].find(i) == atoms_known_by[jproc].end()) {
          err_ss << "(MPI) Inconsistent (bonded) neighbor-list for atoms:\n"
                 << "      " << j << ", " << i
               //<< "  (on procs: " << jproc << ", " << *(iprocs.begin()) << ")"
                 <<"\n";
          error->all(FLERR, err_ss.str().c_str());
        }
        // CHECK 3:
        // ...OR fails to know there is a bond between them.
        if (bonded_pairs_known_by[jproc].find(p) ==
            bonded_pairs_known_by[jproc].end()) {
          err_ss << "(MPI) Inconsistent bond between atoms:\n"
                 << "      " << j << ", " << i
               //<< "  (on procs: " << jproc << ", " << *(iprocs.begin()) << ")"
                 <<"\n";
          error->all(FLERR, err_ss.str().c_str());
        }
      } //loop over the processers which know about atom j

    } //loop over all of the pairs of bonded atoms in the system

  } //if (me == 0)

  // cleanup:
  delete [] bond_atoms_local;
  if (me == 0) {
    delete [] bond_atoms_known_unsorted;
    delete [] atom_counts;
    delete [] atom_displacements;
    delete [] bond_counts_x3;
    delete [] bond_displacements;
  }

} // CheckBondConsistency()

#endif  //#ifdef DEBUG_FIX_BOND_MODIFY


/* ---------------------------------------------------------------------- */

bool FixBondModify::SatisfiesAtomProperties(int i1, int i2) { 
  return (
          (
           ((((atom->type[i1] >= atype1old_lo) &&
              (atom->type[i1] <= atype1old_hi)) || 
             (atype1old_lo > atype1old_hi))
            &&
            ((! atom->q_flag) ||
             (((atom->q[i1] >= qold1_lo) &&
               (atom->q[i1] <= qold1_hi)) ||
              (qold1_lo > qold1_hi))))
           &&
           ((((atom->type[i2] >= atype2old_lo) &&
              (atom->type[i2] <= atype2old_hi)) || 
             (atype2old_lo > atype2old_hi))
            &&
            ((! atom->q_flag) ||
             (((atom->q[i2] >= qold2_lo) &&
               (atom->q[i2] <= qold2_hi)) ||
              (qold2_lo > qold2_hi))))
           )
          && (ignore_tags_flag || (atom->tag[i1] < atom->tag[i2]))
          );
              
}




/* ---------------------------------------------------------------------- */

void FixBondModify::post_integrate()
{
  int i,j,k,m,n,i1,i2,n1,n3;
  double delx,dely,delz,rsq;
  tagint *slist;
  char dbg_msg[MSGLENGTH];

  if ((update->ntimestep - nevery_delay) % nevery) return;

  // check that all procs have needed ghost atoms within ghost cutoff
  // only if neighbor list has changed since last check

  if (lastcheck < neighbor->lastcall) check_ghosts();


  #ifdef DEBUG_FIX_BOND_MODIFY
  fprintf(screen,"checking consistency before timestep %d\n", update->ntimestep);
  CheckAtomConsistency();
  CheckBondConsistency();
  #endif

  // acquire updated ghost atom positions
  // necessary b/c are calling this after integrate, but before Verlet comm

  comm->forward_comm();

  // resize bond partner list and initialize it
  // probability array overlays distsq array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(partner);
    memory->destroy(finalpartner);
    memory->destroy(distsq);
    nmax = atom->nmax;
    memory->create(partner,nmax,"bond/break:partner");
    memory->create(finalpartner,nmax,"bond/break:finalpartner");
    memory->create(distsq,nmax,"bond/break:distsq");
    probability = distsq;
  }

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    partner[i] = 0;
    finalpartner[i] = 0;
    if (prioritize_long_bonds)
      distsq[i] = 0.0;
    else
      distsq[i] = INFINITY;
  }

  // loop over bond list
  // setup possible partner list of bonds to break

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int n_which_bond = -1;
  bool atype_changed = false;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    int btype = bondlist[n][2];
    if (!(mask[i1] & groupbit)) continue;
    if (!(mask[i2] & groupbit)) continue;

    if ((btypeold_lo <= btypeold_hi) &&
          (! ((btype >= btypeold_lo) &&
              (btype <= btypeold_hi))))
        continue;

    // Now consider the properties of the atoms participating in the bond
    if (SatisfiesAtomProperties(i1, i2)) {
      //i1_permute = i1;
      //i2_permute = i2;
    }
    else if (SatisfiesAtomProperties(i2, i1)) {
      //i1_permute = i2;
      //i2_permute = i1;
    }
    else
      continue;

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    rsq = delx*delx + dely*dely + delz*delz;
    if ((bond_rminsq >= 0) && (rsq <= bond_rminsq)) continue;
    if ((bond_rmaxsq >= 0) && (rsq >= bond_rmaxsq)) continue;




    //NOTE: For any given atom, at most one bond can be modified per iteration
    //      If an atom is bonded to multiple atoms which all satisfy the
    //      requirements, then we must choose one of these bonds.
    //      This can be arbitrary.  We try to guess a reasonable choice below:
    if (prioritize_long_bonds) {
      // If there is a minimum length requirement for the bond to be modified
      // then give priority to the longest bonds.  Typically users add
      // these requirements to bonds which are about to be broken.
      // The physical justification for choosing atoms pairs with long bonds
      // is that atoms whose bonds are more stretched are more likely to break
      if (rsq > distsq[i1]) {
        partner[i1] = tag[i2];
        distsq[i1] = rsq;
      }
      if (rsq > distsq[i2]) {
        partner[i2] = tag[i1];
        distsq[i2] = rsq;
      }
    }
    else {
      // In other cases, we (arbitrarily) give priority to the shortest bonds.
      // (The physical justification for this choice is that atoms which are
      //  closer together are more likely to be interacting strongly.)
      if (rsq <= distsq[i1]) {
        partner[i1] = tag[i2];
        distsq[i1] = rsq;
      }
      if (rsq <= distsq[i2]) {
        partner[i2] = tag[i1];
        distsq[i2] = rsq;
      }
    }
    n_which_bond = n;

  } //for (n = 0; n < nbondlist; n++) {

  // reverse comm of partner info

  if (force->newton_bond) comm->reverse_comm_fix(this);

  // each atom now knows its winning partner
  // for prob check, generate random value for each atom with a bond partner
  // forward comm of partner and random value, so ghosts have it

  if (prob_transition < 1.0) {
    for (i = 0; i < nlocal; i++)
      if (partner[i]) probability[i] = random->uniform();
  }

  commflag = 1;
  //comm->forward_comm_fix(this,2);
  comm->forward_comm_fix(this,4);

  // break bonds
  // if both atoms list each other as winning bond partner
  // and probability constraint is satisfied

  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  bool unequal_cutoffs = false;
  nbreak = 0;


  #ifdef ALLOW_BOND_MODIFY_SINGLE_ATOMS
  // First, consider rules which change individual atom types
  // (without regard for neighboring atom types, or nearby bonds)...
  // (Perhaps this should be a separate fix... I leave the feature in for now.)
  if (check_individual_atoms) {
    bool unequal_cutoffs = false;
    for (i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      if ((atypeold_lo <= atypeold_hi) &&
          (! ((atom->type[i] >= atypeold_lo) &&
              (atom->type[i] <= atypeold_hi))))
        continue;
      if (atom->q_flag &&
          (qold_lo < qold_hi) &&
          (! ((atom->q[i] >= qold_lo) &&
              (atom->q[i] <= qold_hi))))
        continue;
      if (random->uniform() > prob_transition) continue;
      #ifndef NDEBUG
      snprintf(dbg_msg, MSGLENGTH,
               "%s[%s], t=%d: atomtype(atomid=" TAGINT_FORMAT ") = %d-->%d, mpirank=%d\n",
               this->style, this->id, update->ntimestep, 
               tag[i], atom->type[i], atypenew, me);
      //error->message(FLERR, dbg_msg);
      fprintf(stderr, dbg_msg);
      #endif
      if (conserve_ke_flag) {
        double sqrt_mass_ratio = sqrt(atom->mass[atom->type[i]] /
                                      atom->mass[atypenew]);
        atom->v[i][0] *= sqrt_mass_ratio;
        atom->v[i][1] *= sqrt_mass_ratio;
        atom->v[i][2] *= sqrt_mass_ratio;
      }
      double **cutsq = force->pair->cutsq;
      if (cutsq[atom->type[i]] != cutsq[atypenew])
        unequal_cutoffs = true;
      atom->type[i] = atypenew;
      if (isfinite(qnew))
        atom->q[i] = qnew;
    } //for (i = 0; i < nlocal; i++) {

    // NOTE: THIS CODE APPEARED IN FixAtomSwap::pre_exchange(). DO WE NEED IT?
    //  if (domain->triclinic) domain->x2lamda(atom->nlocal);
    //  domain->pbc();
    //  comm->exchange();
    //  comm->borders();
    //  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    //  if (modify->n_pre_neighbor) modify->pre_neighbor();
    //  neighbor->build();
    //} else {
    comm->forward_comm_fix(this);  //<---????
    //}

  } //if (check_individual_atoms) {
  #endif //#ifdef ALLOW_BOND_MODIFY_SINGLE_ATOMS



  if (check_bonded_atoms) {
    // Next, consider rules which change bond types or atom-type-pairs
    for (i = 0; i < nlocal; i++) {
      if (partner[i] == 0) continue;
      j = atom->map(partner[i]);
      if (partner[j] != tag[i]) continue;

      // apply probability constraint using RN for atom with smallest ID

      if (prob_transition < 1.0) {
        if (tag[i] < tag[j]) {
          if (probability[i] >= prob_transition) continue;
        } else {
          if (probability[j] >= prob_transition) continue;
        }
      }


      assert(j == atom->map(partner[i]));
      int i1_permute, i2_permute;
      if (SatisfiesAtomProperties(i, j)) {
        i1_permute = i;
        i2_permute = j;
        // Special case: If both orderings satisfy the requirements, then assign
        // the "i1_permute" index to the atom with the lower atom-ID number(tag)
        if (SatisfiesAtomProperties(j, i) && (tag[j] < tag[i])) {
          i1_permute = j;
          i2_permute = i;
        }
      }
      else if (SatisfiesAtomProperties(j, i)) {
        i1_permute = j;
        i2_permute = i;
      }
      if ((atype1new >= 0) || (atype2new >= 0)) {
        #ifndef NDEBUG
        snprintf(dbg_msg, MSGLENGTH,
                 "%s[%s], t=%d: (atom->type[" TAGINT_FORMAT "],atom->type[" TAGINT_FORMAT
                 "]) = (%d,%d)-->(%d,%d), mpirank=%d\n",
                 this->style, this->id, update->ntimestep, 
                 tag[i1_permute], tag[i2_permute],
                 atom->type[i1_permute], atom->type[i2_permute],
                 ((atype1new>=0) ? atype1new : atom->type[i1_permute]),
                 ((atype2new>=0) ? atype2new : atom->type[i2_permute]),
                 me);
        //error->message(FLERR, dbg_msg);
        fprintf(stderr, dbg_msg);
        #endif


        if (conserve_ke_flag) {
          double sqrt_mass_ratio;
          if (atype1new >= 0) {
            sqrt_mass_ratio = sqrt(atom->mass[atom->type[i1_permute]] /
                                     atom->mass[atype1new]);
            atom->v[i1_permute][0] *= sqrt_mass_ratio;
            atom->v[i1_permute][1] *= sqrt_mass_ratio;
            atom->v[i1_permute][2] *= sqrt_mass_ratio;
          }
          if (atype2new >= 0) {
            sqrt_mass_ratio = sqrt(atom->mass[atom->type[i2_permute]] /
                                   atom->mass[atype2new]);
            atom->v[i2_permute][0] *= sqrt_mass_ratio;
            atom->v[i2_permute][1] *= sqrt_mass_ratio;
            atom->v[i2_permute][2] *= sqrt_mass_ratio;
          }
        }


   
        double **cutsq = force->pair->cutsq;
        if (((atype1new >= 0) &&
             (cutsq[atom->type[i1_permute]]!=cutsq[atype1new]))
            ||
            ((atype2new != 0) &&
             (cutsq[atom->type[i2_permute]]!=cutsq[atype2new])))
          unequal_cutoffs = true;

        if (atype1new >= 0) {
          atom->type[i1_permute] = atype1new;
          atype_changed = true;
        }
        if (atype2new >= 0) {
          atom->type[i2_permute] = atype2new;
          atype_changed = true;
        }

      } //if ((atype1new >= 0) || (atype2new >= 0)) {

      if (isfinite(qnew1))
        atom->q[i1_permute] = qnew1;
      if (isfinite(qnew2))
        atom->q[i2_permute] = qnew2;

      // NOTE: THIS CODE APPEARED IN FixAtomSwap::pre_exchange(). DO WE NEED IT?
      //if (unequal_cutoffs) { 
      //  if (domain->triclinic) domain->x2lamda(atom->nlocal);
      //  domain->pbc();
      //  comm->exchange();
      //  comm->borders();
      //  if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      //  if (modify->n_pre_neighbor) modify->pre_neighbor();
      //  neighbor->build();
      //} else {
      comm->forward_comm_fix(this);  // <---- ????
      //}

      
      if (break_bond) {

        //#ifndef NDEBUG
        {
          static const int MSGLENGTH=2048;
          char dbg_msg[MSGLENGTH];
          snprintf(dbg_msg, MSGLENGTH,
                   "%s[%s], t=%d, change_delete_bond(atomid=" TAGINT_FORMAT ",atomid=" TAGINT_FORMAT "), mpirank=%d\n",
                   this->style, this->id, update->ntimestep, 
                   tag[i], tag[j], me);
          error->message(FLERR, dbg_msg);
          //exit(0);
        }
        //#endif

        // delete bond from atom I if I stores it
        // atom J will also do this

        for (m = 0; m < num_bond[i]; m++) {
          if (bond_atom[i][m] == partner[i]) {
            for (k = m; k < num_bond[i]-1; k++) {
              bond_atom[i][k] = bond_atom[i][k+1];
              bond_type[i][k] = bond_type[i][k+1];
            }
            num_bond[i]--;
            break;
          }
        }

        // remove J from special bond list for atom I
        // atom J will also do this, whatever proc it is on

        slist = special[i];
        n1 = nspecial[i][0];
        for (m = 0; m < n1; m++)
          if (slist[m] == partner[i]) break;
        n3 = nspecial[i][2];
        for (; m < n3-1; m++) slist[m] = slist[m+1];
        nspecial[i][0]--;
        nspecial[i][1]--;
        nspecial[i][2]--;

        // store final broken bond partners and count the broken bond once

        finalpartner[i] = tag[j];
        finalpartner[j] = tag[i];
        if (tag[i] < tag[j]) nbreak++;
      } //if (break_bond) {...
      else if (btypenew >= 0) {
        #ifndef NDEBUG
        snprintf(dbg_msg, MSGLENGTH,
                 "%s[%s], t=%d: bond(atomid[" TAGINT_FORMAT "],atomid[" TAGINT_FORMAT
                 "]) = %d --> %d, mpirank=%d\n",
                 this->style, this->id, update->ntimestep, 
                 tag[i1_permute], tag[i2_permute],
                 bondlist[n_which_bond][2], btypenew,
                 me);
        //error->message(FLERR, dbg_msg);
        fprintf(stderr, dbg_msg);
        #endif
        bondlist[n_which_bond][2] = btypenew;
        //CONTINUEHERE:  THIS PROBABLY DOES NOT WORK IN PARALLEL. TEST IT
      } //else if (btypenew != -1) {...

    } //for (n = 0; n < nbondlist; n++) {...
  } //if (check_bonded_atoms) {...

  // tally stats

  MPI_Allreduce(&nbreak,&breakcount,1,MPI_INT,MPI_SUM,world);   //<--DO WE NEED THIS? -A 2819-7-08
  breakcounttotal += breakcount;  //<--DO WE NEED THIS? -A 2819-7-08
  atom->nbonds -= breakcount;

  // trigger reneighboring if any bonds were broken
  // this insures neigh lists will immediately reflect the topology changes
  // done if no bonds broken

  if (breakcount) next_reneighbor = update->ntimestep;
  //if (breakcount || atype_changed) next_reneighbor = update->ntimestep; //<--DO WE NEED THIS? -A 2819-7-08
  //if (breakcount || atype_changed) neighbor->build();  //<--DO WE NEED THIS?
  if (!breakcount) return;

  // communicate final partner and 1-2 special neighbors
  // 1-2 neighs already reflect broken bonds

  commflag = 2;
  comm->forward_comm_fix(this);

  // create list of broken bonds that influence my owned atoms
  //   even if between owned-ghost or ghost-ghost atoms
  // finalpartner is now set for owned and ghost atoms so loop over nall
  // OK if duplicates in broken list due to ghosts duplicating owned atoms
  // check J < 0 to insure a broken bond to unknown atom is included
  //   i.e. bond partner outside of cutoff length

  nbreak = 0;
  for (i = 0; i < nall; i++) {
    if (finalpartner[i] == 0) continue;
    j = atom->map(finalpartner[i]);
    if (j < 0 || tag[i] < tag[j]) {
      if (nbreak == maxbreak) {
        maxbreak += DELTA;
        memory->grow(broken,maxbreak,2,"bond/break:broken");
      }
      broken[nbreak][0] = tag[i];
      broken[nbreak][1] = finalpartner[i];
      nbreak++;
    }

  }

  // update special neigh lists of all atoms affected by any broken bond
  // also remove angles/dihedrals/impropers broken by broken bonds

  update_topology();

  // DEBUG
  // print_bb();

  //#ifdef DEBUG_FIX_BOND_MODIFY
  //fprintf(screen,"checking consistency after timestep %d\n", update->ntimestep);
  //CheckAtomConsistency();
  //CheckBondConsistency();
  //#endif

} //FixBondModify::post_integrate()



/* ----------------------------------------------------------------------
   insure all atoms 2 hops away from owned atoms are in ghost list
   this allows dihedral 1-2-3-4 to be properly deleted
     and special list of 1 to be properly updated
   if I own atom 1, but not 2,3,4, and bond 3-4 is deleted
     then 2,3 will be ghosts and 3 will store 4 as its finalpartner
------------------------------------------------------------------------- */

void FixBondModify::check_ghosts()
{
  int i,j,n;
  tagint *slist;

  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (i = 0; i < nlocal; i++) {
    slist = special[i];
    n = nspecial[i][1];
    for (j = 0; j < n; j++)
      if (atom->map(slist[j]) < 0) flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall)
    error->all(FLERR,"Fix bond/modify needs ghost atoms from further away");
  lastcheck = update->ntimestep;
}

/* ----------------------------------------------------------------------
   double loop over my atoms and broken bonds
   influenced = 1 if atom's topology is affected by any broken bond
     yes if is one of 2 atoms in bond
     yes if both atom IDs appear in atom's special list
     else no
   if influenced:
     check for angles/dihedrals/impropers to break due to specific broken bonds
     rebuild the atom's special list of 1-2,1-3,1-4 neighs
------------------------------------------------------------------------- */

void FixBondModify::update_topology()
{
  int i,j,k,n,influence,influenced,found;
  tagint id1,id2;
  tagint *slist;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int nlocal = atom->nlocal;

  nangles = 0;
  ndihedrals = 0;
  nimpropers = 0;

  //printf("NBREAK %d: ",nbreak);
  //for (i = 0; i < nbreak; i++)
  //  printf(" %d %d,",broken[i][0],broken[i][1]);
  //printf("\n");

  for (i = 0; i < nlocal; i++) {
    influenced = 0;
    slist = special[i];

    for (j = 0; j < nbreak; j++) {
      id1 = broken[j][0];
      id2 = broken[j][1];

      influence = 0;
      if (tag[i] == id1 || tag[i] == id2) influence = 1;
      else {
        n = nspecial[i][2];
        found = 0;
        for (k = 0; k < n; k++)
          if (slist[k] == id1 || slist[k] == id2) found++;
        if (found == 2) influence = 1;
      }
      if (!influence) continue;
      influenced = 1;

      if (angleflag) break_angles(i,id1,id2);
      if (dihedralflag) break_dihedrals(i,id1,id2);
      if (improperflag) break_impropers(i,id1,id2);
    }

    if (influenced) rebuild_special_one(i);
  }

  int newton_bond = force->newton_bond;

  int all;
  if (angleflag) {
    MPI_Allreduce(&nangles,&all,1,MPI_INT,MPI_SUM,world);
    if (!newton_bond) all /= 3;
    atom->nangles -= all;
  }
  if (dihedralflag) {
    MPI_Allreduce(&ndihedrals,&all,1,MPI_INT,MPI_SUM,world);
    if (!newton_bond) all /= 4;
    atom->ndihedrals -= all;
  }
  if (improperflag) {
    MPI_Allreduce(&nimpropers,&all,1,MPI_INT,MPI_SUM,world);
    if (!newton_bond) all /= 4;
    atom->nimpropers -= all;
  }
}

/* ----------------------------------------------------------------------
   re-build special list of atom M
   does not affect 1-2 neighs (already include effects of new bond)
   affects 1-3 and 1-4 neighs due to other atom's augmented 1-2 neighs
------------------------------------------------------------------------- */

void FixBondModify::rebuild_special_one(int m)
{
  int i,j,n,n1,cn1,cn2,cn3;
  tagint *slist;

  tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;

  // existing 1-2 neighs of atom M

  slist = special[m];
  n1 = nspecial[m][0];
  cn1 = 0;
  for (i = 0; i < n1; i++)
    copy[cn1++] = slist[i];

  // new 1-3 neighs of atom M, based on 1-2 neighs of 1-2 neighs
  // exclude self
  // remove duplicates after adding all possible 1-3 neighs

  cn2 = cn1;
  for (i = 0; i < cn1; i++) {
    n = atom->map(copy[i]);
    slist = special[n];
    n1 = nspecial[n][0];
    for (j = 0; j < n1; j++)
      if (slist[j] != tag[m]) copy[cn2++] = slist[j];
  }

  cn2 = dedup(cn1,cn2,copy);

  // new 1-4 neighs of atom M, based on 1-2 neighs of 1-3 neighs
  // exclude self
  // remove duplicates after adding all possible 1-4 neighs

  cn3 = cn2;
  for (i = cn1; i < cn2; i++) {
    n = atom->map(copy[i]);
    slist = special[n];
    n1 = nspecial[n][0];
    for (j = 0; j < n1; j++)
      if (slist[j] != tag[m]) copy[cn3++] = slist[j];
  }

  cn3 = dedup(cn2,cn3,copy);

  // store new special list with atom M

  nspecial[m][0] = cn1;
  nspecial[m][1] = cn2;
  nspecial[m][2] = cn3;
  memcpy(special[m],copy,cn3*sizeof(int));
}

/* ----------------------------------------------------------------------
   break any angles owned by atom M that include atom IDs 1 and 2
   angle is broken if ID1-ID2 is one of 2 bonds in angle (I-J,J-K)
------------------------------------------------------------------------- */

void FixBondModify::break_angles(int m, tagint id1, tagint id2)
{
  int j,found;

  int num_angle = atom->num_angle[m];
  int *angle_type = atom->angle_type[m];
  tagint *angle_atom1 = atom->angle_atom1[m];
  tagint *angle_atom2 = atom->angle_atom2[m];
  tagint *angle_atom3 = atom->angle_atom3[m];

  int i = 0;
  while (i < num_angle) {
    found = 0;
    if (angle_atom1[i] == id1 && angle_atom2[i] == id2) found = 1;
    else if (angle_atom2[i] == id1 && angle_atom3[i] == id2) found = 1;
    else if (angle_atom1[i] == id2 && angle_atom2[i] == id1) found = 1;
    else if (angle_atom2[i] == id2 && angle_atom3[i] == id1) found = 1;
    if (!found) i++;
    else {
      for (j = i; j < num_angle-1; j++) {
        angle_type[j] = angle_type[j+1];
        angle_atom1[j] = angle_atom1[j+1];
        angle_atom2[j] = angle_atom2[j+1];
        angle_atom3[j] = angle_atom3[j+1];
      }
      num_angle--;
      nangles++;
    }
  }

  atom->num_angle[m] = num_angle;
}

/* ----------------------------------------------------------------------
   break any dihedrals owned by atom M that include atom IDs 1 and 2
   dihedral is broken if ID1-ID2 is one of 3 bonds in dihedral (I-J,J-K.K-L)
------------------------------------------------------------------------- */

void FixBondModify::break_dihedrals(int m, tagint id1, tagint id2)
{
  int j,found;

  int num_dihedral = atom->num_dihedral[m];
  int *dihedral_type = atom->dihedral_type[m];
  tagint *dihedral_atom1 = atom->dihedral_atom1[m];
  tagint *dihedral_atom2 = atom->dihedral_atom2[m];
  tagint *dihedral_atom3 = atom->dihedral_atom3[m];
  tagint *dihedral_atom4 = atom->dihedral_atom4[m];

  int i = 0;
  while (i < num_dihedral) {
    found = 0;
    if (dihedral_atom1[i] == id1 && dihedral_atom2[i] == id2) found = 1;
    else if (dihedral_atom2[i] == id1 && dihedral_atom3[i] == id2) found = 1;
    else if (dihedral_atom3[i] == id1 && dihedral_atom4[i] == id2) found = 1;
    else if (dihedral_atom1[i] == id2 && dihedral_atom2[i] == id1) found = 1;
    else if (dihedral_atom2[i] == id2 && dihedral_atom3[i] == id1) found = 1;
    else if (dihedral_atom3[i] == id2 && dihedral_atom4[i] == id1) found = 1;
    if (!found) i++;
    else {
      for (j = i; j < num_dihedral-1; j++) {
        dihedral_type[j] = dihedral_type[j+1];
        dihedral_atom1[j] = dihedral_atom1[j+1];
        dihedral_atom2[j] = dihedral_atom2[j+1];
        dihedral_atom3[j] = dihedral_atom3[j+1];
        dihedral_atom4[j] = dihedral_atom4[j+1];
      }
      num_dihedral--;
      ndihedrals++;
    }
  }

  atom->num_dihedral[m] = num_dihedral;
}

/* ----------------------------------------------------------------------
   break any impropers owned by atom M that include atom IDs 1 and 2
   improper is broken if ID1-ID2 is one of 3 bonds in improper (I-J,I-K,I-L)
------------------------------------------------------------------------- */

void FixBondModify::break_impropers(int m, tagint id1, tagint id2)
{
  int j,found;

  int num_improper = atom->num_improper[m];
  int *improper_type = atom->improper_type[m];
  tagint *improper_atom1 = atom->improper_atom1[m];
  tagint *improper_atom2 = atom->improper_atom2[m];
  tagint *improper_atom3 = atom->improper_atom3[m];
  tagint *improper_atom4 = atom->improper_atom4[m];

  int i = 0;
  while (i < num_improper) {
    found = 0;
    if (improper_atom1[i] == id1 && improper_atom2[i] == id2) found = 1;
    else if (improper_atom1[i] == id1 && improper_atom3[i] == id2) found = 1;
    else if (improper_atom1[i] == id1 && improper_atom4[i] == id2) found = 1;
    else if (improper_atom1[i] == id2 && improper_atom2[i] == id1) found = 1;
    else if (improper_atom1[i] == id2 && improper_atom3[i] == id1) found = 1;
    else if (improper_atom1[i] == id2 && improper_atom4[i] == id1) found = 1;
    if (!found) i++;
    else {
      for (j = i; j < num_improper-1; j++) {
        improper_type[j] = improper_type[j+1];
        improper_atom1[j] = improper_atom1[j+1];
        improper_atom2[j] = improper_atom2[j+1];
        improper_atom3[j] = improper_atom3[j+1];
        improper_atom4[j] = improper_atom4[j+1];
      }
      num_improper--;
      nimpropers++;
    }
  }

  atom->num_improper[m] = num_improper;
}

/* ----------------------------------------------------------------------
   remove all ID duplicates in copy from Nstart:Nstop-1
   compare to all previous values in copy
   return N decremented by any discarded duplicates
------------------------------------------------------------------------- */

int FixBondModify::dedup(int nstart, int nstop, tagint *copy)
{
  int i;

  int m = nstart;
  while (m < nstop) {
    for (i = 0; i < m; i++)
      if (copy[i] == copy[m]) {
        copy[m] = copy[nstop-1];
        nstop--;
        break;
      }
    if (i == m) m++;
  }

  return nstop;
}

/* ---------------------------------------------------------------------- */

void FixBondModify::post_integrate_respa(int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */


int FixBondModify::pack_forward_comm(int n, int *list, double *buf,
                                    int pbc_flag, int *pbc)
{
  int i,j,k,m,ns;
  m = 0;
  int *type = atom->type;
  double *q = atom->q;


  #ifndef NDEBUG
  char dbg_msg[MSGLENGTH];
  snprintf(dbg_msg, MSGLENGTH,
           "pack_forward_comm() invoked, %s[%s], t=%d, commflag=%d, mpirank=%d\n",
           this->style, this->id, update->ntimestep, commflag, me);
  //error->message(FLERR, dbg_msg);
  fprintf(stderr, dbg_msg);
  #endif


  if (atom->q_flag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
      buf[m++] = q[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = type[j];
    }
  }





  if (commflag == 1) {
    //m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(partner[j]).d;
      buf[m++] = probability[j];
    }
    return m;
  }
  else {
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    //m = 0;
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = ubuf(finalpartner[j]).d;
      ns = nspecial[j][0];
      buf[m++] = ubuf(ns).d;
      for (k = 0; k < ns; k++)
      buf[m++] = ubuf(special[j][k]).d;
    }
    return m;
  }
}


/* ---------------------------------------------------------------------- */


void FixBondModify::unpack_forward_comm(int n, int first, double *buf)
{
  int i,j,m,ns,last;

  int *type = atom->type;
  double *q = atom->q;

  #ifndef NDEBUG
  char dbg_msg[MSGLENGTH];
  snprintf(dbg_msg, MSGLENGTH,
           "unpack_forward_comm() invoked, %s[%s], t=%d, mpirank=%d\n",
           this->style, this->id, update->ntimestep, me);
  //error->message(FLERR, dbg_msg);
  fprintf(stderr, dbg_msg);
  #endif

  m = 0;
  last = first + n;






  if (atom->q_flag) {
    for (i = first; i < last; i++) {
      type[i] = static_cast<int> (buf[m++]);
      q[i] = buf[m++];
    }
  } else {
    for (i = first; i < last; i++)
      type[i] = static_cast<int> (buf[m++]);
  }






  if (commflag == 1) {
    //m = 0;
    //last = first + n;
    for (i = first; i < last; i++) {
      partner[i] = (tagint) ubuf(buf[m++]).i;
      probability[i] = buf[m++];
    }

  } else {

    int **nspecial = atom->nspecial;
    tagint **special = atom->special;

    //m = 0;
    //last = first + n;
    for (i = first; i < last; i++) {
      finalpartner[i] = (tagint) ubuf(buf[m++]).i;
      ns = (int) ubuf(buf[m++]).i;
      nspecial[i][0] = ns;
      for (j = 0; j < ns; j++)
        special[i][j] = (tagint) ubuf(buf[m++]).i;
    }
  }
}


/* ---------------------------------------------------------------------- */

int FixBondModify::pack_reverse_comm(int n, int first, double *buf)
{
  int i,j,m,last;

  #ifndef NDEBUG
  char dbg_msg[MSGLENGTH];
  snprintf(dbg_msg, MSGLENGTH,
           "pack_reverse_comm() invoked, %s[%s], t=%d, mpirank=%d\n",
           this->style, this->id, update->ntimestep, me);
  //error->message(FLERR, dbg_msg);
  fprintf(stderr, dbg_msg);
  #endif

  m = 0;
  last = first + n;

  //if (atom->q_flag) {
  //  for (i = first; i < last; i++) {
  //    j = list[i];
  //    buf[m++] = type[j];
  //    buf[m++] = q[j];
  //  }
  //} else {
  //  for (i = first; i < last; i++) {
  //    j = list[i];
  //    buf[m++] = type[j];
  //  }
  //}


  for (i = first; i < last; i++) {
    buf[m++] = ubuf(partner[i]).d;
    buf[m++] = distsq[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixBondModify::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  #ifndef NDEBUG
  char dbg_msg[MSGLENGTH];
  snprintf(dbg_msg, MSGLENGTH,
           "unpack_reverse_comm() invoked, %s[%s], t=%d, mpirank=%d\n",
           this->style, this->id, update->ntimestep, me);
  //error->message(FLERR, dbg_msg);
  fprintf(stderr, dbg_msg);
  #endif

  //if (atom->q_flag) {
  //  for (i = 0; i < n; i++) {
  //    type[i] = static_cast<int> (buf[m++]);
  //    q[i] = buf[m++];
  //  }
  //} else {
  //  for (i = 0; i < n; i++) {
  //    type[i] = static_cast<int> (buf[m++]);
  //}

  for (i = 0; i < n; i++) {
    j = list[i];
    if (((prioritize_long_bonds) && (buf[m+1] > distsq[j])) ||
        ((! prioritize_long_bonds) && (buf[m+1] < distsq[j])))
    {
      partner[j] = (tagint) ubuf(buf[m++]).i;
      distsq[j] = buf[m++];
    } else m += 2;
  }
}

/* ---------------------------------------------------------------------- */

void FixBondModify::print_bb()
{
  for (int i = 0; i < atom->nlocal; i++) {
    printf("TAG " TAGINT_FORMAT ": %d nbonds: ",atom->tag[i],atom->num_bond[i]);
    for (int j = 0; j < atom->num_bond[i]; j++) {
      printf(" %d",atom->bond_atom[i][j]);
    }
    printf("\n");
    printf("TAG " TAGINT_FORMAT ": %d nangles: ",atom->tag[i],atom->num_angle[i]);
    for (int j = 0; j < atom->num_angle[i]; j++) {
      printf(" %d %d %d,",atom->angle_atom1[i][j],
             atom->angle_atom2[i][j],atom->angle_atom3[i][j]);
    }
    printf("\n");
    printf("TAG " TAGINT_FORMAT ": %d ndihedrals: ",atom->tag[i],atom->num_dihedral[i]);
    for (int j = 0; j < atom->num_dihedral[i]; j++) {
      printf(" %d %d %d %d,",atom->dihedral_atom1[i][j],
             atom->dihedral_atom2[i][j],atom->dihedral_atom3[i][j],
             atom->dihedral_atom4[i][j]);
    }
    printf("\n");
    printf("TAG " TAGINT_FORMAT ": %d %d %d nspecial: ",atom->tag[i],
           atom->nspecial[i][0],atom->nspecial[i][1],atom->nspecial[i][2]);
    for (int j = 0; j < atom->nspecial[i][2]; j++) {
      printf(" %d",atom->special[i][j]);
    }
    printf("\n");
  }
}

/* ---------------------------------------------------------------------- */

void FixBondModify::print_copy(const char *str, tagint m,
                              int n1, int n2, int n3, int *v)
{
  printf("%s %i: %d %d %d nspecial: ",str,m,n1,n2,n3);
  for (int j = 0; j < n3; j++) printf(" %d",v[j]);
  printf("\n");
}

/* ---------------------------------------------------------------------- */

double FixBondModify::compute_vector(int n)  //<--DO WE NEED THIS? -A 2819-7-08
{
  if (n == 0) return (double) breakcount;    //<--DO WE NEED THIS?
  return (double) breakcounttotal;           //<--DO WE NEED THIS?
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBondModify::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0;
  bytes += atom->nmax * 2 * sizeof(int);
  bytes += 2*nmax * sizeof(tagint);
  bytes += nmax * sizeof(double);
  return bytes;
}
