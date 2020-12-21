LAMMPS code needed to test these examples
===============================

The examples included in this repository require extra features that
are not included in LAMMPS (as of 2020-12-21).

To run these examples, you must:

1) Obtain the LAMMPS source code.  (https://github.com/lammps/lammps)
2) Copy the .cpp and .h files in this directory into the *src/*
subdirectory distributed with LAMMPS
3) Compile LAMMPS.  Here is one way to do this.
Enter the *src/* subdirectory.  Then enter:
```
make yes-molecule
make yes-mc
make yes-user-reaction
make yes-kspace
make yes-manybody
make yes-rigid
make yes-user-misc
make yes-user-fep
make yes-replica
make yes-asphere
make yes-class2
make yes-replica
make mpi
```
4) Copy the LAMMPS binary ("lmp_mpi" in this case) into somewhere in your $PATH,
or edit your ~/.bashrc or ~/.cshrc files to include the *src/* subdirectory.
If you do not know what a `PATH` environment variable is and are curious,
read [this](http://www.linfo.org/path_env_var.html).

### WARNING: Experimental code

This code is highly experimental.
It will likely be changed significantly in the future.
Alternatively, there is a high chance the features of this code may be
merged into *fix bond/react* (which has many advanced features this code lacks).

I just wanted to make this code available in its current state so that people
can try the examples.  However these examples will probably change significantly
in the future as well.  Eventually, this entire github repository could be
removed, and its contents merged into other projects (like LAMMPS or
moltemplate).

### WARNING: The code here might not be up to date

This code was chosen to work with the examples in their current form as of
2020-12-21.  In 2021, recent versions of this code should be available here:

1) https://github.com/jewettaij/lammps/blob/mca_attempt_2019/src/USER-MISC/fix_bond_modify.cpp
2) https://github.com/jewettaij/lammps/blob/mca_attempt_2019/src/USER-MISC/fix_bond_modify.h
3) https://github.com/jewettaij/lammps/blob/mca_attempt_2019/src/USER-MISC/fix_bond_new.cpp
4) https://github.com/jewettaij/lammps/blob/mca_attempt_2019/src/USER-MISC/fix_bond_new.h

If I ever decide to submit this code to be included with LAMMPS (and if the
developers accept the submission), then the most up-to-date version of this
code will eventually be distributed with LAMMPS in the *src/USER-MISC/*
subdirectory.
