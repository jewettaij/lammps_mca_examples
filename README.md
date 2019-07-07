lammps_mca_examples
===========

## Molecular Cellular Automata Examples for LAMMPS

The simulation of complex chains of chemical reactions that make life
possible cannot be run with traditional molecular dynamics software.
At the coarse grained level such simulations have been out of reach
to scientists without writing thousands of lines of custom code
dedicated to the particular problem they are studying.

This repository contains simple examples of MCA usage in
[LAMMPS](https://lammps.sandia.gov)
MCA is a language to describe the behavior of simple "agents"
which can make decisions in response to collisions
(or other changes in their local environment).
When connected together, such agents can make intelligent decisions
of arbitrary complexity, typically with just a few lines of code.

Videos of MCA simulations were too large to include in this repository,
However some MCA videos are available on youtube.

   https://www.youtube.com/watch?v=0cuEeCdyokU
   https://www.youtube.com/watch?v=QO4LbHGAgxU
   https://www.youtube.com/watch?v=EEbt07vZHew

## WARNING: These examples are not yet ready for public use.

This repository was intended as a private place to back up files
before releasing them formally.
(I an not willing to pay for private repositories on github.)
Please be warned that the LAMMPS MCA code needed by these examples
is currently buggy, and will change significantly
(or may be incorporated into other existing fixes).
Consequently, the syntax of the LAMMPS fix commands
(in the various .odp and .lt files)
will change significantly over time.

## Requirements

These examples require
[LAMMPS](https://lammps.sandia.gov), 
[MOLTEMPLATE](https://moltemplate.org), 
and the code for the new MCA lammps fix.
(This code was intentionally omitted but will be released soon.
 -A 2019-7-07)

## License

Although currently not recommended for use,
these examples are available under the terms of the MIT license.
