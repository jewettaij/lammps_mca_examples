Molecular Cellular Automata Examples for LAMMPS
===========

The simulation of complex chains of chemical reactions
(such as the reactions that make life possible),
cannot be run with traditional molecular dynamics software.
At the coarse grained level, such simulations have been out of reach
to scientists without writing thousands of lines of custom
simulation code dedicated to the particular problem they are studying.

MCA is a language to describe the behavior of simple "agents"
which can make decisions in response to collisions
(or other changes in their local environment).
When connected together, such agents can make intelligent decisions
of arbitrary complexity, typically with just a few lines of MCA code.

This repository contains simple examples of MCA usage in
[LAMMPS](https://lammps.sandia.gov), as well as documentation for their use
in the form of slides for several short talks I gave explaining these examples,
and [LAMMPS code](LAMMPS_code_needed) needed to make these examples work.
*(The slides are available in both PDF and ODP format.
ODP files require [libre-office](https://www.libreoffice.org) to view.)*

Although videos of MCA simulations were too large to include in this repository,
some can be found here:

   https://www.youtube.com/watch?v=0cuEeCdyokU

   https://www.youtube.com/watch?v=QO4LbHGAgxU

   https://www.youtube.com/watch?v=EEbt07vZHew


## WARNING: These examples are not yet ready for public use.

This repository was intended as a private place
to back up files that I eventually plan to release.

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
and the new code for LAMMPS which enables the MCA features.
That code (and compilation instructions) is available
[here](LAMMPS_code_needed).
Again, this code (and these examples) could change or be deleted in the future.

## License

Although not recommended for use, these examples (and code)
are available under the terms of the MIT license.
