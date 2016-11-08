## traj.c

Some C code to analyze molecular dynamics trajectories. Reads/writes lammpstrj 
and xyz files. Focus is on finding the bonding structure. See the executable 
files (and their source codes) as examples how to use:

- overlap checks if there are any overlapping atoms. Much faster than 
  comparing every atom with every other.

- bond determines the bonding structure of a trajectory, consisting of CHNO 
  elements. It outputs some very system specific results (number of biphenyl 
  rings, CO, H2, N2) and total fragments.

- unwrap again determines bonding structure, and "fixes" the molecules broken 
  due to periodic boundary conditions.

Code define 4 types, read the header files and the source for detailed 
information.

- `Crystal` type defined in `crystal.h` holds information about the system: 
  number of atoms, their coordinates and PBC info.
- `CoarseBox` type defined in `boxes.h` divides space in the small boxes and 
  places each atom to the closest box. This is used for determining which atom 
  is close to which one. To use this, one clearly needs a `Crystal`.
- `BondingInfo` type defined in `bonding.h` holds bonding information of the 
  system. One needs a `Crystal` and a `CoarseBox` to use this.
- `Fragments` holds the information about molecules or fragments; that is the 
  atoms that have a bonding connection.

