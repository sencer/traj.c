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
