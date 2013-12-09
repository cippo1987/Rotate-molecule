This Fortran code generate a poscar in cartesian coordinates, and eventually in cartesian poscar (POSCARc) with fractional coordinates of a perfect cubic perovskite structure
with a MA molecule in the center of the cell.
In order to operate the software need an input.
The input structure is simple and it is the following:
#Start f input
LATTICE PARAMETER
ANGLE1 ANGLE2.
#end of input
CAVEAT:
1)To define a rotation 3 angles should be needed. In this case the rotation coaxial to the molecule is ignored.
2)To avoid problem with data precision the format for the lattice is f4.3 (i.e. 4.534) while the one for the angles is  f5.2 (i.e 00.00 or 45.00, writing just 0, 45 or 80.3 leads to errors)

The code is mean to be adaptive for any molecule inside the cage but it should be directly changed inside.
In short the software it is just a rotator subroutine that prints out vasp\aims input.

 
