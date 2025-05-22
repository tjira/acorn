typedef unsigned long ulong;

void coulomb(double *ints, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis);
void kinetic(double *ints, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis);
void nuclear(double *ints, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis);
void overlap(double *ints, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis);
