#include <stdbool.h>

typedef unsigned long ulong;

void coulomb(double *I, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis);
void kinetic(double *I, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis);
void nuclear(double *I, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis);
void overlap(double *I, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis);

void fock(double *F, ulong natoms, const double *anums, const double *coords, ulong nbasis, const double *basis, const double *D, bool G);
