/* Header file for C-compatible variables/functions in energy.F90 */

int natoms;          /* Number of atoms - set on first call to either function below        */
double model_energy; /* Total energy of the system of particles as most recently calculated */
double skin;         /* Verlet neighbour list skin width                                    */
double lj_lrc;       /* Invariant part of long range tail correction (missing N*N/volume)   */

/* Calculates the total energy of the system of particles */
double compute_model_energy(int n, int d, double *pos, int dh2, int dh1, double *hmatrix, int n2, int *species);

/* Calculates the energy of particle i interacting with its neighbours */
double compute_local_real_energy(int i, int n, int d, double *pos, int dh2, int dh1, double *hmatrix, int n2, int *species);

/* Recalculates the Verlet neighbour list. Should be called if atoms have moved sufficiently far */
/* that the originally computed list is invalidated.                                             */

/* May need to expose additional functions here if cell changes (need to recompute image vectors) */
void compute_neighbour_list(int n, int d, double *pos, int dh2, int dh1, double *hmatrix, int n2, int *species);

/* Interfaces to utility functions for testing only */
double fc(int isp, int jpc, double r);
double fs(int isp, double zc);
double tijk(double c, int isp, int jsp, int ksp);
double tijk_fast(double c, int isp, int jsp, int ksp); 
