/**
 * This code was developed by students for a final Master's project
 * As it was developed by students and not professionals, the code may require adjustments / improvements
 * This is an adaptation of an existing code originally written in Gerris
 * Some parts were improved with the use of AI.
 *
 * ===============================================
 * Debris Flow Numerical Simulation — Basilisk C
 * ===============================================

 * Simulates the propagation of a debris flow over a real terrain using:
 *   - DEM loaded (.asc)
 *   - Quadtree adaptive grid
 *   - Saint-Venant shallow-water equations
 *   - Voellmy friction law (basal friction µ + turbulent friction ξ)
 *   - Adaptive Mesh Refinement (AMR)
 *
 * Inputs:
 *   - DEM without rockfall deposit  (without_rockfall.asc)
 *   - Release area / mobile layer   (only_rockfall.asc)
 *   - Difference of DEMs            (DoD.asc)
 *
 * Output:
 *   - Flow depth raster (localDepth) as ESRI ASCII
 *   - Simulation log (simulation_log.txt)
 */


#include "grid/quadtree.h"
#include "saint-venant.h"
#include "output.h"
#include "utils.h"
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>


/* --------------------- PARAMETERS --------------------- */

#ifndef DRY
# define DRY   1e-4   // Minimum flow depth threshold — cells below this are considered dry
#endif

#ifndef GRAV
# define GRAV  9.81   // Gravitational acceleration (m/s²)
#endif

#ifndef LEVEL
# define LEVEL 11      // Maximum refinement level — grid resolution = LENGTH / 2^LEVEL
#endif

#ifndef MU
# define MU    0.155  // 0.12 // Voellmy basal friction coefficient (dimensionless)
#endif

#ifndef XI
# define XI    100.   // 500. // Voellmy turbulent friction coefficient (m/s²)
#endif

#ifndef EPS
# define EPS   1e-6   // Small epsilon to avoid division by zero
#endif

#ifndef END
# define END 600.     // Simulation end time (s)
#endif

/* Bottom left coordinate - Domain origin — should match the raster extent */
#ifndef TX
# define TX 614323.75     // X origin (m) = xllcorner
#endif

#ifndef TY
# define TY 5151476.2499  // Y origin (m) = yllcorner
#endif

#ifndef LENGTH
# define LENGTH 5117.5   // Domain size (m) — for ncols=2048 and cellsize=2.5 = 5 120
#endif                   // Here the domain is slightly smaller (2047x2047) -> 2047*2.5 = 5 117.5

#ifndef PATH
# define PATH "terrain"  // Directory containing input .asc files
#endif


/* ---------------------- VARIABLES --------------------- */

/* saint-venant.h is already providing :
   scalar h[];  - flow depth
   vector u[];  - depth-averaged velocity
   scalar zb[]; - bed topography
*/

scalar RA[];          // Release area — initial flow depth (m)
scalar localDepth[];  // Flow depth corrected for slope angle (m)
scalar Vtotal[];      // Total velocity magnitude (m/s)
scalar DoD[];         // Difference of DEMs — observed erosion/deposition (m)


/* -------------------- RUN METADATA  ------------------- */

/* Needed to measure simulation run time */
static time_t run_t0 = 0;  // Wall-clock start time
static clock_t cpu_t0 = 0; // CPU start time

static char exe_name[256] = "unknown_exe";   // Executable name
static char src_name[256] = "unknown_src";   // Source file name
static char dod_file[1024] = "unknown_dod";  // Path to DoD raster
static char dem_file[1024] = "unknown_dem";  // Path to DEM raster
static char ra_file[1024]  = "unknown_ra";   // Path to release area raster

// Set the volume to 0 
static double Vobs_dep_fixed = 0.;   // Observed deposition volume (m³)
static double Vobs_ero_fixed = 0.;   // Observed erosion volume (m³)  

scalar dzdx[];  // Bed slope in x-direction (dz/dx)
scalar dzdy[];  // Bed slope in y-direction (dz/dy)


/* ---------------------- UTILITIES --------------------- */

static inline double clampd (double v, double lo, double hi)
{
  return v < lo ? lo : (v > hi ? hi : v);
}

/**
 * Reads an ESRI ASCII raster (.asc) and projects it onto the Basilisk grid.
 *
 * - Parses the 6-line header (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA)
 * - Stores values in a temporary 2D array (row 0 = south)
 * - Projects onto Basilisk cells using bilinear interpolation
 *
 * @param s               Target Basilisk scalar field
 * @param filename        Path to the .asc file
 * @param nodata_to_zero  If 1, NODATA values are replaced by 0
 */


static void read_esri_asc_into (scalar s, const char *filename, int nodata_to_zero)
{
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Cannot open %s\n", filename);
    exit(1);
  }

  int ncols = 0, nrows = 0;
  double xll = 0., yll = 0., cellsize = 0., nodata_read = 0.;
  char key[64];

  // Header (6-line)
  for (int k = 0; k < 6; k++) {
    if (fscanf(fp, "%63s", key) != 1) {
      fprintf(stderr, "Header read error in %s\n", filename);
      exit(1);
    }
    if (!strcmp(key, "ncols"))                     fscanf(fp, "%d", &ncols);
    else if (!strcmp(key, "nrows"))                fscanf(fp, "%d", &nrows);
    else if (!strcmp(key, "xllcorner") ||
             !strcmp(key, "xllcenter"))            fscanf(fp, "%lf", &xll);
    else if (!strcmp(key, "yllcorner") ||
             !strcmp(key, "yllcenter"))            fscanf(fp, "%lf", &yll);
    else if (!strcmp(key, "cellsize"))             fscanf(fp, "%lf", &cellsize);
    else if (!strcmp(key, "NODATA_value") ||
             !strcmp(key, "nodata_value"))         fscanf(fp, "%lf", &nodata_read);
    else {
      fprintf(stderr, "Unknown header key '%s' in %s\n", key, filename);
      exit(1);
    }
  }

  if (ncols <= 0 || nrows <= 0 || cellsize <= 0.) {
    fprintf(stderr, "Bad header in %s (ncols=%d nrows=%d cell=%g)\n",
            filename, ncols, nrows, cellsize);
    exit(1);
  }

  fprintf(stderr,
          "%s:\n  ncols=%d nrows=%d\n  xll=%0.12g yll=%0.12g\n  cell=%g nodata=%g\n",
          filename, ncols, nrows, xll, yll, cellsize, nodata_read);

  // Allocate 2D array — row 0 = south (ESRI stores data from north to south, we flip it)
  double **a = (double **) malloc(nrows*sizeof(double *));
  for (int j = 0; j < nrows; j++)
    a[j] = (double *) malloc(ncols*sizeof(double));

  // Read raster values and flip rows (ESRI: north→south, we store: south→north)
  for (int jf = 0; jf < nrows; jf++) {
    int j = nrows - 1 - jf; 
    for (int i = 0; i < ncols; i++) {
      double v;
      if (fscanf(fp, "%lf", &v) != 1) {
        fprintf(stderr, "Data read error in %s at file row=%d col=%d\n",
                filename, jf, i);
        exit(1);
      }
      if (nodata_to_zero && v == nodata_read)
        v = 0.;
      a[j][i] = v;
    }
  }
  fclose(fp);

  // Project the raster onto Basilisk grid using bilinear interpolation
  foreach() {
    double xr = (x - xll)/cellsize - 0.5;
    double yr = (y - yll)/cellsize - 0.5;

    xr = clampd(xr, 0., ncols - 1.001);
    yr = clampd(yr, 0., nrows - 1.001);

    int i0 = (int) floor(xr), j0 = (int) floor(yr);
    int i1 = i0 + 1,        j1 = j0 + 1;
    if (i1 >= ncols) i1 = ncols - 1;
    if (j1 >= nrows) j1 = nrows - 1;

    double tx = xr - i0, ty = yr - j0;

    double v00 = a[j0][i0], v10 = a[j0][i1];
    double v01 = a[j1][i0], v11 = a[j1][i1];

    s[] = (1.-tx)*(1.-ty)*v00 + tx*(1.-ty)*v10 + (1.-tx)*ty*v01 + tx*ty*v11;
  }
  boundary({s});

  for (int j = 0; j < nrows; j++) free(a[j]);
  free(a);
}


/* ----------------- EXPORT ESRI ASCII  ----------------- */
/*
    -Exports a Basilisk scalar field as an ESRI ASCII raster (.asc).
    -Values outside the domain are written as -9999 (NODATA).
*/

  void output_asc (scalar s, const char *fname)
{
  FILE *fp = fopen(fname,"w");
  if (!fp) {
    fprintf(stderr,"Cannot write %s\n",fname);
    return;
  }

  int N = 1 << LEVEL;     // Number of cells per side = 2^LEVEL
  double dxout = L0/N;    // Output cell size (m)

  /* Write ESRI ASCII header */
  fprintf(fp,"ncols %d\n",N);
  fprintf(fp,"nrows %d\n",N);
  fprintf(fp,"xllcorner %.10f\n",X0);
  fprintf(fp,"yllcorner %.10f\n",Y0);
  fprintf(fp,"cellsize %.10f\n",dxout);
  fprintf(fp,"NODATA_value -9999\n");

  /* Write values row by row, north to south (ESRI convention) */
  for (int j = N-1; j >= 0; j--) {
    for (int i = 0; i < N; i++) {

      Point point = locate(X0 + (i + 0.5)*dxout,
                           Y0 + (j + 0.5)*dxout);

      double val = -9999.;

      if (point.level >= 0)
        val = interpolate(s,
                          X0 + (i + 0.5)*dxout,
                          Y0 + (j + 0.5)*dxout);

      fprintf(fp,"%.6f ",val);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);
}

/* -----------------------  MAIN  ----------------------- */

 /*
 * Entry point — sets up the simulation domain and launches the solver.
 *
 * Domain:
 *   - Origin at (TX, TY) in projected coordinates (m) (Bottom left coordinate of the asc)
 *   - Square domain of side LENGTH (m)
 *   - Initial grid at 2^LEVEL cells
 */ 


int main (int argc, char **argv)
{
  size (LENGTH);
  origin (TX, TY);
  init_grid (1 << LEVEL);

  /* Solver settings */
  G   = GRAV;     // gravitational acceleration
  CFL = 0.5;      // courant number (stability criterion) - Courant-Friedrichs-Lewy number
  DT  = 1e-1;     // maximum time step (s)

  /* Store run start times */
  run_t0 = time(NULL);
  cpu_t0 = clock();
  char start_time[32];
  struct tm *tm0 = localtime(&run_t0);
  strftime(start_time, sizeof(start_time), "%H:%M:%S", tm0);
  
  /* Print the inputs */
  fprintf(stderr,
  "\n========================================\n"
  "Simulation started\n"
  "Start time  : %s\n"
  "Code        : %s\n"
  "DoD file    : %s/DoD.asc\n"
  "MU          : %g\n"
  "XI          : %g\n"
  "LEVEL       : %d\n"
  "END         : %g\n"
  "LENGTH      : %g\n"
  "========================================\n\n",
  start_time, argv[0] , PATH, MU, XI, LEVEL, END, LENGTH);


  // executable name
  strncpy(exe_name, argv[0], sizeof(exe_name)-1);
  exe_name[sizeof(exe_name)-1] = '\0';

  // source file (compile-time)
  strncpy(src_name, __FILE__, sizeof(src_name)-1);
  src_name[sizeof(src_name)-1] = '\0';

  run();
}



/* ------------------- INITIALISATION ------------------- */

/**
 * Sets up the simulation at t = 0:
 *   - Loads the DEM, release area and DoD from ESRI ASCII rasters
 *   - Computes bed slopes (dzdx, dzdy) for the friction law
 *   - Sets initial flow depth h = RA and velocity u = 0
 */

event init (t = 0)
{
  char fz[1024], fr[1024], fd[1024];
  snprintf(fz, sizeof(fz), "%s/without_rockfall.asc", PATH);   // DEM without rockfall
  snprintf(fr, sizeof(fr), "%s/only_rockfall.asc", PATH);      // Release area
  snprintf(fd, sizeof(fd), "%s/DoD.asc", PATH);                // Difference of DEMs — loaded for reference only - not used in the simulation (used in post processing)

  /* Store file paths for logging */
  strncpy(dem_file, fz, sizeof(dem_file)-1);
  dem_file[sizeof(dem_file)-1] = '\0';

  strncpy(ra_file, fr, sizeof(ra_file)-1);
  ra_file[sizeof(ra_file)-1] = '\0';

  strncpy(dod_file, fd, sizeof(dod_file)-1);
  dod_file[sizeof(dod_file)-1] = '\0';
  
  /* Load rasters onto Basilisk scalar fields */  
  read_esri_asc_into (zb, fz, 0);      // Bed topography
  read_esri_asc_into (RA, fr, 1);      // Release area 
  read_esri_asc_into (DoD, fd, 0);     // DoD    

  /* Compute bed slopes */
  foreach() {
  dzdx[] = (zb[1,0] - zb[-1,0])/(2.*Delta);
  dzdy[] = (zb[0,1] - zb[0,-1])/(2.*Delta);
  }
  boundary({dzdx, dzdy});

  fprintf(stderr, "OBS fixed volumes: Vdep=%g m3  Vero=%g m3\n",
          Vobs_dep_fixed, Vobs_ero_fixed);

  /* Initial conditions: flow depth = release area, velocity = 0 everywhere */
  foreach() {
    h[]   = (RA[] > DRY ? RA[] : 0.);
    u.x[] = u.y[] = 0.;
  }
  boundary({h, u, zb, RA, DoD});

  /* Print initial volume and field statistics */
  double Vinit = 0.;
  foreach(reduction(+:Vinit))
  Vinit += h[] * sq(Delta);
  fprintf(stderr, "Volume initial h[] : %.4g m3\n", Vinit);
  stats sZ = statsf(zb), sR = statsf(RA), sH = statsf(h);
  fprintf(stderr, "zb: min=%g max=%g avg=%g\n", sZ.min, sZ.max, sZ.sum/sZ.volume);
  fprintf(stderr, "RA: min=%g max=%g avg=%g\n", sR.min, sR.max, sR.sum/sR.volume);
  fprintf(stderr, " h: min=%g max=%g avg=%g\n", sH.min, sH.max, sH.sum/sH.volume);
}


/* ----------------- VOELLMY FRICTION  ------------------- */
 
/** 
 * Applies the Voellmy friction law at each time step.
 *
 * The Voellmy model splits basal resistance into two terms:
 *   - Coulomb friction : µ * N  (solid friction, proportional to normal stress)
 *   - Turbulent drag   : g * u² / ξ  (velocity-dependent, turbulent-like)
 *
 * The friction is applied as a velocity reduction factor F:
 *   u_new = F * u_old
 *
 * where F < 1 accounts for the energy dissipated by friction over dt.
 * Cells below the dry threshold are zeroed out.
 *
 * @param µ   Coulomb friction coefficient (dimensionless)
 * @param ξ   Turbulent friction coefficient (m/s²)
 */

event voellmy (i++)
{
  foreach() {
    if (h[] > DRY) {

      double DX = dzdx[];
      double DY = dzdy[];

      /* Slope geometry */
      double TAN2PHI = DX*DX + DY*DY;          //slope angle
      double SIN2PHI = TAN2PHI/(1. + TAN2PHI); //slope angle
      double COSPHI  = 1./sqrt(1. + TAN2PHI);  //slope angle

      /* Flow velocity and direction */
      double Ux = u.x[], Uy = u.y[];
      double Vel = sqrt(Ux*Ux + Uy*Uy) + EPS; // speed (+ eps avoids division by zero)

      /* Angle between velocity vector and slope direction */
      double TANPSI = -(DX*Ux + DY*Uy)/Vel;
      double COSPSI = 1./sqrt(1. + TANPSI*TANPSI);

      /* Velocity reduction factor — combines slope, Coulomb and turbulent friction */
      double F = Vel /
        (Vel + dt*GRAV*
         ( SIN2PHI*TANPSI
         + MU*COSPHI*COSPSI                    // Coulomb friction term
         + Vel*Vel/(h[]*XI*COSPHI*COSPSI)));   // turbulent friction term

      u.x[] *= F;
      u.y[] *= F;
    }
    else {
      h[] = 0.;
      u.x[] = u.y[] = 0.;
    }
  }
  boundary({h, u});
}


/* --------------------- DIAGNOSTICS -------------------- */

/*  * Computes output fields
 *   localDepth — flow depth perpendicular to the slope surface (m)
 *                h is the vertical depth; on steep terrain it overestimates
 *                the real flow thickness, so we project it onto the slope:
 *                localDepth = h / sqrt(1 + dzdx² + dzdy²)
 *
 *   Vtotal     — total velocity magnitude including the vertical component
 *                induced by the slope (m/s)
*/

event diagnostics (t = END)  // t = END as we want only the final deposit
{                            // If you want an animation : set to t = += 5 for example
  foreach() {

    double DX = dzdx[];
    double DY = dzdy[];

    if (h[] > DRY) {
      Vtotal[] =
        sqrt( sq(u.x[]) + sq(u.y[]) + sq(u.x[]*DX + u.y[]*DY) ) / h[];
      localDepth[] =
        h[]/sqrt(1. + DX*DX + DY*DY);
    }
    else {
      Vtotal[] = 0.;
      localDepth[] = 0.;
    }
  }
  boundary({localDepth});
}


/* ------------------------- AMR ------------------------ */

/**
 * Refines and coarsens the quadtree grid every time step.
 *
 * Uses a wavelet-based error estimator on the flow depth h:
 *   - Cells where h varies significantly are refined (up to LEVEL)
 *   - Cells where h is nearly constant are coarsened (down to 4)
 *
 * This concentrates computational effort where the flow is active
 * and reduces cost in dry or flat regions.
 *
 * Tolerance: 1e-3 m — adjust to control refinement aggressiveness.
 */

event adapt (i += 1)
{
  int minl = (t > 0 ? 1 : 5);

adapt_wavelet ((scalar *){h},
               (double []){1e-3},
               maxlevel = LEVEL,
               minlevel = 4);
}


/* ----------------------- OUTPUTS ---------------------- */

/*
* If you want to output an .asc file at regular time steps: 
* Use the same time steps in the Diagnostics event 
*/

/**
 * Exports localDepth as an ESRI ASCII raster at regular time intervals.
 * Each file is named with the current µ, ξ, refinement level and time,
 * allowing post-processing and animation from the output sequence.
 *
 * Note: localDepth must be computed before export — the slope correction
 * (h / sqrt(1 + dzdx² + dzdy²)) is applied inline here rather than
 * relying on the diagnostics event, so this event is self-contained.
 */

/*
event output_fields (i++)
{
  static double next = 0.;   // Next output time
  char name[256];

  if (t >= next - 1e-12) {

    sprintf(name,
            "results/localDepth-mu%0.3f-xi%0.0f-L%d-t%04.0f.asc", // Choose a path and a name
            (double)mu, (double)xi, LEVEL, next);

    output_asc(localDepth, name);

    fprintf(stderr,"\n========================================\n");
    fprintf(stderr,"asc saved (t = %.2f)\n", t);

    next += 30.;   // next output
  }
}
*/



/* ------------------------- END ------------------------ */

/**
 * Final event — runs once at t = END.
 *
 * 1. Computes and exports the final localDepth as an ESRI ASCII raster
 * 2. Writes a simulation log (simulation_log.txt) with all run metadata:
 *      - Start / end timestamps
 *      - Wall-clock and CPU time
 *      - Input file paths
 *      - Rheological parameters (µ, ξ)
 *      - Simulation settings (LEVEL, END, LENGTH)
 *      - Total number of iterations
 */

event end (t = END)
{
  /* ---- final output ----*/
  char name[256];
  sprintf(name, "results/localDepth-mu%0.3f-xi%0.0f-L%d-t%06.0f.asc",
          (double)MU, (double)XI, LEVEL, t);
  output_asc(localDepth, name);

  fprintf(stderr,"\n========================================\n");
  fprintf(stderr,"Final output saved : %s \n", name);

  /* --- Compute elapsed times --- */
  time_t run_t1 = time(NULL);
  clock_t cpu_t1 = clock();

  double wall_s = difftime(run_t1, run_t0);
  double cpu_s  = (double)(cpu_t1 - cpu_t0)/CLOCKS_PER_SEC;

  /* --- Format  --- */
  char t0_str[64], t1_str[64];
  struct tm tm0_copy, tm1_copy;
  struct tm *tm0, *tm1;

  tm0 = localtime(&run_t0);
  tm0_copy = *tm0;
  tm1 = localtime(&run_t1);
  tm1_copy = *tm1;

  strftime(t0_str, sizeof(t0_str), "%Y-%m-%d %H:%M:%S", &tm0_copy);
  strftime(t1_str, sizeof(t1_str), "%Y-%m-%d %H:%M:%S", &tm1_copy);

    /* --- Append to simulation log --- */
    FILE *ftxt = fopen("simulation_log.txt", "a");
    if (ftxt) {
      fprintf(ftxt,
        "--------------------------------------------------\n"
        "Start            : %s\n"
        "End              : %s\n"
        "Wall (s)         : %.2f\n"
        "CPU  (s)         : %.2f\n"
        "Source           : %s\n"
        "Executable       : %s\n"
        "DEM              : %s\n"
        "RA               : %s\n"
        "DoD              : %s\n"
        "MU               : %g\n"
        "XI               : %g\n"
        "LEVEL            : %d\n"
        "END              : %g\n"
        "LENGTH           : %g\n"
        "Iterations       : %d\n"
        "--------------------------------------------------\n\n",
        t0_str, t1_str, wall_s, cpu_s,
        src_name, exe_name,
        dem_file, ra_file, dod_file,
        (double)MU, (double)XI, LEVEL, (double)END, (double)LENGTH,i
      );
      fclose(ftxt);
    }


  fprintf(stderr, "Simulation finished at t=%g\n", t);
  fprintf(stderr, "Runtime : %.0f s  (%.2f min)\n", wall_s, wall_s/60.);
}