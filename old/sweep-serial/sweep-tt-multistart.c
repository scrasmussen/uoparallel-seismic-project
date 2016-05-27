/********************************************************************************/
/* Given a velocity field v[nx][ny][nz] for a set of points (i,j,k) (where	*/
/* 0 <= i < nx, 0 <= j < ny, 0 <= k < nz) layed out on a grid with delta unit	*/
/* distance, compute the minimum travel time, tt[nx][ny][nz][numstart], for all	*/
/* points to the numstart starting points.  The program is called as follows:	*/
/*										*/
/*	sweep-tt-multistart vfile fsfile startfile       			*/
/*										*/
/* vfile is the velocity field file and has the format:				*/
/*										*/
/*	nx ny nz								*/
/*	v[i][j][k] for every point (i,j,k) in row-major order			*/
/*										*/
/* fsfile is the forward star offset file and has the format:			*/
/*										*/
/*	starsize								*/
/*	oi oj ok for every forward star offset (oi,oj,ok)			*/
/*										*/
/* startfile contains starting points and has the format:			*/
/*										*/
/*	numstart								*/
/*	si sj sk for every starting point					*/
/*										*/
/* The program writes to "output.tt" the following:				*/
/*										*/
/*	nx ny nz								*/
/*	tt[i][j][k] for every point (i,j,k) in row-major order			*/
/*										*/
/* for every starting point.							*/
/* (Note, the program currently exits before this is done.)			*/
/********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define	FSRADIUSMAX	7	/* maximum radius forward star */
#define	FSMAX		818	/* maximum # of points in a forward star */
#define MODELMAX	250	/* maximum model dimension in X,Y,Z */
#define STARTMAX	12	/* maximum starting points */

struct FS {			/* forward start offset */
  int		i, j, k;	/* point coordinates */
  float		d;		/* distance to star center (0,0,0)*/
};

struct MODEL {			/* model point */
  float		v;		/* velocity */
  float		tt[STARTMAX];	/* travel time for starting points */
};

struct START {			/* starting point */
  int		i, j , k;	/* point coordinates */
};

int		changed[STARTMAX];

struct FS	fs[FSMAX];
struct MODEL	model[MODELMAX][MODELMAX][MODELMAX];
struct START	start[STARTMAX];

int sweepXYZ();

int main(int argc, char* argv[]) {
  int		i, j, k, l, m, nx, ny, nz, oi, oj, ok, s;
  int		numradius, starsize, anychange, numstart, numsweeps=0;
  int		fsindex[FSRADIUSMAX];
  float		delta, delay;
  FILE		*vfile, *fsfile, *ttfile, *startfile;

  /**
   * TODO:
   *   1. Initialize MPI library
   *   2. Get number of processes
   *   3. Get rank of current process
   *
   *   Note: Only 4 starting points are used.  Code should work
   *         with 1 or 4 processors (should check).  Usually should
   *         write an MPI code to work with an arbitrary number of
   *         processors unless algorithm doesn't allow it.
   */

  /**
   * TODO: Make sure there aren't any race conditions.  If there are, fix
   *       with barriers.
   */

  /**
   * TODO: Only rank 0 should do I/O
   */

  /* open velocity model file */
  vfile = fopen(argv[1],"r");
  if(vfile == NULL) {
    printf("Cannot open velocity model file: %s\n", argv[1]);
    exit(1);
  }
  printf("Velocity model file: %s\n", argv[1]);

  /* open forward star offset file */
  fsfile = fopen(argv[2],"r");
  if(fsfile == NULL) {
    printf("Cannot open forward star offset file: %s\n", argv[2]);
    exit(1);
  }
  printf("Forward star offset file: %s\n", argv[2]);

  /* open file with starting points */
  startfile = fopen(argv[3],"r");
  if(startfile == NULL) {
    printf("Cannot open starting points file: %s\n", argv[4]);
    exit(1);
  }
  printf("Starting points file: %s\n", argv[3]);

  /* get delta */
  delta = 10.0;
  printf("Delta: %f\n", delta);

  /* read velocity model (modified to read original input for Fortran) */
  nx = 241; ny = 241; nz = 51;
  printf("Velocity model dimensions: %i %i %i\n",nx, ny, nz);
  int ir, jr, kr;   /* read indices */
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
         fscanf(vfile, "%d,%d,%d,%f", &ir, &jr, &kr, &model[i][j][k].v);
         if (ir != i+1 || jr != j+1 || kr != k+1) {
            printf("ERROR: index error reading velocity model "
                   "(%d,%d,%d) (%d,%d,%d)\n", i,j,k, ir,jr,kr);
            exit(1);
         }
      }
    }
  }
  printf("Velocity data read\n");

  /* read forward star offsets */
  starsize = 0;
  fscanf(fsfile, "%i", &starsize);
  printf("Forward star size: %d\n", starsize);

  for (i=0; i<FSRADIUSMAX; i++) {
    fsindex[i] = 0;
  }
  numradius = 0;
  for (i=0; i<starsize; i++) {
    fscanf(fsfile, "%i %i %i", &fs[i].i, &fs[i].j, &fs[i].k);
    fs[i].d = sqrt(fs[i].i*fs[i].i + fs[i].j*fs[i].j + fs[i].k*fs[i].k);
    if ((numradius+1) < fs[i].d) {
      fsindex[numradius] = i;
      numradius++;
    }
    fs[i].d = delta * fs[i].d; 
  }
  printf("Forward star offsets read\n");
  for (i=0; i<FSRADIUSMAX; i++) {
    printf("numradius: %d, fsindex[%d]: %d\n", numradius, i, fsindex[i]);
  }

  /* read starting points */
  fscanf(startfile, "%i", &numstart);
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
	for (s=0; s<numstart; s++) {
	  model[i][j][k].tt[s] = INFINITY;
	}
      }
    }
  }
  for (s=0; s<numstart; s++) {
    fscanf(startfile, "%i %i %i", &start[s].i, &start[s].j, &start[s].k);
    model[start[s].i][start[s].j][start[s].k].tt[s] = 0;
    printf("starting point %d: %d %d %d\n", s, start[s].i, start[s].j, start[s].k);
  }
  printf("Starting points read\n");
  
  /**
   * TODO: Rank 0 should broadcast velocity model and forward star data.
   */

  /**
   * TODO: Rank 0 should assign starting points to other ranks (in this
   *       case, 1 starting point per process).  Other processes should
   *       read its starting point.
   */

  /* sweep until no change in travel times occur */
  anychange = 1;
  while (anychange) {
    /* sweep forward on X */
    numsweeps++;
    anychange = 0;
    printf("sweep %d begin\n", numsweeps);

    /**
     * TODO: A rank should only do its starting point, not all of them.
     *
     * Note: You may want to use also write out the rank of the process
     *       that is reporting progress.
     */
    for (s=0; s<numstart; s++) {
      changed[s] = 0;
      changed[s] += sweepXYZ(nx, ny, nz, s, 0, starsize-1);
      printf(">>> start %d: changed == %d\n", s, changed[s]);
    }
    for (s=0; s<numstart; s++) {
      anychange += changed[s];
    }
    printf("sweep %d finished: anychange = %d\n", numsweeps, anychange);
  }

  /**
   * TODO: Rank 0 should gather results from other processes or it's OK to have
   *       each process write its own output.  Make sure that the filenames
   *       don't collide.  At least think about how rank 0 could gather output
   *       results.
   */

  /* TODO: Remove exit statement so output can complete. */
  exit(0);

  /* print travel times */
  ttfile = fopen("output.tt","w");
  if(ttfile == NULL) {
    printf("Can not open travel time output file: %s\n", "output.tt");
    exit(1);
  }
  fprintf(ttfile, "%d %d %d\n", nx, ny, nz);
  for (s=0; s<numstart; s++) {
    fprintf(ttfile, "starting point: %d\n", s);
    for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
	for (k=0; k<nz; k++) {
	  /* use %g for doubles */
	  fprintf(ttfile, "travel time for (%d,%d,%d): %f %d %d %d\n",
		  i, j, k, model[i][j][k].tt[s], 0, 0, 0);
	}
      }
    }
  }

  /**
   * TODO: Shutdown MPI library
   */

} /* main */


int sweepXYZ(int nx, int ny, int nz, int s, int starstart, int starstop) {
  int	i, j, k, l, oi, oj, ok;
  int	change = 0;
  float	delay = 0.0, tt = 0.0, tto = 0.0;
  
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
	for (l=starstart; l<starstop; l++) {
	  /* find point in forward star based on offsets */
	  oi = i+fs[l].i; oj = j+fs[l].j; ok = k+fs[l].k;
	  /* if (oi,oj,ok) is outside the boundaries, then skip */
	  if ((oi < 0) || (oi > nx-1)
	      || (oj < 0) || (oj > ny-1)
	      || (ok < 0) || (ok > nz-1)) {
	    continue;
	  }
	  /* compute delay from (i,j,k) to (oi,oj,ok) with end point average */
	  delay = fs[l].d * (model[i][j][k].v + model[oi][oj][ok].v) / 2.0;
	  /* update travel times for all starting points */
	  /* if (i,j,k) is starting point, then skip */
	  if ((i == start[s].i) && (j == start[s].j) && (k == start[s].k)) {
	    continue;
	  }
	  tt = model[i][j][k].tt[s];
	  tto = model[oi][oj][ok].tt[s];
	  /* if offset point has infinity travel time, then update */
	  if ((tt == INFINITY) && (tto == INFINITY)) {
	    continue;
	  }
	  if ((tt != INFINITY) && (tto == INFINITY)) {
	    model[oi][oj][ok].tt[s] = delay + tt;
	    change += 1;
	    continue;
	  }
	  if ((tt == INFINITY) && (tto != INFINITY)) {
	    model[i][j][k].tt[s] = delay + tto;
	    change += 1;
	    continue;
	  }
	  if ((tt != INFINITY) && (tto != INFINITY)) {
	    /* if a shorter travel time through (oi,oj,ok), update (i,j,k) */
	    if ((delay + tto) < tt) {
	      model[i][j][k].tt[s] = delay + tto;
	      change += 1;
	    }
	    /* if a shorter travel time through (i,j,k), update (oi,oj,ok) */
	    else if ((delay + tt) < tto) {
	      model[oi][oj][ok].tt[s] = delay + tt;
	      change += 1;
	    }
	  }
	}
      }
    }
  }
  return(change);

} /* end sweepXYZ */ 
