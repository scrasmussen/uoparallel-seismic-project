////////////////////////////////////////////////////////////////////////////////
// sweep-tt-multistart.c - using VELOCITYBOX and FLOATBOX
// vim: set tabstop=2 softtabstop=2 shiftwidth=2 expandtab :
////////////////////////////////////////////////////////////////////////////////

/********************************************************************************/
/* Given a velocity field v[nx][ny][nz] for a set of points (i,j,k) (where	*/
/* 0 <= i < nx, 0 <= j < ny, 0 <= k < nz) layed out on a grid with delta unit	*/
/* distance, compute the minimum travel time, tt[nx][ny][nz][numstart], for all	*/
/* points to the numstart starting points.  The program is called as follows:	*/
/*										*/
/*	sweep-tt-multistart vfile fsfile startfile       			*/
/*										*/
// vfile is the velocity field file and has the .vbox format.
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

#include "../include/iovelocity.h"

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
//struct MODEL	model[MODELMAX][MODELMAX][MODELMAX];
struct START	start[STARTMAX];

struct VELOCITYBOX vbox; // stores JUST velocities
struct FLOATBOX ttboxes[STARTMAX]; // stores JUST travel times, one volume per starting point

int sweepXYZ();

int main(int argc, char* argv[]) {
  int		i, j, k, l, m, nx, ny, nz, oi, oj, ok, s;
  int		numradius, starsize, anychange, numstart, numsweeps=0;
  int		fsindex[FSRADIUSMAX];
  float		delta, delay;
  FILE		*vfile, *fsfile, *ttfile, *startfile;

  const char *velocity_model_file = argv[1];

  /* open velocity model file */
  printf( "Loading velocity model file: %s...", velocity_model_file ); fflush( stdout );
  if( !vboxloadbinary( &vbox, velocity_model_file ) ) {
  //if( !vboxloadtext( &vbox, velocity_model_file ) ) {
    printf( "Cannot open velocity model file: %s\n", velocity_model_file );
    exit(1);
  }
  nx = vbox.box.nx;
  ny = vbox.box.ny;
  nz = vbox.box.nz;
  printf( " done.\n" ); fflush( stdout );
  printf( "Velocity model dimensions: %d x %d x %d\n", nx, ny, nz );

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
  // initialize travel times for all starting points
  for( s = 0; s < numstart; s++ ) {
    // prepare travel time volumes
    boxalloc( &ttboxes[s], nx, ny, nz );
    boxsetall( ttboxes[s], INFINITY );

    // set the starting point to have a travel time of 0
    fscanf( startfile, "%i %i %i", &i, &j, &k );
    boxput( ttboxes[s], i, j, k, 0 );
    printf( "starting point %d: %d %d %d\n", s, i, j, k );
    start[s].i = i; start[s].j = j; start[s].k = k;
  }
  printf("Starting points read\n");
  
  /* sweep until no change in travel times occur */
  anychange = 1;
  while (anychange) {
    /* sweep forward on X */
    numsweeps++;
    anychange = 0;
    printf("sweep %d begin\n", numsweeps);

    for (s=0; s<numstart; s++) {
      changed[s] = 0;
      changed[s] += sweepXYZ(nx, ny, nz, s, 0, starsize-1);
      printf(">>> start %d: changed == %d\n", s, changed[s]);
    }
    for (s=0; s<numstart; s++) {
      anychange += changed[s];
    }
    printf("sweep %d finished: anychange = %d\n", numsweeps, anychange);

    // temporary: break after one sweep, record results, quit
    break;
  }

  /* TODO: Remove exit statement so output can complete. */
  //exit(0);

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
            i, j, k, boxget( ttboxes[s], i, j, k ), 0, 0, 0 );
        }
      }
    }
  }
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
          delay = fs[l].d * (boxget(vbox.box, i, j, k) + boxget(vbox.box, oi, oj, ok)) / 2.0;
          /* update travel times for all starting points */
          /* if (i,j,k) is starting point, then skip */
          if ((i == start[s].i) && (j == start[s].j) && (k == start[s].k)) {
            continue;
          }
          tt = boxget( ttboxes[s], i, j, k );
          tto = boxget( ttboxes[s], oi, oj, ok );
          /* if offset point has infinity travel time, then update */
          if ((tt == INFINITY) && (tto == INFINITY)) {
            continue;
          }
          if ((tt != INFINITY) && (tto == INFINITY)) {
            boxput( ttboxes[s], oi, oj, ok, delay + tt );
            change += 1;
            continue;
          }
          if ((tt == INFINITY) && (tto != INFINITY)) {
            boxput( ttboxes[s], i, j, k, delay + tto );
            change += 1;
            continue;
          }
          if ((tt != INFINITY) && (tto != INFINITY)) {
            /* if a shorter travel time through (oi,oj,ok), update (i,j,k) */
            if ((delay + tto) < tt) {
              boxput( ttboxes[s], i, j, k, delay + tto );
              change += 1;
            }
            /* if a shorter travel time through (i,j,k), update (oi,oj,ok) */
            else if ((delay + tt) < tto) {
              boxput( ttboxes[s], oi, oj, ok, delay + tt );
              change += 1;
            }
          }
        }
      }
    }
  }
  return(change);

} /* end sweepXYZ */ 
