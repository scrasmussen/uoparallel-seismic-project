/********************************************************************************/
/* Given a velocity field v[nx][ny][nz] for a set of points (i,j,k) (where	*/
/* 0 <= i < nx, 0 <= j < ny, 0 <= k < nz) layed out on a grid with delta unit	*/
/* distance, compute the minimum travel time, tt[nx][ny][nz][numstart], for all	*/
/* points to the numstart starting points.  The program is called as follows:	*/
/*										*/
/*	  wave-multistart vfile fsfile startfile delta				*/
/*										*/
/* vfile is the velocity field file and has the format:				*/
/*										*/
/*	  nx ny nz								*/
/*	  v[i][j][k] for every point (i,j,k) in row-major order			*/
/*										*/
/* fsfile is the forward star offset file and has the format:			*/
/*										*/
/*	  starsize								*/
/*	  oi oj ok for every forward star offset (oi,oj,ok)			*/
/*										*/
/* startfile contains starting points and has the format:			*/
/*										*/
/*	  numstart								*/
/*	  si sj sk for every starting point					*/
/*										*/
/* The program writes to "output.tt" the following:				*/
/*										*/
/*	  nx ny nz								*/
/*	  tt[i][j][k] for every point (i,j,k) in row-major order		*/
/*										*/
/* for every starting point.							*/
/* (Note, the program currently exits before this is done.)			*/
/********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#define	FSRADIUSMAX	7	/* maximum radius forward star */
#define	FSMAX		1190	/* maximum # of points in a forward star */
#define MODELMAX	500	/* maximum model dimension in X,Y,Z */
#define STARTMAX	111	/* maximum starting points */

#define fscanf_error(A, B, expected) { \
  int count = (A); \
  if (count == EOF) { \
    if (ferror(B)) { \
      perror("fscanf"); \
    } else { \
      fprintf(stderr, "Error: fscanf reached end of file, no matching characters, no matching failure\n"); \
    } \
    exit (99); \
  } else if (count != expected) { \
    fprintf(stderr, "Error: fscanf successfully matched and assigned %d input items, %d expected\n", count, expected); \
    exit (99); \
  } \
}

struct FS {		/* forward star offset */
  int	i, j, k;	/* point coordinates */
  float	d;		/* distance to star center (0,0,0) */
};

struct MODEL {		/* model point */
  float	v;		/* velocity */
  float	tt[STARTMAX];	/* travel time for starting points */
};

struct START {		/* starting point */
  int	i, j , k;	/* point coordinates */
};

int	nx, ny, nz;		/* model x,y,z dimensions */
int	changed[STARTMAX];	/* change flag table */

struct FS	fs[FSMAX];		/* forward star offset table */
struct MODEL	***model;		/* model */
struct START	start[STARTMAX];	/* starting point table */

// Forward declaration of methods
int waveBottomUp(int s, int starstart, int starstop);
int testconvergence(int s, int starstart, int starstop);
int waveTopDown(int s, int starstart, int starstop);
int starcompute(int i, int j, int k, int s, int starstart, int starstop);

void allocateModel(int nx, int ny, int nz) {
  int x, y;

  model = (struct MODEL***)(malloc(sizeof(struct MODEL**) * nx));
  for (x = 0; x < nx; x++) {
    model[x] = (struct MODEL**)(malloc(sizeof(struct MODEL*) * ny));
    for (y = 0; y < ny; y++) {
      model[x][y] = (struct MODEL*)(malloc(sizeof(struct MODEL) * nz));
    }
  }
} /* allocateModel */

long long timevalDiff(struct timeval *difference, struct timeval *end_time, struct timeval *start_time) {
  struct timeval temp_diff;

  if(difference==NULL) {
    difference=&temp_diff;
  }
  difference->tv_sec =end_time->tv_sec -start_time->tv_sec ;
  difference->tv_usec=end_time->tv_usec-start_time->tv_usec;

  while(difference->tv_usec<0) {
    difference->tv_usec+=1000000;
    difference->tv_sec -=1;
  }
  return 1000000LL * difference->tv_sec + difference->tv_usec;

} /* timevalDiff */

int main(int argc, char* argv[]) {
  int			i, j, k, s;
  int			radius, starsize, anychange, numstart, numwaves, maxthreads;
  int			fsindex[FSRADIUSMAX+1];
  float			delta;
  struct timeval	begintime, endtime, difference;
  FILE			*vfile, *fsfile, *ttfile, *startfile;

  /* open velocity model file */
  vfile = fopen(argv[1],"r");
  if(vfile == NULL) {
    printf("Can not open velocity model file: %s\n", argv[1]);
    exit(1);
  }
  printf("Velocity model file: %s\n", argv[1]);
  
  /* open forward star offset file */
  fsfile = fopen(argv[2],"r");
  if(fsfile == NULL) {
    printf("Can not open forward star offset file: %s\n", argv[2]);
    exit(1);
  }
  printf("Forward star offset file: %s\n", argv[2]);

  /* open file with starting points */
  startfile = fopen(argv[3],"r");
  if(startfile == NULL) {
    printf("Can not open starting points file: %s\n", argv[4]);
    exit(1);
  }
  printf("Starting points file: %s\n", argv[3]);

  /* get delta */
  delta = atof(argv[4]);
  printf("Delta: %f\n", delta);

  /* read velocity model */
  nx = 0; ny = 0; nz = 0;
  fscanf_error(fscanf(vfile, "%i %i %i", &nx, &ny, &nz), vfile, 3);
  printf("Velocity model dimensions: %i %i %i\n",nx, ny, nz);
  allocateModel(nx,ny,nz);
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      for (k=0; k<nz; k++) {
	fscanf_error(fscanf(vfile, "%f", &model[i][j][k].v), vfile, 1);
      }
    }
  }
  printf("Velocity data read\n");

  /* read forward star offsets */
  starsize = 0;
  fscanf_error(fscanf(fsfile, "%i", &starsize), fsfile, 1);
  printf("Forward star size: %d\n", starsize);

  for (i=0; i<FSRADIUSMAX; i++) {
    fsindex[i] = 0;
  }
  radius = 0;
  for (i=0; i<starsize; i++) {
    fscanf_error(fscanf(fsfile, "%i %i %i", &fs[i].i, &fs[i].j, &fs[i].k), fsfile, 3);
    fs[i].d = sqrt(fs[i].i*fs[i].i + fs[i].j*fs[i].j + fs[i].k*fs[i].k);
    /* printf("%i %i %i %f\n", fs[i].i, fs[i].j, fs[i].k, fs[i].d); */
    if (fs[i].d > (radius+1)) {
      fsindex[radius] = i;
      radius++;
    }
    fs[i].d = delta * fs[i].d; 
  }
  fsindex[radius] = starsize;
  printf("Forward star offsets read\n");
  for (i=0; i<FSRADIUSMAX+1; i++) {
    printf("radius: %d, fsindex[%d]: %d\n", radius, i, fsindex[i]);
  }

  /* read starting points */
  fscanf_error(fscanf(startfile, "%i", &numstart), startfile, 1);
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
    fscanf_error(fscanf(startfile, "%i %i %i", &start[s].i, &start[s].j, &start[s].k), startfile, 3);
    model[start[s].i][start[s].j][start[s].k].tt[s] = 0;
    printf("starting point %d: %d %d %d\n", s, start[s].i, start[s].j, start[s].k);
  }
  printf("Starting points read\n");
  
  gettimeofday(&begintime, NULL);

  /* preprocessing with a smaller 3-point forward star */
  printf("wavefront 3-point star preprocess begin\n");
  for (s=0; s<numstart; s++) {
    waveTopDown(s, 0, fsindex[3]);
    waveBottomUp(s, 0, fsindex[3]);
  }
  printf("\nwavefront 3-point star preprocess end\n");

  /* generate waves until no change in travel times occur */
  numwaves = 0;
  anychange = 1;

  maxthreads = omp_get_max_threads();
  printf("maxthreads = %d\n", maxthreads);

  while (anychange) {
    numwaves++;
    anychange = 0;
    printf("wavefront %d begin (top down)\n", numwaves);
    // #pragma omp parallel for private(s)
    for (s=0; s<numstart; s++) {
      changed[s] = 0;
      changed[s] += waveTopDown(s, 0, starsize-1);
      printf("\n>>> start %d: changed == %d\n", s, changed[s]);
    }
    for (s=0; s<numstart; s++) {
      anychange += changed[s];
    }
    printf("wavefront %d finished (top down): anychange = %d\n", numwaves, anychange);

    if (anychange == 0)
      continue;

    numwaves++;
    anychange = 0;
    printf("wavefront %d begin (bottom up)\n", numwaves);
    // #pragma omp parallel for private(s)
    for (s=0; s<numstart; s++) {
      changed[s] = 0;
      changed[s] += waveBottomUp(s, 0, starsize-1);
      printf("\n>>> start %d: changed == %d\n", s, changed[s]);
    }
    for (s=0; s<numstart; s++) {
      anychange += changed[s];
    }
    printf("wavefront %d finished (bottom up): anychange = %d\n", numwaves, anychange);
  }
  
  gettimeofday(&endtime, NULL);
  printf("begin time: %ld.%ld\n", begintime.tv_sec, begintime.tv_usec);
  printf("end time: %ld.%ld\n", endtime.tv_sec, endtime.tv_usec);
  long long milliseconds = timevalDiff(&difference, &endtime, &begintime);
  double total_us = (difference.tv_sec * 1000000) + difference.tv_usec;
  double timeper = (total_us / numwaves) / 1000000.0;
  printf("Duration: %ld.%ld seconds, %f per wave\n", difference.tv_sec, difference.tv_usec, timeper);
  
  for (s=0; s<numstart; s++) {
    if (testconvergence(s, 0, starsize-1) != 1) {
      printf("start point %d: NOT CONVERGED\n", s);
    } else {
      printf("start point %d: CONVERGED\n", s);
    }
  }

  exit(0);	/* do not print out the travel times */

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

} /* main */

/*************************************************/
/* testconvergence() does a check on convergence */
/*************************************************/
int testconvergence(int s, int starstart, int starstop) {
  int	i, j, k, l, oi, oj, ok;
  int	converged = 1;
  float	delay, tt, tto, temptt, tempv;
  
#pragma omp parallel for shared(converged,s,starstart,starstop)
  for (i=0; i<nx; i++) {
	if (!converged) continue;
    for (j=0; j<ny; j++) {
	  if (!converged) break;
      for (k=0; k<nz; k++) {
	    if (!converged) break;
	temptt = model[i][j][k].tt[s];	/* keeps smallest trave time */
	tempv = model[i][j][k].v;	/* velocity at star center point */
	for (l=starstart; l<starstop; l++) {
	  /* find point in forward star based on offsets */
	  oi = i+fs[l].i; oj = j+fs[l].j; ok = k+fs[l].k;
	  /* if (oi,oj,ok) is outside the boundaries, then skip */
	  if ((oi < 0) || (oi > nx-1) ||
	      (oj < 0) || (oj > ny-1) ||
	      (ok < 0) || (ok > nz-1)) {
	    continue;
	  }
	  /* compute delay from (i,j,k) to (oi,oj,ok) with end point averaging */
	  delay = fs[l].d * (tempv + model[oi][oj][ok].v) / 2.0;
	  /* update the travel time */
	  tt = temptt;
	  tto = model[oi][oj][ok].tt[s];
	  /* if either travel times are infinity, not converged */
	  if ((tt == INFINITY) || (tto == INFINITY)) {
	    converged = 0;
	    continue;
	  }
	  /* if a shorter travel time through offset, then update current */
	  if ((delay + tto) < tt) {
	    converged = 0;
	  }
	  /* if a shorter travel time through current, then update offset */
	  else if ((delay + tt) < tto) {
	    converged = 0;
	  }
	}
      }
    }
  }
  return(converged);

} /* end testconvergence */

/***************************************************/
/* starcompute() does the forward star computation */
/***************************************************/
int starcompute(int i, int j, int k, int s, int starstart, int starstop) {
  /* i,j,k     : star center point */
  /* s	       : starting point index */
  /* starstart : start starting index */
  /* starstop  : start starting index */
  int	     l, oi, oj, ok;
  int	     change = 0;
  float	       delay, tt, tto, temptt, tempv;
  
  temptt = model[i][j][k].tt[s]; /* keeps smallest trave time */
  tempv = model[i][j][k].v;	    /* velocity at star center point */
  for (l=starstart; l<starstop; l++) {
    /* find point in forward star based on offsets */
    oi = i+fs[l].i; oj = j+fs[l].j; ok = k+fs[l].k;
    /* if (oi,oj,ok) is outside the boundaries, then skip */
    if ((oi < 0) || (oi > nx-1) ||
	(oj < 0) || (oj > ny-1) ||
	(ok < 0) || (ok > nz-1)) {
      continue;
    }
    /* compute delay from (i,j,k) to (oi,oj,ok) with end point averaging */
    delay = fs[l].d * (tempv + model[oi][oj][ok].v) / 2.0;
    /* update the travel time */
    tt = temptt;
    tto = model[oi][oj][ok].tt[s];
    /* if both travel times are infinity, then skip */
    if ((tt == INFINITY) && (tto == INFINITY)) {
      continue;
    }
    /* if offset point has infinity travel time, then update it to current */
    if ((tt != INFINITY) && (tto == INFINITY)) {
      model[oi][oj][ok].tt[s] = delay + tt;
      change += 1;
      continue;
    }
    /* if current point has infinity travel time, then update it to offset*/
    if ((tt == INFINITY) && (tto != INFINITY)) {
      temptt = delay + tto;
      change += 1;
      continue;
    }
    if ((tt != INFINITY) && (tto != INFINITY)) {
      /* if a shorter travel time through offset, then update current */
      if ((delay + tto) < tt) {
	temptt = delay + tto;
	change += 1;
      }
      /* if a shorter travel time through current, then update offset */
      else if ((delay + tt) < tto) {
	model[oi][oj][ok].tt[s] = delay + tt;
	change += 1;
      }
    }
  }
  model[i][j][k].tt[s] = temptt;
  return(change);

} /* starcompute */

/*********************************************************************/
/* waveTopDown() uses a pseudo 2D wave to generate points in a plane */
/* for processing as the plane moves from the top down.		     */
/*********************************************************************/
int waveTopDown(int s, int starstart, int starstop) {
  int	     i, j, k, r, rmax, starti, startj;
  int	     change = 0;
  
  if (nx > ny) {
    rmax = nx;
  } else {
    rmax = ny;
  }
  starti = start[s].i;
  startj = start[s].j;
#pragma omp parallel for schedule(dynamic,1) private(i,j,k,r) reduction(+:change)
  for (k=nz-1; k>-1; k--) { /* loop through Z planes */
    //int tid = omp_get_thread_num();
    //printf("%d start: %d, k: %d\n", tid, s, k);
    printf("."); fflush(stdout);
    change += starcompute(starti, startj, k, s, starstart, starstop);
    // #pragma omp parallel for schedule(dynamic,1) private(i,j,r) reduction(+:change)
    for (r=1; r<rmax; r++) {
      j = startj-r;
      for (i=(starti-r); i<(starti+r); i++) {
	/* if (i,j,k) is outside the boundaries, then skip */
	if ((i < 0) || (i > nx-1) || (j < 0) || (j > ny-1)) {
	  continue;
	}
	change += starcompute(i, j, k, s, starstart, starstop);
      }
      i = starti+r;
      for (j=(startj-r); j<(startj+r); j++) {
	/* if (i,j,k) is outside the boundaries, then skip */
	if ((i < 0) || (i > nx-1) || (j < 0) || (j > ny-1)) {
	  continue;
	}
	change += starcompute(i, j, k, s, starstart, starstop);
      }
      j = startj+r;
      for (i=(starti+r); i>(starti-r); i--) {
	/* if (i,j,k) is outside the boundaries, then skip */
	if ((i < 0) || (i > nx-1) || (j < 0) || (j > ny-1)) {
	  continue;
	}
	change += starcompute(i, j, k, s, starstart, starstop);
      }
      i = starti-r;
      for (j=(startj+r); j>(startj-r); j--) {
	/* if (i,j,k) is outside the boundaries, then skip */
	if ((i < 0) || (i > nx-1) || (j < 0) || (j > ny-1)) {
	  continue;
	}
	change += starcompute(i, j, k, s, starstart, starstop);
      }
    }
  }
  return(change);

} /* end waveTopDown */

/**********************************************************************/
/* waveBottomUp() uses a pseudo 2D wave to generate points in a plane */
/* for processing as the plane moves from the bottom up.	      */
/**********************************************************************/
int waveBottomUp(int s, int starstart, int starstop) {
  int	     i, j, k, r, rmax, starti, startj;
  int	     change = 0;
  
  if (nx > ny) {
    rmax = nx;
  } else {
    rmax = ny;
  }
  starti = start[s].i;
  startj = start[s].j;
#pragma omp parallel for schedule(dynamic,1) private(i,j,k,r) reduction(+:change)
  for (k=0; k<nz; k++) {
    //int tid = omp_get_thread_num();
    //printf("%d start: %d, k: %d\n", tid, s, k);
    printf("."); fflush(stdout);
    change += starcompute(starti, startj, k, s, starstart, starstop);
    // #pragma omp parallel for schedule(dynamic,1) private(i,j,r) reduction(+:change)
    for (r=1; r<rmax; r++) {
      j = startj-r;
      for (i=(starti-r); i<(starti+r); i++) {
	/* if (i,j,k) is outside the boundaries, then skip */
	if ((i < 0) || (i > nx-1) || (j < 0) || (j > ny-1)) {
	  continue;
	}
	change += starcompute(i, j, k, s, starstart, starstop);
      }
      i = starti+r;
      for (j=(startj-r); j<(startj+r); j++) {
	/* if (i,j,k) is outside the boundaries, then skip */
	if ((i < 0) || (i > nx-1) || (j < 0) || (j > ny-1)) {
	  continue;
	}
	change += starcompute(i, j, k, s, starstart, starstop);
      }
      j = startj+r;
      for (i=(starti+r); i>(starti-r); i--) {
	/* if (i,j,k) is outside the boundaries, then skip */
	if ((i < 0) || (i > nx-1) || (j < 0) || (j > ny-1)) {
	  continue;
	}
	change += starcompute(i, j, k, s, starstart, starstop);
      }
      i = starti-r;
      for (j=(startj+r); j>(startj-r); j--) {
	/* if (i,j,k) is outside the boundaries, then skip */
	if ((i < 0) || (i > nx-1) || (j < 0) || (j > ny-1)) {
	  continue;
	}
	change += starcompute(i, j, k, s, starstart, starstop);
      }
    }
  }
  return(change);
} /* end waveBottomUp */
