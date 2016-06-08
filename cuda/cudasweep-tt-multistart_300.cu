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

#include "iovelocity.h"
#include "timing.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define	FSRADIUSMAX	7	/* maximum radius forward star */
#define	FSMAX		818	/* maximum # of points in a forward star */
#define MODELMAX	250	/* maximum model dimension in X,Y,Z */
#define STARTMAX	4	/* maximum starting points */

#define GRIDX 256
#define GRIDY 256
#define GRIDZ 1
#define BLOCKX 1
#define BLOCKY 1
#define BLOCKZ 64

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
__constant__ struct FS dc_fs[FSMAX];

struct START	start[STARTMAX];
__constant__ struct START	dc_start[STARTMAX];

struct VELOCITYBOX vbox; // stores JUST velocities
__constant__ struct VELOCITYBOX dc_vbox;

struct FLOATBOX ttboxes[STARTMAX]; // stores JUST travel times, one volume per starting point
__constant__ struct FLOATBOX dc_ttboxes[STARTMAX];


void cudaRun(int, int);
__global__ void
cudaWorker(int, int, int, int, int, int, struct FS *, struct START *, struct VELOCITYBOX *, struct FLOATBOX *,long *);
__device__ int 
sweepXYZ(int, int, int, int, int, int, int, int, struct FS *, float *, float *);

int main(int argc, char* argv[]) {
  int		i, j, k, nx, ny, nz, s;
  int		numradius, starsize, numstart;
  int		fsindex[FSRADIUSMAX];
  float		delta;
  FILE		*fsfile, *ttfile, *startfile;

  const char *velocity_model_file = argv[1];

  /* open velocity model file */
  printf( "Loading velocity model file: %s...", velocity_model_file ); fflush( stdout );
  //if( !vboxloadbinary( &vbox, velocity_model_file ) ) {
  if( !vboxloadtext( &vbox, velocity_model_file ) ) {
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
  
  cudaSetDevice(0);
	cudaRun(numstart, starsize);

  // /* print travel times */
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

void cudaRun(
    int numstart, 
    int starsize
)
{
  struct FS     *pd_fs;
  struct START  *pd_start;
  struct VELOCITYBOX  *pd_vbox; // stores JUST velocities
  struct FLOATBOX     *pd_ttboxes; // stores JUST travel times, one volume per starting point
  int i, j, nx = vbox.box.nx, ny = vbox.box.ny, nz = vbox.box.nz;
  cudaError_t err;
  
  //copy fs to device
  cudaMemcpyToSymbol(dc_fs, fs, sizeof(fs));
  size_t fssize = sizeof(struct FS)*FSMAX;
  err = cudaMalloc( (void**)&pd_fs, fssize );
  if(err != cudaSuccess)
    printf("fs malloc error\n");
  err = cudaMemcpy( pd_fs, fs, fssize, cudaMemcpyHostToDevice );
  if(err != cudaSuccess)
    printf("fs copy error: %d\n", (int)fssize);
  printf("1\n");
  
  //copy start points to device
  cudaMemcpyToSymbol(dc_start, start, sizeof(start));
  size_t startsize = sizeof(struct START)*STARTMAX;
  err = cudaMalloc( (void**)&pd_start, startsize );
  if(err != cudaSuccess)
    printf("start malloc error\n");
  err = cudaMemcpy( pd_start, start, startsize, cudaMemcpyHostToDevice );
  if(err != cudaSuccess)
    printf("start copy error\n");
  printf("2\n");
  
  //copy velosity box to device
  size_t vboxsize = sizeof(struct VELOCITYBOX);
  size_t flatbytes = (size_t)nx * ny * nz * sizeof(float);
  float *pd_vboxflat;
  err = cudaMalloc( (void **)&pd_vbox, vboxsize );
  if(err != cudaSuccess)
    printf("vbox malloc error\n");

  err = cudaMalloc( (void **)&pd_vboxflat, flatbytes );
  if(err != cudaSuccess)
    printf("pd_vboxflat malloc error\n");
  
  struct VELOCITYBOX dummyvbox;
  memcpy( &dummyvbox, &vbox, sizeof(struct VELOCITYBOX) );
  dummyvbox.box.flat = pd_vboxflat;
  err = cudaMemcpy( dummyvbox.box.flat, vbox.box.flat, flatbytes, cudaMemcpyHostToDevice );
  if(err != cudaSuccess)
    printf( "pd_vboxflat copy error\n" );
  err = cudaMemcpy( pd_vbox, &dummyvbox, vboxsize, cudaMemcpyHostToDevice );
  if(err != cudaSuccess)
    printf( "vbox copy error\n" );
	cudaMemcpyToSymbol(dc_vbox, &dummyvbox, sizeof(dummyvbox));
  printf( "3\n" );
  
  //copy travel time boxes to device
  size_t boxessize = sizeof(struct FLOATBOX)*STARTMAX;
  err = cudaMalloc( (void **)&pd_ttboxes, boxessize );
  if(err != cudaSuccess)
    printf("boxes malloc error\n");
  
	struct FLOATBOX dummybox[STARTMAX];
  for(i=0; i<STARTMAX; i++){
    float *pd_boxflat;
    err = cudaMalloc( (void **)&pd_boxflat, flatbytes );
    if(err != cudaSuccess)
      printf("pd_boxflat malloc error\n");
    
    memcpy(dummybox+i, ttboxes+i, sizeof(struct FLOATBOX));
    dummybox[i].flat = pd_boxflat;
    
    err = cudaMemcpy( dummybox[i].flat, ttboxes[i].flat, flatbytes, cudaMemcpyHostToDevice );
    if(err != cudaSuccess)
      printf( "boxflat %d copy error\n", i );
  }
	err = cudaMemcpy( pd_ttboxes, dummybox, sizeof(struct FLOATBOX) * STARTMAX, cudaMemcpyHostToDevice );
  if(err != cudaSuccess)
    printf( "box %d copy error\n", i );
	cudaMemcpyToSymbol(dc_ttboxes, dummybox, sizeof(dummybox));
	
  printf("4\n");
  
	
  const int tNum = GRIDX * BLOCKX * GRIDY * BLOCKY * GRIDZ * BLOCKZ ;
  //const int blkNum = GRIDX * GRIDY * GRIDZ;
  //const int blkSize = BLOCKX * BLOCKY * BLOCKZ;
  long *pd_anychange, *anychange;
  double sweepTime = 0, dataTransTime = 0;
  
  err = cudaMalloc(&pd_anychange, sizeof(long) * tNum);
  if(err != cudaSuccess)
    printf( "pd_anychange malloc error\n");
    
  anychange = (long*)malloc(sizeof(long) * tNum);
  printf("5\n");
  
  int nDevices;
  cudaGetDeviceCount(&nDevices);
  printf("device: %d\n", nDevices);
  for (int i = 0; i < nDevices; i++) {
      cudaDeviceProp prop;
      cudaGetDeviceProperties(&prop, i);
      printf("Device Number: %d\n", i);
      printf("  Device name: %s\n", prop.name);
      printf("  Memory Clock Rate (KHz): %d\n",
             prop.memoryClockRate);
      printf("  Memory Bus Width (bits): %d\n",
             prop.memoryBusWidth);
      printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
             2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
    }
  
  
  for(i=0; i<numstart; i++){
    long sweepNum = 0, changeSum = 1;
    while (changeSum) {
      sweepTime = 0; dataTransTime = 0; changeSum = 0;
      sweepNum++;
      err = cudaMemset(pd_anychange, 0, sizeof(long) * tNum);
      if(err != cudaSuccess)
        printf( "pd_anychange memset error\n");
      
			reset_and_start_timer();
      dim3 gridDim(GRIDX,GRIDY,GRIDZ);
      dim3 blockDim(BLOCKX,BLOCKY,BLOCKZ);
      cudaWorker<<<gridDim,blockDim>>>(
        nx, ny, nz, 
        i, 
        0, starsize-1, //Note: change the range to the original starsize only reduce 5ms time.
        dc_fs, 
        dc_start, 
        pd_vbox, 
        pd_ttboxes, 
        pd_anychange
      );
      cudaDeviceSynchronize();
      sweepTime = get_elapsed_msec();
      
      if(err != cudaSuccess)
        printf("  cudaGetLastError() returned %d: %s\n", err, cudaGetErrorString(err));
      
      reset_and_start_timer();
      err = cudaMemcpy( anychange, pd_anychange, sizeof(long) * tNum, cudaMemcpyDeviceToHost );
      if(err != cudaSuccess)
        printf( "anychange copy error: %d\n", err);
      dataTransTime = get_elapsed_msec();
      
      for(j = 0; j < tNum; j++){
				changeSum += anychange[j];
      }
      
      printf(" start point: %d, sweep %d: %d changes, sweep %f, data trans %f\n", 
        i, sweepNum, changeSum, sweepTime, dataTransTime);
    }
  }

	printf("6\n");
	
  for(i=0; i<STARTMAX; i++){
    struct FLOATBOX ttboxbuff;
		 err = cudaMemcpy( &ttboxbuff, pd_ttboxes+i, sizeof(struct FLOATBOX), cudaMemcpyDeviceToHost );
    if(err != cudaSuccess)
      printf( "box %d copy error\n", i );
    
		err = cudaMemcpy( ttboxes[i].flat, ttboxbuff.flat, flatbytes, cudaMemcpyDeviceToHost );
    if(err != cudaSuccess)
      printf( "boxflat %d copy error\n", i );
  }
	printf("7\n");
  
  cudaFree(pd_fs);
  cudaFree(pd_start);
  cudaFree(pd_vbox);
  cudaFree(pd_vboxflat);
  cudaFree(pd_ttboxes);
  cudaFree(dummybox[i].flat);
  cudaFree(pd_anychange);
  free(anychange);
}

__global__ 
void cudaWorker(
    int d_nx, int d_ny, int d_nz,
    int d_s, 
    int d_starstart, int d_starend,
    struct FS *pd_fs,
    struct START *pd_start,
    struct VELOCITYBOX *pd_vbox,
    struct FLOATBOX *pd_ttboxes,
    long *pd_anychange
)
{
  //int d_blktid = threadIdx.z + threadIdx.y * blockDim.z + threadIdx.x * blockDim.z * blockDim.y;
  int d_blkid = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
	int d_glbtid = d_blkid * (blockDim.x * blockDim.y * blockDim.z)
								+ (threadIdx.z * (blockDim.x * blockDim.y))
								+ (threadIdx.y * blockDim.x) + threadIdx.x;
  //int blkSize = blockDim.x * blockDim.y * blockDim.z;
	
	pd_anychange[d_glbtid] = sweepXYZ(
    d_nx, d_ny, d_nz,
    dc_start[d_s].i, dc_start[d_s].j, dc_start[d_s].k, 
    d_starstart, d_starend, 
    dc_fs, 
    pd_vbox->box.flat,
    pd_ttboxes[d_s].flat
  );
}

__device__ int 
sweepXYZ(
    int nx, int ny, int nz, 
    int startx, int starty, int startz,
    int starstart, int starstop,
    struct FS *fs,
    float *vboxflat,
    float *ttboxflat
) 
{
  int	i, j, k, l, oi, oj, ok;
  float	delay = 0.0, tt = 0.0, tto = 0.0, ttd = 0.0, ttod = 0.0;
  int sx = nz * ny;
  int d_blktid = threadIdx.z + threadIdx.y * blockDim.z + threadIdx.x * blockDim.z * blockDim.y;
  __shared__ int change;
  if(d_blktid == 0)
    change = 0;
  __syncthreads();
  
	i = blockIdx.x * blockDim.x + threadIdx.x;
	j = blockIdx.y * blockDim.y + threadIdx.y;
	k = blockIdx.z * blockDim.z + threadIdx.z;
  
	
	if(i >= nx || j >= ny || k >= nz)
		return 0;
	
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
    int iIdx = k+nz*j+i*sx; int oIdx = ok+nz*oj+oi*sx;
		delay = dc_fs[l].d * (vboxflat[iIdx] + vboxflat[oIdx]) / 2.0;
		tt = ttboxflat[iIdx];
		tto = ttboxflat[oIdx];

    /* if a shorter travel time through (oi,oj,ok), update (i,j,k) */ 
    if ((delay + tto) < tt) {
      ttboxflat[iIdx] = delay + tto;
      if(change == 0)
        change = 1;
    }
    /* if a shorter travel time through (i,j,k), update (oi,oj,ok) */
    else if ((delay + tt) < tto) {
      ttboxflat[oIdx] = delay + tt;
      if(change == 0)
        change = 1;
    }
  }
  return(change);

} /* end sweepXYZ */ 