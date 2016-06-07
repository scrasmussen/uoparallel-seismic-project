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
#include <unistd.h>

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
#define DEVNUM 3
const int starSplit[4] = {0, 330, 580, 818};

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
struct VELOCITYBOX vbox; // stores JUST velocities
struct FLOATBOX ttboxes[STARTMAX]; // stores JUST travel times, one volume per starting point

//The cuda wrapper for 
void cudaRun(int, int);
__global__ 
void cudaWorker(
    int d_nx, int d_ny, int d_nz,
    int d_starstart, int d_starend,
    struct FS *pd_fs,
    float *pd_vboxflat,
    float *pd_ttboxflat,
    int *pd_anychange
);
__device__ int 
sweepXYZ(
    int nx, int ny, int nz, 
    int starstart, int starstop,
    struct FS *fs,
    float *vboxflat,
    float *ttboxflat
);
__global__ void 
MergeBoxes(
    float *pd_ttboxflat, 
    float *pd_ttbflatMergeCache, 
    int devNum, 
    int nCells
);

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
  
	cudaRun(numstart, starsize);

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

void cudaRun(
    int numstart, 
    int starsize
)
{
  //constants
  const int blkNum = GRIDX * GRIDY * GRIDZ;
  //const int blkSize = BLOCKX * BLOCKY * BLOCKZ;
  //const int tNum = blkNum * blkSize;

  //host variables
  int anychange[DEVNUM][blkNum];
  int i, j, nx = vbox.box.nx, ny = vbox.box.ny, nz = vbox.box.nz, devIdx = 0, devNum = DEVNUM;
  int nCells = nx * ny * nz;
  size_t flatbytes = (size_t)nCells * sizeof(float);
  cudaError_t err;
  
  //Cuda variables
  dim3 gridDim(GRIDX,GRIDY,GRIDZ);
  dim3 blockDim(BLOCKX,BLOCKY,BLOCKZ);
  float *ppd_vboxflat[DEVNUM];
  float *pppd_ttboxflat[DEVNUM][STARTMAX];
  float *ppd_ttbflatMergeCache[DEVNUM];
  int *ppd_anychange[DEVNUM];
  cudaStream_t streams[DEVNUM];
  
  //allocate device memory
  for(devIdx = 0; devIdx < devNum; devIdx++){
    cudaSetDevice(devIdx);
    err = cudaMalloc((void **)&ppd_vboxflat[devIdx], flatbytes);
    if(err != cudaSuccess)
      printf("ppd_vboxflat malloc error\n");
    err = cudaMalloc((void **)&ppd_anychange[devIdx], sizeof(int) * blkNum);
    if(err != cudaSuccess)
      printf( "ppd_anychange malloc error\n");
    for(i=0; i<STARTMAX; i++){
      err = cudaMalloc((void **)&pppd_ttboxflat[devIdx][i], flatbytes);
      if(err != cudaSuccess)
        printf("pppd_ttboxflat malloc error\n");
    }
    err = cudaMalloc((void **)&ppd_ttbflatMergeCache[devIdx], flatbytes*(devNum-1));
    if(err != cudaSuccess)
      printf( "ppd_ttbflatMergeCache malloc error\n");
  }

  //lock up host memory for async transfer
  cudaHostRegister(fs, sizeof(fs), cudaHostRegisterDefault);
  cudaHostRegister(start, sizeof(start), cudaHostRegisterDefault);
  cudaHostRegister(vbox.box.flat, flatbytes, cudaHostRegisterDefault);
  for(i=0; i<STARTMAX; i++)
    cudaHostRegister(ttboxes[i].flat, flatbytes, cudaHostRegisterDefault);
  
  //async copy memory from host to device
  for(devIdx = 0; devIdx < devNum; devIdx++){
    cudaSetDevice(devIdx);
    cudaStreamCreate(&streams[devIdx]);
    
    //copy fs to device
    err = cudaMemcpyToSymbolAsync(dc_fs, fs, sizeof(fs), 0, cudaMemcpyHostToDevice, streams[devIdx]);
    if(err != cudaSuccess)
      printf("dc_fs copy error\n");
    printf("1\n");
    
    //copy velosity box to device
    err = cudaMemcpyAsync(ppd_vboxflat[devIdx], vbox.box.flat, flatbytes, cudaMemcpyHostToDevice, streams[devIdx]);
    if(err != cudaSuccess)
      printf( "ppd_vboxflat copy error\n" );
    printf( "2\n" );
    
    //copy travel time boxes to device
    for(i=0; i<STARTMAX; i++){
      err = cudaMemcpyAsync(pppd_ttboxflat[devIdx][i], ttboxes[i].flat, flatbytes, cudaMemcpyHostToDevice, streams[devIdx]);
      if(err != cudaSuccess)
        printf( "pppd_ttboxflat %d copy error\n", i );
    }
    printf("3\n");
  }
  cudaStreamSynchronize(0);
  
  //run algorithm
  double tSweep = 0.0, tChangeTrans = 0.0, tSum = 0.0, tMerge = 0.0, tBoxTrans = 0.0, tTotal = 0.0;
  for(i=0; i<numstart; i++){
    int sweepNum = 0, changeSum = 1;
    while (changeSum) {//run until no changes
      changeSum = 0;
      sweepNum++;
      
      //run splited forward stars on different devices
      reset_and_start_timer();
      for(devIdx=0; devIdx<devNum; devIdx++){
        cudaSetDevice(devIdx);
        err = cudaMemset(ppd_anychange[devIdx], 0, sizeof(int) * blkNum);
        if(err != cudaSuccess)
          printf( "ppd_anychange memset error\n");
        
        cudaWorker<<<gridDim,blockDim>>>(
          nx, ny, nz, 
          starSplit[devIdx], starSplit[devIdx+1]-1, //Note: change the range to the original starsize only reduce 5ms time.
          NULL, //kernel is pulling foward stars from __constant__ so pass NULL
          ppd_vboxflat[devIdx], 
          pppd_ttboxflat[devIdx][i],
          ppd_anychange[devIdx]
        );
      }
      cudaStreamSynchronize(0); //sync all devices
      tSweep = get_elapsed_msec();
      
      if(err != cudaSuccess) //check error
        printf("  cudaGetLastError() returned %d: %s\n", err, cudaGetErrorString(err));
      
      //pull back and check changes
      reset_and_start_timer();
      for(devIdx=0; devIdx<devNum; devIdx++){
        cudaSetDevice(devIdx);
        err = cudaMemcpyAsync(anychange[devIdx], ppd_anychange[devIdx], sizeof(int) * blkNum, cudaMemcpyDeviceToHost, streams[devIdx]);
        if(err != cudaSuccess)
          printf("anychange copy error: %d\n", err);
      }
      cudaStreamSynchronize(0); //sync all devices
      tChangeTrans = get_elapsed_msec();
      
      reset_and_start_timer();
      for(devIdx=0; devIdx<devNum; devIdx++)
        for(j = 0; j < blkNum; j++)
          changeSum += anychange[devIdx][j];
      tSum = get_elapsed_msec();
      
      //sync travel time from all devices and merge them
      reset_and_start_timer();
      printf("nCells = %d\n", nCells);
      float *pCacheBasePtr;
      int buffered, devFrom;
      for(devIdx=0; devIdx<devNum; devIdx++){
        cudaSetDevice(devIdx);
        buffered = 0;
        for(devFrom=0; devFrom<devNum; devFrom++){
          if(devFrom != devIdx){
            pCacheBasePtr = ppd_ttbflatMergeCache[devIdx] + nCells*buffered;
            cudaMemcpyPeerAsync(pCacheBasePtr, devIdx, pppd_ttboxflat[devFrom][i], devFrom, flatbytes, streams[devIdx]);
            buffered++;
          }
        }
      }    
      cudaStreamSynchronize(0); //sync all devices
      tBoxTrans = get_elapsed_msec();
      reset_and_start_timer();
      
      for(devIdx=0; devIdx<devNum; devIdx++){
        cudaSetDevice(devIdx);
        MergeBoxes<<<gridDim,blockDim,0,streams[devIdx]>>>(
          pppd_ttboxflat[devIdx][i], 
          ppd_ttbflatMergeCache[devIdx], 
          devNum, 
          nCells
        );
      }
      cudaStreamSynchronize(0);
      tMerge = get_elapsed_msec();
      
      //output statistics
      tTotal = tSweep + tChangeTrans + tSum + tBoxTrans + tMerge;
      printf(" start point: %d, sweep %d: %d changes, sweep %g, change trans %g\n\
sum %g, box trans %g, merg %g, total %g\n", 
        i, sweepNum, changeSum, tSweep, tChangeTrans, tSum, tBoxTrans, tMerge, tTotal);
    }

    devIdx = 0;
    cudaSetDevice(devIdx);
    err = cudaMemcpy(ttboxes[i].flat, pppd_ttboxflat[devIdx][i], flatbytes, cudaMemcpyDeviceToHost);
    if(err != cudaSuccess)
      printf( "pppd_ttboxflat %d copy error\n", i );
  }
	printf("6\n");
  
  for(devIdx=0; devIdx<devNum; devIdx++){
    cudaFree(ppd_vboxflat[devIdx]);
    cudaFree(ppd_anychange[devIdx]);
    for(i=0; i<STARTMAX; i++)
      err = cudaFree(pppd_ttboxflat[devIdx][i]);
    cudaFree(ppd_ttbflatMergeCache[devIdx]);
    cudaStreamDestroy(streams[devIdx]);
  }
  
}

__global__ 
void cudaWorker(
    int d_nx, int d_ny, int d_nz,
    int d_starstart, int d_starend,
    struct FS *pd_fs,
    float *pd_vboxflat,
    float *pd_ttboxflat,
    int *pd_anychange
)
{
  int d_blktid = threadIdx.z + threadIdx.y * blockDim.z + threadIdx.x * blockDim.z * blockDim.y;
  int d_blkid = blockIdx.z + blockIdx.y * gridDim.z + blockIdx.x * gridDim.z * gridDim.y;
	// int d_glbtid = d_blkid * (blockDim.x * blockDim.y * blockDim.z)
								// + (threadIdx.z * (blockDim.x * blockDim.y))
								// + (threadIdx.y * blockDim.x) + threadIdx.x;
  const int blkSize = BLOCKX*BLOCKY*BLOCKZ;
	__shared__ int blkChange[blkSize];
  
	blkChange[d_blktid] = sweepXYZ(
    d_nx, d_ny, d_nz,
    d_starstart, d_starend, 
    dc_fs, 
    pd_vboxflat,
    pd_ttboxflat
  );
  __syncthreads();

  //reduction
  for (int s=blkSize/2; s>0; s>>=1) {
    if (d_blktid < s)
      blkChange[d_blktid] += blkChange[d_blktid + s];
    __syncthreads();
  }
  // write result for this block to global mem
  if (d_blktid == 0) pd_anychange[d_blkid] = blkChange[0];
}

__device__ int 
sweepXYZ(
    int nx, int ny, int nz, 
    int starstart, int starstop,
    struct FS *fs,
    float *vboxflat,
    float *ttboxflat
) 
{
  int	i, j, k, l, oi, oj, ok, iIdx, oIdx;
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
    
		//pre-compute all the needed values
    iIdx = k+nz*j+i*sx; oIdx = ok+nz*oj+oi*sx;
		delay = fs[l].d * (vboxflat[iIdx] + vboxflat[oIdx]) / 2.0;
		tt = ttboxflat[iIdx];
		tto = ttboxflat[oIdx];
    ttd = tt + delay;
    ttod = tto + delay;
    
    //if the difference between two values is greater than delay
    //do value switches using pre-calculated values.
    if(ttd < tto || ttod < tt){
      ttboxflat[iIdx] = fminf(tt, ttod);
      ttboxflat[oIdx] = fminf(tto, ttd);
      if(change == 0)
        change = 1;
    }
  }
  return(change);

} /* end sweepXYZ */ 

__global__ void 
MergeBoxes(
    float *pd_ttboxflat, 
    float *pd_ttbflatMergeCache, 
    int devNum, 
    int nCells)
{
  int d_blkid = blockIdx.z + blockIdx.y * gridDim.z + blockIdx.x * gridDim.z * gridDim.y;
	int d_glbtid = d_blkid * blockDim.x * blockDim.y * blockDim.z
								+ (threadIdx.x * blockDim.z * blockDim.y)
								+ (threadIdx.y * blockDim.z) + threadIdx.z;
  float res = INFINITY;
  res = fminf(res, pd_ttboxflat[d_glbtid]);
  int devFrom;
  for(devFrom=0; devFrom<devNum-1; devFrom++)
    if(d_glbtid < nCells)
      res = fminf(res, pd_ttbflatMergeCache[devFrom*nCells + d_glbtid]);
  pd_ttboxflat[d_glbtid] = res;
}