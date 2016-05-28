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
#include "mpi.h"
#include <omp.h>
#include <sys/time.h>  


#define	FSRADIUSMAX	7	/* maximum radius forward star */
#define	FSMAX		818	/* maximum # of points in a forward star */
#define MODELMAX	250	/* maximum model dimension in X,Y,Z */
#define STARTMAX	4 //12	/* maximum starting points */

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
struct MODEL	model[MODELMAX][MODELMAX][60], model_new[MODELMAX+10][MODELMAX+10][204], transferS1[6][8][255][204],transferR1[6][8][255][204],transferS2[6][255][8][204],transferR2[6][255][8][204];
struct START	start[STARTMAX], start_new;
int startinew,startjnew,stopinew,stopjnew;
int sizeLocal=241,ghostcell=7;




int sweepXYZ(int nx, int ny, int nz, int s, int starstart, int starstop);

int main(int argc, char* argv[]) {
  int		i, j, k, l, m, nx, ny, nz, oi, oj, ok, s;
  int		numradius, starsize, anychange, numstart, numsweeps=0, numOfTasks=16;
  int		fsindex[FSRADIUSMAX];
  float		delta, delay;
  FILE		*vfile, *fsfile, *ttfile, *startfile;
  int           numtasks, taskid, len,tag;
  taskid=0;
  tag=1;
  int mmm;
  struct timeval t1, t2,t3;
  double elapsedTime;
  gettimeofday(&t3, NULL);
  

  int           sizeOfTasks[16][7]= {
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 7, 7, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 7, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 7, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 7, 7, 0},

  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 7, 0, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 7, 0},

  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 7, 0, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 0, 0},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 7, 0},

  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 7, 0, 0, 7},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 0, 7},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 0, 7},
  {(sizeLocal+2*ghostcell), (sizeLocal+2*ghostcell), 204, 0, 0, 7, 7},
  };



  int 	neighbour[16][9]=
  {
  {0, -1, -1, -1, 1, 5, 4, -1,-1},
  {1, -1, -1, -1, 2, 6, 5, 4, 0},
  {2, -1, -1, -1, 3, 7, 6, 5, 1},
  {3, -1, -1, -1, -1, -1, 7, 6, 2},

  {4, -1, 0, 1, 5, 9, 8, -1,-1},
  {5, 0, 1, 2, 6, 10, 9, 8, 4},
  {6, 1, 2, 3, 7, 11, 10, 9, 5},
  {7, 2, 3, -1, -1, -1, 11, 10, 6},

  {8, -1, 4, 5, 9, 13, 12, -1,-1},
  {9, 4, 5, 6, 10, 14, 13, 12, 8},
  {10, 5, 6, 7, 11, 15, 14, 13, 9},
  {11, 6, 7, -1, -1, -1, 15, 14, 10},

  {12, -1, 8, 9, 13, -1, -1, -1,-1},
  {13, 8, 9, 10, 14, -1, -1, -1, 12},
  {14, 9, 10, 11, 15, -1, -1, -1, 13},
  {15, 10, 11, -1, -1, -1, -1, -1, 14}
  };

int maxx=254;
int 	communicationMatrix[128][14];

  for (i=0;i<16;i++) {

     for(j=0;j<8;j++) {
        communicationMatrix[8*i+j][0]=i;
        communicationMatrix[8*i+j][1]=neighbour[i][j+1]; 
        if(communicationMatrix[8*i+j][0]!=communicationMatrix[8*i+j][1] && neighbour[i][j+1]>-1){
           if (j==0) {
              communicationMatrix[8*i+j][2]=ghostcell;
              communicationMatrix[8*i+j][3]=2*ghostcell-1;
              communicationMatrix[8*i+j][4]=ghostcell;
              communicationMatrix[8*i+j][5]=2*ghostcell-1;
              communicationMatrix[8*i+j][6]=0;
              communicationMatrix[8*i+j][7]=203;
              communicationMatrix[8*i+j][8]=maxx-ghostcell+1;
              communicationMatrix[8*i+j][9]=maxx;
              communicationMatrix[8*i+j][10]=maxx-ghostcell+1;
              communicationMatrix[8*i+j][11]=maxx;
              communicationMatrix[8*i+j][12]=0;
              communicationMatrix[8*i+j][13]=203;
           }
           else if (j==1) {
              communicationMatrix[8*i+j][2]=ghostcell;
              communicationMatrix[8*i+j][3]=maxx-ghostcell;
              communicationMatrix[8*i+j][4]=ghostcell;
              communicationMatrix[8*i+j][5]=2*ghostcell-1;
              communicationMatrix[8*i+j][6]=0;
              communicationMatrix[8*i+j][7]=203;
              communicationMatrix[8*i+j][8]=ghostcell;
              communicationMatrix[8*i+j][9]=maxx-ghostcell;
              communicationMatrix[8*i+j][10]=maxx-ghostcell+1;
              communicationMatrix[8*i+j][11]=maxx;
              communicationMatrix[8*i+j][12]=0;
              communicationMatrix[8*i+j][13]=203;
           } 
           else if (j==2) {
              communicationMatrix[8*i+j][2]=maxx-2*ghostcell+1;
              communicationMatrix[8*i+j][3]=maxx-ghostcell;
              communicationMatrix[8*i+j][4]=ghostcell;
              communicationMatrix[8*i+j][5]=2*ghostcell-1;
              communicationMatrix[8*i+j][6]=0;
              communicationMatrix[8*i+j][7]=203;
              communicationMatrix[8*i+j][8]=0;
              communicationMatrix[8*i+j][9]=ghostcell-1;
              communicationMatrix[8*i+j][10]=maxx-ghostcell+1;
              communicationMatrix[8*i+j][11]=maxx;
              communicationMatrix[8*i+j][12]=0;
              communicationMatrix[8*i+j][13]=203;
           }
           else if (j==3) {
              communicationMatrix[8*i+j][2]=maxx-2*ghostcell+1;
              communicationMatrix[8*i+j][3]=maxx-ghostcell;
              communicationMatrix[8*i+j][4]=ghostcell;
              communicationMatrix[8*i+j][5]=maxx-ghostcell;
              communicationMatrix[8*i+j][6]=0;
              communicationMatrix[8*i+j][7]=203;
              communicationMatrix[8*i+j][8]=0;
              communicationMatrix[8*i+j][9]=ghostcell-1;
              communicationMatrix[8*i+j][10]=ghostcell;
              communicationMatrix[8*i+j][11]=maxx-ghostcell;
              communicationMatrix[8*i+j][12]=0;
              communicationMatrix[8*i+j][13]=203;
           } 
           else if (j==4) {
              communicationMatrix[8*i+j][2]=maxx-2*ghostcell+1;
              communicationMatrix[8*i+j][3]=maxx-ghostcell;
              communicationMatrix[8*i+j][4]=maxx-2*ghostcell+1;
              communicationMatrix[8*i+j][5]=maxx-ghostcell;
              communicationMatrix[8*i+j][6]=0;
              communicationMatrix[8*i+j][7]=203;
              communicationMatrix[8*i+j][8]=0;
              communicationMatrix[8*i+j][9]=ghostcell-1;
              communicationMatrix[8*i+j][10]=0;
              communicationMatrix[8*i+j][11]=ghostcell-1;
              communicationMatrix[8*i+j][12]=0;
              communicationMatrix[8*i+j][13]=203;
           }
           else if (j==5) {
              communicationMatrix[8*i+j][2]=ghostcell;
              communicationMatrix[8*i+j][3]=maxx-ghostcell;
              communicationMatrix[8*i+j][4]=maxx-2*ghostcell+1;
              communicationMatrix[8*i+j][5]=maxx-ghostcell;
              communicationMatrix[8*i+j][6]=0;
              communicationMatrix[8*i+j][7]=203;
              communicationMatrix[8*i+j][8]=ghostcell;
              communicationMatrix[8*i+j][9]=maxx-ghostcell;
              communicationMatrix[8*i+j][10]=0;
              communicationMatrix[8*i+j][11]=ghostcell-1;
              communicationMatrix[8*i+j][12]=0;
              communicationMatrix[8*i+j][13]=203;
           } 
           else if (j==6) {
              communicationMatrix[8*i+j][2]=ghostcell;
              communicationMatrix[8*i+j][3]=2*ghostcell-1;
              communicationMatrix[8*i+j][4]=maxx-2*ghostcell+1;
              communicationMatrix[8*i+j][5]=maxx-ghostcell;
              communicationMatrix[8*i+j][6]=0;
              communicationMatrix[8*i+j][7]=203;
              communicationMatrix[8*i+j][8]=maxx-ghostcell+1;
              communicationMatrix[8*i+j][9]=maxx;
              communicationMatrix[8*i+j][10]=0;
              communicationMatrix[8*i+j][11]=ghostcell-1;
              communicationMatrix[8*i+j][12]=0;
              communicationMatrix[8*i+j][13]=203;

           }
           else if (j==7) {
              communicationMatrix[8*i+j][2]=ghostcell;
              communicationMatrix[8*i+j][3]=2*ghostcell-1;
              communicationMatrix[8*i+j][4]=ghostcell;
              communicationMatrix[8*i+j][5]=maxx-ghostcell;
              communicationMatrix[8*i+j][6]=0;
              communicationMatrix[8*i+j][7]=203;
              communicationMatrix[8*i+j][8]=maxx-ghostcell+1;
              communicationMatrix[8*i+j][9]=maxx;
              communicationMatrix[8*i+j][10]=ghostcell;
              communicationMatrix[8*i+j][11]=maxx-ghostcell;
              communicationMatrix[8*i+j][12]=0;
              communicationMatrix[8*i+j][13]=203;
           } 

       }
       else {
              communicationMatrix[8*i+j][2]=-1;
              communicationMatrix[8*i+j][3]=-1;
              communicationMatrix[8*i+j][4]=-1;
              communicationMatrix[8*i+j][5]=-1;
              communicationMatrix[8*i+j][6]=-1;
              communicationMatrix[8*i+j][7]=-1;
              communicationMatrix[8*i+j][8]=-1;
              communicationMatrix[8*i+j][9]=-1;
              communicationMatrix[8*i+j][10]=-1;
              communicationMatrix[8*i+j][11]=-1;
              communicationMatrix[8*i+j][12]=-1;
              communicationMatrix[8*i+j][13]=-1;
       } 

     }
  }

//|| (communicationMatrix[i][5]-communicationMatrix[i][4])!=(communicationMatrix[i][11]-communicationMatrix[i][10]) || (communicationMatrix[i][7]-communicationMatrix[i][13])!=(communicationMatrix[i][9]-communicationMatrix[i][12]) (communicationMatrix[i][3]-communicationMatrix[i][2])!=(communicationMatrix[i][9]-communicationMatrix[i][8])



//printf("\nvalue of i %d \n",i);
//exit(0);

///////////
//6 5 4
//7   3
//0 1 2
//////////

/*  int 	communicationMatrix[16][8]= {
	{0,0,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY},
	{0,1,114,120,0,120,0,50},
	{0,2,0,120,114,120,0,50},
	{0,3,114,120,114,120,0,50},
	{1,0,7,13,0,120,0,50},
	{1,1,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY},
	{1,2,7,13,114,120,0,50},
	{1,3,7,126,114,120,0,50},
	{2,0,0,120,7,13,0,50},
	{2,1,114,120,7,13,0,50},
	{2,2,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY},
	{2,3,114,120,7,126,0,50},
	{3,0,7,13,7,13,0,50},
	{3,1,7,126,7,13,0,50},
	{3,2,7,13,7,126,0,50},
	{3,3,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY,INFINITY}
	}; */

/*
  int 	communicationMatrix[16][14]= {
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,1,(maxx-13),(maxx-7),0,(maxx-7),0,203,0,6,0,(maxx-7),0,203},
	{0,2,0,(maxx-7),(maxx-13),(maxx-7),0,203,0,(maxx-7),0,6,0,203},
	{0,3,(maxx-13),(maxx-7),(maxx-13),(maxx-7),0,203,0,6,0,6,0,203},
	{1,0,7,13,0,(maxx-7),0,203,(maxx-6),maxx,0,(maxx-7),0,203},
	{1,1,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,2,7,13,(maxx-13),(maxx-7),0,203,(maxx-6),maxx,0,6,0,203},
	{1,3,7,maxx,(maxx-13),(maxx-7),0,203,7,maxx,0,6,0,203},
	{2,0,0,(maxx-7),7,13,0,203,0,(maxx-7),(maxx-6),maxx,0,203},
	{2,1,(maxx-13),(maxx-7),7,13,0,203,0,6,(maxx-6),maxx,0,203},
	{2,2,0,0,0,0,0,0,0,0,0,0,0,0},
	{2,3,(maxx-13),(maxx-7),7,maxx,0,203,0,6,7,maxx,0,203},
	{3,0,7,13,7,13,0,203,(maxx-6),maxx,(maxx-6),maxx,0,203},
	{3,1,7,maxx,7,13,0,203,7,maxx,(maxx-6),maxx,0,203},
	{3,2,7,13,7,maxx,0,203,(maxx-6),maxx,7,maxx,0,203},
	{3,3,0,0,0,0,0,0,0,0,0,0,0,0}
	};
*/
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
   MPI_Status stat;
   MPI_Status stats[32];
   MPI_Request reqs[32];
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   MPI_Datatype fstype, modeltype, starttype,oldtypes[2], oldtypes1[1];
   MPI_Aint    offsets[2], extent, offsets1[1];
   int          blockcounts[2], blockcounts1[1];

MPI_Aint intex, floatex;
MPI_Type_extent(MPI_INT, &intex);
MPI_Type_extent(MPI_FLOAT, &floatex);




	offsets[0] = 0;
	oldtypes[0] = MPI_INT;
	blockcounts[0] = 3;
	MPI_Type_extent(MPI_INT, &extent);
	offsets[1] = 3 * extent;
	oldtypes[1] = MPI_FLOAT;
	blockcounts[1] = 1;
   MPI_Type_struct(2, blockcounts, offsets, oldtypes, &fstype);
   MPI_Type_commit(&fstype); 



	//offset1 = 0;
	//oldtype1 = MPI_FLOAT;
	//blockcount1 = STARTMAX+1;
	offsets1[0] = 0;
	oldtypes1[0] = MPI_FLOAT;
	blockcounts1[0] =STARTMAX+1;
   MPI_Type_struct(1, blockcounts1, offsets1, oldtypes1, &modeltype);
   MPI_Type_commit(&modeltype); 

	//offset1 = 0;
	//oldtype1 = MPI_INT;
	//blockcount1 = 3;

	offsets1[0] = 0;
	oldtypes1[0] = MPI_INT;
	blockcounts1[0] = 3;

   MPI_Type_struct(1, blockcounts1, offsets1, oldtypes1,  &starttype);
   MPI_Type_commit(&starttype);

  /**
   * TODO: Make sure there aren't any race conditions.  If there are, fix
   *       with barriers.
   */

  /**
   * TODO: Only rank 0 should do I/O
   */

  /* open velocity model file */
//printf("Cannot open velocity model file: %s\n", argv[1]);
if(taskid==0){


  for (i=0;i<128;i++) {

if ((communicationMatrix[i][7]-communicationMatrix[i][6])!=(communicationMatrix[i][13]-communicationMatrix[i][12]) ){

printf("%d  ,  %d,  %d,  %d,  %d,  %d,  %d,  %d,  %d,  %d,  %d,  %d,  %d,  %d\n", communicationMatrix[i][0],communicationMatrix[i][1],communicationMatrix[i][2],communicationMatrix[i][3],communicationMatrix[i][4],communicationMatrix[i][5],communicationMatrix[i][6],communicationMatrix[i][7],communicationMatrix[i][8],communicationMatrix[i][9],communicationMatrix[i][10],communicationMatrix[i][11],communicationMatrix[i][12],communicationMatrix[i][13]);

}


  }
//exit(0);

printf("\ntask number : %d \n", taskid);
  //vfile = fopen(argv[1],"r");
  vfile = fopen("./mpi/shortest-path/velocity-241-241-51.txt","r");
  if(vfile == NULL) {
    printf("Cannot open velocity model file: %s\n", argv[1]);
    exit(1);
  }
  printf("Velocity model file: %s\n", argv[1]);

  /* open forward star offset file */
  //fsfile = fopen(argv[2],"r");
  fsfile = fopen("./mpi/shortest-path/818-FS.txt","r");
  if(fsfile == NULL) {
    printf("Cannot open forward star offset file: %s\n", argv[2]);
    exit(1);
  }
  printf("Forward star offset file: %s\n", argv[2]);

  /* open file with starting points */
  //startfile = fopen(argv[3],"r");
  startfile = fopen("./mpi/shortest-path/start-4.txt","r");
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
//printf("..... Monil Velocity data read %f", model[241][241][51].v);

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
  //printf("...... Monil Forward star offsets %i %i %i", fs[starsize].i, &fs[starsize].j, &fs[starsize].k);
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

 // commented to stop the initialization of start points
  for (s=0; s<numstart; s++) {
    fscanf(startfile, "%i %i %i", &start[s].i, &start[s].j, &start[s].k);
    //model[start[s].i][start[s].j][start[s].k].tt[s] = 0;
    printf("starting point %d: %d %d %d\n", s, start[s].i, start[s].j, start[s].k);
  }
  printf("Starting points read\n");



//start_new=start[0];
//if( numtasks>1) {
  for (i=1; i<numtasks; i++)  {
     MPI_Send(&start, STARTMAX, starttype, i, tag, MPI_COMM_WORLD);
//MPI_Send(&model, 250*250*250, modeltype, i, tag, MPI_COMM_WORLD);
//MPI_Send(&fs, FSMAX, fstype, i, tag, MPI_COMM_WORLD);
  }
//}

//printf("Taskid : %d: starting point %d: %d %d %d\n", taskid, start[0].i, start[0].j, start[0].k);

} //if ends here
else
{//STARTMAX
printf("\ntask number : %d \n", taskid);
MPI_Recv(&start, STARTMAX, starttype, 0, tag, MPI_COMM_WORLD, &stat);
//MPI_Recv(&model, 250*250*250, modeltype, 0, tag, MPI_COMM_WORLD, &stat);
//MPI_Recv(&fs, FSMAX, fstype, 0, tag, MPI_COMM_WORLD, &stat);
//printf("Taskid : %d: starting point %d: %d %d %d\n", taskid, start[taskid].i, start[taskid].j, start[taskid].k);

}  
MPI_Bcast(fs, FSMAX, fstype, 0, MPI_COMM_WORLD);
MPI_Bcast(model, 250*250*250, modeltype, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);



/*
* This section is for data copy to all tasks
*/
int task=0;
int starti, startj;
numstart=1;
for(task=0;task<numOfTasks;task++){

if (task==taskid){
   //nx = 127+1; ny = 127+1; nz = 51;
    nx = sizeOfTasks[taskid][0]; 
    ny = sizeOfTasks[taskid][1]; 
    nz = sizeOfTasks[taskid][2];
   //printf("Velocity model dimensions: %i %i %i\n",nx, ny, nz);
   //int ir, jr, kr;   /* read indices */
   //starti=sizeOfTasks[taskid][3];
   //startj=sizeOfTasks[taskid][4];
   startinew=sizeOfTasks[taskid][3];  // To restrict the access of empty ghost cell
   startjnew=sizeOfTasks[taskid][4]; // To restrict the access of empty ghost cell
   stopinew=sizeOfTasks[taskid][5];  // To restrict the access of empty ghost cell
   stopjnew=sizeOfTasks[taskid][6];  // To restrict the access of empty ghost cell
   printf("for Task : %d  Velocity model dimensions: %i %i %i \n",taskid, nx, ny, nz);
int m=0;
for (m=0;m<4;m++) {
   for (i=0; i<241; i++) {
     for (j=0; j<241; j++) {
       for (k=0; k<51; k++) {
           model_new[i+7][j+7][k+m*51].v=model[i][j][k].v; 
 	for (s=0; s<numstart; s++) {
 	  model_new[i+7][j+7][k+m*51].tt[s]=model[i][j][k].tt[s];
 	 }
       }
     }
   }
}
} //if ends here

} //for loop for tasks ends here



  anychange = 1;  
  mmm=0;
  s=0;
  tag=0;
  int send=0,receive=0,sender=0, reciever=0,reqnumber=0, count1=0, count2=0;
  int new_changeS[16], new_changeR[16];

if (taskid==1) {
//printf("\n \n \n \n");
//printf("rank: %d velocity %f traveltime %f \n",taskid,model_new[nx-1][ny-1][nz-1].v,model_new[nx-1][ny-1][nz-1].tt[s]);
//printf("rank: %d velocity %f traveltime %f \n",taskid,model_new[242][242][203].v,model_new[242][242][203].tt[s]);

}

//MPI_Barrier(MPI_COMM_WORLD);

//exit(0);

//data preperation mpi part

    reqnumber=0;count1=0;count2=0;
    for (send=0; send<128; send++) {
        sender= communicationMatrix[send][0];
        reciever= communicationMatrix[send][1];

        if (sender==taskid && sender!=reciever && reciever >-1){

         if (communicationMatrix[send][3]-communicationMatrix[send][2] < 9) {
            for (i=communicationMatrix[send][2]; i<communicationMatrix[send][3]+1; i++) {
               for (j=communicationMatrix[send][4]; j<communicationMatrix[send][5]+1; j++) {
                  for (k=communicationMatrix[send][6]; k<communicationMatrix[send][7]+1; k++) {
 	             //for (s=0; s<numstart; s++) {

 	                transferS1[count1][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].tt[s]=model_new[i][j][k].tt[s];
 	                transferS1[count1][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].v=model_new[i][j][k].v;
                  }
               }
             }                        
//printf("sender: %d, reciever: %d, i %d j %d k %d \n",sender,reciever, i, j, k);
            MPI_Isend(&transferS1[count1], 8*255*204, modeltype, reciever, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
            reqnumber++;
            count1++; 	             //}
         }else{
            for (i=communicationMatrix[send][2]; i<communicationMatrix[send][3]+1; i++) {
               for (j=communicationMatrix[send][4]; j<communicationMatrix[send][5]+1; j++) {
                  for (k=communicationMatrix[send][6]; k<communicationMatrix[send][7]+1; k++) {

 	                transferS2[count2][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].tt[s]=model_new[i][j][k].tt[s];
 	                transferS2[count2][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].v=model_new[i][j][k].v;
                        //count2++;                       
                      
                  }
               }
             }
            MPI_Isend(&transferS2[count2], 255*8*204, modeltype, reciever, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
           reqnumber++;
           count2++;  
          }  /// if else finished for comparing the value of x



        }//for comparing sender and reciever
        // for task
             
      } //sender for ends here


count1=0;count2=0;

    for (receive=0; receive<128; receive++) {
        sender= communicationMatrix[receive][0];
        reciever= communicationMatrix[receive][1];
        if (reciever==taskid && sender!=reciever && reciever > -1){

             //MPI_Irecv(&transferR1[sender], 130*130*51, modeltype, sender, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
             //reqnumber++;
        if (communicationMatrix[receive][9]-communicationMatrix[receive][8] < 9) {
             MPI_Irecv(&transferR1[count1], 8*255*204, modeltype, sender, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
             reqnumber++; count1++;
        }else{
             MPI_Irecv(&transferR2[count2], 255*8*204, modeltype, sender, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
             reqnumber++; count2++;             
        }


        }//for comparing sender and reciever
        
             
      } //reciever for ends here


    //now block until requests are complete
    MPI_Waitall(reqnumber, reqs, stats);

reqnumber=0; 
int ii=0; count1=0;count2=0;

    for (receive=0; receive<128; receive++) {
        sender= communicationMatrix[receive][0];
        reciever= communicationMatrix[receive][1];
        if (reciever==taskid && sender!=reciever && reciever >-1){
           //printf("sender: %d, reciever: %d, i %d i %d j %d j %d k %d k %d \n",sender,reciever, communicationMatrix[receive][8], communicationMatrix[receive][9]+1, communicationMatrix[receive][10],communicationMatrix[receive][11]+1,communicationMatrix[receive][12],communicationMatrix[receive][13]+1 );
ii=0;
        if (communicationMatrix[receive][9]-communicationMatrix[receive][8] < 9) {
           for (i=communicationMatrix[receive][8]; i<communicationMatrix[receive][9]+1; i++) {
               for (j=communicationMatrix[receive][10]; j<communicationMatrix[receive][11]+1; j++) {
                  for (k=communicationMatrix[receive][12]; k<communicationMatrix[receive][13]+1; k++) {

                            model_new[i][j][k].tt[s]=transferR1[count1][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].tt[s];
                            model_new[i][j][k].v=transferR1[count1][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].v;

                  }
               }
             } 
         count1++;
         }else{
           for (i=communicationMatrix[receive][8]; i<communicationMatrix[receive][9]+1; i++) {
               for (j=communicationMatrix[receive][10]; j<communicationMatrix[receive][11]+1; j++) {
                  for (k=communicationMatrix[receive][12]; k<communicationMatrix[receive][13]+1; k++) {

                            model_new[i][j][k].tt[s]=transferR2[count2][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].tt[s];     
                            model_new[i][j][k].v=transferR2[count2][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].v;  
                  }
               }
             }
         count2++;          
         } // if else finished here
	             

        }//for comparing sender and reciever
             //printf ( "\n value of i : %d",i);
      } //reciever for ends here

//if (taskid==1) {

//printf("rank: %d velocity %f traveltime %f \n",taskid,model_new[nx-1][ny-1][nz-1].v,model_new[nx-1][ny-1][nz-1].tt[s]);
//printf("rank: %d velocity %f traveltime %f \n",taskid,model_new[242][242][203].v,model_new[242][242][203].tt[s]);

//}
MPI_Barrier(MPI_COMM_WORLD);

//exit(0);

//data preperation done




//memory handling ends here

//sweep activities begin here
s=0;

if (taskid==0) {
 gettimeofday(&t1, NULL); // time counting

model_new[start[s].i][start[s].j][start[s].k].tt[s] = 0;
}



  while (anychange) {
//for (m=0;m<5;m++){
    numsweeps++;
    anychange = 0;
    //for (s=0; s<numstart; s++) {
    for (task=0; task<numOfTasks; task++) {
      if(taskid==task){
        changed[s] = 0;

        changed[s] = sweepXYZ(nx, ny, nz, s, 0, 818-1);
        //printf(">>> start point %d: by Task: %d, changed == %d\n", s,taskid, changed[s]);
        anychange = changed[s];
        new_changeS[taskid]=anychange;
        //if (anychange>0) MPI_Bcast(&anychange, 1, MPI_INT, taskid, MPI_COMM_WORLD);
        //if (numsweeps<4) anychange=1;
        printf("sweep %d finished by task: %d for >>> start point %d: anychange = %d\n",  numsweeps,taskid, s, changed[s]);
        //if (numsweeps<4) anychange=1;



    reqnumber=0;count1=0;count2=0;
    for (send=0; send<128; send++) {
        sender= communicationMatrix[send][0];
        reciever= communicationMatrix[send][1];

        if (sender==taskid && sender!=reciever && reciever >-1){

         if (communicationMatrix[send][3]-communicationMatrix[send][2] < 9) {
            for (i=communicationMatrix[send][2]; i<communicationMatrix[send][3]+1; i++) {
               for (j=communicationMatrix[send][4]; j<communicationMatrix[send][5]+1; j++) {
                  for (k=communicationMatrix[send][6]; k<communicationMatrix[send][7]+1; k++) {
 	             //for (s=0; s<numstart; s++) {

 	                transferS1[count1][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].tt[s]=model_new[i][j][k].tt[s];
 	                //transferS1[count1][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].v=model_new[i][j][k].v;
                  }
               }
             }                        
//printf("sender: %d, reciever: %d, i %d j %d k %d \n",sender,reciever, i, j, k);
            MPI_Isend(&transferS1[count1], 8*255*204, modeltype, reciever, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
            reqnumber++;
            count1++; 	             //}
         }else{
            for (i=communicationMatrix[send][2]; i<communicationMatrix[send][3]+1; i++) {
               for (j=communicationMatrix[send][4]; j<communicationMatrix[send][5]+1; j++) {
                  for (k=communicationMatrix[send][6]; k<communicationMatrix[send][7]+1; k++) {

 	                transferS2[count2][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].tt[s]=model_new[i][j][k].tt[s];
 	                //transferS2[count2][i-communicationMatrix[send][2]][j-communicationMatrix[send][4]][k].v=model_new[i][j][k].v;
                        //count2++;                       
                      
                  }
               }
             }
            MPI_Isend(&transferS2[count2], 255*8*204, modeltype, reciever, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
           reqnumber++;
           count2++;  
          }  /// if else finished for comparing the value of x



        }//for comparing sender and reciever
        // for task
             
      } //sender for ends here


count1=0;count2=0;

    for (receive=0; receive<128; receive++) {
        sender= communicationMatrix[receive][0];
        reciever= communicationMatrix[receive][1];
        if (reciever==taskid && sender!=reciever && reciever > -1){

             //MPI_Irecv(&transferR1[sender], 130*130*51, modeltype, sender, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
             //reqnumber++;
        if (communicationMatrix[receive][9]-communicationMatrix[receive][8] < 9) {
             MPI_Irecv(&transferR1[count1], 8*255*204, modeltype, sender, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
             reqnumber++; count1++;
        }else{
             MPI_Irecv(&transferR2[count2], 255*8*204, modeltype, sender, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
             reqnumber++; count2++;             
        }


        }//for comparing sender and reciever
        
             
      } //reciever for ends here


    //now block until requests are complete
    MPI_Waitall(reqnumber, reqs, stats);
 
   //printf("Taskid %d, anychange %d", taskid,anychange);
//reqnumber=0;
reqnumber=0;

    for(i=0;i<16;i++){
      if (i!=taskid) {
            //MPI_Send(&anychange, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
            MPI_Isend(&new_changeS[taskid], 1, MPI_INT, i, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
           reqnumber++;
         } 
    }
    for(i=0;i<16;i++){
      if (i!=taskid){
       // MPI_Recv(&new_changeR[taskid], 1, MPI_INT, i, tag, MPI_COMM_WORLD);
       MPI_Irecv(&new_changeR[i], 1, MPI_INT, i, tag, MPI_COMM_WORLD, &reqs[reqnumber]);
       reqnumber++;
      }
    }    
//MPI_Recv(&start, STARTMAX, starttype, 0, tag, MPI_COMM_WORLD, &stat); 

MPI_Waitall(reqnumber, reqs, stats);

    for(i=0;i<16;i++){
      if (new_changeR[i]>0){ anychange=new_changeR[i];
      }
    }
    //printf("Taskid %d, anychange %d", taskid,anychange);


reqnumber=0; 
int ii=0; count1=0;count2=0;

    for (receive=0; receive<128; receive++) {
        sender= communicationMatrix[receive][0];
        reciever= communicationMatrix[receive][1];
        if (reciever==taskid && sender!=reciever && reciever >-1){
           //printf("sender: %d, reciever: %d, i %d i %d j %d j %d k %d k %d \n",sender,reciever, communicationMatrix[receive][8], communicationMatrix[receive][9]+1, communicationMatrix[receive][10],communicationMatrix[receive][11]+1,communicationMatrix[receive][12],communicationMatrix[receive][13]+1 );
ii=0;
        if (communicationMatrix[receive][9]-communicationMatrix[receive][8] < 9) {
           for (i=communicationMatrix[receive][8]; i<communicationMatrix[receive][9]+1; i++) {
               for (j=communicationMatrix[receive][10]; j<communicationMatrix[receive][11]+1; j++) {
                  for (k=communicationMatrix[receive][12]; k<communicationMatrix[receive][13]+1; k++) {

                            model_new[i][j][k].tt[s]=transferR1[count1][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].tt[s];
                            //model_new[i][j][k].v=transferR1[count1][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].v;

                  }
               }
             } 
         count1++;
         }else{
           for (i=communicationMatrix[receive][8]; i<communicationMatrix[receive][9]+1; i++) {
               for (j=communicationMatrix[receive][10]; j<communicationMatrix[receive][11]+1; j++) {
                  for (k=communicationMatrix[receive][12]; k<communicationMatrix[receive][13]+1; k++) {

                            model_new[i][j][k].tt[s]=transferR2[count2][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].tt[s];     
                            //model_new[i][j][k].v=transferR2[count2][i-communicationMatrix[receive][8]][j-communicationMatrix[receive][10]][k].v;  
                  }
               }
             }
         count2++;          
         } // if else finished here
	             

        }//for comparing sender and reciever
             //printf ( "\n value of i : %d",i);
      } //reciever for ends here



      } //task comparing if ends here
    } //starting for ends here for ends here


  MPI_Barrier(MPI_COMM_WORLD);
  } //while ends here


MPI_Barrier(MPI_COMM_WORLD);


//time measurement for sweep ends here

if (taskid==0) {
    gettimeofday(&t2, NULL);

    // compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t1.tv_sec);      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    printf( "\n elapsed Second for sweeps only: %f \n", elapsedTime);
}



if (taskid==0) {
  /* print travel times */

  ttfile = fopen("mpi/shortest-path/output.tt","w");
  if(ttfile == NULL) {
    printf(".......................Can not open travel time output file: %s\n", "output.tt");
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
}

MPI_Barrier(MPI_COMM_WORLD);
MPI_Type_free(&modeltype);
MPI_Type_free(&fstype);
MPI_Type_free(&starttype);
MPI_Finalize();


    gettimeofday(&t2, NULL);

    // compute and print the elapsed time in millisec
    elapsedTime = (t2.tv_sec - t3.tv_sec);      // sec to ms
    elapsedTime += (t2.tv_usec - t3.tv_usec) / 1000000.0;   // us to ms
    printf( "\n elapsed Second for whole program: %f \n", elapsedTime);



  /**
   * TODO: Shutdown MPI library
   */

} /* main */


int sweepXYZ(int nx, int ny, int nz, int s, int starstart, int starstop) {
  int	i, j, k, l, oi, oj, ok;
  int	change = 0;
  float	delay = 0.0, tt = 0.0, tto = 0.0;

//  for (i=0; i<nx; i++) {
//    for (j=0; j<ny; j++) {
  for (i=7; i<nx-7; i++) {
    for (j=7; j<ny-7; j++) {
      for (k=0; k<nz; k++) {
	for (l=starstart; l<starstop; l++) {
	  /* find point in forward star based on offsets */
	  oi = i+fs[l].i; oj = j+fs[l].j; ok = k+fs[l].k;
   //startinew=sizeOfTasks[taskid][3];  // To restrict the access of empty ghost cell
   //startjnew=sizeOfTasks[taskid][4]; // To restrict the access of empty ghost cell
   //stopinew=sizeOfTasks[taskid][5];  // To restrict the access of empty ghost cell
   //stopjnew=sizeOfTasks[taskid][6]; 
	  /* if (oi,oj,ok) is outside the boundaries, then skip */
	  if ((oi < 0 + startinew) || (oi > nx-1-stopinew)
	      || (oj < 0 + startjnew) || (oj > ny-1 - stopjnew)
	      || (ok < 0) || (ok > nz-1)) {

	    continue;
	  }
	  /* compute delay from (i,j,k) to (oi,oj,ok) with end point average */
	  delay = fs[l].d * (model_new[i][j][k].v + model_new[oi][oj][ok].v) / 2.0;
	  /* update travel times for all starting points */
	  /* if (i,j,k) is starting point, then skip */

	  if ((i == start[s].i) && (j == start[s].j) && (k == start[s].k)) {

	  	    continue;
	  }
	  tt = model_new[i][j][k].tt[s];
	  tto = model_new[oi][oj][ok].tt[s];
	  /* if offset point has infinity travel time, then update */
	  if ((tt == INFINITY) && (tto == INFINITY)) {

	    continue;
	  }
	  if ((tt != INFINITY) && (tto == INFINITY)) {
	    model_new[oi][oj][ok].tt[s] = delay + tt;
	    change += 1;
	    continue;
	  }
	  if ((tt == INFINITY) && (tto != INFINITY)) {
	    model_new[i][j][k].tt[s] = delay + tto;
	    change += 1;
	    continue;
	  }
	  if ((tt != INFINITY) && (tto != INFINITY)) {
	    /* if a shorter travel time through (oi,oj,ok), update (i,j,k) */
	    if ((delay + tto) < tt) {
	      model_new[i][j][k].tt[s] = delay + tto;
	      change += 1;
	    }
	    /* if a shorter travel time through (i,j,k), update (oi,oj,ok) */
	    else if ((delay + tt) < tto) {
	      model_new[oi][oj][ok].tt[s] = delay + tt;
	      change += 1;
	    }
	  }
	}
      }
    }
  }
  return(change);

} /* end sweepXYZ */ 
