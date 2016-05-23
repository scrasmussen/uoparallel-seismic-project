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
struct START	start[STARTMAX], start_new;





int sweepXYZ(int nx, int ny, int nz, int s, int starstart, int starstop);

int main(int argc, char* argv[]) {
  int		i, j, k, l, m, nx, ny, nz, oi, oj, ok, s;
  int		numradius, starsize, anychange, numstart, numsweeps=0;
  int		fsindex[FSRADIUSMAX];
  float		delta, delay;
  FILE		*vfile, *fsfile, *ttfile, *startfile;
  int           numtasks, taskid, len,tag;
  taskid=0;
  tag=1;
  int mmm;
    struct timeval t1, t2;
    double elapsedTime;
gettimeofday(&t1, NULL);
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
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   MPI_Datatype fstype, modeltype, starttype,oldtypes[2], oldtypes1[1];
   MPI_Aint    offsets[2], extent, offsets1[1];
   int          blockcounts[2], blockcounts1[1];

MPI_Aint intex, floatex;
MPI_Type_extent(MPI_INT, &intex);
MPI_Type_extent(MPI_FLOAT, &floatex);

/*
int blocks[3]={1,1,1,1};
MPI_Datatype types[4]={MPI_INT,MPI_INT,MPI_INT,MPI_FLOAT};
MPI_Aint displacements[4];
displacements[0] = 0;
displacements[1] = intex;
displacements[2] = intex+intex;
displacements[3] = intex+intex+intex;
MPI_Type_struct(4, blocks, displacements, types, &fstype);
MPI_Type_commit(&fstype); 
/*
int blocks1[2]={1,12};
MPI_Datatype types1[2]={MPI_FLOAT,MPI_FLOAT};
MPI_Aint displacements1[2];
//MPI_Aint intex, floatex;
//MPI_Type_extent(MPI_INT, &intex);
//MPI_Type_extent(MPI_FLOAT, &floatex);
displacements1[0] = 0;
displacements1[1] = floatex;

MPI_Type_struct(2, blocks1, displacements1, types1, &modeltype);
MPI_Type_commit(&modeltype); 
/*
int blocks2[3]={1,1,1};
MPI_Datatype types2[3]={MPI_INT,MPI_INT,MPI_INT};
MPI_Aint displacements2[3];
//MPI_Aint intex, floatex;
//MPI_Type_extent(MPI_INT, &intex);
//MPI_Type_extent(MPI_FLOAT, &floatex);
displacements2[0] = 0;
displacements2[1] = intex;
displacements2[2] = intex+ intex;

MPI_Type_struct(3, blocks2, displacements2, types2, &starttype);
MPI_Type_commit(&starttype); */

//MPI_Aint array_of_displaysments[1];
//MPI_Aint address1, address2;

//MPI_Get_address(&fs[0],&address1);
//MPI_Get_address(&fs[0].i,&address2);
//array_of_displaysments[0] = address2 - address1;


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
  for (s=0; s<numstart; s++) {
    fscanf(startfile, "%i %i %i", &start[s].i, &start[s].j, &start[s].k);
    model[start[s].i][start[s].j][start[s].k].tt[s] = 0;
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
printf("Taskid : %d: starting point %d: %d %d %d\n", taskid, start[taskid].i, start[taskid].j, start[taskid].k);

}  
MPI_Bcast(fs, FSMAX, fstype, 0, MPI_COMM_WORLD);
MPI_Bcast(model, 250*250*250, modeltype, 0, MPI_COMM_WORLD);
MPI_Barrier(MPI_COMM_WORLD);

//printf("..... Monil: taskid: %d Velocity data read %f",taskid, model[241][241][51].v);
  //printf("...... Monil: taskid: %d Forward star offsets %i %i %i",taskid, fs[818].i, &fs[818].j, &fs[818].k);
//printf("Taskid : %d: fs value : %f \n", taskid, fs[0].d);
//printf("Taskid : %d: model value : %f \n", taskid, model[0][0][0].tt[taskid]);



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
  
  //mmm=0;
   nx = 241; ny = 241; nz = 51;
   numstart=4;
  while (anychange) {
    /* sweep forward on X */
    numsweeps++;
    anychange = 0;
    //printf("sweep %d begin by task: %d \n", numsweeps, taskid);

    /**
     * TODO: A rank should only do its starting point, not all of them.
     *
     * Note: You may want to use also write out the rank of the process
     *       that is reporting progress.
     */
 
    for (s=0; s<numstart; s++) {
      //changed[taskid] = 0;
      if(taskid==s){
        changed[s] = 0;
        changed[s] += sweepXYZ(nx, ny, nz, s, 0, 818-1);
        printf(">>> start point %d: by Task: %d, changed == %d\n", s,taskid, changed[s]);
        //printf(">>> Task: %d, changed == %d\n", taskid, mmm);
        //}
        //for (s=0; s<1; s++) {
        //anychange += mmm;
        anychange += changed[s];
        printf("sweep %d finished by task: %d : anychange = %d\n",  numsweeps,taskid, anychange);
      }
      else if (numtasks==1){
        changed[s] = 0;
        changed[s] += sweepXYZ(nx, ny, nz, s, 0, 818-1);
        printf(">>> start point %d: by Task: %d, changed == %d\n", s,taskid, changed[s]);
        //printf(">>> Task: %d, changed == %d\n", taskid, mmm);
        //}
        //for (s=0; s<1; s++) {
        //anychange += mmm;
        anychange += changed[s];
        printf("sweep %d finished by task: %d : anychange = %d\n",  numsweeps,taskid, anychange);        
        }
    } //for ends here


  } //while ends here


  /**
   * TODO: Rank 0 should gather results from other processes or it's OK to have
   *       each process write its own output.  Make sure that the filenames
   *       don't collide.  At least think about how rank 0 could gather output
   *       results.
   */

  /* TODO: Remove exit statement so output can complete. */
MPI_Barrier(MPI_COMM_WORLD);

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
    elapsedTime = (t2.tv_sec - t1.tv_sec);      // sec to ms
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000000.0;   // us to ms
    printf( "\n elapsed Second: %f \n", elapsedTime);


  /**
   * TODO: Shutdown MPI library
   */

} /* main */


int sweepXYZ(int nx, int ny, int nz, int s, int starstart, int starstop) {
  int	i, j, k, l, oi, oj, ok;
  int	change = 0;
  float	delay = 0.0, tt = 0.0, tto = 0.0;
//printf("inside function for %d\n",s);
//printf("inside sweep Taskid : %d: fs value : %f \n", s, fs[0].d);
//printf("inside sweep Taskid : %d: model value : %f \n", s, model[0][0][0].tt[s]);
    //printf("inside sweep starting point %d: %d %d %d\n", s, start[s].i, start[s].j, start[s].k);  
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
