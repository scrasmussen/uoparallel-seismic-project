// example_mpihelpers.c - Cody Hatfield
// A working example of all the functions in ../include/mpihelpers.h

#include "../include/mpihelpers.h"

#include <stdio.h> /* for printf and stdout, not necessary otherwise */

int main()
{
    int my_rank = 0;
    int width = 2;
    int i = 0;

    struct POINT3D src_min;
    point3dset(&src_min, 1, 1, 1);
    struct POINT3D src_dims;
    point3dset(&src_dims, 241, 241, 51);
    struct POINT3D min;
    struct POINT3D dims;

    for(i = 0; i < 8; ++i){
        mpicalculatemycoordinates(&min, &dims, i, 8, src_min, src_dims);
        printf("%d min - %d, %d, %d\n", i, min.x, min.y, min.z);
        printf("%d dims - %d, %d, %d\n", i, dims.x, dims.y, dims.z);
    }
    printf("\n");


    for(i = 0; i < 8; ++i){
        printf("%d Neighbor for %d: %d\n", i, my_rank, mpifindneighborrank(my_rank, i, width) );
    }
    printf("\n");

    my_rank = 2;
    for(i = 0; i < 8; ++i){
        printf("%d Neighbor for %d: %d\n", i, my_rank, mpifindneighborrank(my_rank, i, width) );
    }
    printf("\n");

    for(i = 0; i < 8; ++i){
      mpigetsendcoordinates(&min, &dims, i, 7, 1, 1, 125, 125, 51);
      printf("Send coordinates for %d: %d, %d, %d \t %d, %d, %d\n", i, min.x, min.y, min.z, dims.x, dims.y, dims.z);
    }
    printf("\n");

    for(i = 0; i < 8; ++i){
      mpigetreceivecoordinates(&min, &dims, i, 7, 1, 1, 125, 125, 51);
      printf("Receive coordinates for %d: %d, %d, %d \t %d, %d, %d\n", i, min.x, min.y, min.z, dims.x, dims.y, dims.z);
    }
    printf("\n");


}
