// example_mpihelpers.c - Cody Hatfield
// A working example of all the functions in ../include/mpihelpers.h

#include "../include/mpihelpers.h"

#include <stdio.h> /* for printf and stdout, not necessary otherwise */

int main()
{
    int my_rank = 0;
    int width = 2;
    int i = 0;

    for(i = 0; i < 8; ++i){
        printf("%d Neighbor for %d: %d\n", i, my_rank, mpifindneighborrank(my_rank, i, width) );
    }
    printf("\n");

    my_rank = 2;
    for(i = 0; i < 8; ++i){
        printf("%d Neighbor for %d: %d\n", i, my_rank, mpifindneighborrank(my_rank, i, width) );
    }
    printf("\n");

    struct TempPoint tp;
    for(i = 0; i < 8; ++i){
      mpigetsendcoordinates(i, 7, 1, 1, 125, 125, 51, &tp);
      printf("Send coordinates for %d: %d, %d, %d \t %d, %d, %d\n", i, tp.x, tp.y, tp.z, tp.x_size, tp.y_size, tp.z_size);
    }
    printf("\n");

    for(i = 0; i < 8; ++i){
      mpigetreceivecoordinates(i, 7, 1, 1, 125, 125, 51, &tp);
      printf("Receive coordinates for %d: %d, %d, %d \t %d, %d, %d\n", i, tp.x, tp.y, tp.z, tp.x_size, tp.y_size, tp.z_size);
    }
    printf("\n");


}
