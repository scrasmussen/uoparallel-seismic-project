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
}
