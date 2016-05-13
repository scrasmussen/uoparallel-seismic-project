// example_floatbox.c - Atlee Brink
// A working example of all the functions in ../include/floatbox.h

#include "../include/floatbox.h"

#include <stdio.h> /* for printf and stdout, not necessary otherwise */

int main()
{
    struct FLOATBOX box;

    // allocate heap memory for box
    if( !boxalloc( &box, 241, 241, 51 ) ) { /* handle error */ }

    // print box metadata
    boxfprint( stdout, "example: ", "\t", box );

    // put the value 4.567f to coordinates [1][2][3]
    boxput( box, 1, 2, 3, 4.567f );

    // get a value from coordinates [1][2][3] and print it
    printf( "value at (1,2,3): %g\n", boxget( box, 1, 2, 3 ) );

    // free heap memory
    boxfree( &box );
}
