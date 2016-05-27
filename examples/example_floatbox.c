// example_floatbox.c - Atlee Brink
// A working example of all the functions in ../include/floatbox.h

#include "floatbox.h"

#include <stdio.h> /* for printf and stdout, not necessary otherwise */

int main()
{
    struct FLOATBOX box;

    // allocate heap memory for box
    if( !boxalloc( &box, 241, 241, 51 ) ) { /* handle error */ }

    // print box metadata
    boxfprint( stdout, "example: ", "\t", box );

    // put the value 4.567f to coordinates [1][2][3]
    const float putval = 4.567f;
    boxput( box, 1, 2, 3, putval );

    // get a value from coordinates [1][2][3] and print it
    printf( "value at (1,2,3): %g (should be %g)\n", boxget( box, 1, 2, 3 ), putval );

    printf( "boxvolume( box ): %zu\n", boxvolume( box ) );

    const float fillval = 9.876f;
    int bad = 0;
    printf( "calling boxsetall(..)...\n" ); fflush( stdout );
    boxsetall( box, fillval );
    printf( "validating values...\n" ); fflush( stdout );
    for( int x = 0; !bad && x < box.size.x; x++ ) {
        for( int y = 0; !bad && y < box.size.y; y++ ) {
            for( int z = 0; !bad && z < box.size.z; z++ ) {
                if( boxget( box, x, y, z ) != fillval ) {
                    printf( "error: boxsetall(..) or boxget(..) went bad at (%d, %d, %d)!\n",
                        x, y, z );
                    bad = 1;
                }
            }
        }
    }
    if( !bad ) printf( "boxsetall(..) passed\n" );

    // free heap memory
    boxfree( &box );
}
