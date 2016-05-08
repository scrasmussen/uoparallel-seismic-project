#include "floatbox.h"

#include <stdio.h>

int main()
{
    int nx, ny, nz;
    nx = ny = nz = 3192; // make a very large box (flat indices well beyond INT_MAX)

    struct FLOATBOX box;
    if( !boxalloc( &box, nx, ny, nz ) ) {
        printf( "boxalloc failed\n" );
        return 0;
    }

    /*
    for( int x = 0; x < nx; x++ ) {
        for( int y = 0; y < ny; y++ ) {
            for( int z = 0; z < nz; z++ ) {
                boxput( box, x, y, z, (float)x + y + z );
            }
        }
    }


    float sum = 0.f;
    for( int x = 0; x < nx; x++ ) {
        for( int y = 0; y < ny; y++ ) {
            for( int z = 0; z < nz; z++ ) {
                sum += boxget( box, x, y, z );
            }
        }
    }

    printf( "sum: %g\n", sum );
    */

    // explicitly store a key value at the very last index
    box.flat[ 3192l * 3192l * 3192l - 1l ] = 1.234f;

    // check that boxget() indexing works correctly
    float get = boxget( box, nx-1, ny-1, nz-1 );
    printf( "last: %g (should be 1.234)\n", get );

    // show the index computed by boxindex()
    size_t idx = boxindex( box, nx-1, ny-1, nz-1 );
    printf( "index of last element by boxindex(): %zu\n", idx );

    // note: since boxindex() is used for both boxget and boxput, if it
    //       works for one, it will work for the other.

    boxfree( &box );

    return 0;
}
