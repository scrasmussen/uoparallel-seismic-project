#include "floatbox.h"

#include <stdio.h>

int main()
{
    int nx, ny, nz;
    nx = ny = nz = 3192;

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

    box.flat[ 3192l * 3192l * 3192l - 1l ] = 1.234f;

    float get = boxget( box, nx-1, ny-1, nz-1 );
    printf( "last: %g\n", get );

    size_t idx = boxindex( box, nx-1, ny-1, nz-1 );
    printf( "index of last element by boxindex(): %zu\n", idx );

    boxfree( &box );

    return 0;
}
