#include "iovelocity.h"

#include <stdio.h>

int main() {

    struct VELOCITYBOX vbox;

    int x = 200, y = 200, z = 40;

    printf( "reading text format velocities...\n" );
    fflush( stdout );
    if( !vboxloadtext( &vbox, "velocity-241-241-51.txt" ) ) {
        printf( "error: vboxloadtext() failed\n" );
        return 0;
    }

    printf( "text file read! here's what I interpreted:\n" );
    vboxfprint( stdout, "", "  ", vbox );

    printf( "value at: (%d, %d, %d): %g\n",
        x, y, z, boxget( vbox.box, x, y, z ) );

    /*
    printf( "writing binary format velocities...\n" );
    fflush( stdout );
    if( !vboxstorebinary( "velocity-241-241-51.vbox", vbox ) ) {
        printf( "error: vboxstorebinary() failed\n" );
        return 0;
    }

    vboxfree( &vbox );
    */

    struct VELOCITYBOX vboxb;

    if( !vboxloadbinary( &vboxb, "velocity-241-241-51.vbox" ) ) {
        printf( "error: vboxloadbinary() failed\n" );
        return 0;
    }

    printf( "binary file read! here's what I interpreted:\n" );
    vboxfprint( stdout, "", "  ", vboxb );

    printf( "value at: (%d, %d, %d): %g\n",
        x, y, z, boxget( vboxb.box, x, y, z ) );

    if( !vboxstorebinary( "velocity-241-241-51-rewrite.vbox", vboxb ) ) {
        printf( "error: vboxstorebinary() failed\n" );
        return 0;
    }

    vboxfree( &vboxb );

    return 0;
}
