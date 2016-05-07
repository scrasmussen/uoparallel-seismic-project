#include "iovelocity.h"

#include <stdio.h>

int main() {

    struct VELOCITYBOX vbox;

    printf( "reading text format velocities...\n" );
    fflush( stdout );
    if( !vboxloadtext( &vbox, "velocity-241-241-51.txt" ) ) {
        printf( "error: vboxloadtext() failed\n" );
        return 0;
    }

    printf( "file read! here's what I interpreted:\n" );
    vboxfprint( stdout, "", "  ", vbox );

    printf( "writing binary format velocities...\n" );
    fflush( stdout );
    if( !vboxstorebinary( "velocity-241-241-51.vbox", vbox ) ) {
        printf( "error: vboxstorebinary() failed\n" );
        return 0;
    }

    vboxfree( &vbox );

    return 0;
}
