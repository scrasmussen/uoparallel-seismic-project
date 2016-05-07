#include "iovelocity.h"

#include <stdio.h>

int main() {

    struct VELOCITYBOX vbox;

    //if( !vboxloadtext( &vbox, "velocity-241-241-51.txt" ) ) {
    if( !vboxloadtext( &vbox, "velbad.txt" ) ) {
        printf( "error: velocityloadtext() failed\n" );
        return 0;
    }

    printf( "file read! here's what I interpreted:\n" );
    vboxfprint( stdout, "", "  ", vbox );

    return 0;
}
