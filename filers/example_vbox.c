#include "iovelocity.h"
#include <stdio.h> /* for stdout, not necessary otherwise */
int main() {
    struct VELOCITYBOX vbox;
    if( !vboxloadtext( &vbox, "velocity-241-241-51.txt" ) ) { /* handle error */ }
    if( !vboxstorebinary( "velocities.vbox", vbox ) ) { /* handle error */ }
    vboxfree( &vbox );
    if( !vboxloadbinary( &vbox, "velocities.vbox" ) ) { /* handle error */ }
    vboxfprint( stdout, "example: ", "\t", vbox );
    vboxfree( &vbox );
}
