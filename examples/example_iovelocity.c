// example_iovelocity.c - Atlee Brink
// A working example of all the functions in ../include/iovelocity.h,
// and shows how to get and put values from the VELOCITYBOX (through its FLOATBOX).

#include "../include/iovelocity.h" /* also includes ../include/floatbox.h */

#include <stdio.h> /* for stdout, not necessary otherwise */

const char *text_velocity_file = "../docs/velocity-241-241-51.txt";
const char *vbox_output_file = "velocities.vbox";
const int x = 101, y = 105, z = 25; // test coordinates

int main()
{
    struct VELOCITYBOX vbox;

    // load a text-format velocity file
    printf( "loading old velocity model %s...", text_velocity_file ); fflush( stdout );
    if( !vboxloadtext( &vbox, text_velocity_file ) ) { /* handle error */ }
    printf( " done.\n" ); fflush( stdout );

    // print the old vbox metadata
    vboxfprint( stdout, "old: ", "\t", vbox );

    // read one value from the vbox in memory
    float oldval = boxget( vbox.box, x, y, z );
    printf( "value at (%d, %d, %d) = %g\n", x, y, z, oldval ); 
    
    // write one value to the vbox in memory
    float newval = oldval * 3.5f + 1.f;
    printf( "modifying value at (%d, %d, %d) to %g\n", x, y, z, newval );
    boxput( vbox.box, x, y, z, newval );
    
    // write modified vbox to new file using new vbox file format
    printf( "storing modified velocity model %s...", vbox_output_file ); fflush( stdout );
    if( !vboxstorebinary( vbox_output_file, vbox ) ) { /* handle error */ }
    printf( " done.\n" ); fflush( stdout );

    // free heap memory used by the velocity box
    vboxfree( &vbox );

    // load the vbox file we just created
    printf( "loading modified velocity model %s...", vbox_output_file); fflush( stdout );
    if( !vboxloadbinary( &vbox, vbox_output_file ) ) { /* handle error */ }
    printf( " done.\n" ); fflush( stdout );

    // print the new vbox metadata
    vboxfprint( stdout, "example: ", "\t", vbox );

    // show that the modified value is present in the new file
    printf( "value at (%d, %d, %d) = %g\n", x, y, z, boxget( vbox.box, x, y, z ) );

    // free heap memory used by the velocity box
    vboxfree( &vbox );

    printf( "all done.\n" );

    return 0;
}
