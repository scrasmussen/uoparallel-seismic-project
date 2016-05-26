// example_iovelocity.c - Atlee Brink
// A working example of all the functions in ../include/iovelocity.h,
// and shows how to get and put values from the VELOCITYBOX (through its FLOATBOX).

#include "../include/velocityboxfiler.h" /* also includes ../include/floatbox.h */

#include <stdio.h> /* for stdout, not necessary otherwise */

const char *text_velocity_file = "../docs/velocity-241-241-51.txt";
const char *vbox_output_file = "velocities.vbox";
const int x = 101, y = 105, z = 25; // test coordinates

int main()
{
    struct VELOCITYBOX vbox;

////////////////////////////////////////////////////////////////////////////////
// reading text-format file
////////////////////////////////////////////////////////////////////////////////

    // load a text-format velocity file
    printf( "loading old velocity model %s...", text_velocity_file ); fflush( stdout );
    if( !vbfileloadtext( &vbox, text_velocity_file ) ) { /* handle error */ }
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

////////////////////////////////////////////////////////////////////////////////
// writing vbox-format file
////////////////////////////////////////////////////////////////////////////////
    
    // write modified vbox to new file using new vbox file format
    printf( "storing modified velocity model %s...", vbox_output_file ); fflush( stdout );
    if( !vbfilestorebinary( vbox_output_file, vbox ) ) { /* handle error */ }
    printf( " done.\n" ); fflush( stdout );

    // free heap memory used by the velocity box
    vboxfree( &vbox );

////////////////////////////////////////////////////////////////////////////////
// reading whole vbox-format file
////////////////////////////////////////////////////////////////////////////////

    // load the vbox file we just created
    printf( "loading modified velocity model %s...", vbox_output_file); fflush( stdout );
    if( !vbfileloadbinary( &vbox, vbox_output_file ) ) { /* handle error */ }
    printf( " done.\n" ); fflush( stdout );

    // print the new vbox metadata
    vboxfprint( stdout, "new: ", "\t", vbox );

    // show that the modified value is present in the new file
    printf( "value at (%d, %d, %d) = %g\n", x, y, z, boxget( vbox.box, x, y, z ) );

    // free heap memory used by the velocity box
    vboxfree( &vbox );

////////////////////////////////////////////////////////////////////////////////
// reading subset of vbox-format file
////////////////////////////////////////////////////////////////////////////////

    struct VBOXOPENFILE vbfile;

    // open (don't read whole file yet)
    printf( "opening modified velocity model %s...", vbox_output_file); fflush( stdout );
    if( !vbfileopenbinary( &vbfile, vbox_output_file ) ) { /* handle error */ }
    printf( " done.\n" ); fflush( stdout );

    // display info about the file
    printf( "origin: (%d, %d, %d)\n", vbfile.min.x, vbfile.min.y, vbfile.min.z );
    printf( "size: (%d, %d, %d)\n", vbfile.dims.x, vbfile.dims.y, vbfile.dims.z );

    // subset region
    int ox = vbfile.min.x + vbfile.dims.x / 3;
    int oy = vbfile.min.y + vbfile.dims.y / 3;
    int oz = vbfile.min.z + vbfile.dims.z / 3;
    int nx = vbfile.dims.x / 3;
    int ny = vbfile.dims.y / 3;
    int nz = vbfile.dims.z / 3;

    printf( "subset origin: (%d, %d, %d)\n", ox, oy, oz );
    printf( "subset size: (%d, %d, %d)\n", nx, ny, nz );

    // load a subset of the file into a VELOCITYBOX
    if( !vbfileloadbinarysubset( &vbox, ox, oy, oz, nx, ny, nz, vbfile ) ) { /* handle error */ }

    // close the open vbfile
    vbfileclosebinary( &vbfile );

    // display value we changed before
    printf( "in loaded subset: value at (%d, %d, %d) = %g\n",
        x, y, z, boxget( vbox.box, x - ox, y - oy, z - oz ) );

    // free heap memory used by the velocity box
    vboxfree( &vbox );


    printf( "all done.\n" );

    return 0;
}
