////////////////////////////////////////////////////////////////////////////////
// mpihelpers.h - 2016.05.25 - Cody Hatfield
////////////////////////////////////////////////////////////////////////////////
//
// Ghost cell ids:
//   Starting from SW, counter-clockwise
//   6 5 4
//   7   3
//   0 1 2
//
////////////////////////////////////////////////////////////////////////////////

#include "point3d.h"
#include "splitsquare.h"

void
mpicalculatemycoordinates (
    struct POINT3D *min,
    struct POINT3D *dims,
    int my_rank,
    int total_nodes,
    struct POINT3D src_min,
    struct POINT3D src_dims
)
// returns the global coordinates for a particular rank
// assumes z is depth
// assumes total_nodes is power of 2
{
    int num_x = splitsquare_numx(total_nodes);
    int num_y = (total_nodes/num_x);

    int dims_x = (src_dims.x/num_x);
    int dims_y = (src_dims.y/num_y);

    int col_rank = my_rank % num_x;
    int row_rank = my_rank / num_x;

    min->x = (dims_x * col_rank) + src_min.x;
    dims->x = col_rank < (num_x - 1) ? dims_x : src_dims.x - (dims_x * (num_x - 1));
    min->y = (dims_y * row_rank) + src_min.y;
    dims->y = row_rank < (num_y - 1) ? dims_y : src_dims.y - (dims_y * (num_y - 1));
    min->z = src_min.z;
    dims->z = src_dims.z;
}

int
mpifindneighborrank (
    int my_rank,
    int ghost_id,   // see ^Ghost cell ids
    int width
)
// returns mpi rank for a ghost cell id
// < 0 is no neighbor
{
  int north     = (my_rank+width < width*width) ? my_rank+width : -1;
  int south     = (my_rank-width);
  int west      = (my_rank % width == 0) ? -1 : my_rank-1;
  int east      = (my_rank % width == width - 1) ? -1 : my_rank+1;

  switch (ghost_id) {
    case 0: // SW
      return (south >= 0 && west >= 0) ? (south-1) : -1;
      break;
    case 1: // S
      return south;
      break;
    case 2: // SE
      return (south >= 0 && east >= 0) ? (south+1) : -1;
      break;
    case 3: // E
      return east;
      break;
    case 4: // NE
      return (north >= 0 && east >= 0) ? (north+1) : -1;
      break;
    case 5: // N
      return north;
      break;
    case 6: // NW
      return (north >= 0 && west >= 0) ? (north-1) : -1;
      break;
    case 7: // W
      return west;
      break;
    default:
      return -1;
      break;
  }
}



void
mpigetsendcoordinates(
    struct POINT3D *min,
    struct POINT3D *dims,
    int ghost_id,       // see ^Ghost cell ids
    int ghost_size,      // probably want 7
    int min_x,
    int min_y,
    int max_x,
    int max_y,
    int z
)
// returns send coordinates for a given ghost cell and dims
{
  if( min == 0 ) return;
  if( dims == 0 ) return;

  switch (ghost_id) {
    case 0:
      min->x = min_x-ghost_size;
      min->y = min_y-ghost_size;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 1:
      min->x = min_x;
      min->y = min_y-ghost_size;
      dims->x = max_x-min_x;
      dims->y = ghost_size;
      break;
    case 2:
      min->x = max_x;
      min->y = min_y-ghost_size;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 3:
      min->x = max_x;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = max_y-min_y;
      break;
    case 4:
      min->x = max_x;
      min->y = max_y;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 5:
      min->x = min_x;
      min->y = max_y;
      dims->x = max_x-min_x;
      dims->y = ghost_size;
      break;
    case 6:
      min->x = min_x-ghost_size;
      min->y = max_y;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 7:
      min->x = min_x-ghost_size;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = max_y-min_y;
      break;
  }
  min->z = 0;
  dims->z = z;
}


void
mpigetreceivecoordinates(
    struct POINT3D *min,
    struct POINT3D *dims,
    int ghost_id,       // see ^Ghost cell ids
    int ghost_size,      // probably want 7
    int min_x,
    int min_y,
    int max_x,
    int max_y,
    int z
)
// returns receive coordinates for a given ghost cell and dims
{
  if( min == 0 ) return;

  switch (ghost_id) {
    case 0:
      min->x = min_x;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 1:
      min->x = min_x;
      min->y = min_y;
      dims->x = max_x-min_x;
      dims->y = ghost_size;
      break;
    case 2:
      min->x = max_x-ghost_size;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 3:
      min->x = max_x-ghost_size;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = max_y-min_y;
      break;
    case 4:
      min->x = max_x-ghost_size;
      min->y = max_y-ghost_size;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 5:
      min->x = min_x;
      min->y = max_y-ghost_size;
      dims->x = max_x-min_x;
      dims->y = ghost_size;
      break;
    case 6:
      min->x = min_x;
      min->y = max_y-ghost_size;
      dims->x = ghost_size;
      dims->y = ghost_size;
      break;
    case 7:
      min->x = min_x;
      min->y = min_y;
      dims->x = ghost_size;
      dims->y = max_y-min_y;
      break;
  }
  min->z = 0;
  dims->z = z;
}
