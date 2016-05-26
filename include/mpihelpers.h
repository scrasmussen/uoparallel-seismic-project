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

struct TempPoint{
  int x, y, z;
  int x_size, y_size, z_size;
};

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
    int ghost_id,       // see ^Ghost cell ids
    int ghost_size,      // probably want 7
    int min_x,
    int min_y,
    int max_x,
    int max_y,
    int z_size,
    struct TempPoint *tp
)
// returns send coordinates for a given ghost cell and size
{
  if( tp == 0 ) return;

  switch (ghost_id) {
    case 0:
      tp->x = min_x-ghost_size;
      tp->y = min_y-ghost_size;
      tp->x_size = ghost_size;
      tp->y_size = ghost_size;
      break;
    case 1:
      tp->x = min_x;
      tp->y = min_y-ghost_size;
      tp->x_size = max_x-min_x;
      tp->y_size = ghost_size;
      break;
    case 2:
      tp->x = max_x;
      tp->y = min_y-ghost_size;
      tp->x_size = ghost_size;
      tp->y_size = ghost_size;
      break;
    case 3:
      tp->x = max_x;
      tp->y = min_y;
      tp->x_size = ghost_size;
      tp->y_size = max_y-min_y;
      break;
    case 4:
      tp->x = max_x;
      tp->y = max_y;
      tp->x_size = ghost_size;
      tp->y_size = ghost_size;
      break;
    case 5:
      tp->x = min_x;
      tp->y = max_y;
      tp->x_size = max_x-min_x;
      tp->y_size = ghost_size;
      break;
    case 6:
      tp->x = min_x-ghost_size;
      tp->y = max_y;
      tp->x_size = ghost_size;
      tp->y_size = ghost_size;
      break;
    case 7:
      tp->x = min_x-ghost_size;
      tp->y = min_y;
      tp->x_size = ghost_size;
      tp->y_size = max_y-min_y;
      break;
  }
  tp->z = 0;
  tp->z_size = z_size;
}


void
mpigetreceivecoordinates(
    int ghost_id,       // see ^Ghost cell ids
    int ghost_size,      // probably want 7
    int min_x,
    int min_y,
    int max_x,
    int max_y,
    int z_size,
    struct TempPoint *tp
)
// returns receive coordinates for a given ghost cell and size
{
  if( tp == 0 ) return;

  switch (ghost_id) {
    case 0:
      tp->x = min_x;
      tp->y = min_y;
      tp->x_size = ghost_size;
      tp->y_size = ghost_size;
      break;
    case 1:
      tp->x = min_x;
      tp->y = min_y;
      tp->x_size = max_x-min_x;
      tp->y_size = ghost_size;
      break;
    case 2:
      tp->x = max_x-ghost_size;
      tp->y = min_y;
      tp->x_size = ghost_size;
      tp->y_size = ghost_size;
      break;
    case 3:
      tp->x = max_x-ghost_size;
      tp->y = min_y;
      tp->x_size = ghost_size;
      tp->y_size = max_y-min_y;
      break;
    case 4:
      tp->x = max_x-ghost_size;
      tp->y = max_y-ghost_size;
      tp->x_size = ghost_size;
      tp->y_size = ghost_size;
      break;
    case 5:
      tp->x = min_x;
      tp->y = max_y-ghost_size;
      tp->x_size = max_x-min_x;
      tp->y_size = ghost_size;
      break;
    case 6:
      tp->x = min_x;
      tp->y = max_y-ghost_size;
      tp->x_size = ghost_size;
      tp->y_size = ghost_size;
      break;
    case 7:
      tp->x = min_x;
      tp->y = min_y;
      tp->x_size = ghost_size;
      tp->y_size = max_y-min_y;
      break;
  }
  tp->z = 0;
  tp->z_size = z_size;
}
