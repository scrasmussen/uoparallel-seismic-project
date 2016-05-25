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


int
mpifindneighborrank (
    int my_rank,
    int ghost_id,
    int width
)
// returns mpi rank for a ghost cell id
// see ^Ghost cell ids
// - North -> i+d if < d^2 else -1
// - South -> i-d
// - West -> -1 if (i % d == 0) else i-1
// - East -> -1 if (i % d == d - 1) else i+1
// - NE -> North + 1 if (North > 0 and East > 0) else -1
// - NW -> North - 1 if (North > 0 and West > 0) else -1
// - SE -> South + 1 if (South > 0 and East > 0) else -1
// - SW -> South - 1 if (South > 0 and West > 0) else -1
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
