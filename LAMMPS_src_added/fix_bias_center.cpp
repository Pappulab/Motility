// clang-format off
/* ----------------------------------------------------------------------
   Copied this from fix_langevin



/* ----------------------------------------------------------------------
   Contributing authors: Gaurav Chauhan, WUSTL
------------------------------------------------------------------------- */

#include "fix_bias_center.h"
#include "comm.h"
#include "stdio.h"
#include "string.h"
#include "atom.h"
#include "neighbor.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "math.h"
#include "random_mars.h"

using namespace LAMMPS_NS;
using namespace FixConst;



/* ---------------------------------------------------------------------- */

// example command:
// fix active FORCE_MAGNITUDE seed FORCE_TYPE (RTP or ABP) D_R(for ABP)
FixBias::FixBias(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"bias_center") != 0 && narg < 4)
    error->all(FLERR,"Illegal fix active command: not enough args");
  f_bias = utils::numeric(FLERR,arg[3],false,lmp);
  r_center=utils::numeric(FLERR,arg[4],false,lmp);
  
  }



/* ---------------------------------------------------------------------- */

int FixBias::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBias::post_force(int /*vflag*/)
{
  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **f = atom->f;

    for (int n = 0; n < nlocal; n++) {

        f[n][0] += (-1*f_bias*(x[n][0]-r_center));  // 
        f[n][1] += (-1*f_bias*(x[n][1]-r_center));
        f[n][2] += (-1*f_bias*(x[n][2]-r_center));
      }

}
