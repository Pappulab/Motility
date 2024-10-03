// clang-format off
/* ----------------------------------------------------------------------
   Copied this from fix_langevin

   Arguments to the updated command would be:
   // fix active FORCE_MAGNITUDE seed FORCE_TYPE (RTP or ABP) D_R(for ABP)


/* ----------------------------------------------------------------------
   Contributing authors: Gaurav Chauhan, WUSTL
------------------------------------------------------------------------- */

#include "fix_directional.h"
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
// fix directional f_active dir_x dir_y dir_z
FixDirectional::FixDirectional(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),random(nullptr)
{
 if ((strcmp(style,"directional") == 0) && narg >6)
 {
      f_active = utils::numeric(FLERR,arg[3],false,lmp);
      
      //Unit vector dir_x, dir_y, dir_z
      dir_x=utils::numeric(FLERR,arg[4],false,lmp); 
      dir_y=utils::numeric(FLERR,arg[5],false,lmp);
      dir_z=utils::numeric(FLERR,arg[6],false,lmp);
  }
  else 
  {
    error->all(FLERR, "Too many arguments");
  }
  
}


/* ---------------------------------------------------------------------- */

int FixDirectional::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDirectional::post_force(int /*vflag*/)
{
  double u_x, u_y, u_z;

  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;


    for (int n = 0; n < nlocal; n++)
      if(mask[n] & groupbit)
     {

        f[n][0] += (-1*f_active*dir_x);  // 
        f[n][1] += (-1*f_active*dir_y);
        f[n][2] += (-1*f_active*dir_z);
      }
  
}
