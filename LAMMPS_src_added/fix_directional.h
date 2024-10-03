/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(directional,FixDirectional)

#else

#ifndef LMP_FIX_DIRECTIONAL_H
#define LMP_FIX_DIRECTIONAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDirectional : public Fix {
 public:
  FixDirectional(class LAMMPS *, int, char **);
  virtual ~FixDirectional() {}
  int setmask();
  virtual void post_force(int);

 protected:
  class RanMars *random;
  double f_active;
  double dir_x,dir_y,dir_z;
};

}

#endif
#endif
