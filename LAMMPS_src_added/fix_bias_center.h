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

FixStyle(bias_center,FixBias)

#else

#ifndef LMP_FIX_BIAS_CENTER_H
#define LMP_FIX_BIAS_CENTER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBias : public Fix {
 public:
  FixBias(class LAMMPS *, int, char **);
  virtual ~FixBias() {}
  int setmask();
  virtual void post_force(int);

 protected:

  double f_bias;
  double r_center;
};

}

#endif
#endif
