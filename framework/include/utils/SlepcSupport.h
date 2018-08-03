/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef SLEPCSUPPORT_H
#define SLEPCSUPPORT_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_SLEPC

#include "Moose.h"

class EigenProblem;
class InputParameters;

namespace Moose
{
namespace SlepcSupport
{
/**
 * @return InputParameters object containing the SLEPC related parameters
 *
 * The output of this function should be added to the the parameters object of the overarching class
 * @see EigenProblem
 */
InputParameters getSlepcValidParams(InputParameters & params);
InputParameters getSlepcEigenProblemValidParams();
void storeSlepcOptions(FEProblemBase & fe_problem, const InputParameters & params);
void storeSlepcEigenProblemOptions(EigenProblem & eigen_problem, const InputParameters & params);
void slepcSetOptions(EigenProblem & eigen_problem, const InputParameters & params);
void setSlepcEigenSolverTolerances(EigenProblem & eigen_problem, const InputParameters & params);
void setSlepcOutputOptions(EigenProblem & eigen_problem);

PetscErrorCode mooseSlepcEigenFormJacobianA(SNES snes, Vec x, Mat jac, Mat pc, void * ctx);
PetscErrorCode mooseSlepcEigenFormJacobianB(SNES snes, Vec x, Mat jac, Mat pc, void * ctx);
PetscErrorCode mooseSlepcEigenFormFunctionA(SNES snes, Vec x, Vec r, void * ctx);
PetscErrorCode mooseSlepcEigenFormFunctionB(SNES snes, Vec x, Vec r, void * ctx);
} // namespace SlepcSupport
} // namespace moose

#endif // LIBMESH_HAVE_SLEPC

#endif // SLEPCSUPPORT_H
