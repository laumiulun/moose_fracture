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

#ifndef NONLINEARSYSTEM_H
#define NONLINEARSYSTEM_H

#include "NonlinearSystemBase.h"

/**
 * Nonlinear system to be solved
 *
 * It is a part of FEProblemBase ;-)
 */
class NonlinearSystem : public NonlinearSystemBase
{
public:
  NonlinearSystem(FEProblemBase & problem, const std::string & name);
  virtual ~NonlinearSystem();

  virtual void solve() override;

  /**
   * Quit the current solve as soon as possible.
   */
  virtual void stopSolve() override;

  /**
   * Returns the current nonlinear iteration number.  In libmesh, this is
   * updated during the nonlinear solve, so it should be up-to-date.
   */
  virtual unsigned int getCurrentNonlinearIterationNumber() override
  {
    return _transient_sys.get_current_nonlinear_iteration_number();
  }

  virtual void setupFiniteDifferencedPreconditioner() override;

  /**
   * Returns the convergence state
   * @return true if converged, otherwise false
   */
  virtual bool converged() override;

  virtual NumericVector<Number> & RHS() override { return *_transient_sys.rhs; }

  virtual NonlinearSolver<Number> * nonlinearSolver() override
  {
    return _transient_sys.nonlinear_solver.get();
  }

  virtual NumericVector<Number> & solutionOld() override
  {
    return *_transient_sys.old_local_solution;
  }

  virtual NumericVector<Number> & solutionOlder() override
  {
    return *_transient_sys.older_local_solution;
  }

  virtual TransientNonlinearImplicitSystem & sys() { return _transient_sys; }

protected:
  TransientNonlinearImplicitSystem & _transient_sys;

private:
  /**
  * Form preconditioning matrix via a standard finite difference method
  * column-by-column. This method computes both diagonal and off-diagonal
  * entrices regardless of the structure pattern of the Jacobian matrix.
  */
  void setupStandardFiniteDifferencedPreconditioner();

  /**
  * According to the nonzero pattern provided in the matrix, a graph is constructed.
  * A coloring algorithm is applied to the graph. The graph is partitioned into several
  * independent subgraphs (colors), and a finte difference method is applied color-by-color
  * to form a preconditioning matrix. If the number of colors is small, this method is much
  * faster than the standard FD. But there is an issue. If the matrix provided by users does not
  * represent the actual structure of the true Jacobian, the matrix computed via coloring could
  * be wrong or inaccurate. In this case, users should switch to the standard finite difference
  * method.
  */
  void setupColoringFiniteDifferencedPreconditioner();

  bool _use_coloring_finite_difference;
};

#endif /* NONLINEARSYSTEM_H */
