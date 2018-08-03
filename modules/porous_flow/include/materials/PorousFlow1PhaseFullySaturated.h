/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOW1PHASEFULLYSATURATED_H
#define POROUSFLOW1PHASEFULLYSATURATED_H

#include "PorousFlowVariableBase.h"

// Forward Declarations
class PorousFlow1PhaseFullySaturated;

template <>
InputParameters validParams<PorousFlow1PhaseFullySaturated>();

/**
 * Base material designed to calculate fluid phase porepressure and saturation
 * for the single-phase situation assuming full saturation where porepressure
 * is the nonlinear variable.
 */
class PorousFlow1PhaseFullySaturated : public PorousFlowVariableBase
{
public:
  PorousFlow1PhaseFullySaturated(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /**
   * Assemble std::vectors of porepressure, saturation and temperature at the quadpoints
   */
  void buildQpPPSS();

  /// Nodal or quadpoint value of porepressure of the fluid phase
  const VariableValue & _porepressure_var;
  /// Gradient(_porepressure at quadpoints)
  const VariableGradient & _gradp_qp_var;
  /// Moose variable number of the porepressure
  const unsigned int _porepressure_varnum;
  /// the PorousFlow variable number of the porepressure
  const unsigned int _p_var_num;
};

#endif // POROUSFLOW1PHASEFULLYSATURATED_H
