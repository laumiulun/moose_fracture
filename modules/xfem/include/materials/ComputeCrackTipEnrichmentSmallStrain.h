/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTECRACKTIPENRICHMENTSMALLSTRAIN_H
#define COMPUTECRACKTIPENRICHMENTSMALLSTRAIN_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "ComputeStrainBase.h"
#include "Assembly.h"
#include "CrackFrontDefinition.h"
#include "EnrichmentFunctionCalculation.h"

/**
 * ComputeCrackTipEnrichmentSmallStrain calculates the sum of standard strain and enrichement strain
 */
class ComputeCrackTipEnrichmentSmallStrain : public ComputeStrainBase,
                                             public EnrichmentFunctionCalculation
{
public:
  ComputeCrackTipEnrichmentSmallStrain(const InputParameters & parameters);
  virtual ~ComputeCrackTipEnrichmentSmallStrain() {}

protected:
  virtual void computeProperties() override;

  virtual void computeQpProperties() override;

  /// enrichment displacement
  std::vector<Real> _enrich_disp;

  /// gradient of enrichment displacement
  std::vector<RealVectorValue> _grad_enrich_disp;

  /// enrichment displacement variables
  std::vector<std::vector<MooseVariable *>> _enrich_variable;

  /// the current shape functions
  const VariablePhiValue & _phi;

  /// gradient of the shape function
  const VariablePhiGradient & _grad_phi;

private:
  /// enrichment function value
  std::vector<Real> _B;
  /// derivatives of enrichment function respect to global cooridnate
  std::vector<RealVectorValue> _dBX;
  /// derivatives of enrichment function respect to crack front cooridnate
  std::vector<RealVectorValue> _dBx;
  /// enrichment function at node I
  std::vector<std::vector<Real>> _BI;
  /// shape function
  const std::vector<std::vector<Real>> * _fe_phi;
  /// gradient of shape function
  const std::vector<std::vector<RealGradient>> * _fe_dphi;
  NonlinearSystem * _nl;
  const NumericVector<Number> * _sln;
};

#endif // COMPUTECRACKTIPENRICHMENTSMALLSTRAIN_H
