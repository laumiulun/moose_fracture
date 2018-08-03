/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEAXISYMMETRIC1DFINITESTRAIN_H
#define COMPUTEAXISYMMETRIC1DFINITESTRAIN_H

#include "Compute1DFiniteStrain.h"
#include "SubblockIndexProvider.h"

/**
 * ComputeAxisymmetric1DFiniteStrain defines a strain increment for finite strains
 * in an Axisymmetric 1D problem. The COORD_TYPE in the Problem block must be set to RZ.
 */
class ComputeAxisymmetric1DFiniteStrain : public Compute1DFiniteStrain
{
public:
  ComputeAxisymmetric1DFiniteStrain(const InputParameters & parameters);

protected:
  void initialSetup() override;

  /// Computes the current dUy/dy for axisymmetric problems
  Real computeGradDispYY() override;

  /// Computes the old dUy/dy for axisymmetric problems
  Real computeGradDispYYOld() override;

  /// Computes the current dUz/dz for axisymmetric problems, where
  /// \f$ \epsilon_{\theta} = \frac{u_r}{r} \f$
  Real computeGradDispZZ() override;

  /// Computes the old dUz/dz for axisymmetric problems, where
  /// \f$ \epsilon_{\theta-old} = \frac{u_{r-old}}{r_{old}} \f$
  Real computeGradDispZZOld() override;

  /// gets its subblock index for current element
  unsigned int getCurrentSubblockIndex() const
  {
    return _subblock_id_provider ? _subblock_id_provider->getSubblockIndex(*_current_elem) : 0;
  };

  /// the old value of the first component of the displacements vector
  const VariableValue & _disp_old_0;

  const SubblockIndexProvider * _subblock_id_provider;

  bool _has_out_of_plane_strain;
  const VariableValue & _out_of_plane_strain;
  const VariableValue & _out_of_plane_strain_old;

  bool _has_scalar_out_of_plane_strain;
  unsigned int _nscalar_strains;
  std::vector<const VariableValue *> _scalar_out_of_plane_strain;
  std::vector<const VariableValue *> _scalar_out_of_plane_strain_old;
};

#endif // COMPUTEAXISYMMETRIC1DFINITESTRAIN_H
