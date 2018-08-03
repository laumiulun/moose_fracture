/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWTOTALGRAVITATIONALDENSITYFULLYSATURATEDFROMPOROSITY_H
#define POROUSFLOWTOTALGRAVITATIONALDENSITYFULLYSATURATEDFROMPOROSITY_H

#include "PorousFlowTotalGravitationalDensityBase.h"

// Forward Declarations
class PorousFlowTotalGravitationalDensityFullySaturatedFromPorosity;

template <>
InputParameters validParams<PorousFlowTotalGravitationalDensityFullySaturatedFromPorosity>();

/**
 * Material designed to provide the density of the porous medium for the
 * fully-saturated case. Density is calculated as a
 * weighted average of the fluid and solid densities:
 * density = phi * rho_f + (1 - phi) * rho_s
 * where phi is porosity, rho_f is fluid density and rho_s is solid
 * density (assumed constant).
 */
class PorousFlowTotalGravitationalDensityFullySaturatedFromPorosity
    : public PorousFlowTotalGravitationalDensityBase
{
public:
  PorousFlowTotalGravitationalDensityFullySaturatedFromPorosity(const InputParameters & parameters);

protected:
  void initQpStatefulProperties() override;
  void computeQpProperties() override;

  /// Solid density
  const Real _rho_s;

  /// fluid density at qps
  const MaterialProperty<std::vector<Real>> & _rho_f_qp;

  /// porosity at qps
  const MaterialProperty<Real> & _porosity_qp;

  /// d(rho_f)/d(PorousFlow variable)
  const MaterialProperty<std::vector<std::vector<Real>>> & _drho_f_qp_dvar;

  /// d(porosity)/d(PorousFlow variable)
  const MaterialProperty<std::vector<Real>> & _dporosity_qp_dvar;
};

#endif // POROUSFLOWTOTALGRAVITATIONALDENSITYFULLYSATURATEDFROMPOROSITY_H
