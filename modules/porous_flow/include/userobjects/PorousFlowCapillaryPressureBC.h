/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOWCAPILLARYPRESSUREBC_H
#define POROUSFLOWCAPILLARYPRESSUREBC_H

#include "PorousFlowCapillaryPressure.h"

class PorousFlowCapillaryPressureBC;

template <>
InputParameters validParams<PorousFlowCapillaryPressureBC>();

/**
 * Brooks-Corey effective saturation, capillary pressure and relative
 * permeability functions.
 *
 * From Brooks, R. H. and A. T. Corey (1966), Properties of porous media affecting
 * fluid flow, J. Irrig. Drain. Div., 6, 61
 */
class PorousFlowCapillaryPressureBC : public PorousFlowCapillaryPressure
{
public:
  PorousFlowCapillaryPressureBC(const InputParameters & parameters);

  virtual Real capillaryPressureCurve(Real saturation) const override;
  virtual Real dCapillaryPressureCurve(Real saturation) const override;
  virtual Real d2CapillaryPressureCurve(Real saturation) const override;

  virtual Real effectiveSaturation(Real pc) const override;
  virtual Real dEffectiveSaturation(Real pc) const override;
  virtual Real d2EffectiveSaturation(Real pc) const override;

protected:
  /// Brooks-Corey exponent lambda
  const Real _lambda;
  /// Threshold entry pressure
  const Real _pe;
};

#endif // POROUSFLOWCAPILLARYPRESSUREBC_H
