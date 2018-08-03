/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CAPPEDWEAKPLANECOSSERATSTRESSUPDATE_H
#define CAPPEDWEAKPLANECOSSERATSTRESSUPDATE_H

#include "CappedWeakPlaneStressUpdate.h"

class CappedWeakPlaneCosseratStressUpdate;

template <>
InputParameters validParams<CappedWeakPlaneCosseratStressUpdate>();

/**
 * CappedWeakPlaneCosseratStressUpdate performs the return-map
 * algorithm and associated stress updates for plastic
 * models that describe capped weak-plane Cosserat plasticity
 *
 * It assumes various things about the elasticity tensor, viz
 * E(i,i,j,k) = 0 except if k=j
 * E(0,0,i,j) = E(1,1,i,j)
 */
class CappedWeakPlaneCosseratStressUpdate : public CappedWeakPlaneStressUpdate
{
public:
  CappedWeakPlaneCosseratStressUpdate(const InputParameters & parameters);

  /**
   * Does the model require the elasticity tensor to be isotropic?
   */
  bool requiresIsotropicTensor() override { return false; }

protected:
  virtual void consistentTangentOperator(const RankTwoTensor & stress_trial,
                                         Real p_trial,
                                         Real q_trial,
                                         const RankTwoTensor & stress,
                                         Real p,
                                         Real q,
                                         Real gaE,
                                         const yieldAndFlow & smoothed_q,
                                         const RankFourTensor & Eijkl,
                                         bool compute_full_tangent_operator,
                                         RankFourTensor & cto) const override;

  virtual void setStressAfterReturn(const RankTwoTensor & stress_trial,
                                    Real p_ok,
                                    Real q_ok,
                                    Real gaE,
                                    const std::vector<Real> & intnl,
                                    const yieldAndFlow & smoothed_q,
                                    const RankFourTensor & Eijkl,
                                    RankTwoTensor & stress) const override;

  virtual RankTwoTensor dqdstress(const RankTwoTensor & stress) const override;

  virtual RankFourTensor d2qdstress2(const RankTwoTensor & stress) const override;
};

#endif // CAPPEDWEAKPLANECOSSERATSTRESSUPDATE_H
