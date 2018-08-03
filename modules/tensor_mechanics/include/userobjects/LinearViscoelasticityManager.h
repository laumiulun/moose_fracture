/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef LINEARVISCOELASTICITYMANAGER_H
#define LINEARVISCOELASTICITYMANAGER_H

#include "ElementUserObject.h"
#include "LinearViscoelasticityBase.h"
#include "RankTwoTensor.h"

class LinearViscoelasticityManager;

template <>
InputParameters validParams<LinearViscoelasticityManager>();

/**
 * This class manages a LinearViscoelasticityBase object. Its primary purpose
 * is to initialize the internal MaterialProperties contained in the viscoelastic
 * model at the beginning of each time step, and update those properties at the
 * end of each time step.
 *
 * Whenever a LinearViscoelasticityBase object is created, it must be associated to
 * one LinearViscoelasticityManager user object, otherwise the viscoelastic
 * creep strains and properties will not be computed accordingly.
 *
 * See LinearViscoelasticityBase for more information.
 */
class LinearViscoelasticityManager : public ElementUserObject
{
public:
  LinearViscoelasticityManager(const InputParameters & parameters);

protected:
  virtual void initialize() override;
  virtual void execute() override;
  virtual void threadJoin(const UserObject & /*uo*/) override{};
  virtual void finalize() override{};

  std::string _stress_name;
  /*
   * The effective stress used for the update of the viscoelastic strain
   * (typically, the real stress).
   */
  const MaterialProperty<RankTwoTensor> & _stress;

  std::string _creep_strain_name;
  /// Name of the creep strain variable used for the update of the viscoelastic strain
  const MaterialProperty<RankTwoTensor> & _creep_strain;

  std::string _elastic_strain_name;
  /// Name of the elastic strain variable used for the update of the viscoelastic strain
  const MaterialProperty<RankTwoTensor> & _elastic_strain;

  /// Name of the viscoelastic model to update
  std::string _viscoelastic_model_name;
  /// Pointer to the viscoelastic model to update
  std::shared_ptr<LinearViscoelasticityBase> _viscoelastic_model;
};

#endif // STRESSTIMESTEPSETUP_H
