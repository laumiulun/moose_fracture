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

#ifndef PROCESSORIDAUX_H
#define PROCESSORIDAUX_H

#include "AuxKernel.h"

// Forward Declarations
class ProcessorIDAux;

template <>
InputParameters validParams<ProcessorIDAux>();

class ProcessorIDAux : public AuxKernel
{
public:
  ProcessorIDAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
};

#endif // PROCESSORIDAUX_H
