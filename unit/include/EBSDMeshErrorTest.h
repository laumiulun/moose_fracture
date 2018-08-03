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

#ifndef EBSDMESHERRORTEST_H
#define EBSDMESHERRORTEST_H

// CPPUnit includes
#include "gtest/gtest.h"

// Moose includes
#include "EBSDMesh.h"
#include "InputParameters.h"
#include "MooseParsedFunction.h"
#include "MooseUnitApp.h"
#include "AppFactory.h"

class EBSDMeshErrorTest : public ::testing::Test
{
protected:
  void SetUp()
  {
    const char * argv[2] = {"foo", "\0"};
    _app = AppFactory::createAppShared("MooseUnitApp", 1, (char **)argv);
    _factory = &_app->getFactory();
  }

  template <typename T>
  void testParam(unsigned int nparam, const char ** param_list, std::string name)
  {
    for (unsigned int i = 0; i < nparam; ++i)
    {
      // create a unique name
      std::ostringstream oss;
      oss << name << "_" << i;

      // generate input parameter set
      InputParameters params = validParams<EBSDMesh>();
      params.addPrivateParam("_moose_app", _app.get());
      params.set<std::string>("_object_name") = oss.str();

      // set a single parameter
      params.set<T>(param_list[i]) = T(1.0);

      // set filename (is a required param but not used in these tests)
      params.set<FileName>("filename") = "DUMMY";

      try
      {
        // construct mesh object
        std::unique_ptr<EBSDMesh> mesh = libmesh_make_unique<EBSDMesh>(params);
        // TODO: fix and uncomment this - it was missing before.
        // FAIL() << "mesh construction should have failed but didn't";
      }
      catch (const std::exception & e)
      {
        std::string msg(e.what());
        ASSERT_TRUE(
            msg.find("Do not specify mesh geometry information, it is read from the EBSD file.") !=
            std::string::npos)
            << "failed with unexpected error: " << msg;
      }
    }
  }

  std::shared_ptr<MooseApp> _app;
  Factory * _factory;
};

#endif // EBSDMESHERRORTEST_H
