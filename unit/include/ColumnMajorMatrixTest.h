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

#ifndef COLUMNMAJORMATRIXTEST_H
#define COLUMNMAJORMATRIXTEST_H

// CPPUnit includes
#include "gtest/gtest.h"

// Moose includes
#include "ColumnMajorMatrix.h"

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

class ColumnMajorMatrixTest : public ::testing::Test
{
protected:
  /**
   * Initialize the various test objects to the right size.
   */
  ColumnMajorMatrixTest()
    : a(3, 3), t(3, 2), two_mat(2, 2), add(3, 2), add_solution(3, 2), sub(3, 2), sub_solution(3, 2)
  {
  }

  void SetUp()
  {
    // Define commonly used matrices for testing
    a(0, 0) = 1;
    a(1, 0) = 2;
    a(2, 0) = 3;
    a(0, 1) = 4;
    a(1, 1) = 5;
    a(2, 1) = 6;
    a(0, 2) = 7;
    a(1, 2) = 8;
    a(2, 2) = 9;

    t(0, 0) = 1;
    t(1, 0) = 2;
    t(2, 0) = 3;
    t(0, 1) = 4;
    t(1, 1) = 5;
    t(2, 1) = 6;

    two_mat(0, 0) = 1;
    two_mat(1, 0) = 2;
    two_mat(0, 1) = 3;
    two_mat(1, 1) = 4;

    add(0, 0) = 6;
    add(1, 0) = 5;
    add(2, 0) = 4;
    add(0, 1) = 1;
    add(1, 1) = 1;
    add(2, 1) = 1;

    add_solution(0, 0) = 7;
    add_solution(1, 0) = 7;
    add_solution(2, 0) = 7;
    add_solution(0, 1) = 5;
    add_solution(1, 1) = 6;
    add_solution(2, 1) = 7;

    sub(0, 0) = 0;
    sub(1, 0) = 1;
    sub(2, 0) = 2;
    sub(0, 1) = 1;
    sub(1, 1) = 1;
    sub(2, 1) = 1;

    sub_solution(0, 0) = 1;
    sub_solution(1, 0) = 1;
    sub_solution(2, 0) = 1;
    sub_solution(0, 1) = 3;
    sub_solution(1, 1) = 4;
    sub_solution(2, 1) = 5;
  }

  ColumnMajorMatrix a, t, two_mat, add, add_solution, sub, sub_solution;
};

#endif // COLUMNMAJORMATRIXTEST_H
