/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef POLYCRYSTALCIRCLES_H
#define POLYCRYSTALCIRCLES_H

#include <array>
#include "PolycrystalUserObjectBase.h"
#include "DelimitedFileReader.h"

// Forward Declarations
class PolycrystalCircles;

template <>
InputParameters validParams<PolycrystalCircles>();

/**
 * PolycrystalCircles creates a polycrystal made up of circles.
 * The locations and radii of the circles are given either
 * through a user input or by reading a .txt file.
 * The file is expected to have a one-line header labeling the
 * colums 'x y z r'.
**/

class PolycrystalCircles : public PolycrystalUserObjectBase
{
public:
  PolycrystalCircles(const InputParameters & parameters);

  // Required functions from PolycrystalUserObjectBase
  virtual void precomputeGrainStructure() override;
  virtual void getGrainsBasedOnPoint(const Point & point,
                                     std::vector<unsigned int> & grains) const override;
  virtual Real getVariableValue(unsigned int op_index, const Point & p) const override;
  virtual unsigned int getNumGrains() const override { return _grain_num; }

protected:
  enum COLS
  {
    X,
    Y,
    Z,
    R
  }; // Names of columns in text file.
  const bool _columnar_3D;
  unsigned int _grain_num; // Number of crystal grains to create

  std::vector<Point> _centerpoints; // x,y,z coordinates of circle centers
  std::vector<Real> _radii;         // Radius for each circular grain created
};

#endif // POLYCRYSTALCIRCLES_H
