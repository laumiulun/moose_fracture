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

#ifndef MOOSEENUMBASE_H
#define MOOSEENUMBASE_H

// C++ includes
#include <string>
#include <set>
#include <vector>
#include <map>

// MOOSE includes
#include "MooseEnumItem.h"

/**
 * The base class for both the MooseEnum and MultiMooseEnum classes.
 */
class MooseEnumBase
{
public:
  /**
   * Constructor that takes a list of enumeration values, and a
   * separate string to set a default for this instance.
   * @param names - a list of names for this enumeration
   * @param allow_out_of_range - determines whether this enumeration will accept values outside of
   *                             its range of defined values.
   */
  MooseEnumBase(std::string names, bool allow_out_of_range = false);

  /**
   * Copy Constructor for use when creating vectors of MooseEnumBases
   * @param other_enum - The other enumeration to copy state from
   */
  MooseEnumBase(const MooseEnumBase & other_enum);

  /**
   * This class must have a virtual destructor since it has derived classes.
   */
  virtual ~MooseEnumBase() = default;

  /**
   * Deprecates various options in the MOOSE enum. For each deprecated option,
   * you may supply an optional new option that will be used in a message telling
   * the user which new option replaces the old one.
   */
  virtual void deprecate(const std::string & name, const std::string & raw_name = "");

  /**
   * Method for returning a vector of all valid enumeration names for this instance
   * @return a vector of names
   */
  std::vector<std::string> getNames() const;

  /**
   * Method for returning the raw name strings for this instance
   * @return a space separated list of names
   */
  std::string getRawNames() const;

  /**
   * IsValid
   * @return - a Boolean indicating whether this Enumeration has been set
   */
  virtual bool isValid() const = 0;

  /**
   * isOutOfRangeAllowed
   * @return - a Boolean indicating whether enum names out of range are allowed
   */
  bool isOutOfRangeAllowed() const { return _allow_out_of_range; }

  /**
   * Return the complete set of available flags.
   */
  const std::set<MooseEnumItem> & items() { return _items; }

  ///@{
  /**
   * Locate an item.
   */
  std::set<MooseEnumItem>::const_iterator find(const MooseEnumItem & other) const;
  std::set<MooseEnumItem>::const_iterator find(const std::string & name) const;
  std::set<MooseEnumItem>::const_iterator find(int id) const;
  ///@}

  /**
   * Compute the next valid ID.
   */
  int getNextValidID() const;

protected:
  MooseEnumBase();

  ///@{
  /**
   * Methods to add possible enumeration value to the enum.
   *
   * The MooseEnum/MultiMooseEnum are not designed to be modified, with respect to the list
   * of possible values. However, this is not the case for the ExecFlagEnum which is a special
   * type of MultiMooseEnum for managing the "execute_on" flags. These methods are used by
   * ExecFlagEnum to allow users to modify the available execute flags for their object.
   */
  void addEnumerationNames(const std::string & names);
  void addEnumerationName(const std::string & raw_name);
  void addEnumerationName(const std::string & name, const int & value);
  void addEnumerationItem(const MooseEnumItem & item);
  ///@}

  /**
   * Method that must be implemented to check derived class values against the _deprecated_names
   */
  virtual void checkDeprecated() const = 0;

  /**
   * Check and warn deprecated values
   */
  void checkDeprecated(const MooseEnumItem & item) const;

  /// Storage for the assigned items
  std::set<MooseEnumItem> _items;

  /// The map of deprecated names and optional replacements
  std::map<MooseEnumItem, MooseEnumItem> _deprecated_items;

  /// Flag to enable enumeration items not previously defined
  bool _allow_out_of_range;
};

#endif // MOOSEENUMBASE_H
