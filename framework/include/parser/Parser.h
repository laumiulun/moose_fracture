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

#ifndef PARSER_H
#define PARSER_H

// MOOSE includes
#include "ConsoleStreamInterface.h"
#include "MooseTypes.h"
#include "InputParameters.h"
#include "Syntax.h"

#include "hit.h"

#include <vector>
#include <string>

// Forward declarations
class ActionWarehouse;
class SyntaxTree;
class MooseApp;
class Factory;
class ActionFactory;
class GlobalParamsAction;
class JsonSyntaxTree;

inline std::string
errormsg(std::string /*fname*/, hit::Node * /*n*/)
{
  return "";
}

template <typename T, typename... Args>
std::string
errormsg(std::string fname, hit::Node * n, T arg, Args... args)
{
  std::stringstream ss;
  if (n && fname.size() > 0)
    ss << fname << ":" << n->line() << ": ";
  else if (fname.size() > 0)
    ss << fname << ":0: ";
  ss << arg;
  ss << errormsg("", nullptr, args...);
  return ss.str();
}

// Expands ${...} substitution expressions with variable values from the tree.
class ExpandWalker : public hit::Walker
{
public:
  ExpandWalker(std::string fname) : _fname(fname) {}
  virtual void
  walk(const std::string & /*fullpath*/, const std::string & /*nodepath*/, hit::Node * n) override
  {
    auto f = dynamic_cast<hit::Field *>(n);
    auto s = f->val();

    auto start = s.find("${");
    while (start < s.size())
    {
      auto end = s.find("}", start);
      if (end != std::string::npos)
      {
        auto var = s.substr(start + 2, end - (start + 2));
        auto curr = n;
        while ((curr = curr->parent()))
        {
          auto src = curr->find(var);
          if (src && src != n && src->type() == hit::NodeType::Field)
          {
            used.push_back(hit::pathJoin({curr->fullpath(), var}));
            s = s.substr(0, start) + curr->param<std::string>(var) +
                s.substr(end + 1, s.size() - (end + 1));

            if (end + 1 - start == f->val().size())
              f->setVal(s, dynamic_cast<hit::Field *>(curr->find(var))->kind());
            else
              f->setVal(s);

            // move end back to the position of the end of the replacement text - not the replaced
            // text since the former is the one relevant to the string for remaining replacements.
            end = start + curr->param<std::string>(var).size();
            break;
          }
        }

        if (curr == nullptr)
          errors.push_back(
              errormsg(_fname, n, "no variable '", var, "' found for substitution expression"));
      }
      else
        errors.push_back(errormsg(_fname, n, "missing substitution expression terminator '}'"));
      start = s.find("${", end);
    }
  }

  std::vector<std::string> used;
  std::vector<std::string> errors;

private:
  std::string _fname;
};

/**
 * Class for parsing input files. This class utilizes the GetPot library for actually tokenizing and
 * parsing files. It is not currently designed for extensibility. If you wish to build your own
 * parser, please contact the MOOSE team for guidance.
 */
class Parser : public ConsoleStreamInterface, public hit::Walker
{
public:
  enum SyntaxFormatterType
  {
    INPUT_FILE,
    YAML
  };

  Parser(MooseApp & app, ActionWarehouse & action_wh);

  virtual ~Parser();

  /// Retrieve the Syntax associated with the passed Action and task
  std::string getSyntaxByAction(const std::string & action, const std::string & task)
  {
    return _syntax.getSyntaxByAction(action, task);
  }

  /**
   * Return the filename that was parsed
   */
  std::string getFileName(bool stripLeadingPath = true) const;

  /**
   * Parse an input file consisting of hit syntax and setup objects
   * in the MOOSE derived application
   */
  void parse(const std::string & input_filename);

  /**
   * This function attempts to extract values from the input file based on the contents of
   * the passed parameters objects.  It handles a number of various types with dynamic casting
   * including vector types
   */
  void extractParams(const std::string & prefix, InputParameters & p);

  /**
   * Creates a syntax formatter for printing
   */
  void initSyntaxFormatter(SyntaxFormatterType type, bool dump_mode);

  /**
   * Use MOOSE Factories to construct a full parse tree for documentation or echoing input.
   */
  void buildFullTree(const std::string & search_string);

  /**
   * Use MOOSE Factories to construct a parameter tree for documentation or echoing input.
   */
  void buildJsonSyntaxTree(JsonSyntaxTree & tree) const;

  void walk(const std::string & fullpath, const std::string & nodepath, hit::Node * n);

  void errorCheck(const Parallel::Communicator & comm, bool warn_unused, bool err_unused);

  std::vector<std::string> listValidParams(std::string & section_name);

protected:
  /**
   * Helper functions for setting parameters of arbitrary types - bodies are in the .C file
   * since they are called only from this Object
   */
  /// Template method for setting any scalar type parameter read from the input file or command line
  template <typename T, typename Base>
  void setScalarParameter(const std::string & full_name,
                          const std::string & short_name,
                          InputParameters::Parameter<T> * param,
                          bool in_global,
                          GlobalParamsAction * global_block);

  template <typename T, typename UP_T, typename Base>
  void setScalarValueTypeParameter(const std::string & full_name,
                                   const std::string & short_name,
                                   InputParameters::Parameter<T> * param,
                                   bool in_global,
                                   GlobalParamsAction * global_block);

  /// Template method for setting any vector type parameter read from the input file or command line
  template <typename T, typename Base>
  void setVectorParameter(const std::string & full_name,
                          const std::string & short_name,
                          InputParameters::Parameter<std::vector<T>> * param,
                          bool in_global,
                          GlobalParamsAction * global_block);

  /**
   * Sets an input parameter representing a file path using input file data.  The file path is
   * modified to be relative to the directory this application's input file is in.
   */
  template <typename T>
  void setFilePathParam(const std::string & full_name,
                        const std::string & short_name,
                        InputParameters::Parameter<T> * param,
                        InputParameters & params,
                        bool in_global,
                        GlobalParamsAction * global_block);

  /**
   * Sets an input parameter representing a vector of file paths using input file data.  The file
   * paths are modified to be relative to the directory this application's input file is in.
   */
  template <typename T>
  void setVectorFilePathParam(const std::string & full_name,
                              const std::string & short_name,
                              InputParameters::Parameter<std::vector<T>> * param,
                              InputParameters & params,
                              bool in_global,
                              GlobalParamsAction * global_block);
  /**
   * Template method for setting any double indexed type parameter read from the input file or
   * command line.
   */
  template <typename T>
  void setDoubleIndexParameter(const std::string & full_name,
                               const std::string & short_name,
                               InputParameters::Parameter<std::vector<std::vector<T>>> * param,
                               bool in_global,
                               GlobalParamsAction * global_block);

  /**
   * Template method for setting any multivalue "scalar" type parameter read from the input file or
   * command line.  Examples include "Point" and "RealVectorValue".
   */
  template <typename T>
  void setScalarComponentParameter(const std::string & full_name,
                                   const std::string & short_name,
                                   InputParameters::Parameter<T> * param,
                                   bool in_global,
                                   GlobalParamsAction * global_block);

  /**
   * Template method for setting several multivalue "scalar" type parameter read from the input
   * file or command line.  Examples include "Point" and "RealVectorValue".
   */
  template <typename T>
  void setVectorComponentParameter(const std::string & full_name,
                                   const std::string & short_name,
                                   InputParameters::Parameter<std::vector<T>> * param,
                                   bool in_global,
                                   GlobalParamsAction * global_block);

  std::unique_ptr<hit::Node> _cli_root = nullptr;
  std::unique_ptr<hit::Node> _root = nullptr;
  std::vector<std::string> _secs_need_first;

  /// The MooseApp this Parser is part of
  MooseApp & _app;
  /// The Factory associated with that MooseApp
  Factory & _factory;
  /// Action warehouse that will be filled by actions
  ActionWarehouse & _action_wh;
  /// The Factory that builds actions
  ActionFactory & _action_factory;
  /// Reference to an object that defines input file syntax
  Syntax & _syntax;

  /// Object for holding the syntax parse tree
  std::unique_ptr<SyntaxTree> _syntax_formatter;

  /// The input file name that is used for parameter extraction
  std::string _input_filename;

  /// The set of all variables extracted from the input file
  std::set<std::string> _extracted_vars;

  /// Boolean to indicate whether parsing has started (sections have been extracted)
  bool _sections_read;

  /// The current parameter object for which parameters are being extracted
  InputParameters * _current_params;

  /// The current stream object used for capturing errors during extraction
  std::ostringstream * _current_error_stream;

private:
  std::string _errmsg;
  std::string _warnmsg;
  std::string hitCLIFilter(std::string appname, int argc, char * argv[]);
  void walkRaw(std::string fullpath, std::string nodepath, hit::Node * n);
};

#endif // PARSER_H
