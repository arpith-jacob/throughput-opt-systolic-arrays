
//  polyhedron-options.hpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 7 2009
//
//  Read polyhedron options from config file

#ifndef __POLYHEDRON_OPTIONS_H__
#   define __POLYHEDRON_OPTIONS_H__

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/utility.hpp>

using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <cstdlib>
using namespace std;

class PolyhedronOptions
{

 public:

   // constructor
   PolyhedronOptions(string config_file)
   {
     try {
       // parent path of config file
       fs::path config_file_path(config_file, fs::native);

       // group of config options
       po::options_description config("Polyhedron Configuration");
       // variable map
       po::variables_map vm;

       // declare configuration options for polyhedron
       config.add_options()
         ("dimensions", po::value<int>(), "Number of unknowns/dimensions in polyhedron")
         ("parameters", po::value<int>(), "Number of parameters in polyhedron")
         ("parameternames", po::value<string>(), "Parameter names separated by whitespace")
         ("parameterinstantiations", po::value<string>(), "Parameter instantiations separated by whitespace")
         ("pipconstraints", po::value<string>(), "File with pip constraints of polyhedron")
         ("dependencies", po::value<string>(), "File with dependencies of polyhedron")
         ("vertices", po::value<string>(), "File with vertices of polyhedron")
         ;

       // specify configuration file
       ifstream ifs(config_file.c_str());

       // read configuration file
       store(parse_config_file(ifs, config), vm);
       notify(vm);

       // read parameters into variables
       if (vm.count("dimensions")) {
         dimensions = vm["dimensions"].as<int>();
       } else {
         throw "Must specify dimensions of polyhedron";
       }

       if (vm.count("parameters")) {
         parameters = vm["parameters"].as<int>();
       } else {
         throw "Must specify parameters of polyhedron";
       }

       if (vm.count("parameternames")) {
         // split parameter names
         char_separator<char> sep(" \t");
         tokenizer< char_separator<char> > tokens(
                          vm["parameternames"].as<string>(),
                          sep
                        );

         // add it into vector
         for (tokenizer< char_separator<char> >::iterator beg=tokens.begin();
              beg!=tokens.end(); ++beg) {
           parameternames.push_back(*beg);
         }

         // number of parameters must equal number of parameter names
         if (parameters != parameternames.size()) {
           throw "Number of parameters do not match number of parameter names";
         }
       } else {
         throw "Must specify parameter names for polyhedron";
       }

       if (vm.count("parameterinstantiations")) {
         // split parameter instantiations
         char_separator<char> sep(" \t");
         tokenizer< char_separator<char> > tokens(
                          vm["parameterinstantiations"].as<string>(),
                          sep
                        );

         // add it into vector
         for (tokenizer< char_separator<char> >::iterator beg=tokens.begin();
              beg!=tokens.end(); ++beg) {
           parameterinstantiations.push_back( atoi(beg->c_str()) );
         }

         // number of parameters must equal number of parameter instantiations
         if (parameters != parameterinstantiations.size()) {
           throw "Number of parameters do not match number of parameter instantiations";
         }
       } else {
         throw "Must specify parameter instantiations for polyhedron";
       }

       if (vm.count("pipconstraints")) {
         // store pip constraints file.  path is relative to that of the
         // configuration file.
         // note:  branch_path() has changed to parent_path() in the new
         // release of boost.
         pipconstraints = config_file_path.branch_path().string() + "/" +
                          vm["pipconstraints"].as<string>();
       } else {
         throw "Must specify pip constraints of polyhedron";
       }

       if (vm.count("dependencies")) {
         // store dependencies file.  path is relative to that of the
         // configuration file.
         // note:  branch_path() has changed to parent_path() in the new
         // release of boost.
         dependencies = config_file_path.branch_path().string() + "/" +
                        vm["dependencies"].as<string>();
       } else {
         throw "Must specify dependencies of polyhedron";
       }

       if (vm.count("vertices")) {
         // store vertices file.  path is relative to that of the
         // configuration file.
         // note:  branch_path() has changed to parent_path() in the new
         // release of boost.
         vertices = config_file_path.branch_path().string() + "/" +
                    vm["vertices"].as<string>();
       } else {
         throw "Must specify vertices of polyhedron";
       }
     }
     catch(exception &err)
     {
       cerr << "Error parsing options: " << err.what() << endl;
       exit (-1);
     }
     catch(const char *err)
     {
       cerr << "Error parsing options: " << err << endl;
       exit (-1);
     }

     return;
   }

   // destructor
   ~PolyhedronOptions()
   {
   }

   // list of options
   unsigned int dimensions;
   unsigned int parameters;
   vector< string > parameternames;
   vector< int > parameterinstantiations;
   string pipconstraints;
   string dependencies;
   string vertices;

};

#endif // __POLYHEDRON_OPTIONS_H__
