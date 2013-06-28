
//  commandline-options.hpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 7 2009
//
//  Read program options from command line

#ifndef __COMMANDLINE_OPTIONS_H__
#   define __COMMANDLINE_OPTIONS_H__

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <string>
#include <iterator>
using namespace std;

class CommandLineOptions
{

 public:

   // constructor
   CommandLineOptions(int argc, char **argv)
   {
     try {
       // group of config options
       po::options_description commandline("Program Options");
       // variable map
       po::variables_map vm;

       // declare commandline options for program
       commandline.add_options()
         ("help,?", "This help screen")
         ("polyhedron,i", po::value<string>(), "Polyhedron configuration file")
         ("magnitude-bound,m", po::value<int>(), "Upper bound on the magnitude of the projection vector")
         ("pe-inefficiency,n", po::value<int>(), "Upper bound on processor inefficiency: (lambda * u) factor")
         ("pe-pipeline-stages,s", po::value<int>(), "Lower bound on number of processor pipeline stages (Minimum delay on each dependency)")
         ;

       // read command line
       store(po::parse_command_line(argc, argv, commandline), vm);
       notify(vm);

       if (vm.count("help")) {
         cout << commandline << "\n";
         exit(-1);
       }

       // read parameters into variables
       if (vm.count("polyhedron")) {
         polyhedron = vm["polyhedron"].as<string>();
       } else {
         cerr << "Must specify polyhedron configuration file\n";
         throw "Incomplete options";
       }

       // read bound for the magnitude of the projection vector
       if (vm.count("magnitude-bound")) {
         magnitudebound = vm["magnitude-bound"].as<int>();
       } else {
         magnitudebound = 3;
       }

       // read upper bound for the inefficiency of a processor
       if (vm.count("pe-inefficiency")) {
         peinefficiency = vm["pe-inefficiency"].as<int>();
         
         if (peinefficiency < 1 || peinefficiency > 100) {
           throw "Processor inefficiency must be between 1 and 100";
         }
       } else {
         peinefficiency = 100;
       }

       // read lower bound for the number of processor pipeline stages
       if (vm.count("pe-pipeline-stages")) {
         pepipelinestages = vm["pe-pipeline-stages"].as<int>();
         
         if (pepipelinestages < 1 || pepipelinestages > 100) {
           throw "Processor pipeline stages must be between 1 and 100";
         }
       } else {
         pepipelinestages = 1;
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
   ~CommandLineOptions()
   {
   }

   // list of options
   string polyhedron;
   int    magnitudebound;
   int    peinefficiency;
   int    pepipelinestages;

};

#endif // __COMMANDLINE_OPTIONS_H__
