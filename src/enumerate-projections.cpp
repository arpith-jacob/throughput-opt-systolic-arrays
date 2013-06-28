
//  enumerate-projects.cpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 7 2009
//
//  Enumerate list of projection vectors for a polyhedron, with a goal of
//  finding a high throughput vector.

#include <cstdio>

#include <barvinok/barvinok.h>
#include <barvinok/evalue.h>
#include <barvinok/util.h>

// pip includes
//#include <piplib/piplib64.h>

// local includes
#include "polyhedron-options.hpp"
#include "commandline-options.hpp"
#include "projection-solver.hpp"
#include "index-enumerator.hpp"
#include "solutions.hpp"

int main(int argc, char **argv)
{

  //
  // parse command line options
  //
  CommandLineOptions clopt(argc, argv);

  //
  // parse polyhedron configuration file
  //
  cout << "Parsing polyhedron configuration file: " << clopt.polyhedron << endl;
  PolyhedronOptions polyopt(clopt.polyhedron);

  //
  // read polyhedron constraints in pip format from file  
  //
  cout << "Reading pip polyhedron: " << polyopt.pipconstraints << endl;
  ProjectionSolver solver(
              polyopt.dimensions,
              polyopt.parameters,
              &polyopt.parameterinstantiations,
              &polyopt.parameternames,
              clopt.pepipelinestages,
              polyopt.pipconstraints,
              polyopt.dependencies,
              polyopt.vertices
             );

  cout << "Magnitude bound for the projection vector: " << clopt.magnitudebound << endl;
  cout << "Processor inefficiency (lambda * u): " << clopt.peinefficiency << endl;
  cout << "Minimum processor pipeline stages (lambda * d): " << clopt.pepipelinestages << endl;

  //
  // Projection vector solutions
  //
  Solutions projsols;

  // this is the projection vector index
  IndexEnumerator pv(polyopt.dimensions, clopt.magnitudebound);

  int candidates = 0;
  while (!pv.end()) {

    //
    //  Check GCD(projection vector) == 1
    //
    //
    // ignore if gcd != 1 or zero vector
    if (pv.gcd() == 1 && !pv.isOverBound()) {

      //
      // Call ILP solver using throughput ILP for this projection vector
      //  
      ProjectionSolution *ps = solver.findThroughput(pv.index);

      //
      // Call ILP solver using schedule ILP for this projection vector
      //
      solver.findSchedule(ps);

      // compute delays induced by schedule
      solver.computeScheduleNetwork(ps);

      // compute allocation matrix
      solver.computeAllocation(ps);

      // compute size of interconnection network links
      solver.computeInterconnectionNetwork(ps);

      // compute number of PEs in this projection
      solver.countPEs(ps);

      // compute throughput for an instance of the problem
      // parameter instances are given in the options file
      ps->computeInstanceBPP();

      // store this solution
      projsols.push_front(ps);

//      ps->print();
//      cout << endl << endl;

      // count number of candidate projection vectors explored
      candidates++;
    }

    //
    // increment projection vector
    //
    pv.incr();
  }
  
  cout << candidates << " projection vectors explored\n";


  //
  // Sort projection vectors by throughput (for an instance of the parameters),
  // utilization, max network length, sum of network lengths and latency
  //
  projsols.Sort();

  cout << "\n\nPrinting solutions\n";
  projsols.printSolutions(clopt.peinefficiency);

  return 0;
}

