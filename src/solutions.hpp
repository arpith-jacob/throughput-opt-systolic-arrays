
//  solutions.hpp
//
//  Arpith Chacko Jacob
//  jarpith@cse.wustl.edu
//  May 11 2009
//
//  All projection solutions.

#ifndef __SOLUTIONS_H__
#   define __SOLUTIONS_H__

#include <list>
using namespace std;

#include "projection-solution.hpp"

class Solutions : public list<ProjectionSolution *>
{

 public:

   // constructor
   Solutions()
   {
   }

   // destructor
   ~Solutions()
   {
   }

   // static member function, sort helper
   // compare two projection solutions; sort by throughput, utilization and
   // latency
   static bool compare_proj_solns (
                                    ProjectionSolution *second,
                                    ProjectionSolution *first
                                  )
   {
     if (first->instance_bpp < second->instance_bpp)
       return true;
     else if (first->instance_bpp == second->instance_bpp)
       if (first->instance_pe_count > second->instance_pe_count)
         return true;
       else if (first->instance_pe_count == second->instance_pe_count)
         if (first->utilization > second->utilization)
           return true;
         else if (first->utilization == second->utilization)
           if (first->latency > second->latency)
             return true;
           else if (first->latency == second->latency)
             if (first->network_max_length > second->network_max_length)
               return true;
             else if (first->network_max_length == second->network_max_length)
               if (first->network_avg_length > second->network_avg_length)
                 return true;
         
     return false;
   }

   // sort projection solutions
   void Sort()
   {
     this->sort( Solutions::compare_proj_solns );
   }

   // print projection solutions with a unimodular change of basis
   void printSolutions(int peinefficiency)
   {

     list<ProjectionSolution *>::iterator i;

     for (i = begin(); i != end(); i++) {
       if ((*i)->utilization <= peinefficiency) {
         (*i)->print();
         cout << endl;
       }
     }

   }

#if 0
   // print all projection solutions, including those with non unimodular
   // change of basis matrices
   void printSolutions()
   {

     list<ProjectionSolution *>::iterator i;

     // print only projection solutions that have a change of basis that
     // is unimodular
     unsigned int prev_instance_bpp = 0;
     for (i = begin(); i != end(); i++) {
       //if (prev_instance_bpp != (*i)->instance_bpp) {
         (*i)->print();
         cout << endl;

       // prev_instance_bpp =  (*i)->instance_bpp;
       //}
     }

   }
#endif

 private:

};

#endif // __SOLUTIONS_H__
