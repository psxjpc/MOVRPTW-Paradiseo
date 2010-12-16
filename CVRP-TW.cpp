/*
 * Copyright (C) DOLPHIN Project-Team, INRIA Futurs, 2006-2007
 * (C) OPAC Team, LIFL, 2002-2007
 *
 * (c) Antonio LaTorre <atorre@fi.upm.es>, 2007
 *
 * ===================================================================
 * Multi-Objective Support added by:
 * ===================================================================
 * Uninversity of Nottingham (UK), ASAP Research Group
 * (c) Juan Castro-Gutierrez <jpcastrog@gmail.com>, 2010
 * (c) Dario Landa-Silva <dario.landasilva@nottingham.ac.uk>, 2010
 *
 * Universidad de La Laguna (Spain), DEIOC
 * (c) José A. Moreno Pérez <jamoreno@ull.es>, 2010
 *
 * -------------------------------------------------------------------
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 * ParadisEO WebSite : http://paradiseo.gforge.inria.fr
 * Contact: paradiseo-help@lists.gforge.inria.fr
 *
 */

// Defines /////////////

// Output debug data
#define VERBOSE false

////////////////////////

// Miscellaneous includes and declarations
#include <iostream>
#include <time.h>

// eo general include
#include <eo>
#include <moeo>
// the real bounds (not yet in general eo include)
#include "utils/eoRealVectorBounds.h"

// Include here whatever specific files for your representation.
// Basically, this should include at least the following:

/** definition of representation:
  * class moeoVRP MUST derive from EO<FitT> for some fitness
  */
#include "moeoVRP.h"

/** definition of initilizqtion:
  * class moeoVRPInit MUST derive from eoInit<moeoVRP>
  */
#include "moeoVRPInit.h"

/** definition of evaluation:
  * class moeoVRPEvalFunc MUST derive from eoEvalFunc<moeoVRP>
  * and should test for validity before doing any computation
  * see tutorial/Templates/evalFunc.tmpl
  */
#include "moeoVRPEvalFunc.h"

/** definitions of operators: write as many classes as types of operators
  * and include them here. In this simple example,
  * one crossover (2->2) and one mutation (1->1) operators are used
  */
#include "moeoVRPQuadCrossover.h"
#include "moeoVRPMutation.h"

/* And (possibly) your personal statistics */
#include "moeoVRPStat.h"

#include "moeoVRPUtils.h"


/* **********************************************************************************
   ********************************************************************************** */

// Use existing modules to define representation independent routines

// Output (stats, population dumps, ...)
#include "do/make_checkpoint.h"

// Simply call to the algo. Stays there for consistency reasons
// No template for that one
#include "do/make_run.h"

// The instanciating fitnesses
#include <eoScalarFitness.h>

#include "doMakes.h"

// Checks for help demand, and writes the status file
// and make_help; in libutils

void make_help (eoParser& _parser);

/* **********************************************************************************
   ********************************************************************************** */

// Now use all of the above, + representation dependent things
int main (int argc, char* argv [])
{
    // Checking parameters
    if (argc != 8)
    {
        std::cout << "Error!: Insuffient Parameters" << std::endl;
        std::cout << "You must provide:" << std::endl;
        std::cout << " - algorithm: moeoNSGAII, moeoMOGA" << std::endl;
        std::cout << " - file containig the (solomon-like) data of the problem" << std::endl;
        std::cout << " - file containig the distance matrices of the problem" << std::endl;
        std::cout << " - file containig the time matrices of the problem" << std::endl;
        std::cout << " - size of population" << std::endl;
        std::cout << " - number of generations" << std::endl;
        std::cout << " - seed" << std::endl;
        std::cout << "Example "<< argv[0] << " moeoNSGAII solomonData.dat distanceMatrixData.dat timeMatrixData.dat 20 1000 0" << std::endl;
        exit(1);
    }


    try
    {
        // ////////////////////// //
        // User parameter reading //
        // ////////////////////// //

        eoParser parser (argc, argv);

        // Parameter for loading a problem instance
        eoValueParam<std::string> algorithmName (std::string(argv[1]), "algorithmName", "MO Algorithm to tackle this instance");
        parser.processParam (algorithmName, "Problem params");
        std::string algorithm = algorithmName.value ();


        // Parameter for loading a problem instance
        eoValueParam<std::string> instanceParam (std::string(argv[2]), "instance", "Instance to be loaded");
        parser.processParam (instanceParam, "Problem params");
        std::string instance = instanceParam.value ();

        // Parameter for loading the distance matrix
        eoValueParam<std::string> instanceDistanceMatrixParam(std::string(argv[3]), "distanceMatrixFile", "File that contains the distance matrix");
        parser.processParam(instanceDistanceMatrixParam, "Problem params");
        std::string distanceMarixFile = instanceDistanceMatrixParam.value();

        // Parameter for loading the time matrix
        eoValueParam<std::string> instanceTimeMatrixParam(std::string(argv[4]), "timeMatrixFile", "File that contains the time matrix");
        parser.processParam(instanceTimeMatrixParam, "Problem params");
        std::string timeMatrixFile = instanceTimeMatrixParam.value();

        // Parameter that establishes the max size of the fleet
        // Not used in this implementation. The size of the fleet is not fixed
        // eoValueParam<unsigned> fleet(atoi(argv[4]), "sizeOfFleet", "Number of vehicle to be used when creating the route-plan");
        // parser.processParam(fleet, "Problem params");
        // unsigned sizeOfFleet = fleet.value();

        // Parameter to set the size of the population
        eoValueParam<unsigned> popSize(unsigned(atoi(argv[5])), "popSize", "Size of the population");
        parser.processParam(popSize, "Algo params");

        // Parameter to set the max number of generations
        eoValueParam<unsigned> maxGenParam(unsigned(atoi(argv[6])), "maxGen", "Maximum number of generations");
        parser.processParam(maxGenParam, "Algo params");

        // Parameter to set the seed to generate random numbers
        eoValueParam<uint32_t> seed(uint32_t(atoi(argv[7])), "seed", "Seed to generate the stochastic behabiour");
        parser.processParam(seed, "Algo params");

        // //////////////////////////////////// //
        // Crossover and Mutation probabilities //
        // //////////////////////////////////// //

        // First read the individual level parameters
        eoValueParam<double> pCrossParam(double(0.8), "pCrossParam", "Probability of Crossover");
        parser.processParam(pCrossParam);
        // minimum check
        if ( (pCrossParam.value() < 0) || (pCrossParam.value() > 1) )
          throw std::runtime_error("Invalid pCross");

        eoValueParam<double> pMutParam(double(0.2), "pMutParam", "Probability of Mutation");
        parser.processParam(pMutParam);
        // minimum check
        if ( (pMutParam.value() < 0) || (pMutParam.value() > 1) )
          throw std::runtime_error("Invalid pMut");


        // Initialization of a given seed
        rng.reseed (seed.value());

        // Load an instance of the problem
        moeoVRPUtils::load(instance.c_str ());
        moeoVRPUtils::loadDistanceMatrix(distanceMarixFile.c_str());
        moeoVRPUtils::loadTimeMatrix(timeMatrixFile.c_str());

        // ////////////////////////// //
        // Keeps all things allocated //
        // ////////////////////////// //

        eoState state;

        // ///////////////////// //
        // The fitness evaluator //
        // ///////////////////// //

        moeoVRPEvalFunc plainEval;

        // Turn that object into an evaluation counter
        eoEvalFuncCounter<moeoVRP> eval (plainEval);

        // ////////////////////// //
        // A genotype initializer //
        // ////////////////////// //

        moeoVRPInit init;

        // ////////////////////////////////////////////// //
        // Now some representation-independent things     //
        // (no need to modify anything beyond this point) //
        // ////////////////////////////////////////////// //

        // Initialize the population
        eoPop<moeoVRP>& pop = do_make_pop (parser, state, init);

        // Variation operators
        eoGenOp<moeoVRP>& op = do_make_op (parser, state);

        // Stopping criteria
        eoContinue<moeoVRP>& genCont = do_make_continue (parser, state);

        // Checkpoint
        eoCheckPoint<moeoVRP>& checkpoint = do_make_checkpoint (parser, state, eval, genCont);


        // ////////// // //////////
        // Statistics - Not used //
        // ////////// // //////////

        moeoVRPStat myStat;
        checkpoint.add (myStat);

        // This one is probably redundant with the one in make_checkpoint, but w.t.h.
        eoIncrementorParam<unsigned> generationCounter ("Gen.");
        checkpoint.add (generationCounter);

        // Need to get the name of the redDir param (if any)
        std::string dirName = parser.getORcreateParam (std::string ("./Res"), "resDir", "Directory to store DISK outputs", '\0', "Output - Disk").value () + "/";

        // Those need to be pointers because of the if's
        eoStdoutMonitor* myStdOutMonitor;
        eoFileMonitor*   myFileMonitor;

#ifdef HAVE_GNUPLOT
        eoGnuplot1DMonitor* myGnuMonitor;
#endif

        // Now check how you want to output the stat:
        bool printVRPStat = parser.createParam (false, "coutVRPStat", "Prints my stat to screen, one line per generation", '\0', "My application").value ();
        bool fileVRPStat = parser.createParam (false, "fileVRPStat", "Saves my stat to file (in resDir", '\0', "My application").value ();
        bool plotVRPStat = parser.createParam (false, "plotVRPStat", "On-line plots my stat using gnuplot", '\0', "My application").value ();

        // Should we write it on StdOut ?
        if (printVRPStat)
        {

            myStdOutMonitor = new eoStdoutMonitor (false);

            // Don't forget to store the memory in the state
            state.storeFunctor (myStdOutMonitor);

            // And of course to add the monitor to the checkpoint
            checkpoint.add (*myStdOutMonitor);

            // And the different fields to the monitor
            myStdOutMonitor->add (generationCounter);
            myStdOutMonitor->add (eval);
            myStdOutMonitor->add (myStat);

        }

        // First check the directory (and creates it if not exists already):
        if (fileVRPStat || plotVRPStat)
            if (!testDirRes (dirName, true))
                throw std::runtime_error ("Problem with resDir");

        // Should we write it to a file ?
        if (fileVRPStat)
        {

            // The file name is hard-coded - of course you can read
            // a string parameter in the parser if you prefer
            myFileMonitor = new eoFileMonitor (dirName + "myStat.xg");

            // Don't forget to store the memory in the state
            state.storeFunctor (myFileMonitor);

            // And of course to add the monitor to the checkpoint
            checkpoint.add (*myFileMonitor);

            // And the different fields to the monitor
            myFileMonitor->add (generationCounter);
            myFileMonitor->add (eval);
            myFileMonitor->add (myStat);

        }

#ifdef HAVE_GNUPLOT

        // Should we PLOT it on StdOut ? (one dot per generation, incremental plot)
        if (plotVRPStat) {

            myGnuMonitor = new eoGnuplot1DMonitor (dirName + "plot_myStat.xg", minimizing_fitness<RouteSet> ());
            // NOTE: you can send commands to gnuplot at any time with the method
            // myGnuMonitor->gnuplotCommand(string)
            // par exemple, gnuplotCommand("set logscale y")

            // Don't forget to store the memory in the state
            state.storeFunctor (myGnuMonitor);

            // And of course to add the monitor to the checkpoint
            checkpoint.add (*myGnuMonitor);

            // And the different fields to the monitor (X = eval, Y = myStat)
            myGnuMonitor->add (eval);
            myGnuMonitor->add (myStat);

        }

#endif


       // ///////////////////////////// //
        // Construction of the algorithm //
        // ///////////////////////////// //
        // initialization of the population

        /*** the representation-independent things ***/

        // definition of the archive
        moeoUnboundedArchive <moeoVRP> archive;
        moeoArchiveUpdater <moeoVRP> updater(archive, pop);
        checkpoint.add(updater);

        // fitness assignment
        moeoDominanceDepthFitnessAssignment <moeoVRP> fitnessAssignment;

        // diversity preservation
        moeoFrontByFrontCrowdingDiversityAssignment <moeoVRP> diversityAssignment;

        // comparator
        moeoFitnessThenDiversityComparator <moeoVRP> comparator;

        // selection scheme
        moeoDetTournamentSelect <moeoVRP> select(comparator, 2);

        // replacement scheme
        moeoElitistReplacement <moeoVRP> replace(fitnessAssignment, diversityAssignment, comparator);

        // breeder
        eoGeneralBreeder <moeoVRP> breed(select, op);

        // algorithms
        moeoEA < moeoVRP >* algo;
        if (algorithm == "moeoNSGAII")
           algo = new moeoNSGAII <moeoVRP> (checkpoint, eval, op);
        else if (algorithm == "moeoMOGA")
           algo = new moeoMOGA <moeoVRP> (checkpoint, eval, op);
        else
           throw std::runtime_error("[ERROR]: Algorithm not defined or not available. Options are: moeoEasyEA, moeoNSGAII, moeoMOGA");


        /*** Go ! ***/

        // help ?
        make_help(parser);

        // first evalution (for printing)
        apply<moeoVRP>(eval, pop);

        // printing of the initial population
        //std::cout << "Initial Population\n";
        //pop.sortedPrintOn(std::cout);
        //std::cout << std::endl;

        // run the algo
        (*algo)(pop);

        // printing of the final population
        //std::cout << "Final Population\n";
        //pop.sortedPrintOn(std::cout);
        //std::cout << std::endl;

        // printing of the final archive
        std::cout << "Final Archive\n";
        archive.sortedPrintOn(std::cout);

    }
    catch (std::exception& e)
    {
        std::cerr << e.what () << std::endl;
    }

    return 0;


}
