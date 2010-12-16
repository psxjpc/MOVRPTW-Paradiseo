/*
 * --------------------------------------------------------------------------
 *
 *                             Copyright (c) 2010
 *                  Juan Castro-Gutierrez <jpcastrog@gmail.com>      (1)
 *             Dario Landa-Silva <dario.landasilva@nottingham.ac.uk> (1)
 *                  José A. Moreno Pérez <jamoreno@ull.es> (2)
 *           --------------------------------------------------------
 *            (1) University of Nottingham (UK) - ASAP Research Group.
 *            (2) Universidad de La Laguna (Spain) - DEIOC.
 *
 * This program is free software (software libre); you can redistribute
 * it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, you can obtain a copy of the GNU
 * General Public License at:
 *                http://www.gnu.org/copyleft/gpl.html
 * or by writing to:
 *           Free Software Foundation, Inc., 59 Temple Place,
 *                 Suite 330, Boston, MA 02111-1307 USA
 *
 * --------------------------------------------------------------------------
 */

#ifndef DOMAKES_H
#define DOMAKES_H

// eo general include
#include <eo>
#include <moeo>

#include <do/make_continue.h>

/** definition of representation:
  * class moeoVRP MUST derive from EO<FitT> for some fitness
  */
#include "moeoVRP.h"

#include "moeoVRPQuadCrossover.h"
#include "moeoVRPMutation.h"


/**
  * \brief Sets the stop criterion.
  * \param _parser Contains the parameters grabbed by the Paradiseo's parser
  * \param _state Stores the object to set the stopping criterion
  */
eoContinue<moeoVRP> & do_make_continue(eoParser& _parser, eoState& _state)
{
    //////////// Stopping criterion ///////////////////
    // the combined continue - to be filled
    eoCombinedContinue<moeoVRP> *continuator = NULL;

    eoValueParam<unsigned> maxGenParam = _parser.getORcreateParam<unsigned>(unsigned(100), "maxGen", "Maximum number of generations");

    if (VERBOSE)
        std::cout << "maxGenParam = " << maxGenParam.value() << std::endl;

    // for each possible criterion, check if wanted, otherwise do nothing
    if (maxGenParam.value()) // positive: -> define and store
    {
        eoGenContinue<moeoVRP> *genCont = new eoGenContinue<moeoVRP>(maxGenParam.value());
        _state.storeFunctor(genCont);
        // and "add" to combined
        continuator = make_combinedContinue<moeoVRP>(continuator, genCont);

        if (!VERBOSE)
           genCont->verbose = false;
    }

    // now check that there is at least one!
    if (!continuator)
        throw std::runtime_error("You MUST provide a stopping criterion");
    // OK, it's there: store in the eoState
    _state.storeFunctor(continuator);

    // and return
    return *continuator;
}


/**
 * This function builds the operators that will be applied to the moeoVRP
 * \param eoParser& _parser to get user parameters
 * \param eoState& _state to store the memory
 */
eoGenOp<moeoVRP> & do_make_op(eoParser& _parser, eoState& _state)
{

  /////////////////////////////
  // Variation operators
  ////////////////////////////

  // the crossover
  ////////////////

  // One point crossover
  eoQuadOp <moeoVRP> *cross = new moeoVRPOnePointCrossover;
  // store in the state
  _state.storeFunctor(cross);
  // relative rate in the combination
  double cross1Rate = _parser.createParam(0.33, "OnePointCrossRate", "Relative rate for the OnePointCrossover").value();
  // creation of the combined operator with this one
  eoPropCombinedQuadOp<moeoVRP> *propXover = new eoPropCombinedQuadOp<moeoVRP>(*cross, cross1Rate);

  // Edge crossover
  cross = new moeoVRPEdgeCrossover;
  // store in the state
  _state.storeFunctor(cross);
  // relative rate in the combination
  double cross2Rate = _parser.createParam(0.33, "EdgeCrossRate", "Relative rate for the EdgeCrossover").value();
  // add new operator
  propXover -> add(*cross, cross2Rate);

  // Generic crossover
  cross = new moeoVRPGenericCrossover;
  // store in the state
  _state.storeFunctor(cross);
  // relative rate in the combination
  double cross3Rate = _parser.createParam(0.33, "GenericCrossRate", "Relative rate for the GenericCrossover").value();
  // add new operator
  propXover -> add(*cross, cross3Rate);

  // Store in the state
  _state.storeFunctor(propXover);

  // the mutation
  ///////////////

  // Swap mutation
  eoMonOp<moeoVRP> *mut = new moeoVRPSwapMutation;
  _state.storeFunctor(mut);
  // its relative rate in the combination
  double mut1Rate = _parser.createParam(0.25, "swapMutationRate", "Relative rate for shift mutation", 0, "Variation Operators").value();
  // creation of the combined operator with this one
  eoPropCombinedMonOp<moeoVRP> *propMutation = new eoPropCombinedMonOp<moeoVRP>(*mut, mut1Rate);
  _state.storeFunctor(propMutation);

  // Insertion mutation
  mut = new moeoVRPInsertionMutation;
  _state.storeFunctor(mut);
  // its relative rate in the combination
  double mut2Rate = _parser.createParam(0.25, "insertionMutationRate", "Relative rate for exchange mutation", 0, "Variation Operators").value();
  // addition of this one to the combined operator
  propMutation -> add(*mut, mut2Rate);

  // Inversion mutation
  mut = new moeoVRPInversionMutation;
  _state.storeFunctor(mut);
  // its relative rate in the combination
  double mut3Rate = _parser.createParam(0.25, "inversionMutationRate", "Relative rate for exchange mutation", 0, "Variation Operators").value();
  // addition of this one to the combined operator
  propMutation -> add(*mut, mut3Rate);

  // Displacement mutation
  mut = new moeoVRPDisplacementMutation;
  _state.storeFunctor(mut);
  // its relative rate in the combination
  double mut4Rate = _parser.createParam(0.25, "displacementMutationRate", "Relative rate for exchange mutation", 0, "Variation Operators").value();
  // addition of this one to the combined operator
  propMutation -> add(*mut, mut4Rate);


  // end of crossover and mutation definitions
  ////////////////////////////////////////////

  // First read the individual level parameters
  eoValueParam<double>& pCrossParam = _parser.getORcreateParam<double>(double(0.25), "pCrossParam", "Probability of Crossover");
  // minimum check
  if ( (pCrossParam.value() < 0) || (pCrossParam.value() > 1) )
    throw std::runtime_error("Invalid pCross");

  eoValueParam<double>& pMutParam = _parser.getORcreateParam<double>(double(0.35), "pMutParam", "Probability of Mutation");
  // minimum check
  if ( (pMutParam.value() < 0) || (pMutParam.value() > 1) )
    throw std::runtime_error("Invalid pMut");

  if (VERBOSE)
  {
      std::cout << "pCrossParam = " << pCrossParam.value() << std::endl;
      std::cout << "pMutParam = " << pMutParam.value() << std::endl;
  }

  eoSequentialOp<moeoVRP> *op = new eoSequentialOp<moeoVRP>;
  _state.storeFunctor(op);

  op -> add(*propXover, pCrossParam.value());
  op -> add(*propMutation, pMutParam.value());

  // return a reference
  return *op;
}

/**
  * \brief Creates the population
  * \param _parser Contains the parameters grab by the Paradiseo's eoParser
  * \param _state Container to keep objects allocated.
  * \param _init Initialisation of the chromosomes.
  */
template <class EOT>
eoPop<EOT>&  do_make_pop(eoParser & _parser, eoState& _state, eoInit<EOT> & _init)
{
    // Seed - In case of non-set value, we set the current time value
    eoValueParam<uint32_t>& seedParam = _parser.getORcreateParam<uint32_t>(uint32_t(time(0)), "seed", "Seed to generate the stochastic behaviour");

    if (VERBOSE)
        std::cout << "seedParam = " << seedParam.value() << std::endl;

    // Size of the population, if non-set then generate 20 individuals
    eoValueParam<unsigned>& popSize = _parser.getORcreateParam<unsigned>(unsigned(20), "popSize", "Size of the population");

    if (VERBOSE)
        std::cout << "popSize = " << popSize.value() << std::endl;

    // create an empty pop and let the state handle the memory
    eoPop<EOT>& pop = _state.takeOwnership(eoPop<EOT>());

    rng.reseed(seedParam.value());

    if (pop.size() < popSize.value()) // missing some guys
    {
        // Init pop from the randomizer: need to use the append function
        pop.append(popSize.value(), _init);
    }

    // for future stateSave, register the algorithm into the state
    _state.registerObject(_parser);
    _state.registerObject(pop);
    _state.registerObject(rng);

    return pop;
}

#endif // DOMAKES_H
