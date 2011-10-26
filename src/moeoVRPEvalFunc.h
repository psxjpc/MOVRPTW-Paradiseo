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

#ifndef _moeoVRPEvalFunc_h
#define _moeoVRPEvalFunc_h

// General includes
#include <stdexcept>
#include <fstream>

// The base definition of eoEvalFunc

// Utilities for the VRP-TW problem
#include "moeoVRPUtils.h"

// Objective Vector
#include "moeoVRPObjectiveVectorTraits.h"
typedef moeoRealObjectiveVector<moeoVRPObjectiveVectorTraits> moeoVRPObjectiveVector;

/**
  * \class eoVRPEvalFunc eoVRPEvalFunc.h
  * \brief Evaluates an individual of type eoVRP.
  */

class moeoVRPEvalFunc : public moeoEvalFunc<moeoVRP>
{

public:

    /**
      * \brief Constructor: nothing to do here.
      */
    moeoVRPEvalFunc ()
    { }

    /**
      * \brief Destructor: nothing to do here.
      */
    ~moeoVRPEvalFunc ()
    { }

    /**
      * \brief Invokes all evaluators
      * \param _moeo The individual to be evaluated.
      */
    void operator () (moeoVRP& _moeo)
    {
        if (_moeo.decoded ())
        {
            if (_moeo.invalid ())
                std::runtime_error("Error: invalid individual presents decoding information.");
        }
        else
        {
            if (!_moeo.invalid ())
            {
                std::cerr << "Warning: valid individual does not present decoding information." << std::endl;
                std::cerr << "         Proceeding to decode..." << std::endl;
            }
            _moeo.decode ();
        }

        // Setting new objective values
        moeoVRPObjectiveVector objVector;
        objVector[0] = sizeOfFleet(_moeo);
        objVector[1] = travelDistance(_moeo);
        //objVector[2] = travelTime(_moeo);
        objVector[2] = makespan(_moeo);
        objVector[3] = waitingTime(_moeo);
        objVector[4] = delayTime(_moeo);
        _moeo.objectiveVector(objVector);

        // Setting cache values
        _moeo.sizeOfFleet(objVector[0]);
        _moeo.length(objVector[1]);
        _moeo.time(objVector[2]);
        _moeo.waitingTime(objVector[3]);
        _moeo.delayTime(objVector[4]);
    }

    /**
      * \brief Computes the number of routes of the individual.
      * \param _moeo The individual to be evaluated.
      */
    unsigned sizeOfFleet(const moeoVRP& _moeo) const
    {
        return _moeo.routes().size();
    }

    /**
      * \brief Computes the travel distance of the individual.
      * \param _moeo The individual to be evaluated.
      */
    double travelDistance(const moeoVRP& _moeo) const
    {
        double length = 0.0;
        unsigned i = 0, j = 0;
        for (i = 0; i < _moeo.routes().size(); i++)
        {
            // Departing from the depot
            length += moeoVRPUtils::distance(0, _moeo.routes()[i][0]);

            for (j = 1; j < _moeo.routes()[i].size(); j++)
                length += moeoVRPUtils::distance(_moeo.routes()[i][j - 1], _moeo.routes()[i][j]);

            // Returning to the depot
            length += moeoVRPUtils::distance(_moeo.routes()[i][j - 1], 0);

        }
        return length;
    }

    /**
      * \brief Computes the travel time of the individual.
      * \param _moeo The individual to be evaluated.
      */
    double travelTime(const moeoVRP& _moeo) const
    {
        double elapsedTime = 0.0, totalElapsedTime = 0.0;
        double readyTime = 0.0, dueTime = 0.0, serviceTime = 0.0;
        unsigned previous = 0, current = 0;

        for (size_t i = 0; i < _moeo.routes().size(); i++)
        {
            // Departing from the depot
            previous = 0;

            for (size_t j = 0; j < _moeo.routes()[i].size(); j++)
            {
                current = _moeo.routes()[i][j];
                moeoVRPUtils::getTimeWindow (current, readyTime, dueTime, serviceTime);
                elapsedTime = std::max(elapsedTime +  moeoVRPUtils::elapsedTime(previous, current) , readyTime);
                previous = current;
                elapsedTime += serviceTime;
            }

            // Returning to the depot
            current = 0;
            elapsedTime += moeoVRPUtils::elapsedTime(previous, current);

            totalElapsedTime += elapsedTime;
            elapsedTime = 0.0;
        }
        return totalElapsedTime;
    }

    /**
      * \brief Computes the travel time of the individual.
      * \param _moeo The individual to be evaluated.
      */
    double makespan(const moeoVRP& _moeo) const
    {
        double elapsedTime = 0.0, totalElapsedTime = 0.0;
        double readyTime = 0.0, dueTime = 0.0, serviceTime = 0.0;
        unsigned previous = 0, current = 0;

        for (size_t i = 0; i < _moeo.routes().size(); i++)
        {
            // Departing from the depot
            previous = 0;

            for (size_t j = 0; j < _moeo.routes()[i].size(); j++)
            {
                current = _moeo.routes()[i][j];
                moeoVRPUtils::getTimeWindow (current, readyTime, dueTime, serviceTime);
                elapsedTime = std::max(elapsedTime +  moeoVRPUtils::elapsedTime(previous, current) , readyTime);
                previous = current;
                elapsedTime += serviceTime;
            }

            // Returning to the depot
            current = 0;
            elapsedTime += moeoVRPUtils::elapsedTime(previous, current);

            if (elapsedTime > totalElapsedTime)
               totalElapsedTime = elapsedTime;

            elapsedTime = 0.0;
        }
        return totalElapsedTime;
    }


    /**
      * \brief Computes the total waiting time of the individual.
      * \param _moeo The individual to be evaluated.
      */
    double waitingTime(const moeoVRP& _moeo) const
    {
        double elapsedTime = 0.0, waitingTime = 0.0;
        double readyTime = 0.0, dueTime = 0.0, serviceTime = 0.0;
        unsigned previous = 0, current = 0;

        for (size_t i = 0; i < _moeo.routes().size(); i++)
        {
            // Departing from the depot
            previous = 0;

            for (size_t j = 0; j < _moeo.routes()[i].size(); j++)
            {
                current = _moeo.routes()[i][j];
                moeoVRPUtils::getTimeWindow (current, readyTime, dueTime, serviceTime);
                elapsedTime += moeoVRPUtils::elapsedTime(previous, current);
                // Do we have to wait?
                if (readyTime > elapsedTime)
                {
                    waitingTime += readyTime - elapsedTime;
                    elapsedTime = readyTime;
                }
                previous = current;
                elapsedTime += serviceTime;
            }

            // Returning to the depot - Here we don't have to wait anyway
            elapsedTime = 0.0;
        }
        return waitingTime;
    }

    /**
      * \brief Computes the total delay of the individual.
      * \param _moeo The individual to be evaluated.
      */

    double delayTime(const moeoVRP& _moeo) const
    {

        double elapsedTime = 0.0, delay = 0.0;
        double readyTime = 0.0, dueTime = 0.0, serviceTime = 0.0;
        unsigned previous = 0, current = 0;

        for (size_t i = 0; i < _moeo.routes().size(); i++)
        {
            // Departing from the depot
            previous = 0;

            for (size_t j = 0; j < _moeo.routes()[i].size(); j++)
            {
                current = _moeo.routes()[i][j];
                moeoVRPUtils::getTimeWindow (current, readyTime, dueTime, serviceTime);
                elapsedTime = std::max(elapsedTime +  moeoVRPUtils::elapsedTime(previous, current) , readyTime);
                // Are we arriving late?
                if (elapsedTime > dueTime)
                   delay += elapsedTime - dueTime;
                previous = current;
                elapsedTime += serviceTime;
            }

            // Are we arriving on time at the depot?
            current = 0;
            elapsedTime += moeoVRPUtils::elapsedTime(previous, current);
            moeoVRPUtils::getTimeWindow (current, readyTime, dueTime, serviceTime);
            if (elapsedTime > dueTime)
               delay += elapsedTime - dueTime;

            elapsedTime = 0.0;
        }
        return delay;
    }


protected:


private:


};

#endif
