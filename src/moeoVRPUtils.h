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

#ifndef _moeoVRP_h
#define _moeoVRP_h

// For swappings elements within vectors
#include <algorithm>
#include <vector>


// The base definition of eoVector
#include <eoVector.h>

// Utilities for the VRP-TW problem
#include "moeoVRPUtils.h"

#include <algorithm>
#include <limits>

// Infinite <double>
// #define INFd (std::numeric_limits<double>::max)()
#define INFd std::numeric_limits<double>::infinity()
// Infinite <integer>
// #define INFi (std::numeric_limits<unsigned long>::max)()
#define INFi std::numeric_limits<long>::infinity()

#include <string>

// Objective Vector
#include "moeoVRPObjectiveVectorTraits.h"
typedef moeoRealObjectiveVector<moeoVRPObjectiveVectorTraits> moeoVRPObjectiveVector;


// TEMP, delay max 30 min = 1800 segs
#define DELAY_MAX 1800


/**
  * \class moeoVRP moeoVRP.h
  * \brief Defines the getoype used to solve the VRP-TW problem.
  * Objectives:
  * - <Size of the fleet>::int is number of vehicles we need to use to satisfy all costumers.
  * – <Travel Distance>::double is length of the route-plan.
  * – <Travel Time>::double is elapsed time since the first delivery vehicle departs from the depot until the last arrived at the depot.
  * – <Waiting Time>::double is the amount of time that vehicles have to wait at each costumer location.
  * – <Delay Time>::double is the amount of time by which the arrival of the vehicles + service time is retarded respect to the closing time of the costumers.
  */

class moeoVRP: public moeoVector<moeoVRPObjectiveVector, double, double, int>
{

public:

    /**
      * \brief Default constructor: initializes variables to safe values.
      */

    moeoVRP () : mLength (0.0), mTime (0.0), mDelayTime (0.0), mWaitingTime (0.0), mSizeOfFleet (0)
    {

    }


    /**
      * \brief Copy contructor: creates a new individual from a given one.
      * \param _orig The individual used to create the new one.
      */

    moeoVRP (const moeoVRP& _orig)
    {
        operator= (_orig);
    }


    /**
      * \brief Default destructor: nothing to do here.
      */

    ~moeoVRP ()
    {

    }


    /**
      * \brief Performs a copy from the invidual passed as argument.
      * \param _orig The individual to copy from.
      * \return A reference to this.
      */

    moeoVRP& operator= (const moeoVRP& _orig)
    {

        // Sanity check
        if (&_orig != this)
        {

            // Cleans both individual and decoding information
            clean ();

            // We call the assignment operator from the base class
            moeoVector<moeoVRPObjectiveVector, double, double, int>::operator= (_orig);

            // And then copy all our attributes
            mRoutes      = _orig.mRoutes;
            mLength      = _orig.mLength;
            mTime        = _orig.mTime;
            mWaitingTime = _orig.mWaitingTime;
            mDelayTime   = _orig.mDelayTime;
            mSizeOfFleet = _orig.mSizeOfFleet;
        }

        return *this;

    }


    /**
      * \brief Returns a string containing the name of the class.
      * \return The string containing the name of the class.
      */

    virtual std::string className () const
    {
        return "moeoVRP";
    }


    /**
      * \brief Prints the individual to a given stream.
      * \param _os The stream to print to.
      */

    void printOn (std::ostream& _os) const
    {
        // Then the individual itself using the base printing method
        //moeoVector<moeoVRPObjectiveVector, double, double, int>::printOn (_os);
        //_os << std::endl;
        // We will use a custom method to know where each route starts and ends

        // Objective values
        _os << this->objectiveVector() << "\t";

        printRoutes(_os);
        // Using 0 as separator for each route within the route-plan

    }


    /**
      * \brief Prints a detailed version of the individual (decoding information, unsatisfied contraints, etc.) to a given stream.
      * \param _os The stream to print to.
      */

    void printAllOn (std::ostream& _os) const
    {
        // Print the individual itself using the base printing method
        moeoVector<moeoVRPObjectiveVector, double, double, int>::printOn (_os);
        _os << std::endl << std::endl;

        // Check if we have decoding information to print
        if (decoded ())
        {
            // First, we print the decoded routes (stored in mRoutes)
            _os << " => Routes: " << std::endl << std::endl;
            printRoutes (_os);
            _os << std::endl << std::endl;

            if (this->invalid ())
                _os << " => Fitness: INVALID." << std::endl << std::endl;
            else
                _os << " => Fitness: " << this->fitness () << std::endl << std::endl;
        }
        else
            std::cerr << "Warning: 'printAllOn' called but the individual was not already decoded." << std::endl;

    }

    /**
      * \brief Write a route-plan in a file.
      * \param _filename is the name of the file the route-plan will be written in
      */

    void writeRoutePlan (std::string _filename) const
    {
        std::ofstream _file (_filename.c_str());
        if (decoded ())
        {
            _file << "0";
            for (size_t i = 0; i < mRoutes.size(); i++)
            {
               for (unsigned j = 0; j < mRoutes [i].size (); j++)
                  _file << " " << mRoutes [i][j];

               _file << " 0";
            }
        }
        else
            std::cerr << "Warning: 'writeRoutePlan' called but the individual was not already decoded." << std::endl;

    }

    /**
      * \brief Reads an individual from a given stream.
      * \param _is The stream to read from.
      */

    void readFrom (std::istream& _is)
    {
        // Read the individual using the method from the base class
        moeoVector<moeoVRPObjectiveVector, double, double, int>::readFrom (_is);
    }


    /**
      * \brief Returns a reference to the decoded individual.
      * \return A reference to the decoded individual.
      */

    const Routes& routes () const
    {

        if (mRoutes.size () == 0)
            std::cerr << "Warning: This individual has not been already decoded." << std::endl;

        return mRoutes;

    }


    /**
      * \brief Aux. method to print a structure of routes.
      * \param _os The stream to print to.
      */

    void printRoutes (std::ostream& _os) const
    {
        _os << "[ " << mRoutes.size() << " routes ] ";

        for (unsigned i = 0; i < mRoutes.size (); i++)
           printRoute (_os, i);

        _os << "0 ";
    }


    /**
      * \brief Aux. method to print only one route.
      * \param _os The stream to print to.
      * \param _p The route to print.
      */

    void printRoute (std::ostream& _os, unsigned _p) const
    {
        _os << "0 ";

        for (unsigned i = 0; i < mRoutes [_p].size (); i++)
            _os << mRoutes [_p][i] << " ";

    }


    /**
      * \brief Cleans the individual (the vector of clients and also the decoding information).
      * \return True if the operation finishes correctly. False otherwise.
      */

    bool clean ()
    {
        this->clear ();
        mRoutes.clear ();
        mSizeOfFleet = 0;
        mLength = 0.0;
        mTime = 0.0;
        mWaitingTime = 0.0;
        mDelayTime = 0.0;
        return true;
    }


    /**
      * \brief Invalidates the decoding information (usually after crossover or mutation).
      * \return True if the operation finishes correctly. False otherwise.
      */

    bool cleanRoutes ()
    {
        mRoutes.clear ();
        mSizeOfFleet = 0;
        mLength = 0.0;
        mTime = 0.0;
        mWaitingTime = 0.0;
        mDelayTime = 0.0;
        return true;
    }


    /**
      * \brief Has this individual been decoded?
      * \return True if has decoding information. False otherwise.
      */

    bool decoded () const
    {
        if (mRoutes.size () == 0)
            return false;

        return true;
    }


    /**
      * \brief Encodes an individual from a set of routes (usually used within crossover). The chromosome
      will have the following structure:

      [ Route_1 ] 0 [ Route_2 ] 0 [ Route_3 ] 0 ... 0 [ Route_n ]

      * \return True if the operation finishes correctly. False otherwise.
      */

    bool encode (Routes& _routes)
    {
        clean ();

        for (unsigned i = 0; i < _routes.size (); i++)
            for (unsigned j = 0; j < _routes [i].size (); j++)
                this->push_back (_routes [i][j]);

        return true;
    }



    /**
      \brief Decodes an individual for its evaluation.
    * \return True if the operation finishes correctly. False otherwise.
    */
    bool decode ()
    {
        double demand = 0.0, time = 0.0;
        double readyTime = 0.0, dueTime = 0.0, serviceTime = 0.0;

        // Data members re-initilization
        cleanRoutes ();

        // Temp route
        Route route;

        //   We force feasibility (capacity) by splitting overloaded routes
        //   and splitting those in which the time windows is exceeded by a maximum delay
        unsigned currentCustomer, previousCustomer = 0;

        //  Depot closing time
		double depotClosingTime = moeoVRPUtils::clients [0].dueTime, returnTripTime;
		
        for (size_t i = 0; i < this->size(); i++)
        {
            currentCustomer = this->at(i);
            demand += moeoVRPUtils::clients [currentCustomer].demand;

            moeoVRPUtils::getTimeWindow (currentCustomer, readyTime, dueTime, serviceTime);
            time += moeoVRPUtils::elapsedTime (previousCustomer, currentCustomer);
            time = std::max(readyTime, time);

            // To calculate if it is feasible to go back to the depot arriving within its time window after visiting currentCustomer
			returnTripTime = time + moeoVRPUtils::elapsedTime (currentCustomer, 0) + serviceTime;
            if ((demand > vehicleCapacity) || (time > (dueTime + DELAY_MAX)) || returnTripTime  > depotClosingTime)
            {
                this->mRoutes.push_back(route);
                route.clear();
				i--;
				previousCustomer = 0;
                demand = 0.0;
                time = 0.0;
            }
            else
            {
                time += serviceTime;
                route.push_back(currentCustomer);
				previousCustomer = currentCustomer;
            }
        }

        // Last route
        if (route.size() != 0)
           this->mRoutes.push_back(route);

        #ifdef DEBUG
            moeoVRPUtils::safetyCheck(mRoutes, "moeoVRP::decoding");
        #endif

        return true;
    }

    /**
      * \brief Returns the number of vehicles need to satisfy all costumers in this route.
      * \return The total number of vehicles used in this route-plan.
      */
    unsigned sizeOfFleet()
    {
        return mSizeOfFleet;
    }

    /**
      * \brief Sets the size of the fleet
      * \param mSizeOfFleet is the size of the fleet to be set up.
      */
    void sizeOfFleet(unsigned mSizeOfFleet)
    {
        this->mSizeOfFleet = mSizeOfFleet;
    }

    /**
    * \brief Returns the total cost (length) of traveling all the routes.
    * \return The total cost (length) of traveling all the routes.
    */
    double length ()
    {
        return mLength;
    }

    /**
      * \brief Sets the lenght of the current route-plan
      * \param mLength is the value of the length to be set up.
      */
    void length (double mLength)
    {
        this->mLength = mLength;
    }

    /**
    * \brief Returns the total cost (travel time) of traveling all routes.
    * \return The total cost (travel time) of traveling all routes.
    */
    double time()
    {
        return mTime;
    }

    /**
      * \brief Sets the elapsed time of all routes.
      * \param mTime is the elapsed time to be set up.
      */
    void time(double mTime)
    {
        this->mTime = mTime;
    }

    /**
    * \brief Returns the total amount of time vehicles have to wait for costumer to open their time windows.
    * \return The total amount of waiting time at all costumers.
    */
    double waitingTime()
    {
        return mWaitingTime;
    }

    /**
    * \brief Sets the waiting time at all costumers
    * \param mWaitingTime is the total ammount of waiting time
    */
    void waitingTime(double mWaitingTime)
    {
        this->mWaitingTime = mWaitingTime;
    }


    /**
    * \brief Returns the total amount of exceed time serving costumers.
    * \return The total delay time in all routes.
    */
    double delayTime()
    {
        return mDelayTime;
    }

    /**
      * \brief Sets the total delay time.
      * \param mDelayTime is the value of the total delay.
      */
    void delayTime(double mDelayTime)
    {
        this->mDelayTime = mDelayTime;
    }

private:

    std::vector<std::vector<int> > mRoutes;   /**< A set of routes containing the decoding information of the individual. */
    int mSizeOfFleet;                         /**< Cached number of vehicled need to serve all costumer defined by the individual. */
    double mLength;                           /**< Cached cost (distance) of traveling the set of routes defined by the individual. */
    double mTime;                             /**< Cached cost (time) of traveling the set of routes defined by the individual. */
    double mWaitingTime;                      /**< Cached of the amount of time vehicles have to wait to start serving defined by the individual. */
    double mDelayTime;                        /**< Cached of the amount of time by which the arrival is late defined by this individual. */
};

#endif
