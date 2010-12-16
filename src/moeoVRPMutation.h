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

#ifndef moeoVRPMutation_H
#define moeoVRPMutation_H

// General includes
#include <algorithm>

// The base definition of eoMonOp
#include <eoOp.h>

// Safety checks
#include "moeoVRPUtils.h"



/**
  * \class moeoVRPSwapMutation moeoVRPMutation.h.h
  * \brief Implementation of the swap mutation operator.
  */

class moeoVRPSwapMutation: public eoMonOp <moeoVRP>
{

public:

    /**
      * \brief Deafult constructor.
      */

    moeoVRPSwapMutation () {  }


    /**
      * \brief Returns a string containing the name of the class. Used to display statistics.
      * \return The string containing the name of the class.
      */

    std::string className () const
    {
        return "moeoVRPSwapMutation";
    }


    /**
      * \brief It exhanges the positions of two clients within the individual.
      * Clients may or may not be in the same route.
      * \param _genotype The genotype being mutated (it will be probably modified).
      * \return True if the individual has been modified. False otherwise.
      */

    bool operator () (moeoVRP& _genotype)
    {
        int p1 = rng.random (_genotype.size ());
        int p2 = -1;

        do {
            p2 = rng.random (_genotype.size ());
        } while (_genotype [p2] == _genotype [p1]);

        std::swap (_genotype [p1], _genotype [p2]);

        return true;
    }


private:


};



/**
  * \class moeoVRPInsertionMutation moeoVRPMutation.h.h
  * \brief Implementation of the insertion mutation operator.
  */

class moeoVRPInsertionMutation: public eoMonOp <moeoVRP>
{

public:

    /**
      * \brief Deafult constructor.
      */

    moeoVRPInsertionMutation () {  }


    /**
      * \brief Returns a string containing the name of the class. Used to display statistics.
      * \return The string containing the name of the class.
      */

    std::string className () const
    {
        return "moeoVRPInsertionMutation";
    }


    /**
      * \brief It selects and individual, erases it from its original position and inserts it somewhere else.
      * The insertion may or may not be within the same route.
      * \param _genotype The genotype being mutated (it will be probably modified).
      * \return True if the individual has been modified. False otherwise.
      */

    bool operator () (moeoVRP& _genotype)
    {
        int p = -1;

        // Selection of the client to be moved
        do {
            p = rng.random (_genotype.size ());
        } while (_genotype [p] == -1);

        // Temporary copy of the client
        unsigned client = _genotype [p];

        _genotype.erase (_genotype.begin () + p);

        p = rng.random (_genotype.size () - 1);
        _genotype.insert (_genotype.begin () + p, client);

        return true;
    }


private:


};


/**
  * \class moeoVRPInversionMutation moeoVRPMutation.h.h
  * \brief Implementation of the inversion mutation operator.
  */

class moeoVRPInversionMutation: public eoMonOp <moeoVRP>
{

public:

    /**
      * \brief Deafult constructor.
      */

    moeoVRPInversionMutation () {  }


    /**
      * \brief Returns a string containing the name of the class. Used to display statistics.
      * \return The string containing the name of the class.
      */

    std::string className () const
    {
        return "moeoVRPInversionMutation";
    }


    /**
      * \brief It selects two positions in the genotype and inverts the clients between them.
      * Clients may or may not be in the same route.
      * \param _genotype The genotype being mutated (it will be probably modified).
      * \return True if the individual has been modified. False otherwise.
      */

    bool operator () (moeoVRP& _genotype)
    {
        int p1 = rng.random (_genotype.size ());
        int p2 = -1;

        do {
            p2 = rng.random (_genotype.size ());
        } while (_genotype [p2] == _genotype [p1]);

        if (p1 > p2)
            std::swap (p1, p2);

        // Reverse the subroute
        reverse (_genotype.begin () + p1, _genotype.begin () + p2 + 1);

        return false;
    }


private:


};


/**
  * \class moeoVRPDisplacementMutation moeoVRPMutation.h.h
  * \brief Implementation of the displacement mutation operator.
  */

class moeoVRPDisplacementMutation: public eoMonOp <moeoVRP>
{

public:

    /**
      * \brief Deafult constructor.
      */

    moeoVRPDisplacementMutation () {  }


    /**
      * \brief Returns a string containing the name of the class. Used to display statistics.
      * \return The string containing the name of the class.
      */

    std::string className () const
    {
        return "moeoVRPDisplacementMutation";
    }


    /**
      * \brief It selects a set of clients, erases them from their original position and inserts them somewhere else.
      * The selected set of clients may cover different routes.
      * \param _genotype The genotype being mutated (it will be probably modified).
      * \return True if the individual has been modified. False otherwise.
      */

    bool operator () (moeoVRP& _genotype)
    {
        int p1 = rng.random (_genotype.size ());
        int p2 = -1;

        do {
            p2 = rng.random (_genotype.size ());
        } while (_genotype [p2] == _genotype [p1]);

        if (p1 > p2)
            std::swap (p1, p2);

        // Temporary copy of the fragment being moved
        Route route;

        for (unsigned i = p1; i <= p2; i++)
            route.push_back (_genotype [i]);

        _genotype.erase (_genotype.begin () + p1, _genotype.begin () + p2 + 1);

        unsigned p = rng.random ((_genotype.size () > 0) ? _genotype.size () - 1 : 0);
        _genotype.insert (_genotype.begin () + p, route.begin (), route.end ());

        return true;

    }


private:


};


#endif
