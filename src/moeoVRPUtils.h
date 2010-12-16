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

#ifndef moeoVRPUtils_h
#define moeoVRPUtils_h

// General includes
#include <cassert>
#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

/**
  * \def PI
  * Guess you know what this constant represents.
  */

#define PI                   3.14159265

/**
  * \def VEHICLE_CAPACITY
  * Hard-coded parameter for the capacity of the vehicles. This
  * should be parametrized in a config file in a future version.
  */

// #define VEHICLE_CAPACITY   200

/**
  * Max capacity of each vehicle within the fleet
  */

double vehicleCapacity = 200;

/**
  * Max number of usable vehicles
  */

unsigned sizeOfFleet = 0;

typedef std::vector<int> Route;
typedef std::vector< Route > Routes;


/**
  * \namespace moeoVRPUtils
  * \brief A set of structures and utility functions for the VRP-TW problem.
  */

namespace moeoVRPUtils {

/**
* \var typedef struct ClientData ClientDataT.
* \brief Renaming of struct ClientData.
*/

/**
* \struct ClientData
* \brief Information regarding each client in the dataset.
* This structure is intended to be used to store the information of each
* client read from the data file.
*/

typedef struct ClientData
{
    unsigned id;            /**< Client ID number. */
    double   x;             /**< Client's 'x' position in the map. */
    double   y;             /**< Client's 'y' position in the map. */
    double   demand;        /**< Client's demand of delivered product. */
    double   readyTime;     /**< Client's beginning of the time window. */
    double   dueTime;       /**< Client's end of the time window. */
    double   serviceTime;   /**< Client's service time (time needed to serve the product). */

} ClientDataT;


static std::vector <ClientDataT> clients;             /**< Vector to store clients's information. */
static std::vector <std::vector <double> > dist;      /**< Distance matrix. */
static std::vector <std::vector <double> > time;      /**< Time matrix */

/**
   * \brief Computes the distance between two clients.
   * The computed distances will be stored in dist.
   */

void computeDistances ()
{
    unsigned numClients = clients.size ();

    dist.resize (numClients) ;

    for (unsigned i = 0; i < dist.size (); i ++)
        dist [i].resize (numClients);

    // Distances computation
    for (unsigned i = 0; i < dist.size (); i ++)
        for (unsigned j = i + 1 ; j < dist.size (); j ++)
        {
            double distX = clients [i].x - clients [j].x;
            double distY = clients [i].y - clients [j].y;

            dist [i][j] = dist [j][i] = sqrt (distX * distX + distY * distY);

        }

}


/**
   * \brief Returns the time window information of a given client.
   * \param _client The client whose information we want to know.
   * \param _readyTime Return value. The beginning of the client's time window.
   * \param _dueTime Return value. The end of the client's time window.
   * \param _serviceTime Return value. The client's service time.
   */

void getTimeWindow (unsigned _client, double& _readyTime, double& _dueTime, double& _serviceTime)
{
    assert (_client >= 0 && _client < clients.size ());

    _readyTime = clients [_client].readyTime;
    _dueTime = clients [_client].dueTime;
    _serviceTime = clients [_client].serviceTime;

}


/**
   * \brief A function to get the distance between two clients.
   * \param _from The first client.
   * \param _to The second client.
   * \return The distance between _from and _to.
   */

double distance (unsigned _from,  unsigned _to)
{
    if (!(_from >= 0 && _from < clients.size ()) || !(_to   >= 0 && _to   < clients.size ()))
    {
        std::cout << "_from: " << _from << " _to: " << _to << std::endl;
        throw std::runtime_error(std::string("Safety check NOT passed!. Error: Costumer ID out of range \n - from: moeoVRPUtils::distance").c_str());
    }

    return dist [_from][_to];
}

/**
   * \brief A function to get the distance between two clients.
   * \param _from The first client.
   * \param _to The second client.
   * \return The distance between _from and _to.
   */

double elapsedTime (unsigned _from,  unsigned _to)
{
    assert (_from >= 0 && _from < clients.size ());
    assert (_to   >= 0 && _to   < clients.size ());

    return time [_from][_to];
}

/**
   * \brief Computes de polar angle between clients.
   * \param _from The first client.
   * \param _to The second client.
   * \return The polar angle between _from and _to.
   */

float polarAngle (unsigned _from, unsigned _to)
{
    assert (_from >= 0 && _from < clients.size ());
    assert (_to   >= 0 && _to   < clients.size ());

    double angle = atan2 (clients [_from].y - clients [_to].y,
                          clients [_from].x - clients [_to].x);

    // To convert it from radians to degrees
    angle *= 180 / PI;

    if (angle < 0)
        angle *= -1;

    return angle;

}


/**
   * \brief Loads the problem distance matrix data from a given file.
   * \param _fileName The file to load distance matrix data from.
   * \warning No error check is performed!
   */

void loadDistanceMatrix(const char* _filename)
{
    //std::cout << " - Reading Distance Matrix";
    std::ifstream file;
    file.open(_filename);

    // - Resize Distance Matrix -
    dist.resize(clients.size());
    for (size_t i = 0; i < dist.size(); i++)
       dist[i].resize(clients.size());


    for (size_t i = 0; i < dist.size(); i++)
       for (size_t j = 0; j < dist[i].size(); j++)
          file >> dist[i][j];

    //std::cout << "...[OK]" << std::endl;
}

/**
   * \brief Loads the problem time matrix data from a given file.
   * \param _fileName The file to load time matrix data from.
   * \warning No error check is performed!
   */

void loadTimeMatrix(const char* _filename)
{
    //std::cout << " - Reading Time Matrix";
    std::ifstream file;
    file.open(_filename);


    // - Resize Time Matrix -
    time.resize(clients.size());
    for (size_t i = 0; i < time.size(); i++)
       time[i].resize(clients.size());


    for (size_t i = 0; i < time.size(); i++)
       for (size_t j = 0; j < time[i].size(); j++)
          file >> time[i][j];


    //std::cout << "...[OK]" << std::endl;
}

/**
   * \brief Loads the problem data from a given file.
   * \param _fileName The file to load data from.
   * \warning No error check is performed!
   */

void load (const char* _fileName)
{
    std::ifstream f (_fileName);

    if (f)
    {

        /** Reading the header of the file **/
        std::string dummy;
        for (int i = 1; i < 10; i++)
        {
           // The 5th line contains the max allowed number of vehicles and capacity
           if (i == 5)
           {
              f >> sizeOfFleet; // We're not going to use the max allowed num of vehicles.
              f >> vehicleCapacity;
           }
           else
              getline(f, dummy);
        }

        //std::cout << " - The max allowed capacity is " << vehicleCapacity << std::endl;
        //std::cout << " - Size of the fleet is " << sizeOfFleet << std::endl;

        ClientDataT client;
        while (f >> client.id) {


            f >> client.x;
            f >> client.y;
            f >> client.demand;
            f >> client.readyTime;
            f >> client.dueTime;
            f >> client.serviceTime;

            clients.push_back (client);

        }

        f.close ();
    }
    else {

        std::cerr << "Error: the file: " << _fileName << " doesn't exist !!!" << std::endl ;
        exit (1);

    }

}


/**
  * \brief Prints a route to the standard output.
  * \param _route The route to print.
  */

void printRoute (const Route& _route)
{
    std::cout << "[";

    for (unsigned i = 0; i < _route.size (); i++)
    {
        std::cout << _route [i];

        if (i != _route.size () -1)
            std::cout << ", ";
    }

    std::cout << "]";

}


/**
  * \brief Prints a set of routes to the standard output.
  * \param _routes The set of routes to print.
  */

void printRoutes (const Routes& _routes)
{
    for (unsigned i = 0; i < _routes.size (); i++)
    {
        std::cout << "[";

        for (unsigned j = 0; j < _routes [i].size (); j++)
        {
            std::cout << _routes [i][j];

            if (j != _routes [i].size () -1)
                std::cout << ", ";
        }

        if (i == _routes.size () -1)
            std::cout << "]";
        else
            std::cout << "]," << std::endl;
    }

    std::cout << "]" << std::endl;

}


/**
  * \brief Swaps the position of two costumers.
  * \param _routePlan The routePlan in which we want to swap costumers.
  * \param _i Position of the first costumer.
  * \param _j Position of the second costumer.
  */
void swap(Route& _routePlan, unsigned i, unsigned j)
{
    assert(_routePlan.size() > j && _routePlan.size() > i);
    assert(i > 0 && j > 0);

    int temp = _routePlan[i];
    _routePlan[i] = _routePlan[j];
    _routePlan[j] = temp;
}


/** \brief Evaluate the travel distance of a route plan
  * \param _routePlan The route plan to calculate its length
  * \return length of the given route plan.
  */
double travelDistance(const std::vector<int>& _routePlan)
{
    double length = 0.0;
    for (size_t i = 1; i < _routePlan.size(); i++)
       length += dist[_routePlan[i - 1]][_routePlan[i]];
    return length;
}

/**
  * \brief Checks if a routePlan is feasible in terms of capacity.
  * \param _routePlan The route plan in which we want to swap costumers.
  * \return <true> if is feasible and <false> otherwise.
  */
bool feasibleCapacity(const Route& _routePlan)
{
    double currentCapacity = 0.0;

    for (size_t i = 1; i < _routePlan.size(); i++)
    {
       if (_routePlan[i] == 0)
          currentCapacity = 0.0;

       else
       {
          currentCapacity += clients[_routePlan[i]].demand;
          if (currentCapacity > vehicleCapacity)
             return false;
       }
    }
    return true;
}

/**
  * \brief Calculates the time that takes from a costumer<i> to a costumer<j>.
  * \param _totalElapsedTime It's the travel time of a given route up to costumer<i>.
  * \param _numberOfViolations Counts the number of violations or costumer not satisfaied within the given time.
  * \param _i Index of the costumer the delivery vehicle departs from.
  * \param _j Index of the costumer the delivery vehicle arrives at.
  * \return <true> if is feasible (no costumer is served late) and <false> otherwise.
  */
bool elapsedTimeBetweenTwoCostumers(double& _totalElapsedTime, unsigned& _numberOfViolations, unsigned _i, unsigned _j)
{

   // Time it takes to go from <costumer i> to <costumer j>
   _totalElapsedTime += time[_i][_j];

   // Time we have to wait if we arrive before the <costumer j> opens
   //    First, we calculate whether we're going to wait or not.
   if (clients[_j].readyTime > _totalElapsedTime)
      _totalElapsedTime = clients[_j].readyTime;
   else if (clients[_j].dueTime < _totalElapsedTime)
   {
      _numberOfViolations++;
      return false;
   }

   // Finally, we have to sum up the servicetime at the <costumer i>
   _totalElapsedTime += clients[_j].serviceTime;
   return true;
}

/**
  * \brief Checks if a routePlan is feasible in terms of time windows.
  * \param _routePlan The route plan to be checked.
  * \return <true> if is feasible and <false> otherwise.
  */
bool feasibleTimeWindows(const Route& _routePlan)
{
    double _routeElapsedTime = 0;
    unsigned _numberOfViolations = 0;

    for (size_t i = 0; i < _routePlan.size() - 1; i++)
    {
       if (_routePlan[i] == 0)
          _routeElapsedTime = 0;

       if (!elapsedTimeBetweenTwoCostumers(_routeElapsedTime, _numberOfViolations, _routePlan[i], _routePlan[i + 1]))
           return false;
    }

    return true;
}

bool safetyCheck(const Route& _routePlan, std::string methodName = "")
{
    Route costumers(dist.size(), 0);
    unsigned numberOfCostumers = 0;
    unsigned numberOfConsecutivesZeros = 0;

    // Check for invalid values
    for (size_t i = 0; i < _routePlan.size(); i++)
    {
        if (_routePlan[i] != 0)
        {
           if (costumers[_routePlan[i]] != 0)
              throw std::runtime_error(std::string("Safety check NOT passed!. Error: Repeated costumer ID \n - from: " + methodName).c_str());
           else
           {
              costumers[_routePlan[i]]++;
              numberOfCostumers++;
           }
           // Reset the counter for consecutive zeros
           numberOfConsecutivesZeros = 0;
        }
        else
        {
           if (numberOfConsecutivesZeros > 1)
              throw std::runtime_error(std::string("Safety check NOT passed!. Error: More than one consecutive zero - from: " + methodName).c_str());
           numberOfConsecutivesZeros++;
        }


        if (int(_routePlan[i]) != _routePlan[i])
           throw std::runtime_error(std::string("Safety check NOT passed!. Error: Costumer ID is not int \n - from: " + methodName).c_str());

        if (_routePlan[i] > dist.size() || _routePlan[i] < 0)
           throw std::runtime_error(std::string("Safety check NOT passed!. Error: Invalid costumer ID \n - from: " + methodName).c_str());

    }
    // Check if all costumers are served
    if (dist.size() - 1 != numberOfCostumers)
       throw std::runtime_error(std::string("Safety check NOT passed!. Error: Not all costumer are served \n - from: " + methodName).c_str());

    // Safety check passed!
    return 0;
}

bool safetyCheck(const std::vector<std::vector<int> >& _routePlan, std::string methodName = "")
{
    Route costumers(dist.size(), 0);
    unsigned numberOfCostumers = 0;

    // Are all costumers in?
    for (size_t i = 0; i < _routePlan.size(); i++)
        for (size_t j = 0; j < _routePlan[i].size(); j++)
        {
            costumers[_routePlan[i][j]]++;
            numberOfCostumers++;

            if (_routePlan[i][j] > dist.size() + 1 || _routePlan[i][j] < 1)
            {
                std::cout << "Costumer ID: " << _routePlan[i][j] << " for i = " << i << " and j = "<< j << " Limit: " << dist.size() + 1 << std::endl;
                printRoutes(_routePlan);
                throw std::runtime_error(std::string("Safety check NOT passed!. Error: Costumer ID out of range \n - from: " + methodName).c_str());
            }
        }

    // Costumer 0 does not exist, it is the depot.
    for (size_t i = 1; i < costumers.size(); i++)
    {
        if (costumers[i] != 1)
        {
           std::cout << "Costumer: " << i << " is being served " << costumers[i] << " times" << std::endl;
           printRoutes(_routePlan);
           throw std::runtime_error(std::string("Safety check NOT passed!. Error: At least a costumer is not served\n - from: " + methodName).c_str());

        }
    }
    // Safety check passed!
    return 0;
}



}
#endif
