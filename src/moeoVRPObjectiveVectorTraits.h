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

#ifndef CVRPTWOBJECTIVEVECTORTRAITS_H
#define CVRPTWOBJECTIVEVECTORTRAITS_H

#include <moeo>


class moeoVRPObjectiveVectorTraits : public moeoObjectiveVectorTraits
{
   public:
       static bool minimizing (int i) { return true; }
       static bool maximizing (int i) { return false; }
       static unsigned int nObjectives () { return 5; }
};


#endif // CVRPTWOBJECTIVEVECTORTRAITS_H
