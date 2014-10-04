 /**************************************************************************
  *  Copyright S. V. Paulauskas 2014                                       *
  *                                                                        *
  *  This program is free software: you can redistribute it and/or modify  *
  *  it under the terms of the GNU General Public License as published by  *
  *  the Free Software Foundation, version 3.0 License.                    *
  *                                                                        *
  *  This program is distributed in the hope that it will be useful,       *
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
  *  GNU General Public License for more details.                          *
  *                                                                        *
  *  You should have received a copy of the GNU General Public License     *
  *  along with this program.  If not, see <http://www.gnu.org/licenses/>. *
  **************************************************************************
*/

/*! \file VandleTimingFunction.hpp
 *  \brief A class to handle the processing of traces
 *  \author S. V. Paulauskas
 *  \date 03 October 2014
 */
#include <cmath>

#include "VandleTimingFunction.hpp"

double VandleTimingFunction::operator()(double *x, double *par) {
    //double beta = 0.17706971465472;
    //double gamma = 0.162673945899791;

    double phase = par[0];
    double amplitude = par[1];
    double beta = par[2];
    double gamma = par[3];
    double diff = x[0] - phase;

    if(x[0] < phase)
        return(15217.6388888889);

    double val = amplitude * exp(-beta*diff) * (1-exp(-pow(gamma*diff,4.)))
        +15217.6388888889;

    return(val);
}
