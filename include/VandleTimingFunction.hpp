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
#ifndef __VANDLETIMINGFUNCITON__HPP__
#define __VANDLETIMINGFUNCITON__HPP__

class VandleTimingFunction {
public:
    VandleTimingFunction() {};
    ~VandleTimingFunction(){};

    double operator() (double *x, double *p);

    void SetBeta(const double &a){beta_ = a;};
    void SetGamma(const double &a){gamma_ = a;};
private:
    double beta_, gamma_;
};
#endif
