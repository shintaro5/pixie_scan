/** \file FittingAnalyzer.hpp
 * \brief Class to fit functions to waveforms
 */
#ifndef __FITTINGANALYZER_HPP_
#define __FITTINGANALYZER_HPP_

#include "TimingInformation.hpp"
#include "Trace.hpp"
#include "TraceAnalyzer.hpp"

class FittingAnalyzer : public TraceAnalyzer,
			public TimingInformation
{
 public:
    FittingAnalyzer();
    virtual void DeclarePlots(void);
    virtual void Analyze(Trace &, const std::string &, const std::string &);
    virtual ~FittingAnalyzer() {};

 private:
    void LoadMask(void);
    double ApplyMask(const std::vector<double> &waveform,
		     const double &qdc, const double &maxval,
		     const double &sigma);
    double CalcWalk(const double &maxValue, const std::string &type,
		    const std::string &subType);
};
#endif // __FITTINGANALYZER_HPP_
// David is awesome.
