/** \file FittingAnalyzer.hpp
 * \brief Class to fit functions to waveforms
 */
#ifndef __FITTINGANALYZER_HPP_
#define __FITTINGANALYZER_HPP_

#include "SiPmtFastTimingFunction.hpp"
#include "TimingInformation.hpp"
#include "Trace.hpp"
#include "TraceAnalyzer.hpp"
#include "VandleTimingFunction.hpp"

class FittingAnalyzer : public TraceAnalyzer,
			public TimingInformation
{
 public:
    FittingAnalyzer();
    virtual ~FittingAnalyzer();

    virtual void DeclarePlots(void);
    virtual void Analyze(Trace &, const std::string &, const std::string &);
 private:
    VandleTimingFunction *vandleFunc_;
    SiPmtFastTimingFunction *siPmtFunc_;
    double CalcWalk(const double &maxValue, const std::string &type,
		    const std::string &subType);
};
#endif // __FITTINGANALYZER_HPP_
// David is awesome.
