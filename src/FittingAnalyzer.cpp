/** \file FittingAnalyzer.cpp
 * \brief Uses a chi^2 minimization to fit waveforms using ROOT
 *
 * Obtains the phase of a waveform using a Chi^2 fitting algorithm
 * implemented through ROOT
 *
 * \author S. V. Paulauskas
 * \date 22 July 2011
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <vector>

#include <ctime>

#include <TF1.h>
#include <TGraphErrors.h>
#include <TFitResult.h>

#include "DammPlotIds.hpp"
#include "FittingAnalyzer.hpp"
#include "VandleTimingFunction.hpp"

namespace dammIds {
    namespace trace {
        namespace fitting {
            const int DD_TRACES      = 0;
            const int D_CHISQPERDOF  = 1;
            const int D_PHASE        = 2;
            const int DD_AMP         = 3;
            const int D_SAT          = 4;
            const int DD_MAXVSQDCMAX = 5;
            const int DD_MAXVALPOS   = 6;
            const int DD_QDCMASK     = 7;
            const int DD_MAXVSTHRESH = 8;
            const int D_SIGMA        = 9;
        }
    }
}

using namespace std;
using namespace dammIds::trace::fitting;

//********** DeclarePlots **********
void FittingAnalyzer::DeclarePlots(void)
{
    DeclareHistogram2D(DD_TRACES, SB, S7, "traces data");
    DeclareHistogram2D(DD_AMP, SE, SC, "Fit Amplitude");
    DeclareHistogram1D(D_PHASE, SE, "Fit X0");
    DeclareHistogram1D(D_CHISQPERDOF, SE, "Chi^2/dof");
    DeclareHistogram1D(D_SAT, S4, "Saturations");
    DeclareHistogram2D(DD_MAXVSQDCMAX, SB, SC, "Max to Max Value Ratio");
    DeclareHistogram2D(DD_MAXVALPOS, S5, SC, "Max Val vs Pos");
    DeclareHistogram2D(DD_QDCMASK, SE, SC, "Max vs Reduced Chi^2");
    DeclareHistogram2D(DD_MAXVSTHRESH, S7, SC, "Max vs Num Bins Tresh");
    DeclareHistogram1D(D_SIGMA, SE, "Standard Dev Baseline");
}


//********** FittingAnalyzer **********
FittingAnalyzer::FittingAnalyzer() : TraceAnalyzer(OFFSET,RANGE)
{
    name = "FittingAnalyzer";
}


//********** Analyze **********
void FittingAnalyzer::Analyze(Trace &trace, const string &detType,
			      const string &detSubtype)
{
    TraceAnalyzer::Analyze(trace, detType, detSubtype);

    if(trace.HasValue("saturation") || trace.empty()) {
	plot(D_SAT,2);
        EndAnalyze();
     	return;
    }

    const double sigmaBaseline = trace.GetValue("sigmaBaseline");
    const double maxVal = trace.GetValue("maxval");
    const double qdc = trace.GetValue("tqdc");
    const double qdcToMax = trace.GetValue("qdcToMax");
    const unsigned int maxPos = (unsigned int)trace.GetValue("maxpos");
    const vector<double> waveform = trace.GetWaveform();

    if(waveform.size() == 0) {
        EndAnalyze();
        return;
    }

    plot(DD_MAXVSQDCMAX, qdcToMax*100+100, maxVal);
    plot(DD_MAXVALPOS, maxPos, maxVal);
    plot(D_SIGMA, sigmaBaseline*100);

    if(sigmaBaseline > 3.0) {
        EndAnalyze();
        return;
    }

    double beta, gamma;
    if (detType == "vandleSmall") {
        beta  = TimingInformation::GetConstant("betaVandle");
        gamma = TimingInformation::GetConstant("gammaVandle");
    }else if (detSubtype == "beta") {
        beta  = TimingInformation::GetConstant("betaBeta");
        gamma = TimingInformation::GetConstant("gammaBeta");
    }else if(detType == "tvandle") {
        beta  = TimingInformation::GetConstant("betaTvandle");
        gamma = TimingInformation::GetConstant("gammaTvandle");
    }else if(detType == "pulser") {
        beta  = TimingInformation::GetConstant("betaPulser");
        gamma = TimingInformation::GetConstant("gammaPulser");
    }else {
        beta  = TimingInformation::GetConstant("betaDefault");
        gamma = TimingInformation::GetConstant("gammaDefault");
    }

    vector<double> xvals;
    for(unsigned int i = 0; i < waveform.size(); i++)
        xvals.push_back(i);

    VandleTimingFunction *fobj = new VandleTimingFunction();
    TF1 *f = new TF1("f", fobj, 0., 1000., 2, "VandleTimingFunction");
    f->SetParameters(5.0, qdc*0.5);

    TGraphErrors *graph =
        new TGraphErrors(waveform.size(), &(xvals[0]), &(waveform[0]));

    for(unsigned int i = 0; i < waveform.size(); i++)
        graph->SetPointError(i, 0.0, sigmaBaseline);

    TFitResultPtr fitResults = graph->Fit(f, "MENSQR", "", 0.,waveform.size());
    double phase = fitResults->Value(0);

    //cout << "Fit Status : " << fitResults << " " << fitResults->Value(0)
    //     << " " << fitResults->Value(1) << endl;

    delete(f);

    trace.InsertValue("phase", phase+maxPos);
    trace.InsertValue("walk", CalcWalk(maxVal, detType, detSubtype));

    //plot(DD_AMP, fitPars.at(1), maxVal);
    //plot(D_PHASE, fitPars.at(0)*1000+100);
    //plot(D_CHISQPERDOF, chisqPerDof);
    //plot(DD_QDCMASK, chisqPerDof, maxVal);

    EndAnalyze();
} //void FittingAnalyzer::Analyze



//********** WalkCorrection **********
double FittingAnalyzer::CalcWalk(const double &val, const string &type,
				 const string &subType)
{
    if(type == "vandleSmall") {
	if(val < 175)
	    return(1.09099*log(val)-7.76641);
	if( val > 3700)
	    return(0.0);
	else
	    return(-(9.13743e-12)*pow(val,3.) + (1.9485e-7)*pow(val,2.)
		   -0.000163286*val-2.13918);

	//Original Function - RevD
	// double f = 92.7907602830327 * exp(-val/186091.225414275) +
	// 	0.59140785215161 * exp(val/2068.14618331387) -
	// 	95.5388835298589;
    }else if(subType == "beta") {
	return(-(1.07908*log10(val)-8.27739));
    }else
	return(0.0);
}
