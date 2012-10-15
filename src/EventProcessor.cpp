/*! \file EventProcessor.cpp
 *
 * implementation of a generic event processor
 * David Miller, August 2009
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <unistd.h>
#include <sys/times.h>

#include "DetectorDriver.hpp"
#include "EventProcessor.hpp"
#include "RawEvent.hpp"

using namespace std;

extern RawEvent rawev; // to access detector summaries

// correlator is declared in PixieStd.cpp
extern map<string, Place*> correlator;
// Static field initialization
map<string, Place*>& EventProcessor::correlator = correlator;

EventProcessor::EventProcessor() : 
  userTime(0.), systemTime(0.), name("generic"), initDone(false), 
  didProcess(false), histo(0, 0, PlotsRegister::R() )
{
    clocksPerSecond = sysconf(_SC_CLK_TCK);
    ofstream hislog("histograms.txt");
    hislog << "Non empty histograms:" << endl;
    hislog.close();
}

EventProcessor::EventProcessor(int offset, int range) : 
  userTime(0.), systemTime(0.), name("generic"), initDone(false), 
  didProcess(false), histo(offset, range, PlotsRegister::R() ) {
    clocksPerSecond = sysconf(_SC_CLK_TCK);
    ofstream hislog("histograms.txt");
    hislog << "Non empty histograms:" << endl;
    hislog.close();
}

EventProcessor::~EventProcessor() 
{
    if (initDone) {
	// output the time usage
	cout << "processor " << name << " : " 
	     << userTime << " user time, "
	     << systemTime << " system time" << endl;
        ofstream hislog("histograms.txt", ios_base::app);
        hislog << "In " << name << ": " << endl;
        histo.PrintNonEmpty(hislog);
        hislog.close();
    }
}

/** Declare plots */
void EventProcessor::DeclarePlots(void)
{
    // no plots for generic processor
}

/** See if the detectors of interest have any events */
bool EventProcessor::HasEvent(void) const
{
    for (map<string, const DetectorSummary*>::const_iterator it = sumMap.begin();
	 it != sumMap.end(); it++) {
	if (it->second->GetMult() > 0) {
	    return true;
	}
    }
    return false;
}

/** Initialize the processor if the detectors that require it are used in 
 * the analysis
 */
bool EventProcessor::Init(DetectorDriver &driver) 
{
    vector<string> intersect;   
    const set<string> &usedDets = driver.GetUsedDetectors();
    
    set_intersection(associatedTypes.begin(), associatedTypes.end(),
		     usedDets.begin(), usedDets.end(), 
		     back_inserter(intersect) );
    
    if (intersect.empty()) {
	return false;
    }

    // make the corresponding detector summary
    for (vector<string>::const_iterator it = intersect.begin();
	 it != intersect.end(); it++) {
	sumMap.insert( make_pair(*it, rawev.GetSummary(*it)) );
    }

    initDone = true;
    cout << "processor " << name << " initialized operating on " 
	 << intersect.size() << " detector type(s)." << endl;

    // cout << "  adding detector summary " << iSum->first
    //	    << " at address " << &iSum->second << endl;
	    
    return true;
}

/** Process an event. In PreProcess correlator is filled and basic analysis is done.
 * More sophisiticated analysis (which might also depend on other processors) should be
 * done in Process() function. PreProcess will be first called for all Processors and then
 * the Process function will be called.*/
bool EventProcessor::PreProcess(RawEvent &event)
{
    if (!initDone)
        return (didProcess = false);
    return (didProcess = true);
}

/** Process an event. PreProcess function should fill correlation tree and all processors
 * should have basic parameters calculated during PreProccessing.*/
bool EventProcessor::Process(RawEvent &event)
{
    if (!initDone)
        return (didProcess = false);

    // start the process timer
    times(&tmsBegin);

    EndProcess();
    return (didProcess = true);
}

/** Wrap up the processing and update the time spent by this processor */
void EventProcessor::EndProcess(void)
{
    tms tmsEnd;

    times(&tmsEnd);

    userTime += (tmsEnd.tms_utime - tmsBegin.tms_utime) / clocksPerSecond;
    systemTime += (tmsEnd.tms_stime - tmsBegin.tms_stime) / clocksPerSecond;

    // reset the beginning time so multiple calls of EndProcess from
    //   derived classes work properly
    times(&tmsBegin);
}

#ifdef useroot
/** This functions adds the branch to the tree that will be responsible 
 * for holding the data generated by this event processor
 */
bool EventProcessor::AddBranch(TTree *tree)
{
    return false;
}

/** This function is called to fill the appropriate root branch. Note that
 * since ROOT branches can't be empty that the data needs to be properly
 * zeroed for events where no detectors of interest to this processor 
 * triggered.
 */
void EventProcessor::FillBranch(void)
{
    // Do nothing
    //   Add code in this function for analysis to put in ROOT tree only
}
#endif // useroot


