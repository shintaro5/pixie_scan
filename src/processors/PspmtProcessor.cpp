#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <signal.h>
#include <limits.h>

#include "PspmtProcessor.hpp"
#include "DammPlotIds.hpp"
#include "Globals.hpp"
#include "Messenger.hpp"

using namespace std;
using namespace dammIds::pspmt;

namespace dammIds{
    namespace pspmt{
        
        // OFFSET = 700//    
        const int D_RAW1=0;
        const int D_RAW2=1;
        const int D_RAW3=2;
        const int D_RAW4=3;
        const int D_RAWD=4;
        const int D_SUM=5;
        const int DD_POS1_RAW=6;
        const int DD_POS2_RAW=7;
        const int DD_POS1=8;
        const int DD_POS2=9;
        
        const int D_ENERGY_TRACE1=10;
        const int D_ENERGY_TRACE2=11;
        const int D_ENERGY_TRACE3=12;
        const int D_ENERGY_TRACE4=13;
        const int D_ENERGY_TRACED=14;
        const int D_ENERGY_TRACESUM=15;
        const int DD_POS1_RAW_TRACE=16;
        const int DD_POS2_RAW_TRACE=17;
        const int DD_POS1_TRACE=18;
        const int DD_POS2_TRACE=19;
        
        const int D_QDC_TRACE1=20;
        const int D_QDC_TRACE2=21;
        const int D_QDC_TRACE3=22;
        const int D_QDC_TRACE4=23;
        const int D_QDC_TRACED=24;
        
        const int DD_ESLEW=30;
        
        
        const int DD_SINGLE_TRACE=77;
    }
}

void PspmtProcessor::PspmtData::Clear(void) {    
}

PspmtProcessor::PspmtProcessor(void) : EventProcessor(OFFSET, RANGE, "pspmt") {
    associatedTypes.insert("pspmt");
}

void PspmtProcessor::DeclarePlots(void) {
 
  const int posBins      = 32; 
  const int energyBins   = 8192;
  const int traceBins    = 128;
  const int traceBins2   = 512;
  const int Bins         = 2500;
  
    //offset 1900
    DeclareHistogram1D(D_RAW1, energyBins, "Pspmt1 Raw");
    DeclareHistogram1D(D_RAW2, energyBins, "Pspmt2 Raw");
    DeclareHistogram1D(D_RAW3, energyBins, "Pspmt3 Raw");
    DeclareHistogram1D(D_RAW4, energyBins, "Pspmt4 Raw");
    DeclareHistogram1D(D_RAWD, energyBins, "Pspmt Dynode");
    DeclareHistogram1D(D_SUM,  energyBins, "Pspmt Sum");
    DeclareHistogram2D(DD_POS1_RAW, Bins, Bins, "Pspmt Pos1 Raw");
    DeclareHistogram2D(DD_POS2_RAW, Bins, Bins, "Pspmt Pos2 Raw");
    DeclareHistogram2D(DD_POS1, posBins, posBins, "Pspmt Pos1");
    DeclareHistogram2D(DD_POS2, posBins, posBins, "Pspmt Pos2");
    
    // Trace energies 
    // 710-
    DeclareHistogram1D(D_ENERGY_TRACE1, energyBins, "Energy1 from trace");
    DeclareHistogram1D(D_ENERGY_TRACE2, energyBins, "Energy2 from trace");
    DeclareHistogram1D(D_ENERGY_TRACE3, energyBins, "Energy3 from trace");
    DeclareHistogram1D(D_ENERGY_TRACE4, energyBins, "Energy4 from trace");
    DeclareHistogram1D(D_ENERGY_TRACED, energyBins, "EnergyD from trace");
    DeclareHistogram1D(D_ENERGY_TRACESUM,  energyBins, "Pspmt Sum");
    DeclareHistogram2D(DD_POS1_RAW_TRACE, posBins, posBins, "Pspmt pos Raw by Trace1");
    DeclareHistogram2D(DD_POS2_RAW_TRACE, posBins, posBins, "Pspmt pos Raw by Trace2");
    DeclareHistogram2D(DD_POS1_TRACE, posBins, posBins, "Pspmt pos by Trace1");
    DeclareHistogram2D(DD_POS2_TRACE, posBins, posBins, "Pspmt pos by Trace2");
    
    // QDCs
    DeclareHistogram1D(D_QDC_TRACE1, energyBins, "Energy1 from QDC");
    DeclareHistogram1D(D_QDC_TRACE2, energyBins, "Energy2 from QDC");
    DeclareHistogram1D(D_QDC_TRACE3, energyBins, "Energy3 from QDC");
    DeclareHistogram1D(D_QDC_TRACE4, energyBins, "Energy4 from QDC");
    DeclareHistogram1D(D_QDC_TRACED, energyBins, "EnergyD from QDC");
    
    // Trace
    DeclareHistogram2D(DD_SINGLE_TRACE, traceBins, traceBins2,"Single trace");
    
}


bool PspmtProcessor::PreProcess(RawEvent &event){
    if (!EventProcessor::PreProcess(event))
        return false;
    
    static const vector<ChanEvent*> &pspmtEvents = sumMap["pspmt"]->GetList();
    
    data_.Clear();
    
    double q1=0,q2=0,q3=0,q4=0,qd=0;
    double qdc1=0,qdc2=0,qdc3=0,qdc4=0,qdcd=0;
    double tre1=0,tre2=0,tre3=0,tre4=0,tred=0;
    
    double qright=0,qleft=0,qtop=0,qbottom=0,qsum=0;
    double xright=0,xleft=0,ytop=0,ybottom=0;
    
    static int traceNum;
    
    for (vector<ChanEvent*>::const_iterator it = pspmtEvents.begin();
         it != pspmtEvents.end(); it++) {
        
        ChanEvent *chan   = *it;
        string subtype    = chan->GetChanID().GetSubtype();
        int    ch         = chan->GetChanID().GetLocation();
        double calEnergy  = chan->GetCalEnergy();
        double pspmtTime  = chan->GetTime();
        Trace trace       = chan->GetTrace();
        
        double trace_energy;
        double trace_time;
        double baseline;
        double qdc;
        
        if(trace.HasValue("filterEnergy")){
            traceNum++;   	  
            trace_time    = trace.GetValue("filterTime");
            trace_energy  = trace.GetValue("filterEnergy");
            baseline      = trace.DoBaseline(2,20);
            qdc             = trace.DoQDC(5,128);
            
            if(ch==0){
                qdc1 = qdc;
                tre1 = trace_energy;
                plot(D_QDC_TRACE1,qdc1);
                plot(D_ENERGY_TRACE1,tre1);
            }else if(ch==1){
                qdc2 = qdc;
                tre2 = trace_energy; 
                plot(D_QDC_TRACE2,qdc2);
                plot(D_ENERGY_TRACE2,tre2);
            }else if(ch==2){
                qdc3 = qdc;
                tre3 = trace_energy; 
                plot(D_QDC_TRACE3,qdc3);
                plot(D_ENERGY_TRACE3,tre3);
            }else if(ch==3){
                qdc4 = qdc;
                tre4 = trace_energy; 	  
                plot(D_QDC_TRACE4,qdc4);
                plot(D_ENERGY_TRACE4,tre4);
            }else if(ch==4){
                qdcd = qdc;
                tred = trace_energy; 
                plot(D_QDC_TRACED,qdcd);
                plot(D_ENERGY_TRACED,tred);
            }
        }

        if(ch==0){
            q1= calEnergy;
            plot(D_RAW1,q1);
        }else if(ch==1){
            q2= calEnergy;
            plot(D_RAW2,q2);
        }else if(ch==2){
            q3= calEnergy;
            plot(D_RAW3,q3);
        }else if(ch==3){
            q4= calEnergy;
            plot(D_RAW4,q4);
        }else if(ch==4){
            qd= calEnergy;
            plot(D_RAWD,qd);
        }
        
        if(q1>0 && q2>0 && q3>0 && q4>0){
            qtop    = (q1+q2)/2;
            qleft   = (q2+q3)/2;
            qbottom = (q3+q4)/2;
            qright  = (q4+q1)/2;
            
            qsum    = (q1+q2+q3+q4)/2;
            xright  = (qright/qsum)*512+100;
            xleft   = (qleft/qsum)*512+100;
            ytop    = (qtop/qsum)*512+100;
            ybottom = (qbottom/qsum)*512+100;
         
        }
        
        
   
	for(vector<int>::iterator ittr = trace.begin();ittr != trace.end();ittr++)
                plot(DD_SINGLE_TRACE,ittr-trace.begin(),traceNum,*ittr);
   
    } // end of channel event
    
    EndProcess();
    return(true);
}

bool PspmtProcessor::Process(RawEvent &event){
    if (!EventProcessor::Process(event))
        return false;

    EndProcess();
    return(true);
}
