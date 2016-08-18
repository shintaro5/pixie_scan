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
        
      // OFFSET = 1900//    
      const int D_RAW1=0; 
      const int D_RAW2=1; 
      const int D_RAW3=2; 
      const int D_RAW4=3; 
      const int D_RAW5=4; 
      const int D_SUM=5; 
      const int DD_POS1_RAW=6; 
      const int DD_POS2_RAW=7; 
      const int DD_POS1=8; 
      const int DD_POS2=9; 
      
      const int D_TRE1=10; 
      const int D_TRE2=11; 
      const int D_TRE3=12; 
      const int D_TRE4=13; 
      const int D_TRE5=14; 
      const int D_TRESUM=15; 
      const int DD_POS1_TRACE=16; 
      const int DD_POS2_TRACE=17; 
      const int DD_POS1_TRE=18; 
      const int DD_POS2_TRE=19;
      
      const int D_QDC1=20; 
      const int D_QDC2=21; 
      const int D_QDC3=22; 
      const int D_QDC4=23; 
      const int D_QDC5=24; 
      const int DD_POS1_QDC=26; 
      const int DD_POS2_QDC=27;
      
      const int DD_POS_NEW_CHE=30;
      const int DD_POS_NEW_TRE=31;
      const int DD_POS_NEW_QDC=32;
      
      const int DD_DOUBLE_TRACE=77; //!< Double traces
      const int DD_SINGLE_TRACE=78; //!< Single traces
   
      const int D_TEMP1 =90;
      const int D_TEMP2 =91;
      const int D_TEMP3 =92;
      const int D_TEMP4 =93;
   
      const int DD_POS_NEW1=95;
      const int DD_POS_NEW2=96;
      const int DD_EE=97;
      const int DD_EE2=98;
      const int DD_EE3=99;


    }
}

void PspmtProcessor::PspmtData::Clear(void) {    
}

PspmtProcessor::PspmtProcessor(void) : EventProcessor(OFFSET, RANGE, "PspmtProcessor") {
    associatedTypes.insert("pspmt");
}

void PspmtProcessor::DeclarePlots(void) {
    
  const int posBins      = 32; 
  const int energyBins   = 8192;
  const int traceBins    = 248;
  const int traceBins2   = 1024;
  const int Bins         = 2000;
  
  // Raw 700-707
  DeclareHistogram1D(D_RAW1, energyBins, "Pspmt1 Raw");
  DeclareHistogram1D(D_RAW2, energyBins, "Pspmt2 Raw");
  DeclareHistogram1D(D_RAW3, energyBins, "Pspmt3 Raw");
  DeclareHistogram1D(D_RAW4, energyBins, "Pspmt4 Raw");
  DeclareHistogram1D(D_RAW5, energyBins, "Pspmt Dynode");
  DeclareHistogram1D(D_SUM,  energyBins, "Pspmt Sum");
  DeclareHistogram2D(DD_POS1_RAW, Bins, Bins, "Pspmt Pos1 Raw");
  DeclareHistogram2D(DD_POS2_RAW, Bins, Bins, "Pspmt Pos2 Raw");
  DeclareHistogram2D(DD_POS1, posBins, posBins, "Pspmt Pos1");
  DeclareHistogram2D(DD_POS2, posBins, posBins, "Pspmt Pos2");
  
    // From QDC and traces 
    // 710-
  DeclareHistogram1D(D_TRE1, energyBins, "Energy1 from trace");
  DeclareHistogram1D(D_TRE2, energyBins, "Energy2 from trace");
  DeclareHistogram1D(D_TRE3, energyBins, "Energy3 from trace");
  DeclareHistogram1D(D_TRE4, energyBins, "Energy4 from trace");
  DeclareHistogram1D(D_TRE5, energyBins, "EnergyD from trace");
  DeclareHistogram1D(D_TRESUM,  energyBins, "Pspmt Sum");
  DeclareHistogram2D(DD_POS1_TRACE, posBins, posBins, "Pspmt pos Raw by Trace1");
  DeclareHistogram2D(DD_POS2_TRACE, posBins, posBins, "Pspmt pos Raw by Trace2");
  
  // 720- QDC
  DeclareHistogram1D(D_QDC1, energyBins, "Energy1 from QDC");
  DeclareHistogram1D(D_QDC2, energyBins, "Energy2 from QDC");
  DeclareHistogram1D(D_QDC3, energyBins, "Energy3 from QDC");
  DeclareHistogram1D(D_QDC4, energyBins, "Energy4 from QDC");
  DeclareHistogram1D(D_QDC5, energyBins, "EnergyD from QDC");
  DeclareHistogram2D(DD_POS1_QDC, Bins, Bins, "Pspmt pos1 Raw by QDC");
  DeclareHistogram2D(DD_POS2_QDC, Bins, Bins, "Pspmt pos2 Raw by QDC");

  // 730 Histograms for new board
  DeclareHistogram2D(DD_POS_NEW_CHE, Bins, Bins, "Pspmt pos1 Raw(CHE)");
  DeclareHistogram2D(DD_POS_NEW_TRE, Bins, Bins, "Pspmt pos1 Raw(TRE)");
  DeclareHistogram2D(DD_POS_NEW_QDC, Bins, Bins, "Pspmt pos1 Raw(QDC)");
  
  // Trace
  DeclareHistogram2D(DD_DOUBLE_TRACE, traceBins, traceBins2,"Double traces");
  DeclareHistogram2D(DD_SINGLE_TRACE, traceBins, traceBins2,"Single trace");
  
  DeclareHistogram1D(D_TEMP1, energyBins, "Gated QDC1");
  DeclareHistogram1D(D_TEMP2, energyBins, "Gated QDC2");
  DeclareHistogram1D(D_TEMP3, energyBins, "Gated QDC3");
  DeclareHistogram1D(D_TEMP4, energyBins, "Gated QDC4");
  
  DeclareHistogram2D(DD_POS_NEW1, Bins, Bins, "for new board");
  DeclareHistogram2D(DD_POS_NEW2, Bins, Bins, "for new board");
  DeclareHistogram2D(DD_EE, energyBins, energyBins, "for new board dev");
  DeclareHistogram2D(DD_EE2, energyBins, energyBins, "for new board dev2");
  DeclareHistogram2D(DD_EE3, energyBins, energyBins, "qall vs q1");

}


bool PspmtProcessor::PreProcess(RawEvent &event){
    if (!EventProcessor::PreProcess(event))
        return false;
    
    static const vector<ChanEvent*> &pspmtEvents = sumMap["pspmt"]->GetList();

    data_.Clear();
    
    // channel energy, trace energy, qdc
    double che1=0,che2=0,che3=0,che4=0,che5=0;
    double tre1=0,tre2=0,tre3=0,tre4=0,tre5=0; 
    double qdc1=0,qdc2=0,qdc3=0,qdc4=0,qdc5=0;
    
    // parameters for old board
    double xche=0,yche=0;
    double xtre=0,ytre=0;
    double xqdc=0,yqdc=0;
    
    // parameters for new board
    double xnew_che=0,ynew_che=0;
    double xnew_tre=0,ynew_tre=0;
    double xnew_qdc=0,ynew_qdc=0;
    //
    
    static int traceNum;
    
    for (vector<ChanEvent*>::const_iterator it = pspmtEvents.begin();
         it != pspmtEvents.end(); it++) {
        
        ChanEvent *chan   = *it;
        string subtype    = chan->GetChanID().GetSubtype();
        int    ch         = chan->GetChanID().GetLocation();
        double calEnergy  = chan->GetCalEnergy();
	Trace trace       = chan->GetTrace();
        double trace_energy;
        double trace_time;
        double baseline;
        double qdc;
	int maxpos;

	
	// getting and plotting channel energies from 1900 to 1905
	if(ch==0){
	  che1= calEnergy;
	  plot(D_RAW1,che1);
        }else if(ch==1){
	  che2= calEnergy;
	  plot(D_RAW2,che2);
        }else if(ch==2){
	  che3= calEnergy;
	  plot(D_RAW3,che3);
        }else if(ch==3){
	  che4= calEnergy;
	  plot(D_RAW4,che4);
        }else if(ch==4){
	  che5= calEnergy;
	  plot(D_RAW5,che5);
        }

	trace.FindMaxInfo(10,20,150);
	// tentatively always true
	bool t=true;
      	if(trace.GetValue("maxval")>20){
	  //if(trace.HasValue("filterEnergy")){
            traceNum++;   	  
	   

	    
	    trace_energy = trace.GetValue("maxval");
            trace_time    = trace.GetValue("filterTime");
	    baseline      = trace.DoBaseline(2,50);
	    maxpos        = trace.GetValue("maxpos");
	    qdc = trace.DoQDC(50,150);
	    qdc = trace.GetValue("tqdc");
	    cout <<" tre " << trace_energy << " base  " << baseline << " max " << maxpos << " " << qdc <<endl;
	    
	    
	    //trace_energy:GetValue is old method
	    //trace_energy  = trace.GetValue("filterEnergy");
	    //qdc:DoQDC is old method 
	    //qdc             = trace.DoQDC(80,150);
	    
	    //    trace.DoDiscrimination(10,82);
	    
	    //cout << "maxpos " <<  trace.FindMaxInfo(10,20,100) << " getvalue " << maxpos<< " " << trace_energy << endl;
	    // for old data
	    // qdc    = trace.DoQDC(20,100);
	    // for new data
	   
	    //cout << "qdc " << qdc << endl;
	    //cout << "tre "<< trace_energy <<" baseline " << baseline << " qdc  " << qdc << endl;
	    

	  
	    
	    
	      
	    
	    
	    if(ch==0){
	      qdc1 = qdc;
	      tre1 = trace_energy;
	      plot(D_QDC1,qdc1);
	      plot(D_TRE1,tre1);
            }else if(ch==1){
	      qdc2 = qdc;
	      tre2 = trace_energy; 
	      plot(D_QDC2,qdc2);
                plot(D_TRE2,tre2);
            }else if(ch==2){
	      qdc3 = qdc;
	      tre3 = trace_energy; 
	      plot(D_QDC3,qdc3);
	      plot(D_TRE3,tre3);
            }else if(ch==3){
	      qdc4 = qdc;
	      tre4 = trace_energy; 	  
	      plot(D_QDC4,qdc4);
	      plot(D_TRE4,tre4);
            }else if(ch==4){
	      qdc5 = qdc;
	      tre5 = trace_energy; 
	      plot(D_QDC5,qdc5);
	      plot(D_TRE5,tre5);
            }
        }
	
        
	xche = GetPositionX(che1,che2,che3,che4,che5);
	yche = GetPositionY(che1,che2,che3,che4,che5);
	
	xtre = GetPositionX(tre1,tre2,tre3,tre4,tre5);
	ytre = GetPositionY(tre1,tre2,tre3,tre4,tre5);
	
	xqdc = GetPositionX(qdc1,qdc2,qdc3,qdc4,qdc5);
	yqdc = GetPositionY(qdc1,qdc2,qdc3,qdc4,qdc5);

	// New board
	xnew_che = GetPositionXNew(che1,che2,che3,che4,che5);
	ynew_che = GetPositionYNew(che1,che2,che3,che4,che5);
	xnew_tre = GetPositionXNew(tre1,tre2,tre3,tre4,tre5);
	ynew_tre = GetPositionYNew(tre1,tre2,tre3,tre4,tre5);
	xnew_qdc = GetPositionXNew(qdc1,qdc2,qdc3,qdc4,qdc5);
	ynew_qdc = GetPositionYNew(qdc1,qdc2,qdc3,qdc4,qdc5);
	
	if(tre1>50 && tre2>50 && tre3>50 && tre4>50){
	  //	if(che1>40 && che2>40 && che3>40 && che4>40){
	
	  // for old board
	  plot(DD_POS1_RAW,xche,yche);
	  plot(DD_POS1_TRE,xtre,ytre);
	  plot(DD_POS1_QDC,xqdc,yqdc);
	  // for new board
	  plot(DD_POS_NEW_CHE,xnew_che,ynew_che);
	  plot(DD_POS_NEW_TRE,xnew_tre,ynew_tre);
	  plot(DD_POS_NEW_QDC,xnew_qdc,ynew_qdc);
	
	  //	  cout << xnew_che << " " << ynew_che << endl;
	  
	}
	
	plot(DD_EE,che1,che2);
	plot(DD_EE2,che1+che2+1000,che1-che2+1000);
	plot(DD_EE3,che1+che2,che3+che4);
	  
	  for(vector<int>::iterator ittr = trace.begin();ittr != trace.end();ittr++)
	    plot(DD_SINGLE_TRACE,ittr-trace.begin(),traceNum,*ittr);
	  
    } 
    // end of channel event
    
    

    EndProcess();
    return(true);
}

bool PspmtProcessor::Process(RawEvent &event){
    if (!EventProcessor::Process(event))
        return false;

    EndProcess();
    return(true);
}

double PspmtProcessor::GetPositionX(double q1,double q2,double q3, double q4,double q5){
 
  double xleft,xright,qsum,qright,qleft;
  
  qsum   = q1+q2+q3+q4;
  qright = (q4+q1)/2;
  qleft  = (q2+q3)/2;
  
  xleft  = (qright/qsum)*4096+100;
  xright = (qleft/qsum)*4096+100;
  

  return xright;
  //return xleft;
}

double PspmtProcessor::GetPositionY(double q1,double q2,double q3, double q4,double q5){
 
  double ytop,ybottom,qsum,qtop,qbottom;
  
  qsum     = q1+q2+q3+q4;
  qtop     = (q1+q2)/2;
  qbottom  = (q3+q4)/2;
  
  ytop    = (qtop/qsum)*4096+100;
  ybottom = (qbottom/qsum)*4096+100;
  
  return ytop;
  // return ybottom;
}

double PspmtProcessor::GetPositionXNew(double q1,double q2,double q3,double q4,double q5){
  
  double xdiff,xpos,xsum;
  
  xdiff = q1-q2;
  xsum  = q1+q2; 
  xpos  = (xdiff/xsum)*1024+1000;
  
  return xpos;
  
}
double PspmtProcessor::GetPositionYNew(double q1,double q2,double q3,double q4,double q5){
  
  double ydiff,ypos,ysum;
 
  ydiff = q3-q4;
  ysum  = q3+q4;
  ypos  = (ydiff/ysum)*1024+1000;
  
  return ypos;
}
