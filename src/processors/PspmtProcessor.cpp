/** \file PspmtProcessor.cpp
 * \brief Class for development of Pspmt
 *\author S. Go
 *\date August 24, 2016
 */
#include <iostream>
#include <sstream>

#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "GeProcessor.hpp"
#include "GetArguments.hpp"
#include "PspmtProcessor.hpp"

#ifdef useroot
static double x_maxval;
static double y_maxval;

static double x_qdc;
static double y_qdc;

static double x_en;
static double y_en;

static double x_fen;
static double y_fen;
#endif

namespace dammIds{
    namespace pspmt{
       // OFFSET = 1900//    
      const int D_RAW1=0;
      const int D_RAW2=1;
      const int D_RAW3=2;
      const int D_RAW4=3;
      const int D_RAW5=4;
      const int D_SUM=5;
      const int DD_POS_CHE=6;
      
      const int D_ENERGY_TRE1=10;
      const int D_ENERGY_TRE2=11;
      const int D_ENERGY_TRE3=12;
      const int D_ENERGY_TRE4=13;
      const int D_ENERGY_TRE5=14;
      const int D_ENERGY_TRESUM=15;
      const int DD_POS_TRE=16;
      
      const int D_QDC1=20;
      const int D_QDC2=21;
      const int D_QDC3=22;
      const int D_QDC4=23;
      const int D_QDC5=24;
      const int D_QDCSUM=25;
      const int DD_POS_QDC=26;
      
      const int DD_POSNEW_CHE=30;
      const int DD_POSNEW_TRE=31;
      const int DD_POSNEW_QDC=32;
      
      const int DD_SINGLE_TRACE=77;
      
    }
} // namespace dammIds

using namespace std;
using namespace dammIds::pspmt;


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
    DeclareHistogram1D(D_RAW5, energyBins, "Pspmt Dynode");
    DeclareHistogram1D(D_SUM,  energyBins, "Pspmt Sum");
    DeclareHistogram2D(DD_POS_CHE, Bins, Bins, "Pspmt Pos1 Raw");
    
 
    // 710-Trace energies 
    DeclareHistogram1D(D_ENERGY_TRE1, energyBins, "Energy1 from trace");
    DeclareHistogram1D(D_ENERGY_TRE2, energyBins, "Energy2 from trace");
    DeclareHistogram1D(D_ENERGY_TRE3, energyBins, "Energy3 from trace");
    DeclareHistogram1D(D_ENERGY_TRE4, energyBins, "Energy4 from trace");
    DeclareHistogram1D(D_ENERGY_TRE5, energyBins, "Energy5 from trace");
    DeclareHistogram1D(D_ENERGY_TRESUM,  energyBins, "Pspmt Sum");
    DeclareHistogram2D(DD_POS_TRE, Bins, Bins, "Pspmt pos Raw by Trace1");
    
    // 1920- QDCs
    DeclareHistogram1D(D_QDC1, energyBins, "Energy1 from QDC");
    DeclareHistogram1D(D_QDC2, energyBins, "Energy2 from QDC");
    DeclareHistogram1D(D_QDC3, energyBins, "Energy3 from QDC");
    DeclareHistogram1D(D_QDC4, energyBins, "Energy4 from QDC");
    DeclareHistogram1D(D_QDC5, energyBins, "Energy5 from QDC");
    DeclareHistogram2D(DD_POS_QDC, posBins, posBins, "Pspmt pos Raw by QDC");

    // 1930~ for new board 
    DeclareHistogram2D(DD_POSNEW_CHE, Bins, Bins, "Position CHE Newboard");
    DeclareHistogram2D(DD_POSNEW_TRE, Bins, Bins, "Position TRE Newboard");
    DeclareHistogram2D(DD_POSNEW_QDC, Bins, Bins, "Position QDC Newboard");

    // Trace
    DeclareHistogram2D(DD_SINGLE_TRACE, traceBins, traceBins2,"Single trace");
    
}

PspmtProcessor::PspmtProcessor(void) : 
  EventProcessor(OFFSET, RANGE, "PspmtProcessor") {
  SetAssociatedTypes();
  ObtainHisName();
  SetupRootOutput();
  associatedTypes.insert("pspmt");
}
// Destructor to close output files and clean up pointers
PspmtProcessor::~PspmtProcessor(){
#ifdef useroot
  prootfile_->Write();
  prootfile_->Close();
  delete(prootfile_);
#endif
}

///Associates this Experiment Processor with template and ge detector types
void PspmtProcessor::SetAssociatedTypes(void) {
    associatedTypes.insert("pspmt");
}

#ifdef useroot
///Sets up ROOT output file, tree, branches, histograms.
void PspmtProcessor::SetupRootOutput(void) {
    stringstream rootname;
    rootname << fileName_ << ".root";
    prootfile_ = new TFile(rootname.str().c_str(),"RECREATE");
    proottree_ = new TTree("data","");
    proottree_->Branch("x_maxval",&x_maxval,"x_maxval/D");
    proottree_->Branch("y_maxval",&y_maxval,"y_maxval/D");
    proottree_->Branch("x_qdc", &x_qdc, "x_qdc/D");
    proottree_->Branch("y_qdc",&y_qdc,"y_qdc/D");
    proottree_->Branch("x_en", &x_en, "x_en/D");
    proottree_->Branch("y_en",&y_en,"y_en/D");
    proottree_->Branch("x_fen", &x_fen, "x_fen/D");
    proottree_->Branch("y_fen",&y_fen,"y_fen/D");
}
#endif

///Obtains the name of the histogram file passed via command line
void PspmtProcessor::ObtainHisName(void) {
    char hisFileName[32];
    GetArgument(1, hisFileName, 32);
    string temp = hisFileName;
    fileName_ = temp.substr(0, temp.find_first_of(" "));
}

bool PspmtProcessor::PreProcess(RawEvent &event){
    if (!EventProcessor::PreProcess(event))
        return false;
    
    //    static const vector<ChanEvent*> &pspmtEvents = sumMap["pspmt"]->GetList();
    static const vector<ChanEvent*> &pspmtEvents = event.GetSummary("pspmt:pspmt")->GetList();

    //    data_.Clear();
    
    if(pspmtEvents.size()>5){
      cerr << "We had too many pspmt events in the event list," << pspmtEvents.size() << endl;
      EndProcess();
      return false;
    }
    
    double che1=0,che2=0,che3=0,che4=0,che5=0;
    double tre1=0,tre2=0,tre3=0,tre4=0,tre5=0;
    double qdc1=0,qdc2=0,qdc3=0,qdc4=0,qdc5=0;
    
    double xche=0,yche=0;
    double xtre=0,ytre=0;
    double xqdc=0,yqdc=0;
    
    double xnche=0,ynche=0;
    double xntre=0,yntre=0;
    double xnqdc=0,ynqdc=0;
    
    static int traceNum;
    
    for (vector<ChanEvent*>::const_iterator it = pspmtEvents.begin();
         it != pspmtEvents.end(); it++) {
        
        ChanEvent *chan   = *it;
        string subtype    = chan->GetChanID().GetSubtype();
        int    ch         = chan->GetChanID().GetLocation();
        double calEnergy  = chan->GetCalEnergy();
        //double pspmtTime  = chan->GetTime();
        Trace trace       = chan->GetTrace();
        
	Trace trc  = (*it)->GetTrace();
	double qdc = trc.GetValue("tqdc");
	double en  = (*it)->GetEnergy();
	double max = trc.GetValue("maxval");
	double fen = trc.GetValue("filterEnergy0");
	
	double trace_energy;
        double trace_time;
        double baseline;
	

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
     

	if(true){
	  traceNum++;   	  
	  trace_time    = trace.GetValue("filterTime");
	  trace_energy  = trace.GetValue("filterEnergy");
	  baseline      = trace.DoBaseline(2,20);
	  // qdc             = trace.DoQDC(5,128);
	    
	
	
	   if(ch==0){
	      qdc1 = qdc;
	      tre1 = trace_energy;
	      plot(D_ENERGY_TRE1,tre1);
	      plot(D_QDC1,qdc1);
            }else if(ch==1){
	      qdc2 = qdc;
	      tre2 = trace_energy; 
	      plot(D_ENERGY_TRE2,tre2);
	      plot(D_QDC2,qdc2);
            }else if(ch==2){
	      qdc3 = qdc;
	      tre3 = trace_energy; 
	      plot(D_ENERGY_TRE3,tre3);
	      plot(D_QDC3,qdc3);
	    }else if(ch==3){
	      qdc4 = qdc;
	      tre4 = trace_energy; 	  
	      plot(D_ENERGY_TRE4,tre4);
	      plot(D_QDC4,qdc4);
	    }else if(ch==4){
	      qdc5 = qdc;
	      tre5 = trace_energy; 
	      plot(D_ENERGY_TRE5,tre5);
	      plot(D_QDC5,qdc5);
	    }
	   
	   if(che1>0 && che2>0 && che3>0 && che4>0){

	  xche=GetPositionX(che1,che2,che3,che4);
	  yche=GetPositionY(che1,che2,che3,che4);
	  xtre=GetPositionX(tre1,tre2,tre3,tre4);
	  ytre=GetPositionY(tre1,tre2,tre3,tre4);
	  xqdc=GetPositionX(qdc1,qdc2,qdc3,qdc4);
	  yqdc=GetPositionY(qdc1,qdc2,qdc3,qdc4);
	  
	  xnche=GetPositionXNew(che1,che2,che3,che4);
	  ynche=GetPositionYNew(che1,che2,che3,che4);
	  xntre=GetPositionXNew(tre1,tre2,tre3,tre4);
	  yntre=GetPositionYNew(tre1,tre2,tre3,tre4);
	  xnqdc=GetPositionXNew(qdc1,qdc2,qdc3,qdc4);
	  ynqdc=GetPositionYNew(qdc1,qdc2,qdc3,qdc4);
	  	  
	  
	  
	  
	  plot(DD_POS_CHE,xche,yche);
	  plot(DD_POS_TRE,xtre,ytre);
	  plot(DD_POS_QDC,xqdc,yqdc);

	  plot(DD_POSNEW_CHE,xnche,ynche);
	  plot(DD_POSNEW_TRE,xntre,yntre);
	  plot(DD_POSNEW_QDC,xnqdc,ynqdc);
	  
	  
        }
        
	   
	}

	
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

void PspmtProcessor::PspmtData::Clear(void) {    
}


double PspmtProcessor::GetPositionX(double q1,double q2,double q3,double q4){
  
  double qright=0,qleft=0,qsum=0;
  double xright=0,xleft=0;
  
  qsum = q1+q2+q3+q4;
  qright = (q4+q1)/2;
  qleft  = (q2+q3)/2;
  
  xright = (qright/qsum)*512+100;
  xleft  = (qleft/qsum)*512+100;
  
  return xright;
}

double PspmtProcessor::GetPositionY(double q1,double q2,double q3,double q4){
  
  double qtop=0,qbottom=0,qsum=0;
  double ytop=0,ybottom=0;
  
  qsum = q1+q2+q3+q4;
  qtop = (q1+q2)/2;
  qbottom  = (q3+q4)/2;
  
  ytop     = (qtop/qsum)*512+100;
  ybottom  = (qbottom/qsum)*512+100;
  
  return ytop;
}


double PspmtProcessor::GetPositionXNew(double q1,double q2,double q3,double q4){
  double xdiff=0,xsum=0,xposnew=0;
  
  xdiff = q1-q2;
  xsum  = q1+q2;
  xposnew  = 512*xdiff/xsum+512;
  
  return xposnew;
}

double PspmtProcessor::GetPositionYNew(double q1,double q2,double q3,double q4){
  double ydiff=0,ysum=0,yposnew=0;
  
  ydiff = q3-q4;
  ysum  = q3+q4;
  yposnew  = 512*ydiff/ysum+512;
  
  return yposnew;
}
