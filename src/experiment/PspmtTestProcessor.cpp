/** \file PspmtTestProcessor.cpp
 * \brief Example class for experiment specific setups
 *\author S. V. Paulauskas
 *\date August 19, 2016
 */
#include <iostream>

#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "GeProcessor.hpp"
#include "GetArguments.hpp"
#include "PspmtTestProcessor.hpp"

#ifdef useroot
static double tof_;
static double tEnergy;
#endif

namespace dammIds {
    namespace experiment {
        const int D_TSIZE  = 0; //!< Size of the event
        const int D_GEENERGY  = 1; //!< Gamma energy 
        const int DD_TENVSGEN  = 2; //!< Energy vs Gamma Energy 
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::experiment;

void PspmtTestProcessor::DeclarePlots(void) {
    DeclareHistogram1D(D_TSIZE, S3, "Num Template Evts");
    DeclareHistogram1D(D_GEENERGY, SA, "Gamma Energy with Cut");
    DeclareHistogram2D(DD_TENVSGEN, SA, SA, "Template En vs. Ge En");
}

PspmtTestProcessor::PspmtTestProcessor() :
    EventProcessor(OFFSET, RANGE, "PspmtTestProcessor") {
    SetAssociatedTypes();
    ObtainHisName();
    SetupRootOutput();
}

///Destructor to close output files and clean up pointers
PspmtTestProcessor::~PspmtTestProcessor() {
#ifdef useroot
    prootfile_->Write();
    prootfile_->Close();
    delete(prootfile_);
#endif
}

///Associates this Experiment Processor with template and ge detector types
void PspmtTestProcessor::SetAssociatedTypes(void) {
    associatedTypes.insert("pspmt");
}

#ifdef useroot
///Sets up ROOT output file, tree, branches, histograms.
void PspmtTestProcessor::SetupRootOutput(void) {
    stringstream rootname;
    rootname << fileName_ << ".root";
    prootfile_ = new TFile(rootname.str().c_str(),"RECREATE");
    proottree_ = new TTree("data","");
    proottree_->Branch("tof",&tof_,"tof/D");
    proottree_->Branch("ten",&tEnergy,"ten/D");
    ptvsge_ = new TH2D("tvsge","",1000,-100,900,16000,0,16000);
    ptsize_ = new TH1D("tsize","",40,0,40);
}
#endif

///Obtains the name of the histogram file passed via command line
void PspmtTestProcessor::ObtainHisName(void) {
    char hisFileName[32];
    GetArgument(1, hisFileName, 32);
    string temp = hisFileName;
    fileName_ = temp.substr(0, temp.find_first_of(" "));
}

///We do nothing here since we're completely dependent on the resutls of others
bool PspmtTestProcessor::PreProcess(RawEvent &event){
    if (!EventProcessor::PreProcess(event))
        return(false);
    return(true);
}

///Main processing of data of interest
bool PspmtTestProcessor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return(false);

    static const vector<ChanEvent*> &pspmtEvts =
            event.GetSummary("pspmt:pspmt")->GetList();

    ///Plot the size of the template events vector in two ways
    plot(D_TSIZE, pspmtEvts.size());
#ifdef useroot
    ptsize_->Fill(pspmtEvts.size());
#endif

    ///Begin loop over template events
    static int counter = 0;
    for(vector<ChanEvent*>::const_iterator tit = pspmtEvts.begin();
        tit != pspmtEvts.end(); ++tit) {
        Trace trc = (*tit)->GetTrace();
        if(trc.GetValue("tqdc") < 0 && counter < 1024) {
            cout << "PspmtTestProcessor::Process - IDENTIFIED A NEGATIVE QDC : " << endl;
            cout << trc.GetValue("maxpos") << " " << trc.GetValue("tqdc")
                 << " " << trc.GetValue("maxval") << " " << counter << endl;

            for(unsigned int i = 0; i < trc.size(); i++)
                plot(DD_TENVSGEN,i,counter, trc[i]);
            counter++;
        }

    }
    EndProcess();
    return(true);
}
