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
static double x_maxval;
static double y_maxval;

static double x_qdc;
static double y_qdc;

static double x_en;
static double y_en;

static double x_fen;
static double y_fen;
#endif

namespace dammIds {
    namespace experiment {
        const int D_TSIZE  = 0; //!< Size of the event
        const int D_GEENERGY  = 1; //!< Gamma energy 
        const int DD_QDCEN0  = 2; //!< Traces w/ Negative qdc values
        const int DD_POS  = 3; //!< Position given by the qdc values
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::experiment;

void PspmtTestProcessor::DeclarePlots(void) {
    DeclareHistogram1D(D_TSIZE, S3, "Num Template Evts");
    DeclareHistogram1D(D_GEENERGY, SA, "Gamma Energy with Cut");
    DeclareHistogram2D(DD_QDCEN0, SC, SC, "Traces w/ Neg QDC");
    DeclareHistogram2D(DD_POS, SB, SB, "Traces w/ Neg QDC");
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
            event.GetSummary("pspmt:anode")->GetList();

    ///Plot the size of the template events vector in two ways
    plot(D_TSIZE, pspmtEvts.size());

    if(pspmtEvts.size() > 4) {
        cerr << "We had too many pspmt events in the event list. Counted "
             << pspmtEvts.size() << endl;
        EndProcess();
        return false;
    }

    ///Begin loop over template events
    static int counter = 0;
    pair<double,double> px_maxval = make_pair(0,0);
    pair<double,double> py_maxval = make_pair(0,0);
    pair<double,double> px_qdc = make_pair(0,0);
    pair<double,double> py_qdc = make_pair(0,0);
    pair<double,double> px_en = make_pair(0,0);
    pair<double,double> py_en = make_pair(0,0);
    pair<double,double> px_fen = make_pair(0,0);
    pair<double,double> py_fen = make_pair(0,0);

    for(vector<ChanEvent*>::const_iterator tit = pspmtEvts.begin();
        tit != pspmtEvts.end(); ++tit) {

        Trace trc = (*tit)->GetTrace();
        double qdc = trc.GetValue("tqdc");
        double en = (*tit)->GetEnergy();
        double max = trc.GetValue("maxval");
        double fen = trc.GetValue("filterEnergy0");

        //Skip this channel if the qdc was less than 0
        if(qdc < 0)
            continue;

        if((*tit)->GetChanID().GetLocation() == 0) {
            plot(DD_QDCEN0, qdc, en);
        }

        //Set the value for the x and y pairs depending on the tag that we have in the channel.
        if((*tit)->GetChanID().HasTag("x1")) {
            px_maxval.first = max;
            px_qdc.first = qdc;
            px_en.first = en;
            px_fen.first = fen;
        }
        if((*tit)->GetChanID().HasTag("x2")) {
            px_maxval.second = max;
            px_qdc.second = qdc;
            px_en.second = en;
            px_fen.second = fen;
        }
        if((*tit)->GetChanID().HasTag("y1")) {
            py_maxval.first = max;
            py_qdc.first = qdc;
            py_en.first = en;
            py_fen.first = fen;
        }
        if((*tit)->GetChanID().HasTag("y2")) {
            py_maxval.second = max;
            py_qdc.second = qdc;
            py_en.second = en;
            py_fen.second = fen;
        }

        //Output an error message if we found a trace with a NEGATIVE qdc
        if(trc.GetValue("tqdc") < 0 && counter < 1024) {
            cerr << "PspmtTestProcessor::Process - IDENTIFIED A NEGATIVE QDC : " << endl;
            cerr << trc.GetValue("maxpos") << " " << trc.GetValue("tqdc")
                 << " " << trc.GetValue("maxval") << " " << counter << endl;
            for(unsigned int i = 0; i < trc.size(); i++)
                //plot(DD_NEGQDC,i,counter, trc[i]);
            counter++;
        }
    }

    //Calculate the position using QDC
    x_qdc = (px_qdc.first - px_qdc.second) / (px_qdc.first + px_qdc.second);
    y_qdc = (py_qdc.first - py_qdc.second) / (py_qdc.first + py_qdc.second);
    //Calculate the position using the Energy
    x_en = (px_en.first - px_en.second) / (px_en.first + px_en.second);
    y_en = (py_en.first - py_en.second) / (py_en.first + py_en.second);
    //Calculate the position using the max value
    x_maxval = (px_maxval.first - px_maxval.second) / (px_maxval.first + px_maxval.second);
    y_maxval = (py_maxval.first - py_maxval.second) / (py_maxval.first + py_maxval.second);
    //Calculate the position using the trace filter energy
    x_fen = (px_fen.first - px_fen.second) / (px_fen.first + px_fen.second);
    y_fen = (py_fen.first - py_fen.second) / (py_fen.first + py_fen.second);

    proottree_->Fill();

    EndProcess();
    return(true);
}
