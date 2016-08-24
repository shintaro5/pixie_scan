/** \file PspmtProcessor.hpp
 *  \brief A processor to handle pixelated PMTs
 *  \author Shintaro Go
 *  \date August 24, 2016
 */

#ifndef __PSPMTPROCESSOR_HPP__
#define __PSPMTPROCESSOR_HPP__
#include "EventProcessor.hpp"

#ifdef useroot
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#endif

#include "RawEvent.hpp"


///Class to handle processing of position sensitive pmts
class PspmtProcessor : public EventProcessor {
public:
    /** Default Constructor */
    PspmtProcessor(void);
    /** Default Destructor */
    ~PspmtProcessor();
    
    /** Declare the plots used in the analysis */
    virtual void DeclarePlots(void);
    /** Preprocess the VANDLE data
     * \param [in] event : the event to preprocess
     * \return true if successful */
    virtual bool PreProcess(RawEvent &event);
    /** Process the event for VANDLE stuff
     * \param [in] event : the event to process
     * \return Returns true if the processing was successful */
    virtual bool Process(RawEvent &event);
private:
     /** Obtain the name of the histogram file */
    void ObtainHisName(void);
    /** Sets the detectors that are associated with this processor */
    void SetAssociatedTypes(void);
  std::string fileName_; //!< String to hold the file name from command line

#ifdef useroot
  /** Method to setup the ROOT output, tree and histograms */
  void SetupRootOutput(void);
  TFile *prootfile_; //! pointer to root file
  TTree *proottree_; //! pointer to root tree
#endif
  
  ///Structure defining what data we're storing
    struct PspmtData {
	///Clears the data from the processor 
        void Clear(void);
    } data_; //!< instance of structure to hold the data 
};
#endif // __PSPMTPROCESSOR_HPP__
