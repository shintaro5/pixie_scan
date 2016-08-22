/** \file PspmtTestProcessor.hpp
 * \brief Example class for experiment specific setups
 *\author S. V. Paulauskas
 *\date August 19, 2016
 */
#ifndef __PSPMTTESTPROCESSOR_HPP_
#define __PSPMTTESTPROCESSOR_HPP_
#include "EventProcessor.hpp"

#ifdef useroot
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#endif

/// Working template class for experiment processors
class PspmtTestProcessor : public EventProcessor {
public:
    /** Default Constructor */
    PspmtTestProcessor();
    /** Default Destructor */
    ~PspmtTestProcessor();
    /** Declare the plots used in the analysis */
    virtual void DeclarePlots(void);

    /** PreProcess does nothing since this is solely dependent on results
     from other Processors*/
    virtual bool PreProcess(RawEvent &event);

    /** Process the event
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
};
#endif
