/** \file ScintProcessor.cpp
 *
 * implementation for scintillator processor
 */
#include <iostream>
#include <cmath>

#include "DammPlotIds.hpp"
#include "DetectorDriver.hpp"
#include "RawEvent.hpp"
#include "Hen3Processor.hpp"
#include "ScintProcessor.hpp"
#include "Trace.hpp"

using namespace std;
using namespace dammIds::hen3;

namespace dammIds {
    namespace hen3 {
        const int D_MULT_HEN3 = 0;
        const int D_MULT_NEUTRON = 1;
        const int D_ENERGY_HEN3 = 2;
        const int D_ENERGY_NEUTRON = 3;
        const int D_TDIFF_HEN3_BETA = 4;
        const int D_TDIFF_NEUTRON_BETA = 5;

        const int DD_DISTR_HEN3 = 10;
        const int DD_DISTR_NEUTRON = 11;
    }
}

Hen3Processor::Hen3Processor() : EventProcessor(OFFSET, RANGE)
{
    name = "3hen";
    associatedTypes.insert("3hen"); // associate with the scint type
}

void Hen3Processor::DeclarePlots(void)
{
    DeclareHistogram1D(D_MULT_HEN3, S6, "3Hen event multiplicity");
    DeclareHistogram1D(D_MULT_NEUTRON, S6, "3Hen real neutron multiplicity");

    DeclareHistogram1D(D_ENERGY_HEN3, SE, "3Hen raw energy");
    DeclareHistogram1D(D_ENERGY_NEUTRON, SE, "Neutron raw energy");

    DeclareHistogram1D(D_TDIFF_HEN3_BETA, SB, "time diff hen3 - beta + 1000 (10 ns/ch)");
    DeclareHistogram1D(D_TDIFF_NEUTRON_BETA, SB, "time diff neutron - beta + 1000 (10 ns/ch)");

    DeclareHistogram2D(DD_DISTR_HEN3, S5, S5, "3Hen event distribution");
    DeclareHistogram2D(DD_DISTR_NEUTRON, S5, S5, "3Hen real neutron distribution");
}

bool Hen3Processor::PreProcess(RawEvent &event)
{
    if (!EventProcessor::PreProcess(event))
        return false;

    /** Notice that places Hen3 and Neutron (counters) are activated as parents
     * of detectors activated here. The thresholds are set in DetectorDriver.cpp
     * in constructors.
     */
    static const DetectorSummary *hen3Summary = event.GetSummary("3hen");
    for (vector<ChanEvent*>::const_iterator it = hen3Summary->GetList().begin(); 
        it != hen3Summary->GetList().end(); it++) {
            double time = (*it)->GetTime();
            double energy = (*it)->GetEnergy();
            int location = (*it)->GetChanID().GetLocation();
            string place = (*it)->GetChanID().GetPlaceName();
            CorrEventData data(time, true, energy);
            TreeCorrelator::get().places[place]->activate(data);
            stringstream neutron;
            neutron << "Neutron" << location;
            TreeCorrelator::get().places[neutron.str()]->activate(data);
    }

    return true;
}

bool Hen3Processor::Process(RawEvent &event)
{
    if (!EventProcessor::Process(event))
        return false;
  
    static const DetectorSummary *hen3Summary = event.GetSummary("3hen");
    
    int hen3_count = dynamic_cast<PlaceCounter*>(TreeCorrelator::get().places["Hen3"])->getCounter();
    int neutron_count = dynamic_cast<PlaceCounter*>(TreeCorrelator::get().places["Neutrons"])->getCounter();

    plot(D_MULT_HEN3, hen3_count);
    plot(D_MULT_NEUTRON, neutron_count);

    for (vector<ChanEvent*>::const_iterator it = hen3Summary->GetList().begin(); 
        it != hen3Summary->GetList().end(); it++) {
            ChanEvent *chan = *it;
            int location = chan->GetChanID().GetLocation();

            double energy = chan->GetEnergy();
            plot(D_ENERGY_HEN3, energy);

            stringstream neutron_name;
            neutron_name << "Neutron" << location;
            if (TreeCorrelator::get().places[neutron_name.str()]->status())
                plot(D_ENERGY_NEUTRON, energy);

            string place = chan->GetChanID().GetPlaceName();
            if (TreeCorrelator::get().places["Beta"]->status()) {
                double hen3_time = chan->GetTime();
                double beta_time = TreeCorrelator::get().places["Beta"]->last().time;
                const double timeResolution = 10e-9 / pixie::clockInSeconds;
                double dt = int(hen3_time - beta_time + 1000);
                if (dt > SB)
                    dt = SB - 1;
                if (dt < 0)
                    dt = 0;
                if (TreeCorrelator::get().places[neutron_name.str()]->status())
                    plot(D_TDIFF_NEUTRON_BETA, dt / timeResolution);
                plot(D_TDIFF_HEN3_BETA, dt / timeResolution);

            }

            /** These plots show He3 bar location hit
             * picture is as looking along the beam line
             * (having ORIC behind, 3Hen in front)
             */
            int xpos = 0;
            int ypos = 0;
            if( 0 <= location && location <= 7 ) {
                ypos = 11;
                xpos = 15 - 2 * location;
            } else if( 8 <= location && location <= 14 ) {
                ypos = 13;
                xpos = 30 - 2 * location;
            } else if( 15 <= location && location <= 20 ) {
                ypos = 15;
                xpos = 43 - 2 * location;
            } else if( 21 <= location && location <= 23 ) {
                ypos = 17;
                xpos = 52 - 2 * location;
            } else if( 24 <= location && location <= 31 ) {
                ypos = 7;
                xpos = 63 - 2 * location;
            } else if( 32 <= location && location <= 38 ) {
                ypos = 5;
                xpos = 78 - 2 * location;
            } else if( 39 <= location && location <= 44 ) {
                ypos = 3;
                xpos = 91 - 2 * location;
            } else if( 45 <= location && location <= 47 ) {
                ypos = 1;
                xpos = 100 - 2 * location;
            }

            plot(DD_DISTR_HEN3, xpos, ypos);
            if (TreeCorrelator::get().places[neutron_name.str()]->status())
                plot(DD_DISTR_NEUTRON, xpos, ypos);
    }

    EndProcess();
    return true;
}
