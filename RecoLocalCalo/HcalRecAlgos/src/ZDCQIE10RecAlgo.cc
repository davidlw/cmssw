#include <algorithm>
#include <iostream>

#include "RecoLocalCalo/HcalRecAlgos/interface/ZDCQIE10RecAlgo.h"

#include "DataFormats/HcalDigi/interface/QIE10DataFrame.h"
#include "DataFormats/HcalRecHit/interface/HcalSpecialTimes.h"

#include "CalibFormats/HcalObjects/interface/HcalCoder.h"
#include "CalibFormats/CaloObjects/interface/CaloSamples.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"

ZDCQIE10Info ZDCQIE10RecAlgo::reconstruct(const QIE10DataFrame& digi,
                                      const int tsToStart,
                                      const double satCorrFactor,  
                                      const HcalCoder& coder,
                                      const HcalCalibrations& calib) const
{
    // Scrap the trailing edge time for now -- until the front-end
    // FPGA firmware is finalized and the description of the FPGA
    // output becomes available
    static const float timeFalling = HcalSpecialTimes::UNKNOWN_T_NOTDC;

    ZDCQIE10Info result;

    CaloSamples cs;
    coder.adc2fC(digi, cs);
//    const int nRead = cs.size();

    // This branch is intended for use with cosmic runs
    double charge = 0.0, energy = 0.0;
    ZDCQIE10Info::raw_type raw[ZDCQIE10Info::N_RAW_MAX];

    // check saturation
    double correctSaturation[2] = {1.0,1.0};
    if(cs[tsToStart] == 127) correctSaturation[0] = satCorrFactor;
/*
    for (unsigned int ts=0; ts<ZDCQIE10Info::N_RAW_MAX; ++ts)
    {
      const QIE10DataFrame::Sample s(digi[ts]);
      const int capid = s.capid();
std::cout<<ts<<" "<<cs[ts]<<" "<<capid<<" "<<calib.respcorrgain(capid)<<" "<<calib.pedestal(capid)<<std::endl;
    }
*/
    for (int ts=tsToStart-1; ts<tsToStart+1; ++ts)
    {
      const QIE10DataFrame::Sample s(digi[ts]);
      const int capid = s.capid();
      const float q = cs[ts]*correctSaturation[ts-tsToStart+1] - calib.pedestal(capid);
      charge += q;
      energy += q*calib.respcorrgain(capid);

      raw[ts] = s.wideRaw();
    }

    // Timing measurement does not appear to be useful here
    const float timeRising = HcalSpecialTimes::UNKNOWN_T_NOTDC;

    // The following ZDCQIE10Info arguments correspond to SOI
    // not stored in the raw data. Essentially, only charge
    // and energy are meaningful.
    result = ZDCQIE10Info(digi.id(), charge, energy,
                         timeRising, timeFalling,
                         raw);
    return result;
}
