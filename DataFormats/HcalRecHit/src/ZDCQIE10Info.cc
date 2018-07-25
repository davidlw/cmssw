#include <algorithm>
#include "DataFormats/HcalRecHit/interface/ZDCQIE10Info.h"

const unsigned ZDCQIE10Info::N_RAW_MAX;
const ZDCQIE10Info::raw_type ZDCQIE10Info::INVALID_RAW;
 

ZDCQIE10Info::ZDCQIE10Info()
    : charge_(0.f),
      energy_(0.f),
      timeRising_(0.f),
      timeFalling_(-1.f),
      raw_{INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW}
{}

ZDCQIE10Info::ZDCQIE10Info(const HcalDetId& id,
                         const float i_charge, const float i_energy,
                         const float i_timeRising, const float i_timeFalling,
                         const raw_type* rawData)
    : id_(id),
      charge_(i_charge),
      energy_(i_energy),
      timeRising_(i_timeRising),
      timeFalling_(i_timeFalling),
      raw_{INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW, INVALID_RAW}
{}

bool ZDCQIE10Info::isDataframeOK() const
{
    bool hardwareOK = true;

    for (unsigned i=0; i<N_RAW_MAX && hardwareOK; ++i)
    {
        const QIE10DataFrame::Sample s(raw_[i]);
        hardwareOK = s.ok();
    }

    return hardwareOK;
}
