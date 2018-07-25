#include "DataFormats/HcalRecHit/interface/ZDCQIE10RecHit.h"

ZDCQIE10RecHit::ZDCQIE10RecHit()
    : hasInfo_(false)
{
}

ZDCQIE10RecHit::ZDCQIE10RecHit(const HcalDetId& id, const ZDCQIE10Info* first)
    : id_(id), hasInfo_(false)
{
    if (first)
    {
        zdcQIE10Info_ = *first;
        hasInfo_ = true;
    }
}

float ZDCQIE10RecHit::charge() const
{
    float q = 0.f;
    if (hasInfo_)
       q = zdcQIE10Info_.charge();
    return q;
}

float ZDCQIE10RecHit::energy() const
{
    float e = 0.f;
    if (hasInfo_)
       e += zdcQIE10Info_.energy();
    return e;
}

const ZDCQIE10Info* ZDCQIE10RecHit::getZDCQIE10Info(const unsigned index) const
{
    if (hasInfo_)
        return &zdcQIE10Info_;
    else
        return nullptr;
}
