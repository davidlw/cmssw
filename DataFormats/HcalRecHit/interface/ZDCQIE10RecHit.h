#ifndef DATAFORMATS_HCALRECHIT_ZDCQIE10RECHIT_H
#define DATAFORMATS_HCALRECHIT_ZDCQIE10RECHIT_H

#include <utility>

#include "DataFormats/HcalRecHit/interface/ZDCQIE10Info.h"

/** \class ZDCQIE10RecHit
*
* Class to contain the info needed to perform ZDC reconstruction
* using QIE10 chips. Groups the information
* provided by a single PMT.
*/
class ZDCQIE10RecHit
{
public:
    typedef HcalDetId key_type;

    ZDCQIE10RecHit();

    // pointer can be nullptr
    ZDCQIE10RecHit(const HcalDetId& id, const ZDCQIE10Info* first);

    inline HcalDetId id() const {return id_;}

    // Get a pointer to the QIE10 info. nullptr will be returned
    // if the info with the given index does not exist or if the
    // index is out of range.
    const ZDCQIE10Info* getZDCQIE10Info(unsigned index) const;

    // Quantities simply added from both anodes
    float charge() const;
    float energy() const;

private:
    HcalDetId id_;

    ZDCQIE10Info zdcQIE10Info_;
    bool hasInfo_;
};

#endif // DATAFORMATS_HCALRECHIT_ZDCQIE10RECHIT_H
