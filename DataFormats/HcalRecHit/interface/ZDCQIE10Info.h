#ifndef DATAFORMATS_HCALRECHIT_ZDCQIE10INFO_H
#define DATAFORMATS_HCALRECHIT_ZDCQIE10INFO_H

#include <limits>

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/QIE10DataFrame.h"

/** \class ZDCQIE10Info
*
* Class to contain the info needed to perform ZDC reconstruction
* using QIE10 chips. Intended for use
* inside ZDCQIE10RecHit.
*/
class ZDCQIE10Info
{
public:
    typedef HcalDetId key_type;
    typedef QIE10DataFrame::Sample::wide_type raw_type;

    static const unsigned N_RAW_MAX = 10;
    static const raw_type INVALID_RAW = std::numeric_limits<raw_type>::max();

    ZDCQIE10Info();

    ZDCQIE10Info(const HcalDetId& id, float charge, float energy,
                float timeRising, float timeFalling,
                const raw_type* rawData);

    inline HcalDetId id() const {return id_;}

    inline float charge() const {return charge_;}
    inline float energy() const {return energy_;}
    inline float timeRising() const {return timeRising_;}
    inline float timeFalling() const {return timeFalling_;}

    bool isDataframeOK() const;

private:
    HcalDetId id_;

    float charge_;
    float energy_;
    float timeRising_;
    float timeFalling_;
    raw_type raw_[N_RAW_MAX];
};

#endif // DATAFORMATS_HCALRECHIT_ZDCQIE10INFO_H
