#ifndef RecoLocalCalo_HcalRecAlgos_ZDCQIE10RecAlgo_h_
#define RecoLocalCalo_HcalRecAlgos_ZDCQIE10RecAlgo_h_

#include "DataFormats/HcalRecHit/interface/HFQIE10Info.h"

class QIE10DataFrame;
class HcalCoder;
class HcalCalibrations;

class ZDCQIE10RecAlgo
{
public:
    inline explicit ZDCQIE10RecAlgo(const bool sumAllTS) : sumAllTS_(sumAllTS) {}

    inline ~ZDCQIE10RecAlgo() {}

    HFQIE10Info reconstruct(const QIE10DataFrame& digi,
                            int tsToUse,
                            const HcalCoder& coder,
                            const HcalCalibrations& calibs) const;
private:
    bool sumAllTS_;
};

#endif // RecoLocalCalo_HcalRecAlgos_ZDCQIE10RecAlgo_h_
