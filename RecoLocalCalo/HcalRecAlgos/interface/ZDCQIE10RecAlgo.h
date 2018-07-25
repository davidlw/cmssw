#ifndef RecoLocalCalo_HcalRecAlgos_ZDCQIE10RecAlgo_h_
#define RecoLocalCalo_HcalRecAlgos_ZDCQIE10RecAlgo_h_

#include "DataFormats/HcalRecHit/interface/ZDCQIE10Info.h"

class QIE10DataFrame;
class HcalCoder;
class HcalCalibrations;

class ZDCQIE10RecAlgo
{
public:
    inline explicit ZDCQIE10RecAlgo() {}

    inline ~ZDCQIE10RecAlgo() {}

    ZDCQIE10Info reconstruct(const QIE10DataFrame& digi,
                            int tsToStart,
                            double satCorrFactor,
                            const HcalCoder& coder,
                            const HcalCalibrations& calibs) const;
private:
};

#endif // RecoLocalCalo_HcalRecAlgos_ZDCQIE10RecAlgo_h_
