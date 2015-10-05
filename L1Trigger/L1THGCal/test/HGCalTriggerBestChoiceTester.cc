#include <iostream>
#include <string>
#include <vector>


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/L1THGCal/interface/HGCFETriggerDigi.h"
#include "DataFormats/L1THGCal/interface/HGCFETriggerDigiFwd.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCTriggerDetId.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCTriggerDetId.h"

#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerFECodecBase.h"
#include "L1Trigger/L1THGCal/interface/fe_codecs/HGCalBestChoiceCodec.h"

#include <stdlib.h> 


class HGCalTriggerBestChoiceTester : public edm::EDAnalyzer 
{
    public:
        explicit HGCalTriggerBestChoiceTester(const edm::ParameterSet& );
        ~HGCalTriggerBestChoiceTester();

        virtual void beginRun(const edm::Run&, const edm::EventSetup&);
        virtual void analyze(const edm::Event&, const edm::EventSetup&);


    private:
        void fill(const l1t::HGCFETriggerDigiCollection&, const HGCEEDigiCollection&);
        // inputs
        edm::EDGetToken inputee_, inputfh_, inputbh_;
        //
        std::unique_ptr<HGCalTriggerGeometryBase> triggerGeometry_; 
        //std::unique_ptr<HGCalTriggerFECodecBase> codec_;
        HGCalBestChoiceCodec codec_;
        edm::Service<TFileService> fs_;
        // histos
        TH1F* triggerCellsPerModule_;
        TH1F* triggerCellData_;
        TH1F* moduleSum_;

};


/*****************************************************************/
HGCalTriggerBestChoiceTester::HGCalTriggerBestChoiceTester(const edm::ParameterSet& conf):
  inputee_(consumes<HGCEEDigiCollection>(conf.getParameter<edm::InputTag>("eeDigis"))),
  inputfh_(consumes<HGCHEDigiCollection>(conf.getParameter<edm::InputTag>("fhDigis"))), 
  inputbh_(consumes<HGCHEDigiCollection>(conf.getParameter<edm::InputTag>("bhDigis"))),
  codec_(conf.getParameterSet("FECodec"))
/*****************************************************************/
{
    //setup geometry 
    const edm::ParameterSet& geometryConfig = conf.getParameterSet("TriggerGeometry");
    const std::string& trigGeomName = geometryConfig.getParameter<std::string>("TriggerGeometryName");
    HGCalTriggerGeometryBase* geometry = HGCalTriggerGeometryFactory::get()->create(trigGeomName,geometryConfig);
    triggerGeometry_.reset(geometry);

    //setup FE codec
    //const edm::ParameterSet& feCodecConfig = conf.getParameterSet("FECodec");
    //const std::string& feCodecName =feCodecConfig.getParameter<std::string>("CodecName");
    //HGCalTriggerFECodecBase* codec = HGCalTriggerFECodecFactory::get()->create(feCodecName,feCodecConfig);
    //codec_.reset(codec);
    codec_.unSetDataPayload();

    // initialize output trees
    triggerCellsPerModule_ = fs_->make<TH1F>("TriggerCellsPerModule","Number of trigger cells per module", 64, 0., 64.);
    triggerCellData_       = fs_->make<TH1F>("TriggerCellData","Trigger cell values", 500, 0., 500.);
    moduleSum_             = fs_->make<TH1F>("ModuleSum","Trigger cell sum in modules", 5000, 0., 5000.);
}



/*****************************************************************/
HGCalTriggerBestChoiceTester::~HGCalTriggerBestChoiceTester() 
/*****************************************************************/
{
}

/*****************************************************************/
void HGCalTriggerBestChoiceTester::beginRun(const edm::Run& /*run*/, 
                                          const edm::EventSetup& es)
/*****************************************************************/
{
    triggerGeometry_->reset();
    HGCalTriggerGeometryBase::es_info info;
    const std::string& ee_sd_name = triggerGeometry_->eeSDName();
    const std::string& fh_sd_name = triggerGeometry_->fhSDName();
    const std::string& bh_sd_name = triggerGeometry_->bhSDName();
    es.get<IdealGeometryRecord>().get(ee_sd_name,info.geom_ee);
    es.get<IdealGeometryRecord>().get(fh_sd_name,info.geom_fh);
    es.get<IdealGeometryRecord>().get(bh_sd_name,info.geom_bh);
    es.get<IdealGeometryRecord>().get(ee_sd_name,info.topo_ee);
    es.get<IdealGeometryRecord>().get(fh_sd_name,info.topo_fh);
    es.get<IdealGeometryRecord>().get(bh_sd_name,info.topo_bh);
    triggerGeometry_->initialize(info);
}



/*****************************************************************/
void HGCalTriggerBestChoiceTester::analyze(const edm::Event& e, 
			      const edm::EventSetup& es) 
/*****************************************************************/
{
    std::unique_ptr<l1t::HGCFETriggerDigiCollection> fe_coll_ptr( new l1t::HGCFETriggerDigiCollection );

    edm::Handle<HGCEEDigiCollection> ee_digis_h;
    edm::Handle<HGCHEDigiCollection> fh_digis_h, bh_digis_h;

    e.getByToken(inputee_,ee_digis_h);
    e.getByToken(inputfh_,fh_digis_h);
    e.getByToken(inputbh_,bh_digis_h);

    const HGCEEDigiCollection& ee_digis = *ee_digis_h;
    const HGCHEDigiCollection& fh_digis = *fh_digis_h;
    const HGCHEDigiCollection& bh_digis = *bh_digis_h;

    for( const auto& module : triggerGeometry_->modules() ) 
    {    
        fe_coll_ptr->push_back(l1t::HGCFETriggerDigi());
        l1t::HGCFETriggerDigi& digi = fe_coll_ptr->back();
        codec_.setDataPayload(*(module.second),ee_digis,fh_digis,bh_digis);
        codec_.encode(digi);
        digi.setDetId( HGCTriggerDetId(module.first) );
        codec_.unSetDataPayload();
    }

}


/*****************************************************************/
void HGCalTriggerBestChoiceTester::fill(const l1t::HGCFETriggerDigiCollection& fe_digis, const HGCEEDigiCollection& ee_digis)
/*****************************************************************/
{
    for( const auto& module : triggerGeometry_->modules() ) 
    { 
        // Trigger cells
        unsigned nFEDigi = 0;
        unsigned moduleSum = 0;
        for(const auto& fe_digi : fe_digis)
        {
            if(fe_digi.getDetId<HGCTriggerDetId>()==module.first)
            {
                HGCalBestChoiceCodec::data_type data;
                data.reset();
                fe_digi.decode<HGCalBestChoiceCodec, HGCalBestChoiceCodec::data_type>(codec_, data);
                for(const auto& tc : data.payload)
                {
                    if(tc>0)
                    {
                        nFEDigi++;
                        moduleSum += tc;
                        triggerCellData_->Fill(tc);
                    }
                }
            }
        }
        triggerCellsPerModule_->Fill(nFEDigi);
        moduleSum_->Fill(moduleSum);
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCalTriggerBestChoiceTester);
