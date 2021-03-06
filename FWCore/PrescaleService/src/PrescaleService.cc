////////////////////////////////////////////////////////////////////////////////
//
// PrescaleService
// ---------------
//
//            04/25/2008 Philipp Schieferdecker <philipp.schieferdecker@cern.ch>
////////////////////////////////////////////////////////////////////////////////


#include "FWCore/PrescaleService/interface/PrescaleService.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <set>
#include <algorithm>


namespace edm {
  namespace service {

    ////////////////////////////////////////////////////////////////////////////////
    // construction/destruction
    ////////////////////////////////////////////////////////////////////////////////

    //______________________________________________________________________________
    PrescaleService::PrescaleService(ParameterSet const& iPS,ActivityRegistry& iReg)
      : configured_(false)
      , forceDefault_(iPS.getParameter<bool>("forceDefault"))
      , lvl1Labels_(iPS.getParameter<std::vector<std::string> >("lvl1Labels"))
      , nLvl1Index_(lvl1Labels_.size())
      , iLvl1IndexDefault_(findDefaultIndex(iPS.getParameter<std::string>("lvl1DefaultLabel"), lvl1Labels_))
      , vpsetPrescales_(iPS.getParameterSetVector("prescaleTable"))
      , prescaleTable_()
    {
      iReg.watchPostBeginJob(this, &PrescaleService::postBeginJob);
    }
      
    //______________________________________________________________________________
    PrescaleService::~PrescaleService() {
    }

    ////////////////////////////////////////////////////////////////////////////////
    // implementation of member functions
    ////////////////////////////////////////////////////////////////////////////////

    void PrescaleService::reconfigure(ParameterSet const& iPS) {
      vpsetPrescales_.clear();
      prescaleTable_.clear();
      lvl1Labels_ = iPS.getParameter<std::vector<std::string> >("lvl1Labels");
      nLvl1Index_ = lvl1Labels_.size();
      iLvl1IndexDefault_ = findDefaultIndex(iPS.getParameter<std::string>("lvl1DefaultLabel"), lvl1Labels_);
      vpsetPrescales_ = iPS.getParameterSetVector("prescaleTable");
      configure();
    }

    void PrescaleService::postBeginJob() {
      if (!configured_) {
        configure();
      }
    }

    //______________________________________________________________________________
    void PrescaleService::configure()
    {
      configured_ = true;

      ParameterSet prcPS = getProcessParameterSet();
      
      // find all HLTPrescaler modules
      std::set<std::string> prescalerModules;
      std::vector<std::string> allModules=prcPS.getParameter<std::vector<std::string> >("@all_modules");
      for(unsigned int i = 0; i < allModules.size(); ++i) {
        ParameterSet const& pset  = prcPS.getParameterSet(allModules[i]);
        std::string moduleLabel = pset.getParameter<std::string>("@module_label");
        std::string moduleType  = pset.getParameter<std::string>("@module_type");
        if (moduleType == "HLTPrescaler") prescalerModules.insert(moduleLabel);
      }
      
      // find all paths with an HLTPrescaler and check for <=1
      std::set<std::string> prescaledPathSet;
      std::vector<std::string> allPaths = prcPS.getParameter<std::vector<std::string> >("@paths");
      for (unsigned int iP = 0; iP < allPaths.size(); ++iP) {
        std::string pathName = allPaths[iP];
        std::vector<std::string> modules = prcPS.getParameter<std::vector<std::string> >(pathName);
        for (unsigned int iM = 0; iM < modules.size(); ++iM) {
          std::string moduleLabel = modules[iM];
          if (prescalerModules.erase(moduleLabel)>0) {
            std::set<std::string>::const_iterator itPath=prescaledPathSet.find(pathName);
            if (itPath==prescaledPathSet.end()) {
              prescaledPathSet.insert(pathName);
            } else {
              throw cms::Exception("DuplicatePrescaler")
                <<"path '"<<pathName<<"' has more than one HLTPrescaler!";
            }
          }
        }
      }

      std::vector<std::string> prescaledPaths;
      for (unsigned int iVPSet=0; iVPSet < vpsetPrescales_.size(); ++iVPSet) {
        ParameterSet psetPrescales = vpsetPrescales_[iVPSet];
        std::string pathName = psetPrescales.getParameter<std::string>("pathName");
        if (prescaledPathSet.erase(pathName) > 0) {
          std::vector<unsigned int> prescales =
            psetPrescales.getParameter<std::vector<unsigned int> >("prescales");
          if (prescales.size()!=nLvl1Index_) {
            throw cms::Exception("PrescaleTableMismatch")
              << "path '" << pathName << "' has unexpected number of prescales";
          }
          prescaleTable_[pathName] = prescales;
        }
        else {
          throw cms::Exception("PrescaleTableUnknownPath")
            <<"path '"<<pathName<<"' is invalid or does not "
            <<"contain any HLTPrescaler";
        }
      }      
    }

    //______________________________________________________________________________
    unsigned int PrescaleService::getPrescale(std::string const& prescaledPath)
    {
      return getPrescale(iLvl1IndexDefault_, prescaledPath);
    }
    
    //______________________________________________________________________________
    unsigned int PrescaleService::getPrescale(unsigned int lvl1Index,
                                              std::string const& prescaledPath)
    {
      if (forceDefault_)
        lvl1Index = iLvl1IndexDefault_;

      if (lvl1Index >= nLvl1Index_) {
        throw cms::Exception("InvalidLvl1Index")
          <<"lvl1Index '"<<lvl1Index<<"' exceeds number of prescale columns";
      }

      if (!configured_) {
        configure();
      }
      
      PrescaleTable_t::const_iterator it = prescaleTable_.find(prescaledPath);
      return (it == prescaleTable_.end()) ? 1 : it->second[lvl1Index];
    }
    
    //______________________________________________________________________________
    unsigned int PrescaleService::findDefaultIndex(std::string const & label, std::vector<std::string> const & labels) {
      for (unsigned int i = 0; i < labels.size(); ++i) {
        if (labels[i] == label) {
          return i;
        }
      }
      // FIXME add a LogWarning if the default is not found ?
      return 0;
    }
    
    //______________________________________________________________________________
    void PrescaleService::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
      edm::ParameterSetDescription desc;

      std::vector<std::string> defaultVector;
      defaultVector.push_back(std::string("default"));
      desc.add<std::vector<std::string> >("lvl1Labels", defaultVector);

      // This default vector<ParameterSet> will be used when
      // the configuration does not include this parameter and
      // it also gets written into the generated cfi file.
      std::vector<edm::ParameterSet> defaultVPSet;
      edm::ParameterSet pset0;
      pset0.addParameter<std::string>("pathName", std::string("HLTPath"));
      std::vector<unsigned> defaultVectorU;
      defaultVectorU.push_back(1u);
      pset0.addParameter<std::vector<unsigned> >("prescales", defaultVectorU);
      defaultVPSet.push_back(pset0);

      edm::ParameterSetDescription validator;
      validator.add<std::string>("pathName");
      validator.add<std::vector<unsigned int> >("prescales");

      desc.addVPSet("prescaleTable", validator, defaultVPSet);

      desc.add<std::string>("lvl1DefaultLabel", std::string("default"));
      desc.add<bool>       ("forceDefault",     false);

      descriptions.add("PrescaleService", desc);
    }

  } // namespace service
} // namespace edm
