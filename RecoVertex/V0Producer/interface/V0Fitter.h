// -*- C++ -*-
//
// Package:    V0Producer
// Class:      V0Fitter
// 
/**\class V0Fitter V0Fitter.h RecoVertex/V0Producer/interface/V0Fitter.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
// $Id: V0Fitter.h,v 1.24 2010/08/05 22:06:39 wmtan Exp $
//
//

#ifndef RECOVERTEX__V0_FITTER_H
#define RECOVERTEX__V0_FITTER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Math/interface/angle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string>
#include <fstream>


class V0Fitter {
 public:
  V0Fitter(const edm::ParameterSet& theParams,
	   const edm::Event& iEvent, const edm::EventSetup& iSetup);
  ~V0Fitter();

  // Switching to L. Lista's reco::Candidate infrastructure for V0 storage
  const reco::VertexCompositeCandidateCollection& getKshorts() const;
  const reco::VertexCompositeCandidateCollection& getLambdas() const;
  const reco::VertexCompositeCandidateCollection& getXis() const;
  const reco::VertexCompositeCandidateCollection& getOmegas() const;

 private:
  // STL vector of VertexCompositeCandidate that will be filled with VertexCompositeCandidates by fitAll()
  reco::VertexCompositeCandidateCollection theKshorts;
  reco::VertexCompositeCandidateCollection theLambdas;
  reco::VertexCompositeCandidateCollection theXis;
  reco::VertexCompositeCandidateCollection theOmegas;

  // Tracker geometry for discerning hit positions
  const TrackerGeometry* trackerGeom;

  const MagneticField* magField;

  edm::InputTag recoAlg;
  edm::InputTag vtxAlg;
  bool useRefTrax;
  bool storeRefTrax;
  bool doKshorts;
  bool doLambdas;
  bool doXis;
  bool doOmegas;

  /*bool doPostFitCuts;
    bool doTkQualCuts;*/

  // Cuts
  double tkDCACut;
  double tkChi2Cut;
  int    tkNhitsCut;
  double chi2Cut;
  double rVtxCut;
  double rVtxSigCut;
  double lVtxCut;
  double lVtxSigCut;
  double collinCut;
  double xiChi2Cut;
  double xiRVtxCut;
  double xiRVtxSigCut;
  double xiLVtxCut;
  double xiLVtxSigCut;
  double xiCollinCut;
  double kShortMassCut;
  double lambdaMassCut;
  double xiMassCut;
  double omegaMassCut;
  double dauTransImpactSigCut;
  double dauLongImpactSigCut;
  double mPiPiCut;
  double innerHitPosCut;

  std::vector<reco::TrackBase::TrackQuality> qualities;

  edm::InputTag vtxFitter;

  // Helper method that does the actual fitting using the KalmanVertexFitter
  void fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  double findV0MassError(const GlobalPoint &vtxPos, std::vector<reco::TransientTrack> dauTracks);

  // Applies cuts to the VertexCompositeCandidates after they are fitted/created.
  //void applyPostFitCuts();

  // Stuff for debug file output.
  std::ofstream mPiPiMassOut;

  inline void initFileOutput() {
    mPiPiMassOut.open("mPiPi.txt", std::ios::app);
  }
  inline void cleanupFileOutput() {
    mPiPiMassOut.close();
  }
};

#endif
