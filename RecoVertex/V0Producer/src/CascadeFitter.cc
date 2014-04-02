// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeFitter
// 
/**\class CascadeFitter CascadeFitter.cc RecoVertex/V0Producer/src/CascadeFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Brian Drell
//         Created:  Fri May 18 22:57:40 CEST 2007
// $Id: CascadeFitter.cc,v 1.56 2011/12/23 08:13:39 innocent Exp $
//
//

#include "RecoVertex/V0Producer/interface/CascadeFitter.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <typeinfo>
#include <memory>
#include <vector>

// Constants
const ParticleMass pionMass = 0.13957018;
const double pionMassSquared = pionMass*pionMass;
const ParticleMass protonMass = 0.938272013;
const double protonMassSquared = protonMass*protonMass;
const ParticleMass kaonMass = 0.493677;
const double kaonMassSquared = kaonMass*kaonMass;
const ParticleMass kShortMass = 0.497614;
const ParticleMass lambdaMass = 1.115683;
const ParticleMass xiMass = 1.32171;
const ParticleMass omegaMass = 1.67245;
const ParticleMass DdMass = 1.86962;
const ParticleMass DsMass = 1.96847;

// Constructor and (empty) destructor
CascadeFitter::CascadeFitter(const edm::ParameterSet& theParameters,
		   const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using std::string;

  // Get the track reco algorithm from the ParameterSet
  recoAlg = theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm");
  vtxAlg  = theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm");
  v0Alg  = theParameters.getParameter<edm::InputTag>("v0RecoAlgorithm");

  batDauTrkMass = theParameters.getParameter<double>(string("batDauTrkMass"));

  // Second, initialize post-fit cuts
  tkChi2Cut = theParameters.getParameter<double>(string("tkChi2Cut"));
  tkNhitsCut = theParameters.getParameter<int>(string("tkNhitsCut"));
  chi2Cut = theParameters.getParameter<double>(string("vtxChi2Cut"));
  rVtxCut = theParameters.getParameter<double>(string("rVtxCut"));
  rVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance2DCut"));
  lVtxCut = theParameters.getParameter<double>(string("lVtxCut"));
  lVtxSigCut = theParameters.getParameter<double>(string("vtxSignificance3DCut"));
  collinCut = theParameters.getParameter<double>(string("collinearityCut"));
  casMassCut = theParameters.getParameter<double>(string("casMassCut"));
  dauTransImpactSigCut = theParameters.getParameter<double>(string("dauTransImpactSigCut"));
  dauLongImpactSigCut = theParameters.getParameter<double>(string("dauLongImpactSigCut"));
  std::vector<std::string> qual = theParameters.getParameter<std::vector<std::string> >("trackQualities");
  for (unsigned int ndx = 0; ndx < qual.size(); ndx++) {
    qualities.push_back(reco::TrackBase::qualityByName(qual[ndx]));
  }

  fitAll(iEvent, iSetup);
}

CascadeFitter::~CascadeFitter() {
}

// Method containing the algorithm for vertex reconstruction
void CascadeFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  // Create std::vectors for Tracks and TrackRefs (required for
  //  passing to the KalmanVertexFitter)
  std::vector<TrackRef> theTrackRefs;
  std::vector<TransientTrack> theTransTracks;

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;
  ESHandle<TrackerGeometry> trackerGeomHandle;
  ESHandle<GlobalTrackingGeometry> globTkGeomHandle;
  //cout << "Check 0" << endl;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByLabel(recoAlg, theTrackHandle);
  iEvent.getByLabel(vtxAlg,  theVertexHandle);
  iEvent.getByLabel(std::string("offlineBeamSpot"), theBeamSpotHandle);
  if( !theTrackHandle->size() ) return;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeomHandle);
  iSetup.get<GlobalTrackingGeometryRecord>().get(globTkGeomHandle);

  trackerGeom = trackerGeomHandle.product();
  magField = bFieldHandle.product();

  edm::Handle<reco::VertexCompositeCandidateCollection> v0Collection;
  iEvent.getByLabel(v0Alg, v0Collection);

  bool isVtxPV = 0;
  double xVtx=-99999.0;
  double yVtx=-99999.0;
  double zVtx=-99999.0;
  double xVtxError=-999.0;
  double yVtxError=-999.0;
  double zVtxError=-999.0;
  const reco::VertexCollection vtxCollection = *(theVertexHandle.product());
  reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();
  if(vtxCollection.size()>0 && !vtxPrimary->isFake() && vtxPrimary->tracksSize()>=2)
  {
    isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  }
  else {
    isVtxPV = 0;
    xVtx = theBeamSpotHandle->position().x();
    yVtx = theBeamSpotHandle->position().y();
    zVtx = 0.0;
    xVtxError = theBeamSpotHandle->BeamWidthX();
    yVtxError = theBeamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }

  // Fill vectors of TransientTracks and TrackRefs after applying preselection cuts.
  for(unsigned int indx = 0; indx < theTrackHandle->size(); indx++) {
    TrackRef tmpRef( theTrackHandle, indx );
    bool quality_ok = true;
    if (qualities.size()!=0) {
      quality_ok = false;
      for (unsigned int ndx_ = 0; ndx_ < qualities.size(); ndx_++) {
	if (tmpRef->quality(qualities[ndx_])){
	  quality_ok = true;
	  break;          
	}
      }
    }
    if( !quality_ok ) continue;

    if( tmpRef->normalizedChi2() < tkChi2Cut &&
        tmpRef->numberOfValidHits() >= tkNhitsCut ) {
      TransientTrack tmpTk( *tmpRef, &(*bFieldHandle), globTkGeomHandle );

      math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauTransImpactSig = dzvtx/dzerror;
      double dauLongImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > dauTransImpactSigCut && fabs(dauLongImpactSig) > dauLongImpactSigCut ) {
        theTrackRefs.push_back( tmpRef );
        theTransTracks.push_back( tmpTk );
      }
    }
  }

  if(v0Collection->size() > 0)
  {
    for(unsigned it=0; it<v0Collection->size(); ++it){
        
      const reco::VertexCompositeCandidate & theV0 = (*v0Collection)[it];

      int v0PDG = theV0.pdgId();
      double dauPlusMass = 0;
      double dauMinusMass = 0;
      double v0Mass = 0;
      if(v0PDG == 3122)
      {
        v0Mass = lambdaMass;
        dauPlusMass = protonMass;
        dauMinusMass = pionMass;
      }
      if(v0PDG == -3122)
      {
        v0Mass = lambdaMass;
        dauMinusMass = protonMass;
        dauPlusMass = pionMass;
      }
      if(fabs(v0PDG) == 310 )  
      {
        v0Mass = kShortMass;
        dauMinusMass = pionMass;
        dauPlusMass = pionMass;
      } 

      // check V0 mass to be within +- 20MeV...
      float massWindow = 0.020;
      if(theV0.mass() > v0Mass + massWindow || theV0.mass() < v0Mass - massWindow) continue;

      //get daughter tracks from V0 candidate
      vector<RecoChargedCandidate> v0daughters;
      vector<TrackRef> theDaughterTracks;

      //check momentum and add pion first, proton second
      if (theV0.daughter(0)->momentum().mag2() < theV0.daughter(1)->momentum().mag2()){
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theV0.daughter(0))) );
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theV0.daughter(1))) );
       } else {
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theV0.daughter(1))) );
         v0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                                  (theV0.daughter(0))) );
       }
       for(unsigned int j = 0; j < v0daughters.size(); ++j) {
         theDaughterTracks.push_back(v0daughters[j].track());
       }

       vector<TransientTrack> tracksForKalmanFit;
       for (unsigned int ndx = 0; ndx < theDaughterTracks.size(); ndx++) {
           tracksForKalmanFit.push_back(TransientTrack(theDaughterTracks[ndx], &(*bFieldHandle)));
       }

       vector<double> v0DauMasses;
       if (tracksForKalmanFit[0].charge() > 0) {
           v0DauMasses.push_back(dauPlusMass);
           v0DauMasses.push_back(dauMinusMass);
       } else {
           v0DauMasses.push_back(dauMinusMass);
           v0DauMasses.push_back(dauPlusMass);
       }

      // Loop over tracks and vertex good charged track pairs
      for(unsigned int trdx = 0; trdx < theTrackRefs.size(); trdx++) {

         if ( theTrackRefs[trdx].isNull() ) continue;

         math::XYZPoint bestvtx(xVtx,yVtx,zVtx);
         double dzvtx = theTrackRefs[trdx]->dz(bestvtx);
         double dxyvtx = theTrackRefs[trdx]->dxy(bestvtx);
         if(fabs(dzvtx)>50 || fabs(dxyvtx)>50) continue;

         bool match = false;
         // check if the pion is used in any v0 candidate in the collection and flag it
         for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
           if (theTrackRefs[trdx]->charge() == theDaughterTracks[j]->charge() &&
               theTrackRefs[trdx]->momentum() == theDaughterTracks[j]->momentum() ) match = true;
           if (match) break;
         }
         if (match) continue; // Track is already used in making the V0

         // check if pion is in *any* good V0
         // Placeholder
       
         if(v0PDG==3122 && theTrackRefs[trdx]->charge()>0) continue;
         if(v0PDG==-3122 && theTrackRefs[trdx]->charge()<0) continue;

         TransientTrack dau1TT(theDaughterTracks[0], &(*bFieldHandle) );
         TransientTrack dau2TT(theDaughterTracks[1], &(*bFieldHandle) );
         TransientTrack batDauTT(theTrackRefs[trdx], &(*bFieldHandle) );

         if (!dau1TT.isValid()) continue;
         if (!dau2TT.isValid()) continue;
         if (!batDauTT.isValid()) continue;

         //Creating a KinematicParticleFactory
         KinematicParticleFactoryFromTransientTrack pFactory;

         //initial chi2 and ndf before kinematic fits.
         float chi = 0.;
         float ndf = 0.;
         vector<RefCountedKinematicParticle> v0Particles;
         const ParticleMass v0daumass1 = v0DauMasses[0];
         const ParticleMass v0daumass2 = v0DauMasses[1];
         float v0daumass1_sigma = v0DauMasses[0]*1.e-6;
         float v0daumass2_sigma = v0DauMasses[1]*1.e-6;
         v0Particles.push_back(pFactory.particle(dau1TT,v0daumass1,chi,ndf,v0daumass1_sigma));
         v0Particles.push_back(pFactory.particle(dau2TT,v0daumass2,chi,ndf,v0daumass2_sigma));

         KinematicParticleVertexFitter fitter;
         RefCountedKinematicTree v0VertexFitTree;
         v0VertexFitTree = fitter.fit(v0Particles);
         if (!v0VertexFitTree->isValid()) continue;
  
         v0VertexFitTree->movePointerToTheTop();

         RefCountedKinematicParticle v0_vFit = v0VertexFitTree->currentParticle();
         RefCountedKinematicVertex v0_vFit_vertex = v0VertexFitTree->currentDecayVertex();

         v0VertexFitTree->movePointerToTheFirstChild();
         RefCountedKinematicParticle v0Pi1 = v0VertexFitTree->currentParticle();
         v0VertexFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle v0Pi2 = v0VertexFitTree->currentParticle();

         // now apply V0 mass constraint
         KinematicParticleFitter csFitterV0;
         float v0Mass_sigma = v0Mass*1.e-6;
         KinematicConstraint * v0_c = new MassKinematicConstraint(v0Mass,v0Mass_sigma);
         // add mass constraint to the V0 fit to do a constrained fit:  

         v0VertexFitTree->movePointerToTheTop();
         v0VertexFitTree = csFitterV0.fit(v0_c,v0VertexFitTree);
         if (!v0VertexFitTree->isValid()) continue;
         v0VertexFitTree->movePointerToTheTop();
         RefCountedKinematicParticle v0_vFit_withMC = v0VertexFitTree->currentParticle();

         if (!v0_vFit_withMC->currentState().isValid()) continue; 

         vector<RefCountedKinematicParticle> casFitParticles;
  
         float batDauTrkMass_sigma = batDauTrkMass*1.e-6;
         casFitParticles.push_back(pFactory.particle(batDauTT,batDauTrkMass,chi,ndf,batDauTrkMass_sigma));
         casFitParticles.push_back(v0_vFit_withMC);

         //fit cas 
         RefCountedKinematicTree casFitTree = fitter.fit(casFitParticles);
         if (!casFitTree->isValid()) continue;

           casFitTree->movePointerToTheTop();
           RefCountedKinematicParticle casCand = casFitTree->currentParticle();
           if (!casCand->currentState().isValid()) continue;

           RefCountedKinematicVertex casDecayVertex = casFitTree->currentDecayVertex();
           if (!casDecayVertex->vertexIsValid()) continue;

           if ( casDecayVertex->chiSquared()<0 || casDecayVertex->chiSquared()>1000 ) continue;

           float casC2Prob =
              ChiSquaredProbability((double)(casDecayVertex->chiSquared()),(double)(casDecayVertex->degreesOfFreedom()));
           if (casC2Prob < 0.0001) continue;

           if ( casCand->currentState().mass() > 3 ) continue;

           // get children from final Cascade fit
           casFitTree->movePointerToTheFirstChild();
           RefCountedKinematicParticle batDauCand = casFitTree->currentParticle();
           casFitTree->movePointerToTheNextChild();
           RefCountedKinematicParticle v0CandMC = casFitTree->currentParticle();
         
           if(!batDauCand->currentState().isValid() || !v0CandMC->currentState().isValid()) continue;

           // get batchlor pion and V0 parameters from Cascade fit
           KinematicParameters batDauKP = batDauCand->currentState().kinematicParameters();
           KinematicParameters v0CandKP = v0CandMC->currentState().kinematicParameters();

           GlobalVector casTotalP = GlobalVector (casCand->currentState().globalMomentum().x(),
                                                 casCand->currentState().globalMomentum().y(),
                                                 casCand->currentState().globalMomentum().z());

           GlobalVector batDauTotalP = GlobalVector(batDauKP.momentum().x(),batDauKP.momentum().y(),batDauKP.momentum().z());
           GlobalVector v0TotalP = GlobalVector(v0CandKP.momentum().x(),v0CandKP.momentum().y(),v0CandKP.momentum().z());

           double batDauTotalE = sqrt( batDauTotalP.mag2() + batDauTrkMass*batDauTrkMass );
           double v0TotalE = sqrt( v0TotalP.mag2() + v0Mass*v0Mass );
           double casTotalE = batDauTotalE + v0TotalE;

           const Particle::LorentzVector casP4(casTotalP.x(), casTotalP.y(), casTotalP.z(), casTotalE);

           Particle::Point casVtx((*casDecayVertex).position().x(), (*casDecayVertex).position().y(), (*casDecayVertex).position().z());
           std::vector<double> casVtxEVec;
           casVtxEVec.push_back( v0_vFit_vertex->error().cxx() );
           casVtxEVec.push_back( v0_vFit_vertex->error().cyx() );
           casVtxEVec.push_back( v0_vFit_vertex->error().cyy() );
           casVtxEVec.push_back( v0_vFit_vertex->error().czx() );
           casVtxEVec.push_back( v0_vFit_vertex->error().czy() );
           casVtxEVec.push_back( v0_vFit_vertex->error().czz() );
           SMatrixSym3D casVtxCovMatrix(casVtxEVec.begin(), casVtxEVec.end());
           const Vertex::CovarianceMatrix casVtxCov(casVtxCovMatrix);
           double casVtxChi2(casDecayVertex->chiSquared());
           double casVtxNdof(casDecayVertex->degreesOfFreedom());
           double casNormalizedChi2 = casVtxChi2/casVtxNdof;

           double casRVtxMag = 99999.0;
           double casLVtxMag = 99999.0;
           double casSigmaRvtxMag = 999.0;
           double casSigmaLvtxMag = 999.0;
           double casV0Angle = -100.0;

           GlobalVector casV0LineOfFlight = GlobalVector (casVtx.x() - xVtx,
                                                         casVtx.y() - yVtx,
                                                         casVtx.z() - zVtx);

           SMatrixSym3D casTotalCov;
           if(isVtxPV) casTotalCov = casVtxCovMatrix + vtxPrimary->covariance();
           else casTotalCov = casVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

           SVector3 distanceVector3D(casV0LineOfFlight.x(), casV0LineOfFlight.y(), casV0LineOfFlight.z());
           SVector3 distanceVector2D(casV0LineOfFlight.x(), casV0LineOfFlight.y(), 0.0);

           casV0Angle = angle(casV0LineOfFlight.x(), casV0LineOfFlight.y(), casV0LineOfFlight.z(),
                             casTotalP.x(), casTotalP.y(), casTotalP.z());
           casLVtxMag = casV0LineOfFlight.mag();
           casRVtxMag = casV0LineOfFlight.perp();
           casSigmaLvtxMag = sqrt(ROOT::Math::Similarity(casTotalCov, distanceVector3D)) / casLVtxMag;
           casSigmaRvtxMag = sqrt(ROOT::Math::Similarity(casTotalCov, distanceVector2D)) / casRVtxMag;

           if( casNormalizedChi2 > chi2Cut ||
               casRVtxMag < rVtxCut ||
               casRVtxMag / casSigmaRvtxMag < rVtxSigCut ||
               casLVtxMag < lVtxCut ||
               casLVtxMag / casSigmaLvtxMag < lVtxSigCut ||
               cos(casV0Angle) < collinCut
           ) continue;

           // Cuts finished, Create the VertexCompositeCandidate object for Cascade that will be stored in the Event
           VertexCompositeCandidate* theCas = 0;
           theCas = new VertexCompositeCandidate(0, casP4, casVtx, casVtxCov, casVtxChi2, casVtxNdof);

           RecoChargedCandidate dauCand(theTrackRefs[trdx]->charge(), Particle::LorentzVector(batDauTotalP.x(),
                                                   batDauTotalP.y(), batDauTotalP.z(),
                                                   batDauTotalE), casVtx);
           dauCand.setTrack(theTrackRefs[trdx]);

           AddFourMomenta addp4;
           // Store the daughter Candidates in the VertexCompositeCandidates 
           // if they pass mass cuts
           theCas->addDaughter(theV0);
           theCas->addDaughter(dauCand);
           if(v0PDG == 3122 && fabs(batDauTrkMass-0.14)<0.1) theCas->setPdgId(3312);
           if(v0PDG == -3122 && fabs(batDauTrkMass-0.14)<0.1) theCas->setPdgId(-3312);
           if(v0PDG == 3122 && fabs(batDauTrkMass-0.5)<0.1) theCas->setPdgId(3334);
           if(v0PDG == -3122 && fabs(batDauTrkMass-0.5)<0.1) theCas->setPdgId(-3334);
           if(fabs(v0PDG) == 310 && fabs(batDauTrkMass-0.14)<0.1 && theTrackRefs[trdx]->charge()>0) theCas->setPdgId(411); 
           if(fabs(v0PDG) == 310 && fabs(batDauTrkMass-0.14)<0.1 && theTrackRefs[trdx]->charge()<0) theCas->setPdgId(-411);  
           if(fabs(v0PDG) == 310 && fabs(batDauTrkMass-0.5)<0.1 && theTrackRefs[trdx]->charge()>0) theCas->setPdgId(431);  
           if(fabs(v0PDG) == 310 && fabs(batDauTrkMass-0.5)<0.1 && theTrackRefs[trdx]->charge()<0) theCas->setPdgId(-431);
           double theCasMass = GetCasMass(theCas->pdgId());
           addp4.set( *theCas );
           if( theCas->mass() < theCasMass + casMassCut &&
               theCas->mass() > theCasMass - casMassCut ) theCass.push_back( *theCas );
           if(theCas) delete theCas;
           theCas = 0;
        }
      }
    }
}

// Get methods

const reco::VertexCompositeCandidateCollection& CascadeFitter::getCass() const {
  return theCass;
}

// Experimental
double CascadeFitter::GetCasMass(int pdgcode) { 

  double mass = 0;

  if(fabs(pdgcode) == 3312) mass=xiMass;
  if(fabs(pdgcode) == 3334) mass=omegaMass;
  if(fabs(pdgcode) == 411)  mass=DdMass;
  if(fabs(pdgcode) == 431)  mass=DsMass;

  return mass;
}
