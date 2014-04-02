import FWCore.ParameterSet.Config as cms

generalCascadeCandidates = cms.EDProducer("CascadeProducer",
                                     
    # InputTag that tells which TrackCollection to use for vertexing
    trackRecoAlgorithm = cms.InputTag('generalTracks'),
    vertexRecoAlgorithm = cms.InputTag('offlinePrimaryVertices'),
    v0RecoAlgorithm = cms.InputTag('generalV0Candidates','Lambda'),

    # set to true, uses tracks refit by the KVF for V0Candidate kinematics
    #  NOTE: useSmoothing is automatically set to FALSE
    #  if using the AdaptiveVertexFitter (which is NOT recommended)
    useSmoothing = cms.bool(True),
                                     
    # Select tracks using TrackBase::TrackQuality.
    # Select ALL tracks by leaving this vstring empty, which
    #   is equivalent to using 'loose'
    #trackQualities = cms.vstring('highPurity', 'goodIterative'),
    trackQualities = cms.vstring('loose'),
                                     
    # The next parameters are cut values.
    # All distances are in cm, all energies in GeV, as usual.

    # --Track quality/compatibility cuts--
    #   Normalized track Chi2 <
    tkChi2Cut = cms.double(5.0),
    #   Number of valid hits on track >=
    tkNhitsCut = cms.int32(0),
    #   Track impact parameter significance >
    dauTransImpactSigCut = cms.double(0.),
    dauLongImpactSigCut = cms.double(0.),

    batDauTrkMass = cms.double(0.139570),

    # --V0 Vertex cuts--
    #   Vertex chi2 < 
    vtxChi2Cut = cms.double(25.0),
    #   Lambda collinearity cut
    #   (UNUSED)
    collinearityCut = cms.double(-2),
    #   Vertex radius cut >
    #   (UNUSED)
    rVtxCut = cms.double(0.0),
    #   V0 decay length from primary cut >
    #   (UNUSED)
    lVtxCut = cms.double(0.0),
    #   Radial vertex significance >
    vtxSignificance2DCut = cms.double(0.0),
    #   3D vertex significance using primary vertex
    #   (UNUSED)
    vtxSignificance3DCut = cms.double(0.0),
    
    #   V0 mass window, Candidate mass must be within these values of
    #     the PDG mass to be stored in the collection
    casMassCut = cms.double(0.30),
)


