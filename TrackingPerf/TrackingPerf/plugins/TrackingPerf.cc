// system include files
#include <memory>

// user include files

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"


#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"


#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"


#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"


#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"


/*#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"*/
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"



#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
#include <TLorentzVector.h>
//includes for track-sim association

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/interface/TrackingParticleIP.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"



#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoTracker/FinalTrackSelectors/plugins/getBestVertex.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" 

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"


#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



using reco::TrackCollection;
using TrackingParticleRefKeyToIndex = std::unordered_map<reco::RecoToSimCollection::index_type, size_t>;


class TrackingPerf : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TrackingPerf(const edm::ParameterSet&);
      ~TrackingPerf();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      
      // ----------member data ---------------------------
  
  void clearVariables();
  edm::InputTag tkTraj_;



  //// GETTOKENDECLARATIONS 


  //track ( and event ) information 
  const edm::EDGetTokenT<edm::View<reco::Track> > trackToken_;
  const edm::EDGetTokenT<edm::View<reco::Track> > trackSrc_;
  const edm::EDGetTokenT<std::vector<Trajectory> > trajSrc_;
  const edm::EDGetTokenT<TrajTrackAssociationCollection> trajTrackAssociationSrc_;
  const edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> trackAssociatorToken_;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;

  const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken_;
  edm::EDGetTokenT<TrackingParticleRefVector> trackingParticleRefToken_;

 
  //jet information 
  const edm::EDGetTokenT<edm::View<reco::Jet> > jetToken_;
  //const edm::EDGetTokenT<edm::View<pat::Jet> > viewJetToken_;
  const edm::EDGetTokenT<pat::JetCollection> jetPuppiToken_;
  const edm::EDGetTokenT<pat::JetCollection> ak8jetToken_;
  const edm::EDGetTokenT<pat::JetCollection> ak10jetToken_;
  //const edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;

  //met information 
  //const edm::EDGetTokenT<std::vector<pat::MET> > metTokenAOD_;
  //const edm::EDGetTokenT<pat::METCollection> metTokenPuppi_;
  //const edm::EDGetTokenT<pat::METCollection> metTokenNoHF_;
  //const edm::EDGetTokenT<pat::METCollection> metTokenMINIAOD_;
  const edm::EDGetTokenT<reco::CaloMETCollection> metToken_;
  


  //gen information
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  const edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
  const edm::EDGetTokenT<reco::GenJetCollection> genTTXJetsToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  const edm::EDGetTokenT<LHEEventProduct> LHEEventProductToken_;

  //pat information
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;

  std::string parametersDefinerName_;


  //electrons
  //const edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronToken_;
  //const edm::EDGetTokenT<pat::ElectronCollection> electronPATToken_;

  //muons 
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;



  //additional 
  bool useCluster_;
  TTree *smalltree;


  edm::Service<TFileService> fs;

  std::string ttrhbuilder_; 

  edm::ESHandle<MagneticField> bField;




  /////// VARIABLES FOR FILLING THE TREE 
  //-----------------------
  //fill the tree per track
  //std::vector<int> tree_track_nclusters;
  std::vector< double > tree_track_charge;
  std::vector< double > tree_track_pt;
  std::vector< double > tree_track_outerPt; 
  std::vector< double > tree_track_eta;
  std::vector< double > tree_track_phi;
  std::vector<int>      tree_track_nhits;
  
  std::vector<int>      tree_track_numberOfValidPixelHits;
  std::vector<int>      tree_track_numberOfValidStripHits;
  std::vector<int>      tree_track_numberOfValidStripTIBHits;
  std::vector<int>      tree_track_numberOfValidStripTIDHits;
  std::vector<int>      tree_track_numberOfValidStripTOBHits;
  std::vector<int>      tree_track_numberOfValidStripTECHits;
  std::vector<int>      tree_track_numberOfValidPixelBarrelHits;
  std::vector<int>      tree_track_numberOfValidPixelEndcapHits;
  std::vector<int>      tree_track_hasValidHitInPixelLayer;
  std::vector<int>      tree_track_trackerLayersWithMeasurement;
  std::vector<int>      tree_track_pixelLayersWithMeasurement;
  std::vector<int>      tree_track_stripTECLayersWithMeasurement ;
  std::vector<int>      tree_track_stripTIBLayersWithMeasurement;
  std::vector<int>      tree_track_stripTIDLayersWithMeasurement;
  std::vector<int>      tree_track_stripTOBLayersWithMeasurement;
  
  
  
  std::vector<double >  tree_track_NChi2;
  //std::vector<double >  tree_track_Quality;
  std::vector<double >  tree_track_isHighQuality;
  std::vector<double >  tree_track_isLoose;
  std::vector<double >  tree_track_isTight;
  std::vector< double>  tree_track_dxy; // dxy with respect to beam spot position
  std::vector< double>  tree_track_dxyError;
  std::vector< double>  tree_track_dz;
  std::vector< double>  tree_track_dzError ;
  std::vector<int>      tree_track_numberOfLostHits;
  std::vector<int>      tree_track_numberOfValidHits;
  std::vector<unsigned int>     tree_track_originalAlgo; // definition as comments at the end of the file, http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_1_3/doc/html/d8/df2/classreco_1_1TrackBase.html#aca7611bd1a33d535cefc72b6e497ece8
  std::vector<unsigned int>     tree_track_algo; 
  std::vector<unsigned short>   tree_track_stopReason;
  std::vector<int>              tree_track_nSimHits   ; 
  std::vector<bool>             tree_track_isSimMatched;
  std::vector<int>              tree_track_nPixel;
  std::vector<int>              tree_track_nStrip;
  std::vector< double >         tree_track_vx;
  std::vector< double >         tree_track_vy; 
  std::vector< double >         tree_track_vz;
  std::vector<double>           tree_track_firsthit_X; 
  std::vector<double>           tree_track_firsthit_Y; 
  std::vector<double>           tree_track_firsthit_Z;
  std::vector<double>           tree_track_firsthit_phi; 
                  
  
  std::vector< double >         tree_track_simtrack_charge;	              
  std::vector< double >         tree_track_simtrack_pt;	              
  std::vector< double >         tree_track_simtrack_eta	;	      
  std::vector< double >         tree_track_simtrack_phi	;	      
  std::vector<bool>             tree_track_simtrack_longLived	 ;
 // std::vector<int>              tree_track_simtrack_matchedHit	;
  std::vector<int>              tree_track_simtrack_pdgId;		
  std::vector<int>              tree_track_simtrack_numberOfTrackerHits  ;
  std::vector<int>              tree_track_simtrack_numberOfTrackerLayers;
  std::vector<double>           tree_track_simtrack_mass  ;	              
  std::vector<int>              tree_track_simtrack_status;	              

  std::vector<double>           tree_track_genVertexPos_X;                
  std::vector<double>           tree_track_genVertexPos_Y;                
  std::vector<double>           tree_track_genVertexPos_Z;                
  
                 
  std::vector<int>              tree_track_recoVertex_idx; 
  std::vector<int>              tree_track_recoJet_idx; 
  
  int    runNumber,eventNumber,lumiBlock; 
  
  std::vector<int>              tree_track_idxSiClusterFirst  ;
  std::vector<int>              tree_track_idxPixClusterFirst ;
  std::vector<int>              tree_track_idxSiClusterLast   ;
  std::vector<int>              tree_track_idxPixClusterLast  ;




  
  
  



  //then fill information of clusters
  //attached to the tracks  ----
  
  //strip cluster info  ----
  std::vector< int >     tree_SiCluster_subDet; 
  std::vector< int >     tree_SiCluster_PetalSide;
  std::vector< int >     tree_SiCluster_LayerNbr; 
  std::vector< int >     tree_SiCluster_WheelSide; 
  std::vector< double >  tree_SiCluster_charge;
  std::vector< double >  tree_SiCluster_SoverN;
  std::vector< double >  tree_SiCluster_noise;
  std::vector< int >     tree_SiCluster_width;
  std::vector< double >  tree_SiCluster_barycenter;
  std::vector< int >     tree_SiCluster_detID;
  
   // strip cluster position ----
  std::vector< float >   tree_SiCluster_locX;
  std::vector< float >   tree_SiCluster_locY;
  std::vector< double >  tree_SiCluster_tsosx;
  std::vector< double >  tree_SiCluster_tsosy;
  std::vector< float >   tree_SiCluster_globX;
  std::vector< float >   tree_SiCluster_globY;
  std::vector< float >   tree_SiCluster_globZ;
  std::vector< float >   tree_SiCluster_tsosglobX;
  std::vector< float >   tree_SiCluster_tsosglobY;
  std::vector< float >   tree_SiCluster_tsosglobZ;
  
  // pixel cluster info----
  std::vector< int >     tree_PixCluster_LayerNbr; 
  std::vector< double >  tree_PixCluster_charge;
  std::vector< int >     tree_PixCluster_isBPix;
  std::vector< bool >    tree_PixCluster_detID;
 
  
  
  // pixel cluster position ----
  std::vector< float >  tree_PixCluster_locX;
  std::vector< float >  tree_PixCluster_locY;
  std::vector< float >  tree_PixCluster_globX;
  std::vector< float >  tree_PixCluster_globY;
  std::vector< float >  tree_PixCluster_globZ;
  std::vector< float >  tree_PixCluster_tsosglobX;
  std::vector< float >  tree_PixCluster_tsosglobY;
  std::vector< float >  tree_PixCluster_tsosglobZ;
  std::vector< float >  tree_PixCluster_tsosX;
  std::vector< float >  tree_PixCluster_tsosY;
  
  
  
  // strips infos ----
  std::vector< int >    tree_Strips_stripCharges;
  std::vector< float > 	tree_Strips_stripGains ;
  std::vector< float > 	tree_Strips_stripNoises;
  std::vector< bool >	tree_Strips_stripQualitiesBad ;

      
  //--------------------------------
  // Tracking Particle infos ------- 
  //--------------------------------
   
  std::vector< double >         tree_simtrack_simtrack_charge;	              
  std::vector< double >         tree_simtrack_simtrack_pt;	              
  std::vector< double >         tree_simtrack_simtrack_eta;	      
  std::vector< double >         tree_simtrack_simtrack_phi;	      
  std::vector<bool>             tree_simtrack_simtrack_longLived;
  //std::vector<int>              tree_simtrack_simtrack_matchedHit;
  std::vector<int>              tree_simtrack_simtrack_pdgId;		
  std::vector<int>              tree_simtrack_simtrack_numberOfTrackerHits;
  std::vector<int>              tree_simtrack_simtrack_numberOfTrackerLayers;
  std::vector<double>           tree_simtrack_simtrack_mass;	              
  std::vector<int>              tree_simtrack_simtrack_status;	              

  std::vector<double>           tree_simtrack_genVertexPos_X;                
  std::vector<double>           tree_simtrack_genVertexPos_Y;                
  std::vector<double>           tree_simtrack_genVertexPos_Z;        
  std::vector<bool>             tree_simtrack_isRecoMatched;
  std::vector<double>           tree_simtrack_pca_dxy;
  std::vector<double>           tree_simtrack_pca_dz   ; 
  std::vector<std::vector<int> > tree_simtrack_trkIdx;  
  
  
  //--------------------------------
  // Beam spot infos ------- 
  //--------------------------------
  float	      tree_bs_PosX;	   
  float	      tree_bs_PosY;	   
  float	      tree_bs_PosZ; 
  
  //--------------------------------
  // vertex infos ------- 
  //--------------------------------
  
  // for the main vertex
  
  std::vector<float> tree_vtx_PosX;	
  std::vector<float> tree_vtx_PosY;	
  std::vector<float> tree_vtx_PosZ; 
  
  std::vector<float> tree_vtx_PosXError;	
  std::vector<float> tree_vtx_PosYError;	
  std::vector<float> tree_vtx_PosZError; 
  
  
 //--------------------------------
 // jet infos ------- 
 //--------------------------------
 
 std::vector<float> tree_jet_pt;
 std::vector<float> tree_jet_eta;
 std::vector<float> tree_jet_phi;
 std::vector<float> tree_jet_vx;
 std::vector<float> tree_jet_vy;
 std::vector<float> tree_jet_vz;


 //--------------------------------
 // met infos ------- 
 //--------------------------------


  

  //--------------------------------
 // gen infos ------- 
 //-------------------------------- 

  std::vector< double > tree_genParticle_charge;
  std::vector< double > tree_genParticle_pt;
  std::vector< double > tree_genParticle_eta;
  std::vector< double > tree_genParticle_phi;
  std::vector< double > tree_genParticle_pdgId;
  std::vector< double > tree_genParticle_vx;
  std::vector< double > tree_genParticle_vy;
  std::vector< double > tree_genParticle_vz;
  std::vector< double > tree_genParticle_mass;
  std::vector< int >    tree_genParticle_statusCode; 
  std::vector< int > tree_genParticle_mother;



 //--------------------------------
 // gen jet infos ------- 
 //--------------------------------


  std::vector<float> tree_genJet_pt;
  std::vector<float> tree_genJet_eta;
  std::vector<float> tree_genJet_phi;
  std::vector<float> tree_genJet_vx;
  std::vector<float> tree_genJet_vy;
  std::vector<float> tree_genJet_vz;
  std::vector<float> tree_genJet_mass; 
  std::vector<float> tree_genJet_energy; 
  std::vector<float> tree_genJet_status; 
  std::vector<float> tree_genJet_pdgId; 


 //--------------------------------
 // gen event info ------- 
 //--------------------------------


 //--------------------------------
 // lhe event infos ------- 
 //--------------------------------

 //--------------------------------
 // PF infos ------- 
 //--------------------------------


 //--------------------------------
 // electrons infos ------- 
 //--------------------------------


  std::vector<float> tree_electron_pt;
  std::vector<float> tree_electron_eta;
  std::vector<float> tree_electron_phi;
  std::vector<float> tree_electron_vx;
  std::vector<float> tree_electron_vy;
  std::vector<float> tree_electron_vz;
  std::vector<float> tree_electron_energy; 




 //--------------------------------
 // muons infos ------- 
 //--------------------------------

  std::vector<float> tree_muon_pt;
  std::vector<float> tree_muon_eta;
  std::vector<float> tree_muon_phi;
  std::vector<float> tree_muon_vx;
  std::vector<float> tree_muon_vy;
  std::vector<float> tree_muon_vz;
  std::vector<float> tree_muon_energy;

  
  
};



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrackingPerf::TrackingPerf(const edm::ParameterSet& iConfig):
	trackToken_(consumes<edm::View<reco::Track> >(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
	trackSrc_( consumes<edm::View<reco::Track> >( iConfig.getParameter<edm::InputTag>("trackLabel") )),
	trackAssociatorToken_(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociator"))),
	beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpot"))),
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertices"))),
	jetToken_(consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("jetInput"))),
  //metToken_(consumes<reco::CaloMETCollection> (iConfig.getParameter<edm::InputTag>("metInput"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetInput"))),
  genEventInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoInput"))),
  LHEEventProductToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEEventProductInput"))),
  pfcandsToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfcands"))),
  parametersDefinerName_(iConfig.getUntrackedParameter<std::string>("parametersDefiner")),
  //electronPATToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronPATInput"))),
  //electronToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electronInput"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonInput"))),
  useCluster_(iConfig.getUntrackedParameter<bool>("useCluster"))
{


    usesResource("TFileService");
   
   //trackSrc_ = consumes<reco::TrackCollection>(src_);
   ttrhbuilder_ = iConfig.getParameter<std::string>("TTRHBuilder");
   
   
   smalltree = fs->make<TTree>("ttree", "ttree");
   
   
   
   
    //-----------------------
   //fill the tree per track
   //smalltree->Branch("tree_track_nclusters",         &tree_track_nclusters);        
   smalltree->Branch("tree_track_pt",		     &tree_track_pt);
   smalltree->Branch("tree_track_outerPt",           &tree_track_outerPt ); 		      
   smalltree->Branch("tree_track_eta",  	     &tree_track_eta ); 	      
   smalltree->Branch("tree_track_phi",  	     &tree_track_phi ); 	      
   smalltree->Branch("tree_track_nhits",	     &tree_track_nhits);	      
   smalltree->Branch("tree_track_NChi2",	     &tree_track_NChi2);	      
   //smalltree->Branch("tree_track_Quality",	     &tree_track_Quality );
   smalltree->Branch("tree_track_isHighQuality",     &tree_track_isHighQuality);
   smalltree->Branch("tree_track_isLoose",	     &tree_track_isLoose);
   smalltree->Branch("tree_track_isTight",	     &tree_track_isTight);        
   smalltree->Branch("tree_track_dxy",  	     &tree_track_dxy );        
   smalltree->Branch("tree_track_dxyError",	     &tree_track_dxyError);	    
   smalltree->Branch("tree_track_dz",		     &tree_track_dz);		
   smalltree->Branch("tree_track_dzError",	     &tree_track_dzError  );	    
   smalltree->Branch("tree_track_numberOfLostHits",  &tree_track_numberOfLostHits); 
   smalltree->Branch("tree_track_numberOfValidHits", &tree_track_numberOfValidHits); 
   smalltree->Branch("tree_track_originalAlgo",      &tree_track_originalAlgo);      
   smalltree->Branch("tree_track_algo", 	     &tree_track_algo); 
   smalltree->Branch("tree_track_stopReason",	     &tree_track_stopReason);   
   //smalltree->Branch("tree_track_nSimHits",	     &tree_track_nSimHits   ); 
   smalltree->Branch("tree_track_isSimMatched",      &tree_track_isSimMatched	); 
   
   
   smalltree->Branch("tree_track_numberOfValidPixelHits",                &tree_track_numberOfValidPixelHits);
   smalltree->Branch("tree_track_numberOfValidStripHits",                &tree_track_numberOfValidStripHits);
   smalltree->Branch("tree_track_numberOfValidStripTIBHits",             &tree_track_numberOfValidStripTIBHits);
   smalltree->Branch("tree_track_numberOfValidStripTIDHits",             &tree_track_numberOfValidStripTIDHits);
   smalltree->Branch("tree_track_numberOfValidStripTOBHits",             &tree_track_numberOfValidStripTOBHits);
   smalltree->Branch("tree_track_numberOfValidStripTECHits",             &tree_track_numberOfValidStripTECHits);
   smalltree->Branch("tree_track_numberOfValidPixelBarrelHits",          &tree_track_numberOfValidPixelBarrelHits);
   smalltree->Branch("tree_track_numberOfValidPixelEndcapHits",          &tree_track_numberOfValidPixelEndcapHits);
   smalltree->Branch("tree_track_hasValidHitInPixelLayer",               &tree_track_hasValidHitInPixelLayer);
   smalltree->Branch("tree_track_trackerLayersWithMeasurement",          &tree_track_trackerLayersWithMeasurement);
   smalltree->Branch("tree_track_pixelLayersWithMeasurement",            &tree_track_pixelLayersWithMeasurement);
   smalltree->Branch("tree_track_stripTECLayersWithMeasurement" ,        &tree_track_stripTECLayersWithMeasurement);
   smalltree->Branch("tree_track_stripTIBLayersWithMeasurement",         &tree_track_stripTIBLayersWithMeasurement);
   smalltree->Branch("tree_track_stripTIDLayersWithMeasurement",         &tree_track_stripTIDLayersWithMeasurement);
   smalltree->Branch("tree_track_stripTOBLayersWithMeasurement",         &tree_track_stripTOBLayersWithMeasurement);
   
   
   
    
   smalltree->Branch("tree_track_nPixel",			  &tree_track_nPixel );
   smalltree->Branch("tree_track_nStrip",			  &tree_track_nStrip );
   smalltree->Branch("tree_track_vx",                             &tree_track_vx ); 
   smalltree->Branch("tree_track_vy",                             &tree_track_vy ); 
   smalltree->Branch("tree_track_vz",                             &tree_track_vz ); 
   smalltree->Branch("tree_track_firsthit_X",                       &tree_track_firsthit_Z); 
   smalltree->Branch("tree_track_firsthit_Y",                       &tree_track_firsthit_Y);
   smalltree->Branch("tree_track_firsthit_Z",                       &tree_track_firsthit_Z);
   smalltree->Branch("tree_track_firsthit_phi",                     &tree_track_firsthit_phi); 
   
   smalltree->Branch("tree_track_simtrack_charge",		  &tree_track_simtrack_charge );		      
   smalltree->Branch("tree_track_simtrack_pt",  		  &tree_track_simtrack_pt );	      
   smalltree->Branch("tree_track_simtrack_eta", 		  &tree_track_simtrack_eta  );  	      
   smalltree->Branch("tree_track_simtrack_phi", 		  &tree_track_simtrack_phi  );  	      
   smalltree->Branch("tree_track_simtrack_longLived",		  &tree_track_simtrack_longLived     );
   //smalltree->Branch("tree_track_simtrack_matchedHit",  	  &tree_track_simtrack_matchedHit   );
   smalltree->Branch("tree_track_simtrack_pdgId",		  &tree_track_simtrack_pdgId ); 		
   smalltree->Branch("tree_track_simtrack_numberOfTrackerHits",   &tree_track_simtrack_numberOfTrackerHits   );
   smalltree->Branch("tree_track_simtrack_numberOfTrackerLayers", &tree_track_simtrack_numberOfTrackerLayers );
   smalltree->Branch("tree_track_simtrack_mass",		  &tree_track_simtrack_mass   );		      
   smalltree->Branch("tree_track_simtrack_status",		  &tree_track_simtrack_status );		      
 
   smalltree->Branch("tree_track_genVertexPos_X", &tree_track_genVertexPos_X ); 		  
   smalltree->Branch("tree_track_genVertexPos_Y", &tree_track_genVertexPos_Y ); 		  
   smalltree->Branch("tree_track_genVertexPos_Z", &tree_track_genVertexPos_Z ); 		  
   smalltree->Branch("tree_track_recoVertex_idx", &tree_track_recoVertex_idx); 
   smalltree->Branch("tree_track_recoJet_idx",    &tree_track_recoJet_idx); 
   
   
   if(useCluster_){
     smalltree->Branch("tree_track_idxSiClusterFirst"   , &tree_track_idxSiClusterFirst);
     smalltree->Branch("tree_track_idxPixClusterFirst"  , &tree_track_idxPixClusterFirst);
     smalltree->Branch("tree_track_idxSiClusterLast"    , &tree_track_idxSiClusterLast);
     smalltree->Branch("tree_track_idxPixClusterLast"   , &tree_track_idxPixClusterLast);
     smalltree->Branch("tree_SiCluster_charge",	  &tree_SiCluster_charge);	 
     smalltree->Branch("tree_SiCluster_SoverN",	  &tree_SiCluster_SoverN);	
     smalltree->Branch("tree_SiCluster_noise",	  &tree_SiCluster_noise );	
     smalltree->Branch("tree_SiCluster_width",	  &tree_SiCluster_width );	
     smalltree->Branch("tree_SiCluster_barycenter", &tree_SiCluster_barycenter);  
     smalltree->Branch("tree_SiCluster_detID",	  &tree_SiCluster_detID );    
   }
    
   //then fill information of clusters
   //attached to the tracks
   
   smalltree->Branch("tree_SiCluster_subDet",	     &tree_SiCluster_subDet);	   
   smalltree->Branch("tree_SiCluster_PetalSide",     &tree_SiCluster_PetalSide );    
   smalltree->Branch("tree_SiCluster_LayerNbr",      &tree_SiCluster_LayerNbr);      
   smalltree->Branch("tree_SiCluster_WheelSide",     &tree_SiCluster_WheelSide );    
  
   
   
    // cluster position ----
   /*smalltree->Branch("tree_SiCluster_locX",	   &tree_SiCluster_locX );	
   smalltree->Branch("tree_SiCluster_locY",	 &tree_SiCluster_locY );    
   smalltree->Branch("tree_SiCluster_tsosx",	  &tree_SiCluster_tsosx);	
   smalltree->Branch("tree_SiCluster_tsosy",	  &tree_SiCluster_tsosy); */	   
   
   //smalltree->Branch("tree_SiCluster_globX",	 &tree_SiCluster_globX );     
   //smalltree->Branch("tree_SiCluster_globY",	 &tree_SiCluster_globY );    
   //smalltree->Branch("tree_SiCluster_globZ",	 &tree_SiCluster_globZ );     
   
   /*smalltree->Branch("tree_SiCluster_tsosglobX", &tree_SiCluster_tsosglobX);  
   smalltree->Branch("tree_SiCluster_tsosglobY", &tree_SiCluster_tsosglobY);  
   smalltree->Branch("tree_SiCluster_tsosglobZ", &tree_SiCluster_tsosglobZ);  */
   
   
    // pixcluster position ----
   //smalltree->Branch("tree_PixCluster_locX",      &tree_PixCluster_locX );	  
   //smalltree->Branch("tree_PixCluster_locY",      &tree_PixCluster_locY );	  
   
   //smalltree->Branch("tree_PixCluster_globX",	  &tree_PixCluster_globX );	
   //smalltree->Branch("tree_PixCluster_globY",	  &tree_PixCluster_globY );    
   //smalltree->Branch("tree_PixCluster_globZ",	  &tree_PixCluster_globZ );	
  
   /*smalltree->Branch("tree_PixCluster_tsosglobX", &tree_PixCluster_tsosglobX);  
   smalltree->Branch("tree_PixCluster_tsosglobY", &tree_PixCluster_tsosglobY);  
   smalltree->Branch("tree_PixCluster_tsosglobZ", &tree_PixCluster_tsosglobZ);   
   smalltree->Branch("tree_PixCluster_tsosX",	  &tree_PixCluster_tsosX);  
   smalltree->Branch("tree_PixCluster_tsosY",	  &tree_PixCluster_tsosY); */
   
   
   smalltree->Branch("tree_PixCluster_isBPix",   &tree_PixCluster_isBPix); 
   smalltree->Branch("tree_PixCluster_LayerNbr", &tree_PixCluster_LayerNbr); 
   if(useCluster_) smalltree->Branch("tree_PixCluster_charge",   &tree_PixCluster_charge);
   smalltree->Branch("tree_PixCluster_detID",	 &tree_PixCluster_detID);
   
   smalltree->Branch("runNumber",  &runNumber,  "runNumber/I");
   smalltree->Branch("eventNumber",&eventNumber,"eventNumber/I");
   smalltree->Branch("lumiBlock"  ,&lumiBlock,  "lumiBlock/I");
   smalltree->Branch("tree_bs_PosX", &tree_bs_PosX,  "tree_bs_PosX/F"  );
   smalltree->Branch("tree_bs_PosY", &tree_bs_PosY,  "tree_bs_PosY/F" );
   smalltree->Branch("tree_bs_PosZ", &tree_bs_PosZ,  "tree_bs_PosZ/F"  );
  
   smalltree->Branch("tree_vtx_PosX", &tree_vtx_PosX);        
   smalltree->Branch("tree_vtx_PosY", &tree_vtx_PosY);        
   smalltree->Branch("tree_vtx_PosZ", &tree_vtx_PosZ); 
   
   smalltree->Branch("tree_vtx_PosXError", &tree_vtx_PosXError);        
   smalltree->Branch("tree_vtx_PosYError", &tree_vtx_PosYError);        
   smalltree->Branch("tree_vtx_PosZError", &tree_vtx_PosZError); 
  
  
  //  int runNumber,eventNumber,lumiBlock;
  //
  
  
    // Tracking particle info ----
  
   
   smalltree->Branch("tree_simtrack_simtrack_charge",		     &tree_simtrack_simtrack_charge );  		      
   smalltree->Branch("tree_simtrack_simtrack_pt",		     &tree_simtrack_simtrack_pt );	      
   smalltree->Branch("tree_simtrack_simtrack_eta",		     &tree_simtrack_simtrack_eta  );	      
   smalltree->Branch("tree_simtrack_simtrack_phi",		     &tree_simtrack_simtrack_phi  );	      
   smalltree->Branch("tree_simtrack_simtrack_longLived",	     &tree_simtrack_simtrack_longLived     );
   //smalltree->Branch("tree_simtrack_simtrack_matchedHit",	     &tree_simtrack_simtrack_matchedHit   );
   smalltree->Branch("tree_simtrack_simtrack_pdgId",		     &tree_simtrack_simtrack_pdgId );		
   //smalltree->Branch("tree_simtrack_simtrack_numberOfTrackerHits",   &tree_simtrack_simtrack_numberOfTrackerHits   );
   //smalltree->Branch("tree_simtrack_simtrack_numberOfTrackerLayers", &tree_simtrack_simtrack_numberOfTrackerLayers );
   smalltree->Branch("tree_simtrack_simtrack_mass",		     &tree_simtrack_simtrack_mass   );  		      
   smalltree->Branch("tree_simtrack_simtrack_status",		     &tree_simtrack_simtrack_status );  		      
 
   smalltree->Branch("tree_simtrack_genVertexPos_X", &tree_simtrack_genVertexPos_X );		  
   smalltree->Branch("tree_simtrack_genVertexPos_Y", &tree_simtrack_genVertexPos_Y );		  
   smalltree->Branch("tree_simtrack_genVertexPos_Z", &tree_simtrack_genVertexPos_Z );		  
   smalltree->Branch("tree_simtrack_isRecoMatched",  &tree_simtrack_isRecoMatched  );
   
   
   smalltree->Branch("tree_simtrack_pca_dxy",&tree_simtrack_pca_dxy);
   smalltree->Branch("tree_simtrack_pca_dz", &tree_simtrack_pca_dz);
   smalltree->Branch("tree_simtrack_trkIdx", &tree_simtrack_trkIdx); 
   
   
   smalltree->Branch("tree_jet_pt"  , &tree_jet_pt);
   smalltree->Branch("tree_jet_eta" , &tree_jet_eta);
   smalltree->Branch("tree_jet_phi" , &tree_jet_phi);
   smalltree->Branch("tree_jet_vx"  , &tree_jet_vx);
   smalltree->Branch("tree_jet_vy" , &tree_jet_vy);
   smalltree->Branch("tree_jet_vz" , &tree_jet_vz);


   ////gen info 
   smalltree->Branch("tree_genParticle_pt"  , &tree_genParticle_pt);
   smalltree->Branch("tree_genParticle_eta" , &tree_genParticle_eta);
   smalltree->Branch("tree_genParticle_phi" , &tree_genParticle_phi);
   smalltree->Branch("tree_genParticle_charge" , &tree_genParticle_charge);
   smalltree->Branch("tree_genParticle_pdgId" , &tree_genParticle_pdgId);
   smalltree->Branch("tree_genParticle_vx"  , &tree_genParticle_vx);
   smalltree->Branch("tree_genParticle_vy" , &tree_genParticle_vy);
   smalltree->Branch("tree_genParticle_vz" , &tree_genParticle_vz);
   smalltree->Branch("tree_genParticle_mass" , &tree_genParticle_mass);
   smalltree->Branch("tree_genParticle_statusCode", &tree_genParticle_statusCode); 
   smalltree->Branch("tree_genParticle_mother" , &tree_genParticle_mother);

   //genJet info 
   smalltree->Branch("tree_genJet_pt"  , &tree_genJet_pt);
   smalltree->Branch("tree_genJet_eta" , &tree_genJet_eta);
   smalltree->Branch("tree_genJet_phi" , &tree_genJet_phi);
   smalltree->Branch("tree_genJet_vx"  , &tree_genJet_vx);
   smalltree->Branch("tree_genJet_vy" , &tree_genJet_vy);
   smalltree->Branch("tree_genJet_vz" , &tree_genJet_vz);
   smalltree->Branch("tree_genJet_mass", &tree_genJet_mass); 
   smalltree->Branch("tree_genJet_energy", &tree_genJet_energy); 
   smalltree->Branch("tree_genJet_status", &tree_genJet_status); 
   smalltree->Branch("tree_genJet_pdgId", &tree_genJet_pdgId); 



    //electrons info 
   smalltree->Branch("tree_electron_pt"  , &tree_electron_pt);
   smalltree->Branch("tree_electron_eta" , &tree_electron_eta);
   smalltree->Branch("tree_electron_phi" , &tree_electron_phi);
   smalltree->Branch("tree_electron_vx"  , &tree_electron_vx);
   smalltree->Branch("tree_electron_vy" , &tree_electron_vy);
   smalltree->Branch("tree_electron_vz" , &tree_electron_vz); 
   smalltree->Branch("tree_electron_energy", &tree_electron_energy); 



    //muons info 
   smalltree->Branch("tree_muon_pt"  , &tree_muon_pt);
   smalltree->Branch("tree_muon_eta" , &tree_muon_eta);
   smalltree->Branch("tree_muon_phi" , &tree_muon_phi);
   smalltree->Branch("tree_muon_vx"  , &tree_muon_vx);
   smalltree->Branch("tree_muon_vy" , &tree_muon_vy);
   smalltree->Branch("tree_muon_vz" , &tree_muon_vz); 
   smalltree->Branch("tree_muon_energy", &tree_muon_energy);



   runNumber = 0;
   eventNumber = 0;
   lumiBlock = 0;
   tree_bs_PosX= 0;
   tree_bs_PosY= 0;
   tree_bs_PosZ= 0;


   const bool tpRef = iConfig.getUntrackedParameter<bool>("trackingParticlesRef");
   const auto tpTag = iConfig.getUntrackedParameter<edm::InputTag>("trackingParticles");
   if(tpRef) {
     trackingParticleRefToken_ = consumes<TrackingParticleRefVector>(tpTag);
   }
   else {
     trackingParticleToken_ = consumes<TrackingParticleCollection>(tpTag);
   }
   
   
}



TrackingPerf::~TrackingPerf()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackingPerf::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   
   // 
   
   
   clearVariables();
 
   
   using namespace edm;
   using namespace reco;
   using namespace std;



   runNumber = iEvent.id().run();
   std::cout << "runNumber = " << runNumber << std::endl;
   eventNumber = iEvent.id().event();
   std::cout << "eventNumber = "<< eventNumber <<std::endl; 
   lumiBlock = iEvent.luminosityBlock();
  



    //// HANDLES /////


    // do we really need both here ???
   edm::Handle<  edm::View<reco::Track>  > TracksForRes;
   iEvent.getByToken(trackSrc_, TracksForRes);

   edm::Handle<edm::View<reco::Track> > tracksHandle;
   iEvent.getByToken(trackToken_, tracksHandle);
   const edm::View<reco::Track>& tracks = *tracksHandle;

   
   /*edm::ESHandle<TrackerGeometry> TG;
   iSetup.get<TrackerDigiGeometryRecord>().get(TG);
   const TrackerGeometry* theTrackerGeometry = TG.product();*/
    
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get( "TransientTrackBuilder", theB );
  
   edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
   iSetup.get<TransientRecHitRecord>().get(ttrhbuilder_,theTrackerRecHitBuilder);

   Handle<reco::BeamSpot> recoBeamSpotHandle;
   iEvent.getByToken(beamSpotToken_, recoBeamSpotHandle);


   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vertexToken_, vertices);

   edm::Handle<edm::View<reco::Jet> > jets;
   iEvent.getByToken(jetToken_,jets); 

   //edm::EDGetTokenT<std::vector<pat::MET> > metToken_;
   //iEvent.getByToken(metToken_, met);

    //not sure??
    //edm::Handle<edm::View<pat::Jet> > view_jets;
   //iEvent.getByToken(viewJetToken_,view_jets); 

   edm::Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genParticlesToken_,genParticles);

   edm::Handle<reco::GenJetCollection> genJets;
   iEvent.getByToken(genJetToken_,genJets);

   edm::Handle<GenEventInfoProduct> genEventInfo;
   iEvent.getByToken(genEventInfoToken_,genEventInfo);

   edm::Handle<LHEEventProduct> lheEventProduct;
   iEvent.getByToken(LHEEventProductToken_,lheEventProduct);

   edm::Handle<pat::PackedCandidateCollection> pfcands;
   iEvent.getByToken(pfcandsToken_,pfcands);

   //edm::Handle<edm::View<reco::GsfElectron> > electrons;
   //iEvent.getByToken(electronToken_, electrons);

    //edm::Handle<pat::ElectronCollection> electronsPAT;
    //iEvent.getByToken(electronPATToken_,electronsPAT);

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_,muons);

   //geometry & magnetif field 
   edm::ESHandle<TrackerGeometry> pDD;
   iSetup.get<TrackerDigiGeometryRecord>().get(pDD);

   iSetup.get<IdealMagneticFieldRecord>().get(bField); 
 
   edm::ESHandle<TrackerTopology> tTopoHandle;
   iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
   const TrackerTopology* const tTopo = tTopoHandle.product();


   //association between tracks 
   edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
   iEvent.getByToken(trackAssociatorToken_, theAssociator);
   const reco::TrackToTrackingParticleAssociator& associatorByHits = *theAssociator;


   //tracking Particles info and matching to reco 

   TrackingParticleRefVector tmpTP;
   const TrackingParticleRefVector *tmpTPptr = nullptr;
   edm::Handle<TrackingParticleCollection>  TPCollectionH;
   edm::Handle<TrackingParticleRefVector>  TPCollectionHRefVector;
 
   if(!trackingParticleToken_.isUninitialized()) {
     iEvent.getByToken(trackingParticleToken_, TPCollectionH);
     for(size_t i=0, size=TPCollectionH->size(); i<size; ++i) {
       tmpTP.push_back(TrackingParticleRef(TPCollectionH, i));
     }
     tmpTPptr = &tmpTP;
   }
   else {
     iEvent.getByToken(trackingParticleRefToken_, TPCollectionHRefVector);
     tmpTPptr = TPCollectionHRefVector.product();
   }
   const TrackingParticleRefVector& tpCollection = *tmpTPptr;
 
   
   // Fill mapping from Ref::key() to index
   TrackingParticleRefKeyToIndex tpKeyToIndex;
   for(size_t i=0; i<tpCollection.size(); ++i) {
     tpKeyToIndex[tpCollection[i].key()] = i;
   }
   
   
   edm::ESHandle<ParametersDefinerForTP> parametersDefinerH;
   iSetup.get<TrackAssociatorRecord>().get(parametersDefinerName_, parametersDefinerH);
   const ParametersDefinerForTP *parametersDefiner = parametersDefinerH.product();




    /// Some output numbers 

   //int nElectrons = 0; 
   //int nMuons = 0; 

   //// START FILLING THE INFORMATION FOR THE OBJECT 




   //////////////////////////////////
   //////////////////////////////////
   ////////      BS     /////////////
   //////////////////////////////////
   //////////////////////////////////

   BeamSpot const & bs = *recoBeamSpotHandle;

   tree_bs_PosX = bs.x0();
   tree_bs_PosY = bs.y0();
   tree_bs_PosZ = bs.z0();




   //////////////////////////////////
   //////////////////////////////////
   ////////   Tracks  /////////////
   //////////////////////////////////
   //////////////////////////////////
   ///Tracks and their association to other objects 

   edm::RefToBaseVector<reco::Track> trackRefs;
   size_t idxTrack = 0;



   std::map<size_t , int > trackToVextexMap;
   std::map<size_t , int > trackToJetMap;
   std::map<size_t , int > jetToTrackMap;

   ///// TRACK ASSOCIATION TO VERTICES AND JETS 


   int nRecoTracks = tracks.size(); 
   for(edm::View<Track>::size_type i=0; i<tracks.size(); ++i) 
   {
     //---------------------------
     //minimum selection on tracks
     if(tracks[i].pt() < 1 || fabs(tracks[i].eta()) > 2.4) continue;
     //double trackpT = tracks[i].pt();
     trackRefs.push_back(tracks.refAt(i));
    
     //---------------------------
     //loop on vertex to do trac-vertex association
     //---------------------------
     int idxvtx = 0;
     int thematchidx = -1;
     float dzmin = std::numeric_limits<float>::max();
     for(auto const & vertex : *vertices) {
       //---------------------------
       //loop on vertex to do track-vertex association;
       float dz = std::abs(tracks[i].dz( vertex.position()));
       if(dz < dzmin){
         dzmin = dz;
         thematchidx = idxvtx;
       }
       //---------------------------
       idxvtx++;
     }

     trackToVextexMap[idxTrack] = thematchidx;





     //---------------------------
     //loop on jets to do track-jet association
     //---------------------------
 
     int idxJet=0;
      bool found_match = false;
      for(unsigned int ij=0;ij<jets->size();ij++){
    
        const Jet& jet = jets->at(ij);
        TLorentzVector jet4m, track4m;
        jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
        track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
        if( jet4m.DeltaR(track4m) < 0.4){
          found_match = true;
          break;
        }
        else idxJet++;
      }
      if(found_match) trackToJetMap[idxTrack] = idxJet;
      else trackToJetMap[idxTrack] = -1;
      //---------------------------
      //     
      idxTrack++;
 

    } 



   //////////////////////////////////
   //////////////////////////////////
   ////////   vertices  /////////////
   //////////////////////////////////
   //////////////////////////////////


   int nPV = vertices->size(); 
   for(auto const & vertex : *vertices) {
     
     
     tree_vtx_PosX.push_back(vertex.x());
     tree_vtx_PosY.push_back(vertex.y());
     tree_vtx_PosZ.push_back(vertex.z());
     
     tree_vtx_PosXError.push_back(vertex.xError());
     tree_vtx_PosYError.push_back(vertex.yError());
     tree_vtx_PosZError.push_back(vertex.zError());

     nPV++; 
     
     
   }


   //////////////////////////////////
   //////////////////////////////////
   ////////   jets ///  /////////////
   //////////////////////////////////
   //////////////////////////////////
   
   int nJet = jets->size();
   //cout << "number of jets "<<nJet<<endl;  
   
   for(int ij=0;ij<nJet;ij++){
   
     const Jet& jet = jets->at(ij);
     tree_jet_pt.push_back(jet.pt());
     tree_jet_eta.push_back(jet.eta());
     tree_jet_phi.push_back(jet.phi());
     tree_jet_vx.push_back(jet.vx());
     tree_jet_vy.push_back(jet.vy());
     tree_jet_vz.push_back(jet.vz());
   }


   //////////////////////////////////
   //////////////////////////////////
   ////////   MET   /////////////////
   //////////////////////////////////
   //////////////////////////////////



    //////////////////////////////////
   //////////////////////////////////
   //////// gen Information /////////
   //////////////////////////////////
   //////////////////////////////////


   ///GEN PARTICLES 


   int nGenParticles = genParticles->size(); 
   //cout << "number of gen particles "<<nGenParticles<<endl; 
 
   for (auto const & genParticle : *genParticles)
   {
     tree_genParticle_pt.push_back(genParticle.pt());
     tree_genParticle_eta.push_back(genParticle.eta());
     tree_genParticle_phi.push_back(genParticle.phi());
     tree_genParticle_charge.push_back(genParticle.charge());
     tree_genParticle_pdgId.push_back(genParticle.pdgId());
     tree_genParticle_vx.push_back(genParticle.vx());
     tree_genParticle_vy.push_back(genParticle.vy());
     tree_genParticle_vz.push_back(genParticle.vz());
     
 
     tree_genParticle_mass.push_back(genParticle.mass());
     tree_genParticle_statusCode.push_back(genParticle.status()); 

     ///FIX ME : maybe add the pdg id of mother 
     //const Candidate * mom = genParticle.mother();
     //tree_genParticle_mother.push_back(mom); 
   
    }



    ///GEN JETS 

    int nGenJets = genJets->size();
    //cout << "number of gen Jets "<<nGenJets<<endl; 
 
    for (auto const & genJet : *genJets)
    {
      tree_genJet_pt.push_back(genJet.pt());
      tree_genJet_eta.push_back(genJet.eta());
      tree_genJet_phi.push_back(genJet.phi());
      tree_genJet_mass.push_back(genJet.mass());
      tree_genJet_energy.push_back(genJet.energy());
      tree_genJet_status.push_back(genJet.status());
      tree_genJet_pdgId.push_back(genJet.pdgId());
    }


    /// Electrons 

  /*  int nElectrons = electrons->size();


    for (auto const & electron : *electrons)
    {
      tree_electron_pt.push_back(electron.pt());
      tree_electron_eta.push_back(electron.eta());
      tree_electron_phi.push_back(electron.phi());
      tree_electron_vx.push_back(electron.vx());
      tree_electron_vy.push_back(electron.vy());
      tree_electron_vz.push_back(electron.vz());
      tree_electron_energy.push_back(electron.energy());
    }*/


    /// Muons 
    int nMuons = muons->size();


    for (auto const & muon : *muons)
    {
      tree_muon_pt.push_back(muon.pt());
      tree_muon_eta.push_back(muon.eta());
      tree_muon_phi.push_back(muon.phi());
      tree_muon_vx.push_back(muon.vx());
      tree_muon_vy.push_back(muon.vy());
      tree_muon_vz.push_back(muon.vz());
      tree_muon_energy.push_back(muon.energy());
    }




    ///OUTPUT MESSAGES 

    std::cout << "---------------------------"  << std::endl;
    std::cout << "number of reco tracks      "  << nRecoTracks <<std::endl;
    std::cout << "number of primary vertices "  << nPV <<std::endl;
    std::cout << "number of jets             "  << nJet <<std::endl;
    std::cout << "number of gen particles    "  << nGenParticles <<std::endl;
    std::cout << "number of gen jets         "  << nGenJets <<std::endl;
    //std::cout << "number of electrons        "  << nElectrons <<std::endl; 
    std::cout << "number of muons            "  << nMuons <<std::endl; 






    /////////////////
    //prepare assocoaio to tracks by hit
   reco::RecoToSimCollection recSimColl = associatorByHits.associateRecoToSim(trackRefs, tpCollection);
   

   int nTracks = 0; 
   int nUnmatchTrack_fromPU = 0;
   int nUnmatchTrack_fromPV = 0;
   int nUnmatchTrack= 0;
   int nPUTrack= 0;

   
   




   //---------------
   //loops on tracks 
   //---------------
   
   int idx_sicluster_first  = 0;
   int idx_pixcluster_first = 0;
   
   for(size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack) 
    {
     
   
     const auto& itTrack = trackRefs[iTrack];
     //------------------------
     //general track properties
     //------------------------
     //std::cout << "  " << itTrack->pt()  << std::endl;
     tree_track_charge.push_back(itTrack->charge());
     tree_track_pt.push_back(itTrack->pt());
     tree_track_eta.push_back(itTrack->eta());
     tree_track_phi.push_back(itTrack->phi());
     tree_track_nhits.push_back(itTrack->numberOfValidHits());
     tree_track_NChi2.push_back(itTrack->normalizedChi2());
     tree_track_outerPt.push_back(itTrack->outerPt()); 
     tree_track_vx.push_back(itTrack->vx()); 
     tree_track_vy.push_back(itTrack->vy()); 
     tree_track_vz.push_back(itTrack->vz()); 
     tree_track_firsthit_X.push_back(itTrack->innerPosition().X());  
     tree_track_firsthit_Y.push_back(itTrack->innerPosition().Y()); 
     tree_track_firsthit_Z.push_back(itTrack->innerPosition().Z()); 
     tree_track_firsthit_phi.push_back(itTrack->innerPosition().phi()); 


     if( itTrack->quality(reco::TrackBase::highPurity) ){tree_track_isHighQuality.push_back(true);} 
     	else {tree_track_isHighQuality.push_back(false);}
     if( itTrack->quality(reco::TrackBase::loose) )     {tree_track_isLoose.push_back(true);}	
     	else {tree_track_isLoose.push_back(false);}  
     if( itTrack->quality(reco::TrackBase::tight))      {tree_track_isTight.push_back(true);}
     	else {tree_track_isTight.push_back(false);}	  
     
     tree_track_dxy.push_back(  	    itTrack->dxy(bs.position()));
     tree_track_dxyError.push_back(	    itTrack->dxyError());
     tree_track_dz.push_back(		    itTrack->dz(bs.position()));
     tree_track_dzError.push_back(	    itTrack->dzError());
     tree_track_numberOfLostHits.push_back( itTrack->numberOfLostHits());
     tree_track_numberOfValidHits.push_back(itTrack->numberOfValidHits());
  
  
  
  
     tree_track_originalAlgo.push_back(itTrack->originalAlgo());
     tree_track_algo.push_back(itTrack->algo());
     tree_track_stopReason.push_back(itTrack->stopReason());


     //--------------------------------
     //general hit properties of tracks
     //--------------------------------
     
     const reco::HitPattern& hp = itTrack->hitPattern();
     tree_track_nPixel   .push_back(hp.numberOfValidPixelHits());
     tree_track_nStrip   .push_back(hp.numberOfValidStripHits());
      
     tree_track_numberOfValidPixelHits.push_back(hp.numberOfValidPixelHits());
     tree_track_numberOfValidStripHits.push_back(hp.numberOfValidStripHits());
     tree_track_numberOfValidStripTIBHits.push_back(hp.numberOfValidStripTIBHits());
     tree_track_numberOfValidStripTIDHits.push_back(hp.numberOfValidStripTIDHits());
     tree_track_numberOfValidStripTOBHits.push_back(hp.numberOfValidStripTOBHits());
     tree_track_numberOfValidStripTECHits.push_back(hp.numberOfValidStripTECHits());
     tree_track_numberOfValidPixelBarrelHits.push_back(hp.numberOfValidPixelBarrelHits());
     tree_track_numberOfValidPixelEndcapHits.push_back(hp.numberOfValidPixelEndcapHits());
     tree_track_trackerLayersWithMeasurement.push_back(hp.trackerLayersWithMeasurement());
     tree_track_pixelLayersWithMeasurement.push_back(hp.pixelLayersWithMeasurement());
     tree_track_stripTECLayersWithMeasurement.push_back(hp.stripTECLayersWithMeasurement() );
     tree_track_stripTIBLayersWithMeasurement.push_back(hp.stripTIBLayersWithMeasurement());
     tree_track_stripTIDLayersWithMeasurement.push_back(hp.stripTIDLayersWithMeasurement());
     tree_track_stripTOBLayersWithMeasurement.push_back(hp.stripTOBLayersWithMeasurement());
     
     int hitPixelLayer = 0;
     if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) )  hitPixelLayer += 1;
     if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) )  hitPixelLayer += 10;
     if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) )  hitPixelLayer += 100;
     if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) )  hitPixelLayer += 1000;
     
     if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) )  hitPixelLayer += 2;
     if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) )  hitPixelLayer += 20;
     if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) )  hitPixelLayer += 200;
     
     
     
     tree_track_hasValidHitInPixelLayer.push_back(hitPixelLayer);

     //----------------------------
     //matching to simulated tracks
     //----------------------------
     
     int nSimHits = 0;
     bool isSimMatched = false;
     std::vector<int> tpIdx;
     std::vector<float> sharedFraction;
     std::vector<float> tpChi2;
     
     //initialized values for trackingParticle
     
     double genVertexPos_X = -1000;
     double genVertexPos_Y = -1000;
     double genVertexPos_Z = -1000;
     
     double  simtrack_charge		    =-1000;
     double  simtrack_pt		    =-1;
     double  simtrack_eta		    =-1000;
     double  simtrack_phi		    =-1000;
     bool    simtrack_longLived 	    = false;
     //int     simtrack_matchedHit	    = 0;
     int     simtrack_pdgId		    = -1000;
     int     simtrack_numberOfTrackerHits   = 0;
     int     simtrack_numberOfTrackerLayers = 0;
     double  simtrack_mass = 0;
     int     simtrack_status= -1000;
     
     
     auto foundTPs = recSimColl.find(itTrack);
     if (foundTPs != recSimColl.end()) 
     {
         //if (!foundTPs->val.empty()) {
        isSimMatched = true;
        TrackingParticleRef tpr = foundTPs->val[0].first;
        
        nSimHits = tpr->numberOfTrackerHits();
         
   
	    simtrack_charge	      = tpr->charge();  	  
	    simtrack_pt		      = tpr->pt();		  
	    simtrack_eta		      = tpr->eta();	  
	    simtrack_phi		      = tpr->phi();	  
	    simtrack_longLived	      = tpr->longLived();   
	    //simtrack_matchedHit	      = tpr->matchedHit();    
	    simtrack_pdgId 		      = tpr->pdgId();	 
	    simtrack_numberOfTrackerHits	      = tpr->numberOfTrackerHits(); 
	    simtrack_numberOfTrackerLayers       = tpr->numberOfTrackerLayers();
	    simtrack_mass  		      = tpr->mass();
	    simtrack_status		      = tpr->status();
         
   
        //determiine x,y,z position of the genVertex which produced the associated simtrack
	    genVertexPos_X 		 = tpr->vx();
	    genVertexPos_Y 		 = tpr->vy();
	    genVertexPos_Z 		 = tpr->vz();



     }  


     tree_track_nSimHits		      .push_back(nSimHits); 
     tree_track_isSimMatched		      .push_back(isSimMatched);
     

     nTracks++; 

     if( itTrack->pt() < 1 ) std::cout << " there is a pt problem " << itTrack->pt() << std::endl;
     if( !isSimMatched &&  trackToVextexMap[iTrack] != 0 ) nUnmatchTrack_fromPU++;
     if( !isSimMatched &&  trackToVextexMap[iTrack] == 0 ) nUnmatchTrack_fromPV++;
     if( !isSimMatched                                   ) nUnmatchTrack++;
     if( trackToVextexMap[iTrack] != 0                   ) nPUTrack++;




     tree_track_recoVertex_idx.push_back(trackToVextexMap[iTrack]);
     tree_track_recoJet_idx.push_back(trackToJetMap[iTrack]);
     
     tree_track_simtrack_charge 	      .push_back(simtrack_charge);	     
     tree_track_simtrack_pt		      .push_back(simtrack_pt);  		     
     tree_track_simtrack_eta		      .push_back(simtrack_eta); 	     
     tree_track_simtrack_phi		      .push_back(simtrack_phi); 		     
     tree_track_simtrack_longLived	      .push_back(simtrack_longLived);		
     //tree_track_simtrack_matchedHit	      .push_back(simtrack_matchedHit);         
     tree_track_simtrack_pdgId  	      .push_back(simtrack_pdgId);      
     tree_track_simtrack_numberOfTrackerHits  .push_back(simtrack_numberOfTrackerHits); 
     tree_track_simtrack_numberOfTrackerLayers.push_back(simtrack_numberOfTrackerLayers); 
     tree_track_simtrack_mass		      .push_back(simtrack_mass        );	
     tree_track_simtrack_status 	      .push_back(simtrack_status      );	
     
     tree_track_genVertexPos_X  	      .push_back(genVertexPos_X);	     
     tree_track_genVertexPos_Y  	      .push_back(genVertexPos_Y);	     
     tree_track_genVertexPos_Z  	      .push_back(genVertexPos_Z);	


     //----------------------------
     //look a recHit and check info
     //----------------------------
     
     tree_track_idxSiClusterFirst .  push_back(idx_sicluster_first );
     tree_track_idxPixClusterFirst.  push_back(idx_pixcluster_first);
     
     
     for(auto i=itTrack->recHitsBegin(); i!=itTrack->recHitsEnd(); i++) 
     {
       
       TransientTrackingRecHit::RecHitPointer hit=(*theTrackerRecHitBuilder).build(&**i ); 
       
       DetId hitId = hit->geographicalId();
       
       if(hitId.det() != DetId::Tracker) continue;
       if (!(hit->isValid()) )  	 continue;
       
       int subDet = hitId.subdetId();
       
       
       //const BaseTrackerRecHit* bhit = dynamic_cast<const BaseTrackerRecHit*>(&*hit);
       
       //const auto& clusterRef = bhit->firstClusterRef();
       const SiStripRecHit2D* SiStriphit2D = dynamic_cast<const SiStripRecHit2D*>((*hit).hit());
       
       const SiStripRecHit1D* SiStriphit1D = dynamic_cast<const SiStripRecHit1D*>((*hit).hit());
       
       
       if(  SiStriphit2D!=0 || SiStriphit1D!=0 )
       {  
           //const TransientTrackingRecHit::ConstRecHitPointer theTTrechit = //(*itm).recHit();
	        const SiStripCluster* si_cluster = 0;
    
	        if (SiStriphit2D!=0){
	          si_cluster = &*(SiStriphit2D->cluster());
	        }
	        if (SiStriphit1D!=0){
	          si_cluster = &*(SiStriphit1D->cluster());
	        }
    
	        SiStripClusterInfo *clusterInfo = 0;
	        if(useCluster_)  clusterInfo = new SiStripClusterInfo( *si_cluster, iSetup,  hitId ); 
    
    
	        int Cluster_WheelSide = 0;
	        int Cluster_detID = hitId;
	        int Cluster_PetalSide = 0;
	        int Cluster_subDet = -1;
	        int Cluster_LayerNbr = -1;
	        //determine subdte id
	        if(subDet == SiStripDetId::TIB){
	          Cluster_subDet = 0;
	          Cluster_LayerNbr = tTopo->tobLayer(hitId.rawId());
	        } 
	        if(subDet == SiStripDetId::TOB){
	          Cluster_subDet = 1;
	          Cluster_LayerNbr =tTopo->tibLayer( hitId.rawId());
	        } 
	        if(subDet == SiStripDetId::TID){
	          Cluster_subDet = 2;
	          Cluster_WheelSide =tTopo->tidSide(hitId.rawId());
	          Cluster_LayerNbr =tTopo->tidWheel(hitId.rawId());
	        } 
	        if(subDet == SiStripDetId::TEC){
	          Cluster_subDet = 3; 
	          Cluster_WheelSide = tTopo->tecSide(hitId.rawId());
	          Cluster_LayerNbr =tTopo->tecWheel( hitId.rawId());
	          if(tTopo->tecIsFrontPetal(	   hitId.rawId()))  Cluster_PetalSide =  1;
	          else						    Cluster_PetalSide = -1;
    
	        }


            tree_SiCluster_WheelSide.push_back(Cluster_WheelSide);
	        tree_SiCluster_detID.push_back(	 Cluster_detID);
	        tree_SiCluster_PetalSide.push_back(Cluster_PetalSide);
	        tree_SiCluster_subDet.push_back(   Cluster_subDet);
	        tree_SiCluster_LayerNbr.push_back( Cluster_LayerNbr);


            if(useCluster_){
	        tree_SiCluster_SoverN  .push_back( (*clusterInfo).signalOverNoise());
	        tree_SiCluster_noise	  .push_back( (*clusterInfo).noiseRescaledByGain());
	        tree_SiCluster_charge  .push_back( (*clusterInfo).charge());
	        tree_SiCluster_width	  .push_back( si_cluster->amplitudes().size());
	        tree_SiCluster_barycenter   .push_back( si_cluster->barycenter());
	        }
	   
	        /*tree_SiCluster_globX .push_back( hit->globalPosition().x());
	        tree_SiCluster_globY .push_back( hit->globalPosition().y());
	        tree_SiCluster_globZ .push_back( hit->globalPosition().z());*/
	 
	        tree_SiCluster_globX .push_back( 0 );
	        tree_SiCluster_globY .push_back( 0 );
	        tree_SiCluster_globZ .push_back( 0 );
	 
	 
	        idx_sicluster_first++;
	 
        } 


        const SiPixelRecHit* pixelHits= dynamic_cast<const SiPixelRecHit*>((*hit).hit());
        if(  pixelHits != 0  )  {
            const SiPixelCluster* pix_cluster = 0;
            if(pixelHits != 0){
                if(useCluster_) pix_cluster = &*(pixelHits->cluster());
	            DetId clusterDetId(hitId);
                //	const PixelGeomDetUnit * pixeldet = (const PixelGeomDetUnit*) theTrackerGeometry->idToDetUnit(clusterDetId);
	            int PixCluster_LayerNbr = -1;


                if(subDet ==  PixelSubdetector::PixelBarrel){
	             tree_PixCluster_isBPix.push_back(true);
	             PXBDetId pdetId = PXBDetId(hitId);
	             PixCluster_LayerNbr  = pdetId.layer(); 
                }



                if(subDet ==  PixelSubdetector::PixelEndcap){
                 tree_PixCluster_isBPix.push_back(false);
                 PXFDetId pdetId = PXFDetId(hitId);
                 PixCluster_LayerNbr  = pdetId.disk(); 
	            }



                tree_PixCluster_LayerNbr.push_back(PixCluster_LayerNbr);  
	   
	            if(useCluster_) tree_PixCluster_charge.push_back(pix_cluster->charge());


	            //std::cout <<   pixelHits->globalPosition().x()<< std::endl ;

	            tree_PixCluster_detID.push_back(hitId);
	            /*tree_PixCluster_globX .push_back( pixelHits->globalPosition().x());
	            tree_PixCluster_globY .push_back( pixelHits->globalPosition().y());
	            tree_PixCluster_globZ .push_back( pixelHits->globalPosition().z());*/

	            tree_PixCluster_globX .push_back( 0 );
	            tree_PixCluster_globY .push_back( 0 );
	            tree_PixCluster_globZ .push_back( 0 );



	            idx_pixcluster_first++;
	        }
            
        }

     } //end loop on rechits  

     tree_track_idxSiClusterLast  .push_back(idx_sicluster_first);
     tree_track_idxPixClusterLast .push_back(idx_pixcluster_first);

    }//end loop on tracks




    std::cout << "---------------------"  << std::endl;
    std::cout << "nUnmatchTrack_fromPU " << nUnmatchTrack_fromPU << std::endl;
    std::cout << "nUnmatchTrack_fromPV " << nUnmatchTrack_fromPV << std::endl;
    std::cout << "nUnmatchTrack        " << nUnmatchTrack<< std::endl;
    std::cout << "nTrack               " << tree_track_charge.size() << std::endl;
    std::cout << "nPUTrack             " << nPUTrack<< std::endl;
    std::cout << "nNonPUTrack          " << tree_track_charge.size()- nPUTrack<< std::endl;

    //------------------------------
   //loops on tracking particles 
   //------------------------------
   
   
   
   reco::SimToRecoCollection simRecColl = associatorByHits.associateSimToReco(trackRefs, tpCollection);
   
   for(const TrackingParticleRef& tp: tpCollection) {
     
     tree_simtrack_simtrack_charge  		      .push_back( tp->charge());	 
     tree_simtrack_simtrack_pt		      .push_back( tp->pt());		 
     tree_simtrack_simtrack_eta		      .push_back( tp->eta());	 
     tree_simtrack_simtrack_phi		      .push_back( tp->phi());	 
     tree_simtrack_simtrack_longLived	      .push_back( tp->longLived());   
     //tree_simtrack_simtrack_matchedHit	      .push_back( tp->matchedHit());	
     tree_simtrack_simtrack_pdgId		      .push_back( tp->pdgId());    
     tree_simtrack_simtrack_numberOfTrackerHits   .push_back( tp->numberOfTrackerHits()); 
     tree_simtrack_simtrack_numberOfTrackerLayers	      .push_back( tp->numberOfTrackerLayers());
     tree_simtrack_simtrack_mass			      .push_back( tp->mass());
     tree_simtrack_simtrack_status  		      .push_back( tp->status());
     
     tree_simtrack_genVertexPos_X			 .push_back( tp->vx());
     tree_simtrack_genVertexPos_Y			 .push_back( tp->vy());
     tree_simtrack_genVertexPos_Z			 .push_back( tp->vz());
     
     bool isRecoMatched = false;
     
     std::vector<int> tkIdx;
     auto foundTracks = simRecColl.find(tp);
     if(foundTracks != simRecColl.end()) {
       isRecoMatched = true;
       for(const auto trackQuality: foundTracks->val) {
	 
	 tkIdx.push_back(trackQuality.first.key());
       }
       
       
     }
     tree_simtrack_isRecoMatched.push_back(isRecoMatched);
     //Calcualte the impact parameters w.r.t. PCA
     tree_simtrack_trkIdx.push_back(tkIdx);
     
     TrackingParticle::Vector momentum = parametersDefiner->momentum(iEvent,iSetup,tp);
     TrackingParticle::Point vertex = parametersDefiner->vertex(iEvent,iSetup,tp);
     auto dxySim = TrackingParticleIP::dxy(vertex, momentum, bs.position());
     auto dzSim  = TrackingParticleIP::dz(vertex, momentum, bs.position());
     
     tree_simtrack_pca_dxy      .push_back(dxySim);
     tree_simtrack_pca_dz	    .push_back(dzSim);
    }



    smalltree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
TrackingPerf::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrackingPerf::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackingPerf::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrackingPerf);




void TrackingPerf::clearVariables() {
 //-----------------------
  //fill the tree per track
   //tree_track_nclusters.clear();
   tree_track_charge.clear();
   tree_track_pt.clear();
   tree_track_outerPt.clear(); 
   tree_track_eta.clear();
   tree_track_phi.clear();
   tree_track_nhits.clear();
   tree_track_NChi2.clear();
   //tree_track_Quality.clear();
   tree_track_isHighQuality.clear();
   tree_track_isLoose.clear();
   tree_track_isTight.clear();
   
   tree_track_dxy.clear(); 
   tree_track_dxyError.clear();
   tree_track_dz.clear();
   tree_track_dzError.clear(); 
   tree_track_numberOfLostHits.clear();
   tree_track_numberOfValidHits.clear();
   tree_track_originalAlgo.clear(); 
   tree_track_algo.clear(); 
   tree_track_stopReason.clear();
   tree_track_nSimHits.clear(); 
   tree_track_isSimMatched.clear();
   tree_track_nPixel.clear();
   tree_track_nStrip.clear();
   tree_track_vx.clear(); 
   tree_track_vy.clear(); 
   tree_track_vz.clear(); 
   tree_track_firsthit_X.clear(); 
   tree_track_firsthit_Y.clear(); 
   tree_track_firsthit_Z.clear(); 
   tree_track_firsthit_phi.clear(); 
   
   
   
   tree_PixCluster_isBPix.clear();
   tree_PixCluster_LayerNbr.clear();
   tree_PixCluster_detID.clear();
   
   
   tree_track_nSimHits.clear();			 
   tree_track_isSimMatched.clear();		 
   tree_track_simtrack_charge.clear();		 
   tree_track_simtrack_pt.clear();  		 
   tree_track_simtrack_eta.clear();		 
   tree_track_simtrack_phi.clear(); 		 
   tree_track_simtrack_longLived.clear();		 
   //tree_track_simtrack_matchedHit .clear(); 	 
   tree_track_simtrack_pdgId.clear();		 
   tree_track_simtrack_numberOfTrackerHits.clear();
   tree_track_simtrack_numberOfTrackerLayers.clear();
   tree_track_simtrack_mass.clear();		 
   tree_track_simtrack_status.clear();		 

   tree_track_genVertexPos_X.clear();		 
   tree_track_genVertexPos_Y.clear();		 
   tree_track_genVertexPos_Z.clear();

   tree_track_recoVertex_idx.clear();
   tree_track_recoJet_idx.clear();
   tree_track_idxSiClusterFirst.clear();
   tree_track_idxPixClusterFirst.clear();
   tree_track_idxSiClusterLast.clear() ;
   tree_track_idxPixClusterLast.clear();


   tree_track_numberOfValidPixelHits.clear();
   tree_track_numberOfValidStripHits.clear();
   tree_track_numberOfValidStripTIBHits.clear();
   tree_track_numberOfValidStripTIDHits.clear();
   tree_track_numberOfValidStripTOBHits.clear();
   tree_track_numberOfValidStripTECHits.clear();
   tree_track_numberOfValidPixelBarrelHits.clear();
   tree_track_numberOfValidPixelEndcapHits.clear();
   tree_track_hasValidHitInPixelLayer.clear();
   tree_track_trackerLayersWithMeasurement.clear();
   tree_track_pixelLayersWithMeasurement.clear();
   tree_track_stripTECLayersWithMeasurement .clear();
   tree_track_stripTIBLayersWithMeasurement.clear();
   tree_track_stripTIDLayersWithMeasurement.clear();
   tree_track_stripTOBLayersWithMeasurement.clear();


   tree_SiCluster_subDet.clear();
   tree_SiCluster_PetalSide.clear();
   tree_SiCluster_LayerNbr.clear();
   tree_SiCluster_WheelSide.clear();
   tree_SiCluster_charge.clear();
   tree_SiCluster_SoverN.clear();
   tree_SiCluster_noise.clear();
   tree_SiCluster_width.clear();
   tree_SiCluster_barycenter.clear();
   tree_SiCluster_detID.clear();

   tree_SiCluster_locX.clear();
   tree_SiCluster_locY.clear();
   tree_SiCluster_tsosx.clear();
   tree_SiCluster_tsosy.clear();
   tree_SiCluster_globX.clear();
   tree_SiCluster_globY.clear();
   tree_SiCluster_globZ.clear();
   tree_SiCluster_tsosglobX.clear();
   tree_SiCluster_tsosglobY.clear();
   tree_SiCluster_tsosglobZ.clear();

   tree_PixCluster_LayerNbr.clear();
   tree_PixCluster_charge.clear();
   tree_PixCluster_isBPix.clear();
   tree_PixCluster_detID.clear();

   tree_PixCluster_locX.clear();
   tree_PixCluster_locY.clear();
   tree_PixCluster_globX.clear();
   tree_PixCluster_globY.clear();
   tree_PixCluster_globZ.clear();
   tree_PixCluster_tsosglobX.clear();
   tree_PixCluster_tsosglobY.clear();
   tree_PixCluster_tsosglobZ.clear();
   tree_PixCluster_tsosX.clear();
   tree_PixCluster_tsosY.clear();

   tree_Strips_stripCharges.clear();
   tree_Strips_stripGains .clear();
   tree_Strips_stripNoises.clear();
   tree_Strips_stripQualitiesBad .clear();

   tree_simtrack_simtrack_charge.clear();
   tree_simtrack_simtrack_pt.clear();
   tree_simtrack_simtrack_eta.clear();
   tree_simtrack_simtrack_phi.clear();
   tree_simtrack_simtrack_longLived.clear();




   tree_simtrack_simtrack_pdgId.clear();
   tree_simtrack_simtrack_numberOfTrackerHits.clear();
   tree_simtrack_simtrack_numberOfTrackerLayers.clear();
   tree_simtrack_simtrack_mass.clear();
   tree_simtrack_simtrack_status.clear();

   tree_simtrack_genVertexPos_X.clear();
   tree_simtrack_genVertexPos_Y.clear();
   tree_simtrack_genVertexPos_Z .clear();


   tree_track_recoVertex_idx.clear();
   tree_track_recoJet_idx.clear();


   tree_simtrack_isRecoMatched.clear();
   tree_simtrack_pca_dxy.clear();
   tree_simtrack_pca_dz .clear();
   tree_simtrack_trkIdx  .clear();

   tree_vtx_PosX.clear();
   tree_vtx_PosY.clear();
   tree_vtx_PosZ.clear();


   tree_vtx_PosXError.clear();
   tree_vtx_PosYError.clear();
   tree_vtx_PosZError.clear();

   tree_jet_pt.clear();
   tree_jet_eta.clear();
   tree_jet_phi.clear();
   tree_jet_vx.clear();
   tree_jet_vy.clear();
   tree_jet_vz.clear();


   tree_genParticle_pt.clear();  		 
   tree_genParticle_eta.clear();		 
   tree_genParticle_phi.clear(); 
   tree_genParticle_charge.clear();  		 
   tree_genParticle_pdgId.clear();		 
   tree_genParticle_vx.clear();
   tree_genParticle_vx.clear();
   tree_genParticle_vy.clear();
   tree_genParticle_vz.clear();
   tree_genParticle_mass.clear();
   tree_genParticle_statusCode.clear(); 
   tree_genParticle_mother.clear(); 


   tree_genJet_pt.clear();
   tree_genJet_eta.clear();
   tree_genJet_phi.clear();
   tree_genJet_vx.clear();
   tree_genJet_vy.clear();
   tree_genJet_vz.clear();
   tree_genJet_mass.clear();
   tree_genJet_energy.clear();
   tree_genJet_status.clear();
   tree_genJet_pdgId.clear();



   tree_electron_pt.clear();
   tree_electron_eta.clear(); 
   tree_electron_phi.clear();
   tree_electron_vx.clear();
   tree_electron_vy.clear();
   tree_electron_vz.clear();
   tree_electron_energy.clear();


   tree_muon_pt.clear();
   tree_muon_eta.clear();
   tree_muon_phi.clear();
   tree_muon_vx.clear();
   tree_muon_vy.clear();
   tree_muon_vz.clear();
   tree_muon_energy.clear();








} 























     














