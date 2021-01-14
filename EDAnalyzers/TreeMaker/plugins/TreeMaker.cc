// -*- C++ -*-
//
// Package:    EDAnalyzers/TreeMaker
// Class:      TreeMaker
//
/**\class TreeMaker TreeMaker.cc EDAnalyzers/TreeMaker/plugins/TreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//         Created:  Sat, 11 May 2019 13:14:55 GMT
//
//

// system include files
#include <memory>

// user include files

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDefs.h"
#include "DataFormats/Common/interface/MapOfVectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/FWLite/interface/ESHandle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "EDAnalyzers/TreeMaker/interface/Common.h"
#include "EDAnalyzers/TreeMaker/interface/Constants.h"
#include "EDAnalyzers/TreeMaker/interface/TreeOutputInfo.h"

#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <Compression.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVectorD.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

double HGCal_minEta = 1.479;
double HGCal_maxEta = 3.1;

double el_minPt = 10;    //15;
double el_maxPt = 99999; //30;

double ph_minPt = 10;    //15;
double ph_maxPt = 99999; //30;

double _largeVal = 999999999;

class TreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
    explicit TreeMaker(const edm::ParameterSet &);
    ~TreeMaker();

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
    virtual void endJob() override;

    hgcal::RecHitTools recHitTools;

    int minLayer;
    int maxLayer;

    TreeOutputInfo::TreeOutput *treeOutput;
    void WriteParticleKinematicsToTree(
        reco::GenParticle &part,
        std::vector<CLHEP::HepLorentzVector> &,
        std::vector<CLHEP::HepLorentzVector> &);

    // My stuff //
    bool debug;
    bool isGunSample;

    bool storeSimHit;
    bool storeRecHit;

    // GenEventInfoProduct //
    edm::EDGetTokenT<GenEventInfoProduct> tok_generator;

    // Gen particles //
    edm::EDGetTokenT<std::vector<reco::GenParticle>> tok_genParticle;

    // Pileup //
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> tok_pileup;

    // Rho //
    edm::EDGetTokenT<double> tok_rho;

    // SimHits //
    edm::EDGetTokenT<std::vector<PCaloHit>> tok_HGCEESimHit;

    // SimClusters //
    edm::EDGetTokenT<std::vector<SimCluster>> tok_simCluster;

    // RecHits //
    edm::EDGetTokenT<std::vector<reco::PFRecHit>> tok_PFRecHit;
    edm::EDGetTokenT<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> tok_HGCEERecHit;
    edm::EDGetTokenT<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> tok_HGCHEFRecHit;
    edm::EDGetTokenT<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> tok_HGCHEBRecHit;

    // Calo particles //
    edm::EDGetTokenT<std::vector<CaloParticle>> tok_caloParticle;
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
TreeMaker::TreeMaker(const edm::ParameterSet &iConfig)
{
    usesResource("TFileService");
    edm::Service<TFileService> fs;

    // Compression
    //fs->file().SetCompressionAlgorithm(ROOT::kLZMA);
    //fs->file().SetCompressionLevel(8);

    treeOutput = new TreeOutputInfo::TreeOutput("tree", fs);

    minLayer = +9999;
    maxLayer = -9999;

    // My stuff //
    debug = iConfig.getParameter<bool>("debug");
    isGunSample = iConfig.getParameter<bool>("isGunSample");

    storeSimHit = iConfig.getParameter<bool>("storeSimHit");
    storeRecHit = iConfig.getParameter<bool>("storeRecHit");

    // GenEventInfoProduct //
    tok_generator = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("label_generator"));

    // Gen particles //
    tok_genParticle = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("label_genParticle"));

    // Pileup //
    tok_pileup = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("label_pileup"));

    // Rho //
    tok_rho = consumes<double>(iConfig.getParameter<edm::InputTag>("label_rho"));

    // SimHits //
    tok_HGCEESimHit = consumes<std::vector<PCaloHit>>(iConfig.getParameter<edm::InputTag>("label_HGCEESimHit"));

    // RecHits //
    tok_PFRecHit = consumes<std::vector<reco::PFRecHit>>(iConfig.getParameter<edm::InputTag>("label_PFRecHit"));

    tok_HGCEERecHit = consumes<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>>(iConfig.getParameter<edm::InputTag>("label_HGCEERecHit"));
    tok_HGCHEFRecHit = consumes<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>>(iConfig.getParameter<edm::InputTag>("label_HGCHEFRecHit"));
    tok_HGCHEBRecHit = consumes<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>>(iConfig.getParameter<edm::InputTag>("label_HGCHEBRecHit"));

    // SimClusters //
    tok_simCluster = consumes<std::vector<SimCluster>>(iConfig.getParameter<edm::InputTag>("label_simCluster"));

    // Calo particles //
    tok_caloParticle = consumes<std::vector<CaloParticle>>(iConfig.getParameter<edm::InputTag>("label_caloParticle"));
}

TreeMaker::~TreeMaker()
{
    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
    delete treeOutput;
}

//
// member functions
//
void TreeMaker::WriteParticleKinematicsToTree(
    reco::GenParticle &part,
    std::vector<CLHEP::HepLorentzVector> &v_genEl_4mom,
    std::vector<CLHEP::HepLorentzVector> &v_genPh_4mom)
{
    int pdgId = part.pdgId();
    int status = part.status();

    // Check for eleectrons and (hard/promt) photons, otherwise skip
    bool validHardEl = abs(pdgId) == 11 && ((isGunSample && status == 1) ||
                                            (!isGunSample && part.isHardProcess()));
    bool validHardPh = abs(pdgId) == 22 && ((isGunSample && status == 1) ||
                                            (!isGunSample && part.isHardProcess()));
    bool validPromtPh = abs(pdgId) == 22 && Common::isPromptPhoton(part);
    if (!(validHardEl || validHardPh || validPromtPh))
    {
        return;
    }

    double &maxPt = (abs(pdgId) == 11) ? el_maxPt : ph_maxPt;
    double &minPt = (abs(pdgId) == 11) ? el_minPt : ph_minPt;

    // Check if the particle is in the HGCal, otherwise skip
    bool ptetaCut = fabs(part.eta()) > HGCal_minEta && fabs(part.eta()) < HGCal_maxEta && part.pt() > minPt && part.pt() < maxPt;
    if (!ptetaCut)
    {
        return;
    }

    printf("Gen part found: pdgId %d E %0.2f, pT %0.2f, eta %+0.2f, pz %+0.2f \n",
           pdgId, part.energy(), part.pt(), part.eta(), part.pz());

    std::vector<CLHEP::HepLorentzVector> &v_gen_4mom_ref = (abs(pdgId) == 11) ? v_genEl_4mom : v_genPh_4mom;

    CLHEP::HepLorentzVector gen_4mom;
    gen_4mom.setT(part.energy());
    gen_4mom.setX(part.px());
    gen_4mom.setY(part.py());
    gen_4mom.setZ(part.pz());

    v_gen_4mom_ref.push_back(gen_4mom);

    //Attach the properties of the patricle to the relevant vector in the tree.
    if ((pdgId == 11))
    {
        treeOutput->v_genEl_E.push_back(gen_4mom.e());
        treeOutput->v_genEl_px.push_back(gen_4mom.px());
        treeOutput->v_genEl_py.push_back(gen_4mom.py());
        treeOutput->v_genEl_pz.push_back(gen_4mom.pz());
        treeOutput->v_genEl_pT.push_back(gen_4mom.perp());
        treeOutput->v_genEl_eta.push_back(gen_4mom.eta());
        treeOutput->v_genEl_phi.push_back(gen_4mom.phi());

        treeOutput->genEl_n++;
    }
    else if (pdgId == 22)
    {
        treeOutput->v_genPh_E.push_back(gen_4mom.e());
        treeOutput->v_genPh_px.push_back(gen_4mom.px());
        treeOutput->v_genPh_py.push_back(gen_4mom.py());
        treeOutput->v_genPh_pz.push_back(gen_4mom.pz());
        treeOutput->v_genPh_pT.push_back(gen_4mom.perp());
        treeOutput->v_genPh_eta.push_back(gen_4mom.eta());
        treeOutput->v_genPh_phi.push_back(gen_4mom.phi());

        treeOutput->genPh_n++;
    }
}

// ------------ method called for each event  ------------
void TreeMaker::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    using namespace edm;

    long long eventNumber = iEvent.id().event();
    printf("Event %llu \n", eventNumber);

    treeOutput->clear();

    // Load the geometry
    edm::ESHandle<CaloGeometry> geom;
    iSetup.get<CaloGeometryRecord>().get(geom);
    recHitTools.setGeometry(*(geom.product()));

    //////////////////// Run info ////////////////////
    treeOutput->runNumber = iEvent.id().run();
    treeOutput->eventNumber = iEvent.id().event();
    treeOutput->luminosityNumber = iEvent.id().luminosityBlock();
    treeOutput->bunchCrossingNumber = iEvent.bunchCrossing();

    // HGCal Topology
    edm::ESHandle<HGCalTopology> handle_topo_HGCalEE;
    iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", handle_topo_HGCalEE);

    if (!handle_topo_HGCalEE.isValid())
    {
        printf("Error: Invalid HGCalEE topology. \n");
        exit(EXIT_FAILURE);
    }

    // const auto &topo_HGCalEE = *handle_topo_HGCalEE;

    //////////////////// GenEventInfoProduct ////////////////////
    edm::Handle<GenEventInfoProduct> generatorHandle;
    iEvent.getByToken(tok_generator, generatorHandle);
    GenEventInfoProduct generator = *generatorHandle;

    printf("[%llu] Gen. evt. wt. %0.4g \n", eventNumber, generator.weight());

    //////////////////// Gen particle ////////////////////
    edm::Handle<std::vector<reco::GenParticle>> v_genParticle;
    iEvent.getByToken(tok_genParticle, v_genParticle);

    std::vector<CLHEP::HepLorentzVector> v_genEl_4mom;
    std::vector<CLHEP::HepLorentzVector> v_genPh_4mom;

    // Iterate over the generated Particles, filter for valid ones in the HGCal and write the kinematics to the tree
    for (int iPart = 0; iPart < (int)v_genParticle->size(); iPart++)
    {
        reco::GenParticle part = v_genParticle->at(iPart);
        WriteParticleKinematicsToTree(part, v_genEl_4mom, v_genPh_4mom);
    }

    // Pileup
    edm::Handle<std::vector<PileupSummaryInfo>> pileUps_reco;
    iEvent.getByToken(tok_pileup, pileUps_reco);
    treeOutput->pileup_n = Common::getPileup(pileUps_reco);

    // Rho
    edm::Handle<double> handle_rho;
    iEvent.getByToken(tok_rho, handle_rho);
    double rho = *handle_rho;

    treeOutput->rho = rho;

    // SimHit dictionary
    edm::Handle<std::vector<PCaloHit>> v_HGCEESimHit;
    iEvent.getByToken(tok_HGCEESimHit, v_HGCEESimHit);

    std::map<DetId, const PCaloHit *> m_simHit;

    int nSimHit = v_HGCEESimHit->size();
    for (int iSimHit = 0; iSimHit < nSimHit; iSimHit++)
    {
        const PCaloHit *simHit = &(v_HGCEESimHit->at(iSimHit));

        DetId detId(simHit->id());

        m_simHit[detId] = simHit;
    }

    // RecHit dictionary
    edm::Handle<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> v_HGCEERecHit;
    iEvent.getByToken(tok_HGCEERecHit, v_HGCEERecHit);

    edm::Handle<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> v_HGCHEFRecHit;
    iEvent.getByToken(tok_HGCHEFRecHit, v_HGCHEFRecHit);

    edm::Handle<edm::SortedCollection<HGCRecHit, edm::StrictWeakOrdering<HGCRecHit>>> v_HGCHEBRecHit;
    iEvent.getByToken(tok_HGCHEBRecHit, v_HGCHEBRecHit);

    std::map<DetId, const HGCRecHit *> m_recHit;
    std::map<DetId, const HGCRecHit *> m_HGCEERecHit;

    int nHGCEERecHit = v_HGCEERecHit->size();

    for (int iRecHit = 0; iRecHit < nHGCEERecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCEERecHit)[iRecHit];

        m_recHit[recHit->id()] = recHit;
        m_HGCEERecHit[recHit->id()] = recHit;
    }

    //
    int nHGCHEFRecHit = v_HGCHEFRecHit->size();

    for (int iRecHit = 0; iRecHit < nHGCHEFRecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCHEFRecHit)[iRecHit];

        m_recHit[recHit->id()] = recHit;
    }

    //
    int nHGCHEBRecHit = v_HGCHEBRecHit->size();

    for (int iRecHit = 0; iRecHit < nHGCHEBRecHit; iRecHit++)
    {
        const HGCRecHit *recHit = &(*v_HGCHEBRecHit)[iRecHit];

        m_recHit[recHit->id()] = recHit;
    }

    // Sim clusters
    edm::Handle<std::vector<SimCluster>> h_simCluster;
    iEvent.getByToken(tok_simCluster, h_simCluster);
    const std::vector<SimCluster> &v_simCluster = *h_simCluster;

    // <DetId, SimCluster index>
    std::map<DetId, int> m_simClusterHit;

    int iSimClus = 0;

    for (auto const &simClus : v_simCluster)
    {
        auto v_hnf = simClus.hits_and_fractions();

        for (auto &hnf : v_hnf)
        {
            DetId detId(hnf.first);

            m_simClusterHit[detId] = iSimClus;
        }

        iSimClus++;
    }

    // Calo particles
    edm::Handle<std::vector<CaloParticle>> v_caloParticle;
    iEvent.getByToken(tok_caloParticle, v_caloParticle);

    // <DetId, CaloParticle index>
    std::map<DetId, int> m_caloPart_simClustHit;

    int nCaloPart = v_caloParticle->size();

    for (int iCaloPart = 0; iCaloPart < nCaloPart; iCaloPart++)
    {
        CaloParticle caloPart = v_caloParticle->at(iCaloPart);

        treeOutput->v_caloParticle_E.push_back(caloPart.energy());
        treeOutput->v_caloParticle_px.push_back(caloPart.px());
        treeOutput->v_caloParticle_py.push_back(caloPart.py());
        treeOutput->v_caloParticle_pz.push_back(caloPart.pz());

        treeOutput->v_caloParticle_pT.push_back(caloPart.pt());
        treeOutput->v_caloParticle_eta.push_back(caloPart.eta());
        treeOutput->v_caloParticle_phi.push_back(caloPart.phi());

        treeOutput->v_caloParticle_pdgid.push_back(caloPart.pdgId());

        treeOutput->caloParticle_n++;

        //auto v_simCluster = caloPart.simClusters().refVector();

        //int nSimCluster = v_simCluster.size();

        int iSimCluster = 0;

        //for(int iSimCluster = 0; iSimCluster < nSimCluster; iSimCluster++)
        for (auto const &simCluster : caloPart.simClusters())
        {
            //SimCluster simCluster = v_simCluster.at(iSimCluster);

            std::vector<std::pair<uint32_t, float>> v_hit = simCluster.get()->hits_and_fractions();

            //printf(
            //    "[%llu] "
            //    "CaloParticle %d/%d: "
            //    "SimCluster %d/%d: "
            //    "nSimHit %d "
            //    "nHGCEERecHit %d "
            //    "h&f size %d "
            //    "\n",
            //    eventNumber,
            //    iCaloPart+1, nCaloPart,
            //    iSimCluster+1, (int) caloPart.simClusters().size(),
            //    simCluster.get()->numberOfSimHits(),
            //    simCluster.get()->numberOfRecHits(),
            //    (int) v_hit.size()
            //);

            int nHit = v_hit.size();

            for (int iHit = 0; iHit < nHit; iHit++)
            {
                DetId detId(v_hit.at(iHit).first);

                m_caloPart_simClustHit[detId] = iCaloPart;
            }

            iSimCluster++;
        }
    }

    // SimHits
    if (storeSimHit)
    {
        for (int iSimHit = 0; iSimHit < nSimHit; iSimHit++)
        {
            auto simHit = v_HGCEESimHit->at(iSimHit);

            DetId detId(simHit.id());

            int layer = recHitTools.getLayer(simHit.id()) - 1; // Start from 0
            int zside = recHitTools.zside(simHit.id());

            auto position = recHitTools.getPosition(simHit.id());

            //CLHEP::Hep3Vector recHit_3vec(position.x(), position.y(), position.z());
            //
            //printf(
            //    "[%llu] "
            //    "SimHit %d/%d: "
            //    "layer %d, "
            //    "E %0.2e, "
            //    "ET %0.2e (%0.2e), "
            //    "x %0.2f, y %0.2f, z %0.2f,"
            //    "\n",
            //    eventNumber,
            //    iSimHit+1, nSimHit,
            //    layer,
            //    simHit.energy(),
            //    recHitTools.getPt(position, simHit.energy()), simHit.energy()*sin(recHit_3vec.theta()),
            //    position.x(), position.y(), position.z()
            //);

            //if(detId.layer() < minLayer)
            //{
            //    minLayer = detId.layer();
            //}
            //
            //if(detId.layer() > maxLayer)
            //{
            //    maxLayer = detId.layer();
            //}

            bool isCaloParticleMatched = (m_caloPart_simClustHit.find(detId) != m_caloPart_simClustHit.end());

            treeOutput->v_simHit_E.push_back(simHit.energy());

            treeOutput->v_simHit_x.push_back(position.x());
            treeOutput->v_simHit_y.push_back(position.y());
            treeOutput->v_simHit_z.push_back(position.z());

            treeOutput->v_simHit_eta.push_back(position.eta());
            treeOutput->v_simHit_phi.push_back(position.phi());

            treeOutput->v_simHit_ET.push_back(recHitTools.getPt(position, simHit.energy()));

            treeOutput->v_simHit_layer.push_back(layer + 1);
            treeOutput->v_simHit_zside.push_back(zside);
            treeOutput->v_simHit_isCaloParticleMatched.push_back(isCaloParticleMatched);

            int matchedSimClusIdx = (m_simClusterHit.find(detId) != m_simClusterHit.end()) ? m_simClusterHit[detId] : -1;
            treeOutput->v_simHit_matchedSimClusIndex.push_back(matchedSimClusIdx);

            treeOutput->simHit_n++;
        }
    }

    // RecHits
    if (storeRecHit)
    {
        std::vector<int> v_recHit_matchedSimHitIndex = Common::associateRecToSimHit(v_HGCEESimHit, v_HGCEERecHit);

        nHGCEERecHit = v_HGCEERecHit->size();

        for (int iRecHit = 0; iRecHit < nHGCEERecHit; iRecHit++)
        {
            auto recHit = (*v_HGCEERecHit)[iRecHit];

            int layer = recHitTools.getLayer(recHit.id()) - 1; // Start from 0
            int zside = recHitTools.zside(recHit.id());

            auto position = recHitTools.getPosition(recHit.id());

            //CLHEP::Hep3Vector recHit_3vec(position.x(), position.y(), position.z());
            //
            //printf(
            //    "[%llu] "
            //    "RecHit %d/%d: "
            //    "layer %d "
            //    "E %0.2f "
            //    "\n",
            //    eventNumber,
            //    iRecHit+1, nHGCEERecHit,
            //    layer,
            //    recHit.energy()
            //);

            //if(detId.layer() < minLayer)
            //{
            //    minLayer = detId.layer();
            //}
            //
            //if(detId.layer() > maxLayer)
            //{
            //    maxLayer = detId.layer();
            //}

            bool isCaloParticleMatched = (m_caloPart_simClustHit.find(recHit.id()) != m_caloPart_simClustHit.end());

            treeOutput->v_recHit_E.push_back(recHit.energy());

            treeOutput->v_recHit_x.push_back(position.x());
            treeOutput->v_recHit_y.push_back(position.y());
            treeOutput->v_recHit_z.push_back(position.z());

            treeOutput->v_recHit_eta.push_back(position.eta());
            treeOutput->v_recHit_phi.push_back(position.phi());

            treeOutput->v_recHit_ET.push_back(recHitTools.getPt(position, recHit.energy()));

            treeOutput->v_recHit_layer.push_back(layer + 1);
            treeOutput->v_recHit_zside.push_back(zside);

            treeOutput->v_recHit_isCaloParticleMatched.push_back(isCaloParticleMatched);

            treeOutput->v_recHit_matchedSimHitIndex.push_back(v_recHit_matchedSimHitIndex.at(iRecHit));

            int matchedSimClusIdx = (m_simClusterHit.find(recHit.id()) != m_simClusterHit.end()) ? m_simClusterHit[recHit.id()] : -1;
            treeOutput->v_recHit_matchedSimClusIndex.push_back(matchedSimClusIdx);

            //HGCalTopology::DecodedDetId decodetDetId = topo_HGCalEE.decode(recHit.id());
            //
            //treeOutput->v_recHit_iType.push_back(decodetDetId.iType);
            //treeOutput->v_recHit_iCell1.push_back(decodetDetId.iCell1);
            //treeOutput->v_recHit_iCell2.push_back(decodetDetId.iCell2);

            treeOutput->v_recHit_SiThickness.push_back(recHitTools.getSiThickness(recHit.id()));

            treeOutput->recHit_n++;
        }
    }

    // PFRecHits
    //edm::Handle <std::vector <reco::PFRecHit> > v_PFRecHit;
    //iEvent.getByToken(tok_PFRecHit, v_PFRecHit);
    //
    //int nPFRecHit = v_PFRecHit->size();
    ////printf("[%llu] # of PFRecHits: %d \n", eventNumber, nPFRecHit);
    //
    //for(int iPFRecHit = 0; iPFRecHit < nPFRecHit; iPFRecHit++)
    //{
    //    auto recHit = v_PFRecHit->at(iPFRecHit);
    //
    //    //HGCEEDetId detId(recHit.id());
    //    HGCalDetId detId(recHit.detId());
    //
    //    // int layer = detId.layer();
    //    int layer = recHitTools.getLayer(recHit.detId()) - 1; // Start from 0
    //
    //    //int zside = detId.zside();
    //    int zside = recHitTools.zside(recHit.detId());
    //
    //    printf(
    //        "[%llu] "
    //        "RecHit %d/%d: "
    //        "layer %d (%d), "
    //        "zside %+d, "
    //        "E %0.2f, "
    //        "\n",
    //        eventNumber,
    //        iPFRecHit+1, nPFRecHit,
    //        detId.layer(), recHit.layer(),
    //        zside,
    //        recHit.energy()
    //    );
    //
    //    //if(detId.layer() < minLayer)
    //    //{
    //    //    minLayer = detId.layer();
    //    //}
    //    //
    //    //if(detId.layer() > maxLayer)
    //    //{
    //    //    maxLayer = detId.layer();
    //    //}
    //
    //    //if(zside > 0)
    //    //{
    //    //    treeOutput->v_recHit_HGCalEEPlayer_totE.at(layer) += recHit.energy();
    //    //}
    //    //
    //    //else
    //    //{
    //    //    treeOutput->v_recHit_HGCalEEMlayer_totE.at(layer) += recHit.energy();
    //    //}
    //}

    // Fill tree
    treeOutput->fill();

    //#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    //ESHandle<SetupData> pSetup;
    //iSetup.get<SetupRecord>().get(pSetup);
    //#endif

    printf("\n\n");

    fflush(stdout);
    fflush(stderr);
}

// ------------ method called once each job just before starting event loop  ------------
void TreeMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void TreeMaker::endJob()
{

    printf("minLayer = %d, maxLayer = %d \n", minLayer, maxLayer);

    fflush(stdout);
    fflush(stderr);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TreeMaker::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
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
DEFINE_FWK_MODULE(TreeMaker);
