#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "OneTrack.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "Definition.h"
#include "MUV1Geometry.hh"
#include "MUV2Geometry.hh"
#include "TRecoVCandidate.hh"
#include <algorithm>

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

/// \class OneTrack
/// \Brief
/// Short description of your Analyzer
/// \EndBrief
///
/// \Detailed
/// Detailed description of your Analyzer\n\n
/// For examples of working Analyzer you can have a look at the examples in the Examples/ directory:\n
/// LKrPhotonMC \n
/// Pi0Reconstruction \n
/// SkimmingNoStrawMuonsTracks \n
/// Or the framework analyzers that can be found in the Analyzers/ directories: \n
/// CedarMCTester \n
/// VertexCDA \n
/// \n
/// All the following classes are available for you to use. Check their documentation for more information:\n
/// NA62Analysis::manip::TermManip \n
/// NA62Analysis::Analyzer \n
/// NA62Analysis::CounterHandler \n
/// NA62Analysis::DetectorAcceptance \n
/// NA62Analysis::EventFraction \n
/// NA62Analysis::MCSimple \n
/// NA62Analysis::NeuralNetwork \n
/// NA62Analysis::ParticleInterface \n
/// NA62Analysis::ParticleTree \n
/// NA62Analysis::StringBalancedTable \n
/// NA62Analysis::StringTable \n
/// NA62Analysis::UserMethods \n
/// NA62Analysis::Verbose \n
///
/// You might also be interested in checking the documentation for the following classes. However you should not
/// in general have to create instances of these. If necessary a pointer to the existing instance is usually
/// available or provided by specific methods.\n
/// NA62Analysis::Core::IOHisto \n
/// NA62Analysis::Core::IOTree \n
/// NA62Analysis::Core::IOHandler \n
/// NA62Analysis::Core::HistoHandler \n
/// NA62Analysis::Core::HistoHandler::Iterator \n
///
/// \EndDetailed

OneTrack::OneTrack(Core::BaseAnalysis *ba) : Analyzer(ba, "OneTrack")
{
    /// \MemberDescr
    /// \param ba : parent BaseAnalysis
    ///
    /// Specify the trees you want to use and the event class corresponding\n
    /// Don't try to load MCTruth tree (RUN_0 or Event). Use the MCTruthEvent in Process function instead. Problems when opening twice the same tree.\n
    /// Example with RecoEvent\n
    ///     \code
    ///             RequestTree("GigaTracker", new TRecoGigaTrackerEvent);
    /// \endcode
    /// Example with MC Event\n
    ///     \code
    ///             RequestTree("GigaTracker", new TGigaTrackerEvent);
    /// \endcode
    /// Example with generic tree\n
    ///     \code
    ///             RequestTree<MyClass>("MyTree", "BranchName", "MyClass", new MyClass);
    /// \endcode
    //// \n
    /// Initialize DetectorAcceptance if needed\n
    /// use of global instance\n
    ///     \code
    ///             fDetectorAcceptanceInstance = GetDetectorAcceptanceInstance();
    /// \endcode
    /// use of local instance\n
    ///     \code
    ///             fDetectorAcceptanceInstance = new DetectorAcceptance("./NA62.root");
    /// \endcode

    RequestTree("LKr",new TRecoLKrEvent);
    RequestTree("Spectrometer",new TRecoSpectrometerEvent);
    RequestTree("MUV1",new TRecoMUV1Event);
    RequestTree("MUV2",new TRecoMUV2Event);
    RequestTree("MUV3",new TRecoMUV3Event);
    RequestTree("RICH",new TRecoRICHEvent);
    RequestTree("CHOD",new TRecoCHODEvent);
    RequestTree("Cedar",new TRecoCedarEvent);




}

void OneTrack::InitOutput(){
    /// \MemberDescr
    /// Register the output variables of the analyzer.\n
    /// Call: \n
    ///     \code
    ///     RegisterOutput("outputName", &variableName)
    /// \endcode
    /// for each variable that should be in the output of the Analyzer\n
    /// The name of the analyzer will be prepended to the outputName (to avoid collisions with other analyzers)\n
    /// variableName should be the name of a variable declared in the definition of the class\n
    /// \n
    /// Call one of: \n
    ///     \code
    ///     AddParam("paramName", &variableName, defaultValue);
    /// \endcode
    /// for each parameter of the analyzer. These parameters can be set when starting the FW from the command line with the -p option.\n
    /// paramName is the name of the parameter in the command line\n
    /// variableName is the name of the variable that should be declared in the definition of the class\n
    /// defaultValue is the default value if not specified in the command line\n
    /// The allowed types for parameters are the following: bool, int, long, float, double, char, string, TString\n
    /// \n
    /// To create a new TTree in the output file, call: \n
    ///     \code
    ///     void OpenNewTree("TTreeName", "TTreeTitle");
    /// \endcode
    /// TTreeName is the name of the TTree (will be used to refer to this TTree later)\n
    /// TTreeTitle is the title of the TTree\n
    /// \n
    /// To add a branch to the newly created TTree, call: \n
    ///     \code
    ///     void AddBranch<VariableType>("TTreeName", "BranchName", &pointer);
    /// \endcode
    /// VariableType is the type of the variable for this branch (fundamental data type or class)\n
    /// TTreeName is the name of the TTree to add this branch\n
    /// BranchName is the name of the branch\n
    /// pointer is a pointer to the variable (should be declared in the header file)\n
    /// \n
    /// To create a standard TTree containing KineParts (for candidates) use\n
    ///     \code
    ///     CreateStandardTree("TTreeName", "TTreeTitle");
    /// \endcode
    /// \EndMemberDescr
}

void OneTrack::InitHist(){

    //Run/burst related histograms
    BookHisto(new TH1I("BurstID","Burst  ID;BurstID",4000,0,4000));

    BookHisto(new TH1I("STRAW_NCandidates", "Number of STRAW candidates;NCandidates" ,50,0,50));
    BookHisto(new TH1I("LKr_NCandidates " , "Number of LKr candidates;NCandidates" ,50,0,50));
    BookHisto(new TH1I("MUV1_NCandidates" , "Number of MUV1 candidates;NCandidates" ,50,0,50));
    BookHisto(new TH1I("MUV2_NCandidates" , "Number of MUV2 candidates;NCandidates" ,50,0,50));
    BookHisto(new TH1I("MUV3_NCandidates" , "Number of MUV3 candidates;NCandidates" ,50,0,50));
    BookHisto(new TH1I("RICH_NCandidates" , "Number of RICH candidates;NCandidates" ,50,0,50));
    BookHisto(new TH1I("CHOD_NCandidates" , "Number of CHOD candidates;NCandidates" ,50,0,50));

    //STRAW
    BookHisto(new TH1I("STRAW_Nchambers", "Number of chambers per candidate in STRAW; Nchambers", 20, 0, 20));
    BookHisto(new TH1F("STRAW_Chi2", "Track Chi2 from the STRAW ", 200, 0., 200.));

    //CHOD
    BookHisto(new TH2F("CHOD_cda_x_vs_y", "Distance beteen extrapolated track and position in the CHOD  x vs y;x[mm];y[mm]", 300, -300, 300., 300, -300., 300.));
    BookHisto(new TH1F("CHOD_nearest_track_dtrkcl", "Distance beteen extrapolated track and position in the CHOD for the closest track;CHOD_trkd [mm] ", 150, 0., 300.));
    BookHisto(new TH2F("CHOD_nearest_track_x_vs_y", "CHOD candidate x vs y ; x[mm];y[mm]", 520, -1300., 1300., 520, -1300., 1300.));
    BookHisto(new TH2F("CHOD_extrap_x_vs_y", "CHOD extrapolated track x vs y ; x[mm];y[mm]", 260, -1300., 1300., 260, -1300., 1300.));
}

void OneTrack::DefineMCSimple(){

    /// \MemberDescr
    /// Setup of fMCSimple. You must specify the generated MC particles you want.\n
    /// Add particles you want to recover from fMCSimple\n
    ///     \code
    ///     int particleID = fMCSimple.AddParticle(parentID, pdgCode)
    ///     \endcode
    /// parentID :  0=no parent (=beam particle)\n
    ///     ...\n
    /// Example : you want to retrieve the kaon from the beam, the pi0 an pi+ from the beam kaon and the 2 photons coming from the previous pi0 decay :\n
    ///     \code
    ///     int kaonID = fMCSimple.AddParticle(0, 321) //Ask beam kaon (sequence ID=1)
    ///     fMCSimple.AddParticle(kaonID, 211) //Ask pi+ from previous kaon (sequence ID=2)
    ///     int pi0ID = fMCSimple.AddParticle(kaonID, 111) //Ask pi0 from previous kaon (sequence ID=3)
    ///     fMCSimple.AddParticle(pi0ID, 22) //Ask first gamma from previous pi0 (sequence ID=4)
    ///     fMCSimple.AddParticle(pi0ID, 22) //Ask second gamma from previous pi0 (sequence ID=4)
    ///     \endcode
    ///
    /// @see ROOT TDatabasePDG for a list of PDG codes and particle naming convention
    /// \EndMemberDescr
}

void OneTrack::StartOfRunUser(){
    /// \MemberDescr
    /// This method is called at the beginning of the processing (corresponding to a start of run in the normal NA62 data taking)\n
    /// Do here your start of run processing if any
    /// \EndMemberDescr
}

void OneTrack::StartOfBurstUser(){
    /// \MemberDescr
    /// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the beginning of the first file\n
    /// Do here your start/end of burst processing if any
    /// \EndMemberDescr
}

void OneTrack::Process(int iEvent){

    //if(fMCSimple.fStatus == MCSimple::kMissing){printIncompleteMCWarning(iEvent);return;}
    //if(fMCSimple.fStatus == MCSimple::kEmpty){printNoMCWarning();return;}

    TRecoLKrEvent *LKrEvent = (TRecoLKrEvent*)GetEvent("LKr");
    TRecoSpectrometerEvent *SpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");
    TRecoMUV1Event *MUV1Event = (TRecoMUV1Event*)GetEvent("MUV1");
    TRecoMUV2Event *MUV2Event = (TRecoMUV2Event*)GetEvent("MUV2");
    TRecoMUV3Event *MUV3Event = (TRecoMUV3Event*)GetEvent("MUV3");
    TRecoRICHEvent *RICHEvent = (TRecoRICHEvent*)GetEvent("RICH");
    TRecoCHODEvent *CHODEvent = (TRecoCHODEvent*)GetEvent("CHOD");
    TRecoCedarEvent *CedarEvent = (TRecoCedarEvent*)GetEvent("Cedar");
    //Time Offset for all the detectors differences (ATM using only CHOD as reference)
    double LKrOffset   = +115; //Old 112.3 Use GetClusterTime
    double CedarOffset = 0; //old 0
    double MUV1Offset  = -5;//new -5 old -8.5
    double MUV2Offset  = -5;//new -5 old -20.4
    double MUV3Offset  = -11.35;//old -11.35

    //Cuts used on the timedifferences for the Kmu2 selection
    double LKrOffsetCut  = 100;
    double RICHOffsetCut = 100;
    double MUV3OffsetCut = 100;
    double MUV1OffsetCut = 100;
    double MUV2OffsetCut = 100;
    double CedarOffsetCut= 100;

    double mev2gev = 0.001;
    double gev2mev = 1000.;
    double mm2cm = 0.1;
    double cm2mm = 10.;

    //CUTComment:: Only one candidate in the STRAW
    FillHisto("STRAW_NCandidates", SpectrometerEvent->GetNCandidates());
    if(SpectrometerEvent->GetNCandidates() != 1 )return;

    TRecoSpectrometerCandidate* Track;
    Track = ((TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0));


    //Comment:: All necessary STRAW variables are defined here
    double STRAW_P    = Track->GetMomentum();
    double STRAW_Pbf  = Track->GetMomentumBeforeFit();
    double STRAW_bdxdz= Track->GetSlopeXBeforeMagnet();
    double STRAW_bdydz= Track->GetSlopeYBeforeMagnet();
    double STRAW_dxdz = Track->GetSlopeXAfterMagnet();
    double STRAW_dydz = Track->GetSlopeYAfterMagnet();
    int    STRAW_NC   = Track->GetNChambers();
    double STRAW_chi2 = Track->GetChi2();

    //CUTComment:: Check if the the chi2 of the track fitter in the straw is >20
    // and if there are at least 3 chambers fired
    FillHisto("STRAW_Chi2", STRAW_chi2);
    FillHisto("STRAW_Nchambers", STRAW_NC);


    //CUTComment:: Cut Stage 1
    //1.Check if track is positive (K+ beam in NA62)
    //2.Chi2 coming from STRAW has to be smaller than 20
    //3. Number of STRAW chambers fired has to be >= 3
    //4. At least one Cedar candidate
    //5. Cedar sectors > 5
    if(Track->GetCharge() != 1){ return;}

    if( STRAW_chi2 > 20 || STRAW_NC   < 3 ){ return; }


    if( CedarEvent->GetNCandidates() == 0 ){return;}
    for(int iCedarCand=0; iCedarCand < CedarEvent->GetNCandidates(); iCedarCand++){

        TRecoCedarCandidate*  CedarCandidate = ((TRecoCedarCandidate*)CedarEvent->GetCandidate(iCedarCand));
        if(CedarCandidate->GetNSectors() < 5 ){return;}

    }



    //Getting the position vectors and the slopes of the tracks
    //before the magnet @ DCH1
    TVector3 PositionBefore;
    TVector3 SlopesBefore;
    SlopesBefore.SetX(STRAW_bdxdz);
    SlopesBefore.SetY(STRAW_bdydz);
    SlopesBefore.SetZ(1.);
    PositionBefore = Track->GetPositionBeforeMagnet();

    //Building track momentum from STRAW slopes before magnet
    //and the position @ DCH1
    double norm = 1./sqrt(pow(SlopesBefore.X(),2) + pow(SlopesBefore.Y(),2) + 1.  );
    TVector3 TrackP;
    TrackP.SetXYZ(norm*SlopesBefore.X()*STRAW_P, norm*SlopesBefore.Y()*STRAW_P, norm*STRAW_P );


    //Beam momentum and position on Trim5 that will be used
    //for the making of the vertex
    double beam_norm = 1./sqrt(XAngle*XAngle + 1. );
    TVector3 BeamP;
    TVector3 BeamTrim5Pos;
    BeamP.SetXYZ(beam_norm*KEnergy*XAngle*1000,0,beam_norm*KEnergy*1000.);
    BeamTrim5Pos.SetXYZ(0., 0., Ztrim*1000.);

    //Calculating the intersection point between the kaon and
    //the track and getting the cda using the VertexCDA routine by Giuseppe
    double cda = 0;
    TVector3 Vertex;
    Vertex = VertexCDA(PositionBefore, TrackP, BeamTrim5Pos, BeamP, cda );

    //CUTComment:: Cut Stage 1
    //6.Closest Approached Distance > 40.
    //7.Zvtx to be incide of the fiducial volume of NA62 detector (105 - 180 m)
    if(cda > 40.){return;}
    if( Vertex.Z() < 105000 || Vertex.Z() > 180000){return;}

    //Computing the NuNubar 3vector (Missing mass)
    TVector3 NuNubar;
    NuNubar = BeamP - TrackP;


    //Making the 4Momentum of NuNuBar
    //4Momentum of the track (assuming muon mass for Kmu2 selection)
    TLorentzVector PiP;
    PiP.SetVect(TrackP);
    PiP.SetE( sqrt(pow(TrackP.Mag(),2) + pow(PiPlMass*1000,2)));

    //4Momentum of the Kaon (NO GTK yet)
    TLorentzVector KP;
    KP.SetVect(BeamP);
    KP.SetE(sqrt(pow(BeamP.Mag(),2) + pow(KMass*1000,2)));

    //4Momentum of the NuNubar system
    TLorentzVector NuNubarP;
    NuNubarP = KP - PiP;

    //Angle between the track and the Kaon
    double theta;
    theta = TrackP.Dot(BeamP)/(TrackP.Mag()*BeamP.Mag());


    //Getting the position vectors and the slopes of the tracks
    //after the magnet @ DCH4
    TVector3 SlopesAfter;
    TVector3 PositionAfter;
    SlopesAfter.SetX(STRAW_dxdz);
    SlopesAfter.SetY(STRAW_dydz);
    SlopesAfter.SetZ(1.);
    PositionAfter = Track->GetPositionAfterMagnet();

    //Extrapolating the track to the other detectors
    TVector3 RICH_extrap = PositionAfter + ( ZRICHStart*1000 - PositionAfter.Z() )*SlopesAfter;
    TVector3 CHOD_extrap = PositionAfter + ( ZCHODStart*1000 - PositionAfter.Z() )*SlopesAfter;
    TVector3 MUV1_extrap = PositionAfter + ( ZMUV1Start*1000 - PositionAfter.Z() )*SlopesAfter;
    TVector3 MUV2_extrap = PositionAfter + ( ZMUV2Start*1000 - PositionAfter.Z() )*SlopesAfter;
    TVector3 MUV3_extrap = PositionAfter + ( ZMUV3Start*1000 - PositionAfter.Z() )*SlopesAfter;
    TVector3 LKr_extrap  = PositionAfter + ( ZLKrStart *1000 - PositionAfter.Z() )*SlopesAfter;

    double CHODdtrkcl_min;
    double LKrdtrkcl_min;
    double MUV1dtrkcl_min;
    double MUV2dtrkcl_min;
    double MUV3dtrkcl_min;
    vector <double> CEDAR_STRAW_tdiff;
    int CHODClosestTrackIndex = FindClosestCluster(CHODEvent, CHOD_extrap, "CHOD" , CHODdtrkcl_min);
    int LKrTrackClusterIndex  = FindClosestCluster(LKrEvent , LKr_extrap , "LKr"  , LKrdtrkcl_min);
    int MUV1TrackClusterIndex = FindClosestCluster(MUV1Event, MUV1_extrap, "MUV1" , MUV1dtrkcl_min);
    int MUV2TrackClusterIndex = FindClosestCluster(MUV2Event, MUV2_extrap, "MUV2" , MUV2dtrkcl_min);
    int MUV3TrackClusterIndex = FindClosestCluster(MUV3Event, MUV3_extrap, "MUV3" , MUV3dtrkcl_min);



    //Energy scale correction and non-linearity correction for the LKr taken from Giuseppe
    Double_t fEScale = 1.03;
    for(int iLKrCand=0; iLKrCand<LKrEvent->GetNCandidates(); iLKrCand++){
        TRecoLKrCandidate* LKrCluster = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(iLKrCand));

        // ZS non linearity
        Double_t ue = LKrCluster->GetClusterEnergy();
        Double_t ce = ue;
        if (LKrCluster->GetNCells()>9) {
            if (ue<22) ce = ue/(0.7666+0.0573489*log(ue));
            if (ue>=22 && ue<65) ce = ue/(0.828962+0.0369797*log(ue));
            if (ue>=65) ce = ue/(0.828962+0.0369797*log(65));
        }

        LKrCluster->SetClusterEnergy(ce*fEScale);

    }

    //map that will contain the gamma`s`
    map<string, int> gammas;

    for(int iLKrCand=0; iLKrCand<LKrEvent->GetNCandidates(); iLKrCand++){

        TRecoLKrCandidate* LKrCluster = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(iLKrCand));
        TRecoLKrCandidate* AssociatedCluster = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex));

        TVector3 LkrPos;
        LkrPos.SetX ( LKrCluster->GetClusterX() );
        LkrPos.SetY ( LKrCluster->GetClusterY() );

        TVector3 LkrNtPos;
        LkrNtPos.SetX ( AssociatedCluster->GetClusterX() );
        LkrNtPos.SetY ( AssociatedCluster->GetClusterY() );

        double LKrEcluster   = 1000*LKrCluster->GetClusterEnergy(); //  [MeV]
        double LKrEseed      = 1000*LKrCluster->GetClusterSeedEnergy(); //  [MeV]
        double LKrE77        = 1000*LKrCluster->GetCluster77Energy(); //  [MeV]
        int    LKrNcells     = LKrCluster->GetNCells();
        double ClusterTime   = LKrCluster->GetClusterTime();
        double ClusterNtTime = AssociatedCluster->GetClusterTime();

        double LKrClusterEnergy = 1000*LKrCluster->GetClusterEnergy(); //  [MeV]

        if(iLKrCand == LKrTrackClusterIndex){

            //CUTComment:: Choose MIP Pion
            //if(LKrClusterEnergy > 800 ){return;}
            //if(LKrClusterEnergy < 300 ){return;}
            //if(LKrCluster->GetNCells()>5 ) { return;} //  [MeV]
            //CUTComment:: Only high energy deposition in the LKr
            if(LKrClusterEnergy < 1500 ) { return;} //  [MeV]

        }

        if(iLKrCand != LKrTrackClusterIndex){
            double Cluster_distance =sqrt(pow(LkrPos.X()*10. - LkrNtPos.X()*10, 2 ) + pow(LkrPos.Y()*10. - LkrNtPos.Y()*10., 2 ) );
            double Cluster_tdiff = ClusterNtTime - ClusterTime;

            if(fabs(Cluster_tdiff) < 10 && Cluster_distance > 200){
                //cout << "Icl == " << iLKrCand << "Itrk == " << LKrTrackClusterIndex << "Cluster_tdiff == " << Cluster_tdiff << " Dcl == " << Cluster_distance << endl;
                //cout << gammas.size() << endl;
                if( gammas.size() == 0 ){

                    gammas["g1"] = iLKrCand;
                    //cout << "G1 Before == " << gammas["g1"] << "iLKrcand" << iLKrCand << endl;
                } else if( gammas.size()==1 ){
                    gammas["g2"] = iLKrCand;
                    //cout << "G1 == " << gammas["g1"] << "G2 == " << gammas["g2"] << "iLkrcand == " << iLKrCand << endl;
                } else if( gammas.size()==2 ){
                    gammas["g3"] = iLKrCand;
                    cout << "G1 == " << gammas["g1"] << "G2 == " << gammas["g2"]<< "G3 == " << gammas["g3"] <<  endl;
                }
            }

        }

    }
    if(gammas.size() < 2 ){ return;}

    if(CHODClosestTrackIndex < 0. ){return;}

    TVector2 CD_CHODPos   = ((TRecoCHODCandidate*)CHODEvent->GetCandidate(CHODClosestTrackIndex))->GetHitPosition();
    FillHisto("CHOD_cda_x_vs_y", CD_CHODPos.X()*10. - CHOD_extrap.X() , CD_CHODPos.Y()*10. - CHOD_extrap.Y());
    FillHisto("CHOD_nearest_track_dtrkcl", CHODdtrkcl_min);
    FillHisto("CHOD_nearest_track_x_vs_y", CD_CHODPos.X()*10.,CD_CHODPos.Y()*10. );
    FillHisto("CHOD_extrap_x_vs_y", CHOD_extrap.X(), CHOD_extrap.Y() );



}

void OneTrack::PostProcess(){
    /// \MemberDescr
    /// This function is called after an event has been processed by all analyzers. It could be used to free some memory allocated
    /// during the Process.
    /// \EndMemberDescr

}

void OneTrack::EndOfBurstUser(){
    /// \MemberDescr
    /// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the end of the last file\n
    /// Do here your start/end of burst processing if any
    /// \EndMemberDescr
}

void OneTrack::EndOfRunUser(){

/// \MemberDescr
    /// This method is called at the end of the processing (corresponding to a end of run in the normal NA62 data taking)\n
    /// Do here your end of run processing if any\n
    /// \n
    /// You can also use this space to save plots in your output file.\n
    /// If you want to save them all, just call\n
    /// \code
    ///     SaveAllPlots();
    /// \endcode
    /// Or you can just save the ones you want with\n
    /// \code
    ///     histogram->Write()\n
    ///             fHisto.Get...("histoname")->Write();
    /// \endcode
    /// \n
    /// To run over a set of histograms you can use Iterators (HistoHandler::IteratorTH1,
    /// HistoHandler::IteratorTH2, HistoHandler::IteratorTGraph). You can use it to run over
    /// all the histograms or only a subset of histogram by using one of the two forms of
    /// GetIterator...  (replace ... by TH1, TH2 or TGraph)\n
    /// \code
    ///     GetIterator...()
    /// \endcode
    /// will get an Iterator running over all histograms of this type while
    /// \code
    ///     GetIterator...("baseName")
    /// \endcode
    /// will get an Iterator running only over the histograms of this type whose name starts
    /// with baseName.\n
    /// For more details and examples on how to use the Iterator after getting it, please refer
    /// to the HistoHandler::Iterator documentation.\n
    /// Although this is described here, Iterators can be used anywhere after the
    /// histograms have been booked.
    /// \EndMemberDescr

    SaveAllPlots();
}

void OneTrack::DrawPlot(){
    /// \MemberDescr
    /// This method is called at the end of processing to draw plots when the -g option is used.\n
    /// If you want to draw all the plots, just call\n
    /// \code
    ///     DrawAllPlots();
    /// \endcode
    /// Or get the pointer to the histogram with\n
    /// \code
    ///     fHisto.GetTH1("histoName");// for TH1
    ///     fHisto.GetTH2("histoName");// for TH2
    ///     fHisto.GetGraph("graphName");// for TGraph and TGraphAsymmErrors
    ///     fHisto.GetHisto("histoName");// for TH1 or TH2 (returns a TH1 pointer)
    /// \endcode
    /// and manipulate it as usual (TCanvas, Draw, ...)\n
    /// \EndMemberDescr
}

TVector3 OneTrack::VertexCDA(TVector3 pos1, TVector3 p1, TVector3 pos2, TVector3 p2, Double_t &cda)
{
    // Vertex function using CDA method.

    TVector3 d = pos1 - pos2;
    double p12 = p1.Dot(p2);
    double det = p12*p12 - p1.Mag2() * p2.Mag2();
    if (!det) return TVector3(-9999,-9999,-9999);
    double t1 = (p2.Mag2()*d.Dot(p1) - p1.Dot(p2)*d.Dot(p2)) / det;
    double t2 = (p1.Dot(p2)*d.Dot(p1) - p1.Mag2()*d.Dot(p2)) / det;
    TVector3 q1 = pos1 + t1*p1;
    TVector3 q2 = pos2 + t2*p2;
    TVector3 vertex = 0.5*(q1 + q2);
    cda = (q1 - q2).Mag();
    return vertex;

}

int OneTrack::FindClosestCluster( TRecoVEvent* Event, TVector3 Extrap_track, string detector_type, double& minimum){

    std::vector<double> dtrkcl;
    //std::vector<int>    clindex;
    std::vector<double>::iterator VIterator;
    int position= -1;
    for (int iCand=0; iCand < Event->GetNCandidates(); iCand++){
        if(detector_type == "LKr"){
            TRecoLKrCandidate* Cluster = ((TRecoLKrCandidate*)Event->GetCandidate(iCand));
            double clusterx = Cluster->GetClusterX();
            double clustery = Cluster->GetClusterY();
            double distance = sqrt(pow(clusterx*10. - Extrap_track.X(),2 ) + pow(clustery*10. - Extrap_track.Y(),2 ) ) ;
            dtrkcl.push_back(distance);
        }
        if(detector_type == "MUV1"){
            TRecoMUV1Candidate* Cluster = ((TRecoMUV1Candidate*)Event->GetCandidate(iCand));
            //Old Reco
            //double clusterx = Cluster->GetPosition().X()+ 60;
            //double clustery = Cluster->GetPosition().Y()+ 60;
            //New Reco
            double clusterx = Cluster->GetPosition().X();
            double clustery = Cluster->GetPosition().Y();
            double distance = sqrt(pow(clusterx - Extrap_track.X(),2 ) + pow(clustery - Extrap_track.Y(),2 ) ) ;
            dtrkcl.push_back(distance);
        }
        if(detector_type == "MUV2"){
            TRecoMUV2Candidate* Cluster = ((TRecoMUV2Candidate*)Event->GetCandidate(iCand));
            double clusterx = Cluster->GetPosition().X();
            double clustery = Cluster->GetPosition().Y();
            double distance = sqrt(pow(clusterx - Extrap_track.X(),2 ) + pow(clustery - Extrap_track.Y(),2 ) ) ;
            dtrkcl.push_back(distance);
        }
        if(detector_type == "MUV3"){
            TRecoMUV3Candidate* Cluster = ((TRecoMUV3Candidate*)Event->GetCandidate(iCand));
            double clusterx = Cluster->GetX();
            double clustery = Cluster->GetY();
            double distance = sqrt(pow(clusterx - Extrap_track.X(),2 ) + pow(clustery - Extrap_track.Y(),2 ) ) ;
            dtrkcl.push_back(distance);
        }
        if(detector_type == "CHOD"){
            TRecoCHODCandidate* Cluster = ((TRecoCHODCandidate*)Event->GetCandidate(iCand));
            double clusterx = Cluster->GetHitPosition().X();
            double clustery = Cluster->GetHitPosition().Y();
            double distance = sqrt(pow(clusterx*10. - Extrap_track.X(),2 ) + pow(clustery*10. - Extrap_track.Y(),2 ) ) ;
            dtrkcl.push_back(distance);
        }

    }
    //Debugging!!!
    //cout << " Event -----" << detector_type.data() << endl;
    //for ( VIterator = dtrkcl.begin() ; VIterator != dtrkcl.end(); VIterator++ ){
    //    cout << *VIterator << endl;
    //}
    //
    if(dtrkcl.size() !=0){
        minimum= *min_element(dtrkcl.begin(), dtrkcl.end());
        position = distance(dtrkcl.begin(), min_element(dtrkcl.begin(), dtrkcl.end()));
        dtrkcl.clear();
        //    cout << " Minimum = " << minimum << endl;
        //    cout << "Index of minimum = " << position << endl;

        //cout << "WOWW == " <<  clindex(min_element(dtrkcl.begin(), dtrkcl.end())) << endl;
        //cout << " Event OVER-----" << detector_type.data() << endl;
        //cout << minimum << endl;

    }
    return position;

}
