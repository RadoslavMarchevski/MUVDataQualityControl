#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "Kmu2.hh"
#include "Definition.h"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "MUV1Geometry.hh"
#include "FADCEvent.hh"
#include "TDigiVEvent.hh"
#include "TMUV1Digi.hh"
#include "TMUV2Digi.hh"
#include "MUV2Geometry.hh"
#include "TRecoVCandidate.hh"
#include <algorithm>


using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;

/// \class Kmu2
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

Kmu2::Kmu2(Core::BaseAnalysis *ba) : Analyzer(ba, "Kmu2")
{
    /// \MemberDescr
    /// \param ba : parent BaseAnalysis
    ///
    /// Specify the trees you want to use and the event class corresponding\n
    /// Don't try to load MCTruth tree (RUN_0 or Event). Use the MCTruthEvent in Process function instead. Problems when opening twice the same tree.\n
    /// Example with RecoEvent\n
    ///     \code
    ///             RequestTree("GigaTracker", new TRecoGigaTrackerEvent);
    ///             RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Reco");
    ///             RequestTree("GigaTracker", new TRecoGigaTrackerEvent, "Digis");
    /// \endcode
    /// Example with MC Event\n
    ///     \code
    ///             RequestTree("GigaTracker", new TGigaTrackerEvent);
    /// \endcode
    /// Example with generic tree\n
    ///     \code
    ///             RequestTree<MyClass>("MyTree", "BranchName", "MyClass", new MyClass);
    ///             RequestTree("MyTree", "BranchName", "MyClass", new MyClass);
    /// \endcode
    /// Requesting Trigger data\n
    ///     \code
    ///             RequestL0Data();
    ///             RequestL1Data();
    ///             RequestL2Data();
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
    RequestTree("MUV2",new TRecoMUV2Event, "Reco");
    RequestTree("MUV2",new FADCEvent,"Digis");
    RequestTree("MUV3",new TRecoMUV3Event);
    RequestTree("RICH",new TRecoRICHEvent);
    RequestTree("CHOD",new TRecoCHODEvent);
    RequestTree("Cedar",new TRecoCedarEvent);
    //RequestL0Data();

}

void Kmu2::InitOutput(){
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

void Kmu2::InitHist(){
    /// \MemberDescr
    /// Book and Initialize histograms in this function.\n
    /// Same function to Book TH1, TH2, TGraph and TGraphAsymmErrors (anything derived from TH1 or TGraph)
    ///     \code
    ///     BookHisto(histogram*)
    /// \endcode
    /// If isAutotUpdate is true, this histogram will be drawn and updated regularly during the processing (default=false).\n
    /// The refresh interval can be set with (default=10):
    /// \code
    ///     SetUpdateInterval(interval)
    /// \endcode
    /// Defining plots as AutoUpdate and setting the interval can also be done at runtime with a configuration file.\n
    /// \n
    /// Example of booking an histogram: \n
    ///     \code
    ///     BookHisto(new TH2I("PartEnergy", "Energy as a function of particle", 0, 0, 0, Bins, MinEnergy, MaxEnergy));
    /// \endcode
    /// Booking of counters and creation of EventFraction can be done here with\n
    ///     \code
    ///             BookCounter(name)
    ///             NewEventFraction(name)
    /// \endcode
    /// Example\n
    ///     \code
    ///     BookCounter("Total");
    ///     BookCounter("PassCuts");
    ///     NewEventFraction("Cuts");
    /// \endcode
    /// Add the counters to the EventFraction\n
    ///     \code
    ///     AddCounterToEventFraction(EventFractionName, CounterName)
    /// \endcode
    /// Example\n
    ///     \code
    ///     AddCounterToEventFraction("Cuts", "Total");
    ///     AddCounterToEventFraction("Cuts", "PassCuts");
    /// \endcode
    /// Then define which counter represents the sample size\n
    ///     \code
    ///     DefineSampleSizeCounter(EventFractionName, CounterName)
    /// \endcode
    /// Example\n
    ///     \code
    ///     DefineSampleSizeCounter("Cuts", "Total");
    /// \endcode
    /// You can check the content of the input file (directory, histograms, generic key) in any directory with the following methods
    /// \code
    ///     GetListOfKeys("CedarMonitoring"); // CedarMonitoring directory
    ///     GetListOfHisto(""); // Top directory
    ///     GetListOfTH1("CedarMonitoring");
    ///     GetListOfTH2("CedarMonitoring");
    ///     GetListOfTGraph("CedarMonitoring");
    ///     GetListOfDirs("");
    /// \endcode
    /// The first one returns a vector of IOHandler::keyPair containing the name and class name of the object
    /// \code
    ///     vector<IOHandler::keyPair> keys = GetListOfKeys("CedarMonitoring");
    ///     cout << keys[0].name << " " << keys[0].className << endl;
    /// \endcode
    /// The others returns a vector of TString containing the name of the objects.\n
    /// You can retrieve histograms from the input ROOT file (Anything derived from TH1) with\n
    ///     \code
    ///     RequestHistogram("TDirectoryName", "HistogramName", appendOnNewFile);
    /// \endcode
    /// appendOnNewFile is a boolean. If set to true, each time a new file is opened the content
    /// of the histogram will be appended to the content of the previous one. If set to false, the content
    /// of the histogram is replaced each time a new file is opened.
    /// \EndMemberDescr

    //Run/burst related histograms
    //testing MUV candidates

    BookHisto(new TH1I("MUV123", "Track having MUV1,2,3 associated clusters; MUV1+2+3", 5, 0, 5));
    BookHisto(new TH1I("MUV13", "Track having MUV1,3 associated clusters; MUV1+3", 5, 0, 5));
    BookHisto(new TH1I("MUV23", "Track having MUV2,3 associated clusters; MUV2+3", 5, 0, 5));
    BookHisto(new TH1I("MUV3Only", "Track having MUV3 associated cluster; MUV3", 5, 0, 5));

    BookHisto(new TH1I("BurstID","Burst  ID;BurstID",4000,0,4000), "BurstInfo");
    BookHisto(new TH2I("BurstID_vs_MUV123","Burst  ID for MUV1&MUV2&MUV3 events;BurstID"  ,4000,0,4000, 5, 0 , 5));
    BookHisto(new TH2I("BurstID_vs_MUV23","Burst  ID for MUV3&MUV2 !& MUV1 events;BurstID",4000,0,4000, 5, 0 , 5));
    BookHisto(new TH2I("BurstID_vs_MUV13","Burst  ID for MUV3&MUV1 !& MUV2 events;BurstID",4000,0,4000, 5, 0 , 5));
    BookHisto(new TH2I("BurstID_vs_MUV3","Burst  ID for MUV3 only events;BurstID",4000,0,4000, 5,0,5));

    //MUV3 hit plots for the inefficient events
    BookHisto(new TH1I("MUV3Hit_NoMUV1", "MUV3 Hitmap for the MUV1 inefficient events", 200, 0, 200));
    BookHisto(new TH1I("MUV3Hit_NoMUV2", "MUV3 Hitmap for the MUV2 inefficient events", 200, 0, 200));
    BookHisto(new TH1I("MUV3Hit_NoMUV12", "MUV3 Hitmap for the MUV2 and MUV1 inefficient events", 200, 0, 200));
    BookHisto(new TH1I("MUV3Hit_GoodEvent", "MUV3 Hitmap for the fully efficient events", 200, 0, 200));

    //General Kinematic histograms
    BookHisto(new TH1F("TrackP", "STRAW Momentum ; Track_P[MeV]", 100, 0., 100000.));
    BookHisto(new TH1F("TrackPfit_TrackP", "GetMomentum() - GetMomentumBeforeFit() ; Track_P[MeV] - Track_Ppat[MeV]", 100, -50000., 50000.));
    BookHisto(new TH1F("BeamP", "Beam Momentum ; Beam_P[MeV]", 100, 0., 100000.));
    BookHisto(new TH1F("MM2", "Missing mass squared; (P_{K} - P_{#pi} )^2 [GeV^2]",200, -0.2,0.2));
    BookHisto(new TH2F("Track_P_vs_MM2", "Track Momentum vs Missing mass squared;P_{track} [GeV/c]; M_{miss}^2 [GeV^2/c^2]",100, 0., 100., 200, -0.2,0.2));
    BookHisto(new TH2F("Track_P_vs_Theta", " Missing mass squared vs Angle between kaon and #pi; P_{track} [GeV/c];#theta_{K#pi} [rad]",100, 0., 100., 200, 0.,0.02));

    //STRAW
    BookHisto(new TH1I("STRAW_Nchambers", "STRAW number of chambers per candidate; Nchambers", 20, 0, 20));
    BookHisto(new TH1F("TrackChi2", "STRAW Chi2", 200, 0., 200.));


    //MUV1
    BookHisto(new TH1I("MUV1_Ncandidates", "MUV1 number of candidates", 50, 0, 50));
    BookHisto(new TH1I("MUV1_Nhits", "MUV1 number of hits for the associated track cluster", 50, 0, 50));
    BookHisto(new TH2F("MUV1xvsy", "x vs y position in MUV1", 44, -1320., 1320.,44, -1320., 1320.));
    BookHisto(new TH1F("MUV1_RICH_timediff", "Time difference between MUV1 and RICH", 100, -50, 50.));
    BookHisto(new TH1F("MUV1_SeedEnergy", " MUV1 Energy of the two most energetic Horizontal+Vertical channels", 1000, 0, 10000.));
    BookHisto(new TH1F("MUV1_ClusterEnergy", " MUV1 cluster energy", 1000, 0, 10000.));
    BookHisto(new TH1F("MUV1_Eseed_over_Ecl", "MUV1 E_{seed}/E_{cluster}; E_{seed}/E_{cluster}", 120, 0.,1.2));
    BookHisto(new TH2F("MUV1_SeedEnergy_xvsy", " MUV1 SeedEnergy Horizontal + Vertical channels", 1000, 0, 10000., 1000, 0, 1000.));
    BookHisto(new TH2F("MUV1_SeedEnergy_vs_ClusterEnergy", " MUV1 SeedEnergy vs ClusterEnergy", 1000, 0, 10000., 1000, 0, 10000.));
    BookHisto(new TH2F("MUV1_cda_x_vs_y", "Difference between extraplated track and position given by MUV1;#Delta_x [mm];#Delta_y [mm]", 300, -300., 300., 300, -300., 300.));
    BookHisto(new TH1F("MUV1_trk_dist", "Distance beteen extrapolated track and position in the MUV1;MUV1_trkd [mm] ", 1500, 0., 3000.));


    //MUV2
    BookHisto(new TH1I("MUV2_Ncandidates", "MUV2 number of candidates", 50, 0, 50));
    BookHisto(new TH1I("MUV2_Nhits", "MUV2 number of hits for the associated track cluster", 50, 0, 50));
    BookHisto(new TH2F("MUV2xvsy", "x vs y position in MUV2", 22, -1320., 1320.,22, -1320., 1320.));


    BookHisto(new TH2F("MUV2_cda_x_vs_y", "Difference between extraplated track and position given by MUV2;#Delta_x [mm];#Delta_y [mm]", 300, -300., 300., 300, -300., 300.));
    BookHisto(new TH1F("MUV2_trk_dist", "Distance beteen extrapolated track and position in the MUV2;MUV2_trkd [mm] ", 1500, 0., 3000.));
    BookHisto(new TH1F("MUV2_ClusterEnergy", " MUV2 cluster energy", 1000, 0, 10000.));
    BookHisto(new TH1F("MUV2_SeedEnergy", " MUV2 Energy of the two most energetic Horizontal and Vertical channels", 1000, 0, 10000.));
    BookHisto(new TH1F("MUV2_Eseed_over_Ecl", "MUV2 E_{seed}/E_{cluster}; E_{seed}/E_{cluster}", 120, 0.,1.2));
    BookHisto(new TH2F("MUV2_SeedEnergy_xvsy", " MUV2 SeedEnergy Horizontal vs Vertical", 1000, 0, 10000., 1000, 0., 10000.));
    BookHisto(new TH2F("MUV2_SeedEnergy_vs_ClusterEnergy", " MUV2 SeedEnergy vs ClusterEnergy", 1000, 0, 10000., 1000, 0, 10000.));

    //MUV3
    BookHisto(new TH1I("MUV3_Ncandidates", "MUV3 number of candidates", 50, 0, 50));
    BookHisto(new TH2F("MUV3_cda_x_vs_y", "Difference between extraplated track and position given by MUV3;#Delta_x [mm];#Delta_y [mm]", 1000, -1000., 1000., 1000, -1000., 1000.));
    //BookHisto(new TH1F("MUV3_time", "MUV3 cluster Time", 100, -50, 50.));

    //CHOD
    BookHisto(new TH1I("CHOD_Ncandidates", "CHOD number of candidates", 20, 0, 20));
    BookHisto(new TH2F("CHOD_x_vs_y", "CHOD hitposition x vs y ; x[mm];y[mm]", 26, -1300., 1300., 26, -1300., 1300.));
    BookHisto(new TH1F("CHOD_trk_dist", "Distance beteen extrapolated track and position in the CHOD;CHOD_trkd [mm] ", 1500, 0., 3000.));
    BookHisto(new TH1F("CHOD_nt_timediff" , " Time difference between the associated hit and the others ; CHOD_{associated cl} - CHOD_{secondary cl} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("CHOD_nt_dtrk", "Distance beteen extrapolated track and position in the CHOD for the closest track;CHOD_trkd [mm] ", 150, 0., 300.));
    BookHisto(new TH2F("CHOD_cda_x_vs_y", "Distance beteen extrapolated track and position in the CHOD  x vs y;x[mm];y[mm]", 300, -300, 300., 300, -300., 300.));


    //RICH
    BookHisto(new TH1I("RICH_Ncandidates", "RICH number of candidates", 50, 0, 50));
    BookHisto(new TH1F("RICHRadius", "RICH radius", 500, 0., 500.));
    BookHisto(new TH1F("RICHMass", "Mass computed with RICH", 1000, 0., 1.));
    BookHisto(new TH2F("RICHMvsP", "Momentum vs RICH mass", 100, 0., 100., 1000, 0., 1.));
    BookHisto(new TH1F("RICHAngle", "RICH angle", 200, 0., 1.));
    BookHisto(new TH2F("RICHRvsP", "Energy as a function of particle", 100, 0., 100000., 500, 0., 500.));
    BookHisto(new TH1F("RICH_STRAW_dxdzdiff", "Slope difference between x slopes of RICH and STRAW", 200, -0.5, 0.5));
    BookHisto(new TH1F("RICH_STRAW_dydzdiff", "Slope difference between y slopes of RICH and STRAW", 200, -0.5, 0.5));
    BookHisto(new TH2F("RICH_dxdz_vs_dydz", "Slopes taken from the RICH ;dxdz[rad];dydz[rad]", 400, -0.1, 0.1, 400, -0.1, 0.1));
    BookHisto(new TH2F("RICH_x_vs_y", "Position of extrapolated track to the RICH front mirror;x[mm];y[mm]", 220, -1100, 1100., 220, -1100., 1100.));
    BookHisto(new TH2F("RICH_cda_x_vs_y", "Distance beteen extrapolated track and position in the RICH  x vs y", 220, -1100, 1100., 220, -1100., 1100.));

    //LKr
    BookHisto(new TH1I("LKr_Ncandidates", "LKr number of candidates;Ncandidates", 20, 0., 20.));
    BookHisto(new TH2F("LKr_x_vs_y", "LKr hitposition x vs y ;x[mm];y[mm]", 130 , -1300., 1300., 130, -1300., 1300.));
    BookHisto(new TH2F("LKr_cda_x_vs_y", "Distance beteen extrapolated track and position in the LKr  x vs y;x[mm];y[mm]", 300, -300, 300., 300, -300., 300.));
    BookHisto(new TH2F("LKr_Ecl_vs_NCell", "LKr cluster energy vs number of cells ; E_{cluster} [MeV]; Number of Cells", 1000, 0., 10000., 150, 0., 150.));
    BookHisto(new TH1F("LKr_EoP", "LKr E/p", 120, 0., 1.2));
    BookHisto(new TH1F("LKr_Ecl", "LKr cluster energy", 1000, 0., 10000.));
    BookHisto(new TH1F("LKr_Eseed_over_Ecl", "LKr E_{seed}/E_{cluster}; E_{seed}/E_{cluster}", 120, 0.,1.2));
    BookHisto(new TH2F("LKr_Es_Ecl_vs_E77_Ecl", "LKr clusters ; E_{seed}/E_{cluster}; 1 - E_{77}/E_{cluster}", 100, 0.,1., 100, 0.,1.));
    BookHisto(new TH1F("LKr_nearest_track_DDeadCell", "Distance to nearest dead cell in the LKr for the closest track;LKr_DDeadcell [cm] ", 500, 0., 250.));
    BookHisto(new TH1F("LKr_nt_timediff" , " Time difference between the associated LKr cluster and the others ; LKr_{associated cl} - LKr_{secondary cl} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("LKr_nt_dtrk", "Distance beteen the associated cluster in the LKr and the secondary clusters ;ClusterDistance [mm] ", 250, 0., 500.));
    //CEDAR
    BookHisto(new TH1I("CEDAR_Ncandidates", "CEDAR number of candidates; Ncandidates", 20, 0, 20));

    //Other
    BookHisto(new TH1F("Vertex_Z", " Z vertex ; Zvtx [mm]", 500, 0, 500000));
    BookHisto(new TH1F("Vertex_Y", " Y vertex ; Yvtx [mm]", 1000, -500., 500));
    BookHisto(new TH1F("Vertex_X", " X vertex ; Xvtx [mm]", 1000, -500., 500));
    BookHisto(new TH1F("Vertex_cda", " Closest distance approached between the kaon and the track ; cda [mm]", 300, 0, 300));

    BookHisto(new TH2F("STRAW1_x_vs_y", "STRAW hitposition @ Chamber 1 x vs y ; x[mm];y[mm]", 110, -1100., 1100., 110, -1100., 1100.));
    BookHisto(new TH2F("STRAW4_x_vs_y", "STRAW hitposition @ Chamber 4 x vs y ; x[mm];y[mm]", 110, -1100., 1100., 110, -1100., 1100.));

    //Information  about the nearest cluster for all detectors
    BookHisto(new TH1F("CHOD_nearest_track_dtrkcl", "Distance beteen extrapolated track and position in the CHOD for the closest track;CHOD_trkd [mm] ", 150, 0., 300.));
    BookHisto(new TH2F("CHOD_nearest_track_x_vs_y", "CHOD candidate x vs y ; x[mm];y[mm]", 520, -1300., 1300., 520, -1300., 1300.));
    BookHisto(new TH2F("CHOD_extrap_x_vs_y", "CHOD extrapolated track x vs y ; x[mm];y[mm]", 260, -1300., 1300., 260, -1300., 1300.));
    BookHisto(new TH1F("LKr_nearest_track_dtrkcl", "Distance beteen extrapolated track and cluster position in the LKr for the closest track;LKr_trkd [mm] ", 150, 0., 300.));
    BookHisto(new TH2F("LKr_nearest_track_x_vs_y", "LKr candidate x vs y  ;x[mm];y[mm]", 520 , -1300., 1300., 520, -1300., 1300.));
    BookHisto(new TH2F("LKr_extrap_x_vs_y", "LKr extrapolated track x vs y ;x[mm];y[mm]", 260 , -1300., 1300., 260, -1300., 1300.));
    BookHisto(new TH1F("MUV1_nearest_track_dtrkcl", "Distance beteen extrapolated track and cluster position in the MUV1 for the closest track;MUV1_trkd [mm] ", 150, 0., 300.));
    BookHisto(new TH2F("MUV1_extrap_x_vs_y", "Extrapolated x vs y position for the associated track in MUV1;x[mm];y[mm]", 260, -1320., 1320., 260, -1320., 1320.));
    BookHisto(new TH2F("MUV1_nearest_track_x_vs_y", "MUV1  x vs y candidate position of the associated track;x[mm];y[mm]", 260, -1320., 1320., 260, -1320., 1320.));
    BookHisto(new TH2I("MUV1_nearest_track_VvsH", "Associated track cluster Vertical channel ID (x)  vs cluster Horizontal channel ID (y) in MUV1", 44, 0, 44.,44, 0., 44.));
    BookHisto(new TH1F("MUV1_nearest_track_cluster_charge", "Charge of the cluster associated with the track at MUV1;MUV1_Q[fC]  ", 1000, 0., 10000.));
    BookHisto(new TH2F("MUV1_near_charge_vs_dtrkcl", " Cluster charge vs distance for the associated track at MUV1;MUV1_Q[fC];MUV1_dtrkcl  ", 1000, 0., 10000., 100., 0., 100.));
    BookHisto(new TH1F("MUV2_nearest_track_dtrkcl", "Distance beteen extrapolated track and cluster position in the MUV2 for the closest track;MUV2_trkd [mm] ", 150, 0., 300.));
    BookHisto(new TH2F("MUV2_extrap_x_vs_y", "Extrapolated x vs y position in MUV2;x[mm];y[mm]", 260, -1320., 1320., 260, -1320., 1320.));
    BookHisto(new TH2F("MUV2_nearest_track_x_vs_y", "Associated track cluster x vs cluster y position in MUV2", 440 , -1320., 1320., 440, -1320., 1320.));
    BookHisto(new TH2I("MUV2_nearest_track_VvsH", "Associated track cluster Vertical channel ID (x)  vs cluster Horizontal channel ID (y) in MUV2", 22, 0, 22.,22, 0., 22.));
    BookHisto(new TH1F("MUV2_nearest_track_cluster_charge", "Charge of the cluster associated with the track at MUV2;MUV2_Q[fC]  ", 1000, 0., 10000.));
    BookHisto(new TH2F("MUV2_near_charge_vs_dtrkcl", " Cluster charge vs distance for the associated track at MUV2;MUV2_Q[fC];MUV2_dtrkcl  ", 1000, 0., 10000., 200., 0., 200.));
    BookHisto(new TH1F("MUV3_nearest_track_dtrkcl", "Distance beteen extrapolated track and cluster position in the MUV3 for the closest track;MUV3_trkd [mm] ", 60, 0., 4000.));
    BookHisto(new TH2F("MUV3_extrap_x_vs_y", " MUV3 extrapolated x vs y position of the associated track ;x[mm];y[mm]", 260, -1320., 1320., 260, -1320., 1320.));
    BookHisto(new TH2F("MUV3_nearest_track_x_vs_y", " MUV3 x vs y candidate position of the associated track ;x[mm];y[mm]", 440, -1320., 1320., 440, -1320., 1320.));

    //Time differences using CHOD as the reference detector
    BookHisto(new TH1F("RICH_timediff" , " RICH_{time} - CHOD_{time} ; RICH_{time} - CHOD_{time} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("LKr_timediff"  , " LKr_{time} - CHOD_{time}  ; LKr_{time}  - CHOD_{time} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("MUV1_timediff" , " MUV1_{time} - CHOD_{time} ; MUV1_{time} - CHOD_{time} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("MUV2_timediff" , " MUV2_{time} - CHOD_{time} ; MUV2_{time} - CHOD_{time} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("MUV3_timediff" , " MUV3_{time} - CHOD_{time} ; MUV3_{time} - CHOD_{time} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("CEDAR_timediff"," CEDAR_{time} - CHOD_{time} ; CEDAR_{time}- CHOD_{time} [ns]", 400, -100, 100.));

    //Timedifferences between detectors using MUV3 as reference detector
    BookHisto(new TH1F("MUV3_MUV2timediff" , " MUV3_{time} - MUV2_{time} ; MUV3_{time} - MUV2_{time} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("MUV3_MUV1timediff" , " MUV3_{time} - MUV1_{time} ; MUV3_{time} - MUV1_{time} [ns]", 100, -50, 50.));

    BookHisto(new TH1I("MUV1_Ncandidates_MUV13", "MUV1 number of candidates", 50, 0, 50));
    BookHisto(new TH1I("MUV2_Ncandidates_MUV13", "MUV2 number of candidates", 50, 0, 50));
    BookHisto(new TH1I("MUV3_Ncandidates_MUV13", "MUV3 number of candidates", 50, 0, 50));
    BookHisto(new TH1I("MUV1_Ncandidates_MUV23", "MUV1 number of candidates", 50, 0, 50));
    BookHisto(new TH1I("MUV2_Ncandidates_MUV23", "MUV2 number of candidates", 50, 0, 50));
    BookHisto(new TH1I("MUV3_Ncandidates_MUV23", "MUV3 number of candidates", 50, 0, 50));
    BookHisto(new TH1I("MUV1_Ncandidates_MUV3", "MUV1 number of candidates", 50, 0, 50) );
    BookHisto(new TH1I("MUV2_Ncandidates_MUV3", "MUV2 number of candidates", 50, 0, 50) );
    BookHisto(new TH1I("MUV3_Ncandidates_MUV3", "MUV3 number of candidates", 50, 0, 50) );
    BookHisto(new TH2I("BadMUV2_HitMap_MUV13", " Hitmap for the MUV1&MUV3 !&MUV2 events", 23, 0, 22, 23, 0, 22 ));
    BookHisto(new TH2I("BadMUV1_HitMap_MUV13", " Hitmap for the MUV1&MUV3 !&MUV2 events", 45, 0, 44, 45, 0, 44) );
    BookHisto(new TH2I("BadMUV2_HitMap_MUV23", " Hitmap for the MUV2&MUV3 !&MUV1 events", 23, 0, 22, 23, 0, 22 ));
    BookHisto(new TH2I("BadMUV1_HitMap_MUV23", " Hitmap for the MUV2&MUV3 !&MUV1 events", 45, 0, 44, 45, 0, 44 ));
    BookHisto(new TH2I("BadMUV2_HitMap_MUV123", " Hitmap for the MUV2&MUV3&MUV1 events", 23, 0, 22, 23, 0, 22 ) );
    BookHisto(new TH2I("BadMUV1_HitMap_MUV123", " Hitmap for the MUV2&MUV3&MUV1 events", 45, 0, 44, 45, 0, 44)  );
    BookHisto(new TH2I("BadMUV1_HitMap_MUV3", " Hitmap for the MUV3 only events", 45, 0, 44, 45, 0, 44));
    BookHisto(new TH2I("BadMUV2_HitMap_MUV3", " Hitmap for the MUV3 only events", 23, 0, 22, 23, 0, 22));
    BookHisto(new TH1F("TrackPfit_TrackP_MUV13", "GetMomentum() - GetMomentumBeforeFit() ; Track_P[MeV] - Track_Ppat[MeV]", 100, -50000., 50000.));
    BookHisto(new TH1F("TrackPfit_TrackP_MUV23", "GetMomentum() - GetMomentumBeforeFit() ; Track_P[MeV] - Track_Ppat[MeV]", 100, -50000., 50000.));
    BookHisto(new TH1F("TrackPfit_TrackP_MUV3", "GetMomentum() - GetMomentumBeforeFit() ; Track_P[MeV] - Track_Ppat[MeV]", 100, -50000., 50000.) );
    BookHisto(new TH1F("TrackP_MUV13", "STRAW Momentum ; Track_P[MeV]", 100, 0., 100000.));
    BookHisto(new TH1F("TrackP_MUV23", "STRAW Momentum ; Track_P[MeV]", 100, 0., 100000.));
    BookHisto(new TH1F("TrackP_MUV3", "STRAW Momentum ; Track_P[MeV]", 100, 0., 100000.) );

    BookHisto(new TH1F("MUV3_LKr_tdiff_MUV3" , " MUV3_{time} - LKr_{time} for MUV3 only events ; MUV3_{time} - LKr_{time} [ns]", 100, -50, 50.)  );
    BookHisto(new TH1F("MUV3_LKr_tdiff_MUV13" , " MUV3_{time} - LKr_{time} for MUV3 only events ; MUV3_{time} - LKr_{time} [ns]", 100, -50, 50.) );
    BookHisto(new TH1F("MUV3_LKr_tdiff_MUV123" , " MUV3_{time} - LKr_{time} for MUV3 only events ; MUV3_{time} - LKr_{time} [ns]", 100, -50, 50.));
    BookHisto(new TH1F("MUV3_LKr_tdiff_MUV23" , " MUV3_{time} - LKr_{time} for MUV3 only events ; MUV3_{time} - LKr_{time} [ns]", 100, -50, 50.) );

    BookHisto(new TH1F("MUV3_nearest_track_dtrkcl_MUV123", "Distance beteen extrapolated track and cluster position in the MUV3 for the closest track;MUV3_trkd [mm] ", 60, 0., 4000.));
    BookHisto(new TH1F("MUV3_nearest_track_dtrkcl_MUV13", "Distance beteen extrapolated track and cluster position in the MUV3 for the closest track;MUV3_trkd [mm] ", 60, 0., 4000.) );
    BookHisto(new TH1F("MUV3_nearest_track_dtrkcl_MUV23", "Distance beteen extrapolated track and cluster position in the MUV3 for the closest track;MUV3_trkd [mm] ", 60, 0., 4000.) );
    BookHisto(new TH1F("MUV3_nearest_track_dtrkcl_MUV3", "Distance beteen extrapolated track and cluster position in the MUV3 for the closest track;MUV3_trkd [mm] ", 60, 0., 4000.)  );
    BookHisto(new TH1F("MUV1_nt_SW", " Shower width for the associated cluster in MUV1;MUV1_SW [mm] ", 500, 0., 500.));
    BookHisto(new TH2F("MUV1_TrP_SW"," Shower width vs Track P; Track_P[MeV];MUV1_SW[mm]", 100, 0, 100000., 500, 0., 500. ));
    BookHisto(new TH1F("MUV1_nt_timediff" , " Time difference between the associated cluster and the others ; MUV1_{associated cl} - MUV1_{secondary cl} [ns]", 100, -50, 50.));
    BookHisto(new TH2I("MUV1_sc_HitMap", " Hitmap for the secondary clusters in time with the associated cluster in MUV1", 45, 0, 44, 45, 0, 44 ));
    BookHisto(new TH1I("MUV1_hz_nt_chdiff", " Difference in horizontal channels between the cluster and the secondary clusters in MUV1", 44, -22, 22));
    BookHisto(new TH1I("MUV1_vt_nt_chdiff", " Difference in vertical channels between the cluster and the secondary clusters in MUV1", 44, -22, 22));
    BookHisto(new TH2I("MUV1_hz_vt_chdiff", " Channel difference between the secondary clusters and the associated cluster in MUV1;Horizontal Channel;Vertical Channel", 44, -22, 22, 44, -22, 22 ));
    BookHisto(new TH2F("MUV1_charge_vs_SW", " Cluster charge vs SW at MUV1;MUV1_Q[fC];MUV1_SW[mm]  ", 1000, 0., 10000., 500., 0., 500.));

    BookHisto(new TH2F("MUV1_PvsQ", " Charge vs Momentum of the cluster associated with a track at MUV1; TrackP[MeV];MUV1_Q[pC]  ", 100, 0., 100000., 1000., 0., 10000.));
    BookHisto(new TH2F("MUV2_PvsQ", " Charge vs Momentum of the cluster associated with a track at MUV2; TrackP[MeV];MUV2_Q[pC]  ", 100, 0., 100000., 1000., 0., 10000.));

    BookHisto(new TH1F("MUV1_zero_distance_timediff" , " Time difference between the associated and sc, which are very close to eachother  ; MUV1_{associated cl} - MUV1_{secondary cl} [ns]", 100, -50, 50.));


    BookHisto(new TH1F("MUV2_nt_SW", " Shower width for the associated cluster in MUV2;MUV2_SW [mm] ", 500, 0., 500.));
    BookHisto(new TH2F("MUV2_TrP_SW"," Shower width vs Track P; Track_P[MeV];MUV2_SW[mm]", 100, 0, 100000., 500, 0., 500. ));
    BookHisto(new TH1F("MUV2_nt_timediff" , " Time difference between the associated cluster and the others ; MUV2_{associated cl} - MUV2_{secondary cl} [ns]", 100, -50, 50.));
    BookHisto(new TH2I("MUV2_sc_HitMap", " Hitmap for the secondary clusters in time with the associated cluster", 23, 0, 22, 23, 0, 22 ));
    BookHisto(new TH1I("MUV2_hz_nt_chdiff", " Difference in horizontal channels between the cluster and the secondary clusters in MUV2", 22, -11, 11));
    BookHisto(new TH1I("MUV2_vt_nt_chdiff", " Difference in vertical channels between the cluster and the secondary clusters in MUV2", 22, -11, 11));
    BookHisto(new TH2I("MUV2_hz_vt_chdiff", " Channel difference between the secondary clusters and the associated cluster in MUV2;Horizontal Channel;Vertical Channel", 22, -11, 11, 22, -11, 11 ));
    BookHisto(new TH2F("MUV2_charge_vs_SW", " Cluster charge vs SW at MUV2;MUV2_Q[fC];MUV2_SW[mm]  ", 1000, 0., 10000., 500., 0., 500.));
    BookHisto(new TH1F("MUV2_zero_distance_timediff" , " Time difference between the associated and sc, which are very close to eachother  ; MUV2_{associated cl} - MUV2_{secondary cl} [ns]", 100, -50, 50.));

    //Hits histograms for checking for readout failure

    BookHisto(new TH1I("Nhits123_MUV2", "MUV2 number of hits for  MUV1+2+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits123_MUV1", "MUV1 number of hits for  MUV1+2+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits123_LKr",  "LKr  number of hits for  MUV1+2+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits0C13_MUV2", "MUV2 number of hits for events without reconstructed cluster for MUV1+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits0C13_MUV1", "MUV1 number of hits for events without reconstructed cluster for MUV1+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits0C13_LKr",  "LKr  number of hits for events without reconstructed cluster for MUV1+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits0C23_MUV2", "MUV2 number of hits for events without reconstructed cluster for MUV2+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits0C23_MUV1", "MUV1 number of hits for events without reconstructed cluster for MUV2+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits0C23_LKr",  "LKr  number of hits for events without reconstructed cluster for MUV2+3 events", 50, 0, 50));

    BookHisto(new TH1I("Nhits0C23_MUV1_BB", "MUV1  number of hits for events without reconstructed cluster for the inefficient bursts MUV2+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits0C13_MUV2_BB", "MUV2  number of hits for events without reconstructed cluster for the inefficient bursts MUV1+3 events", 50, 0, 50));
    BookHisto(new TH1I("Nhits0C3_MUV2_BB",  "MUV2  number of hits for events without reconstructed cluster for the inefficient bursts MUV3 events", 50, 0, 50)  );
    BookHisto(new TH1I("Nhits0C3_MUV1_BB",  "MUV1  number of hits for events without reconstructed cluster for the inefficient bursts MUV3 events", 50, 0, 50)  );

    BookHisto(new TH1I("0C_ChID_MUV2",  "Channel ID of the hits from the inefficient MUV2 events (MUV1+3)", 200, 100, 300) );
    BookHisto(new TH1I("0C_VChID_MUV2",  "Vertical Channel ID of the hits from the inefficient MUV2 events (MUV1+3)", 22, 1, 23) );
    BookHisto(new TH1I("0C_HChID_MUV2",  "Horizontal Channel ID of the hits from the inefficient MUV2 events (MUV1+3)", 22, 1, 23) );
    BookHisto(new TH1I("0C_ChID_MUV1",  "Channel ID of the hits from the inefficient MUV1 events (MUV2+3)", 200, 100, 300) );
    BookHisto(new TH1I("0C_HChID_MUV1",  "Horizontal Channel ID of the hits from the inefficient MUV1 events (MUV2+3)", 44, 1, 45) );
    BookHisto(new TH1I("0C_VChID_MUV1",  "Vertical Channel ID of the hits from the inefficient MUV1 events (MUV2+3)", 44, 1, 45) );
    BookHisto(new TH1I("0C_VM1_CHOD_t",  "Vertical Channel Time - CHOD time of the hits from the inefficient MUV1 events (MUV2+3)",  800, -400, 400) );
    BookHisto(new TH1I("0C_HM1_CHOD_t",  "Horizontal Channel Time - CHOD time of the hits from the inefficient MUV1 events (MUV2+3)",800, -400, 400));
    BookHisto(new TH1I("0C_VM2_CHOD_t",  "Vertical Channel Time - CHOD time of the hits from the inefficient MUV2 events (MUV1+3)",  800, -400, 400));
    BookHisto(new TH1I("0C_HM2_CHOD_t",  "Horizontal Channel Time - CHOD time of the hits from the inefficient MUV2 events (MUV1+3)",800, -400, 400));

    BookHisto(new TH1I("0C_HChID_diff_M1",  " Horizontal ChannelID_{Hits} - ChannelID_{extrap} for inefficient MUV1 events (MUV2+3)", 80, -40, 40) );
    BookHisto(new TH1I("0C_VChID_diff_M1",  " Vertical ChannelID_{Hits} - ChannelID_{extrap} for inefficient MUV1 events (MUV2+3)", 80, -40, 40)   );
    BookHisto(new TH1I("0C_HChID_diff_M2",  " Horizontal ChannelID_{Hits} - ChannelID_{extrap} for inefficient MUV2 events (MUV1+3)", 40, -20, 20) );
    BookHisto(new TH1I("0C_VChID_diff_M2",  " Vertical ChannelID_{Hits} - ChannelID_{extrap} for inefficient MUV2 events (MUV1+3)", 40, -20, 20)   );

    //Quality checks for Gia`s reconstruction
    BookHisto(new TH1I("Quality",  " MUV1 fQuality variable: 0 - true cluster 1 - time-charge information ambiguous 2 - wrongly reconstructed", 5, 0, 5));
    BookHisto(new TH2F("Q0_nearest_track_x_vs_y", "X vs Y position from MUV1 for fQuality = 0;x[mm];y[mm]", 436, -1308., 1308., 436, -1308., 1308.));
    BookHisto(new TH2F("Q1_nearest_track_x_vs_y", "X vs Y position from MUV1 for fQuality = 1;x[mm];y[mm]", 436, -1308., 1308., 436, -1308., 1308.));
    BookHisto(new TH2F("Q2_nearest_track_x_vs_y", "X vs Y position from MUV1 for fQuality = 2;x[mm];y[mm]", 436, -1308., 1308., 436, -1308., 1308.));

    //25ns difference hits
    BookHisto(new TH1I("ChannelID_25ns_away_M1",  "ChannelID for the hits that are 25 ns away (MUV2+3)", 200, 100, 300));
    BookHisto(new TH1I("ChannelID_25ns_away_M2",  "ChannelID for the hits that are 25 ns away (MUV1+3)", 200, 100, 300));


    //Saving hits
    BookHisto(new TH1I("MUV13_Vsaved_M2" , "No reco cluster in MUV2, hit at the extrapolated vertical strip (MUV1+3 events)", 22, 1, 23) );
    BookHisto(new TH1I("MUV13_Hsaved_M2" , "No reco cluster in MUV2, hit at the extrapolated horizontal strip (MUV1+3 events)", 22, 1, 23) );
    BookHisto(new TH1I("RecVHits_M2" , "Number of hits in the Vertical MUV2 channels that are found in the extrapolated strips (MUV1+3 events)", 20, 0, 20));
    BookHisto(new TH1I("RecHHits_M2" , "Number of hits in the Horizontal MUV2 channels that are found in the extrapolated strips (MUV1+3 events)", 20, 0, 20) );
    BookHisto(new TH1I("RecHits_M2" , " Events that have at least one hit in the Vertical and Horizontal channel (0 reconstructed clusters) in MUV2 (MUV1+3 events)", 20, 0, 20) );
    BookHisto(new TH1I("MUV23_Vsaved_M1" , "No reco cluster in MUV1, hit at the extrapolated vertical strip (MUV2+3 events)", 44, 1, 45) );
    BookHisto(new TH1I("MUV23_Hsaved_M1" , "No reco cluster in MUV1, hit at the extrapolated horizontal strip (MUV2+3 events)", 44, 1, 45) );
    BookHisto(new TH1I("RecVHits_M1" , "Number of hits in the Vertical MUV1 channels that are found in the extrapolated strips (MUV2+3 events)", 20, 0, 20) );
    BookHisto(new TH1I("RecHHits_M1" , "Number of hits in the Horizontal MUV1 channels that are found in the extrapolated strips (MUV2+3 events)", 20, 0, 20) );
    BookHisto(new TH1I("RecHits_M1" , "Events that have at least one hit in the Vertical and Horizontal channel (0 reconstructed clusters) MUV1 (MUV2+3 events)", 20, 0, 20) );


}

void Kmu2::DefineMCSimple(){
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

void Kmu2::StartOfRunUser(){
    /// \MemberDescr
    /// This method is called at the beginning of the processing (corresponding to a start of run in the normal NA62 data taking)\n
    /// Do here your start of run processing if any
    /// \EndMemberDescr
}

void Kmu2::StartOfBurstUser(){
    /// \MemberDescr
    /// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the beginning of the first file\n
    /// Do here your start/end of burst processing if any
    /// \EndMemberDescr
}

void Kmu2::Process(int iEvent){
    /// \MemberDescr
    /// \param iEvent : Event number
    ///
    /// Main process method. Called on each event. Write you analysis here.\n
    /// You can retrieve MC particles from the fMCSimple set with (returns a vector<KinePart*>)\n
    /// \code
    ///     fMCSimple[ particleName""]
    ///     fMCSimple[pdgID]
    /// \endcode
    /// Example\n
    /// \code
    ///     fMCSimple["K+"][index]; //for the kaon
    ///     fMCSimple["pi+"][index]; //for the positive pion
    ///     fMCSimple["gamma"][index]; //for the photon
    /// \endcode
    /// The number in the brackets is the index of the particle (if you asked for two photons in the set, you can ask fMCSimple["gamma"][0] for the first one and fMCSimple["gamma"][1] for the second)\n
    /// \n
    /// If you need a property of a particle, you can make a call to fParticleInterface (instance of the ParticleInterface class).\n
    ///     This class has two methods FindParticle that will return a TParticlePDG with the required particle. You can search by pdgID or by name.\n
    ///     This class also provide two methods to switch between particle name and pdgID if necessary.\n
    ///     Example\n
    /// \code
    ///     double kaonMass = fParticleInterface->FindParticle(321).Mass();
    ///     double pi0Lifetime = fParticleInterface->FindParticle("pi0").Lifetime();
    /// \endcode
    /// You can retrieve the events from the trees with\n
    /// \code
    ///     (eventClass*)GetEvent("detectorName");
    ///     (eventClass*)GetEvent("detectorName", "Digis");
    /// \endcode
    /// You can retrieve data from generic TTrees with\n
    /// \code
    ///     GetObject<MyClass>("treeName");
    /// \endcode
    /// You can retrieve full MC events if available ( GetWithMC() ) with\n
    /// \code
    ///     GetMCEvent();
    ///     GetMCEvent("Digis");
    /// \endcode
    /// You can retrieve RawHeader if available ( GetWithRawHeader() ) with\n
    /// \code
    ///     GetRawHeader();
    ///     GetRawHeader("Digis");
    /// \endcode
    /// You can retrieve Trigger data if requested with\n
    /// \code
    ///     GetL0Data();
    ///     GetL1Data();
    ///     GetL2Data();
    /// \endcode
    /// You can retrieve the histograms you booked (for drawing, changing, filling, ...) with\n
    /// \code
    ///     fHisto.GetTH1("histoName");// for TH1
    ///     fHisto.GetTH2("histoName");// for TH2
    ///     fHisto.GetTGraph("graphName");// for TGraph and TGraphAsymmErrors
    ///     fHisto.GetHisto("histoName");// for TH1 or TH2 (returns a TH1 pointer)
    /// \endcode
    /// To fill the histograms you can use\n
    /// \code
    ///     FillHisto("histoName", values)
    /// \endcode
    /// where values are the same parameters as if you call histogram->Fill(values) (x,y,weight,...)\n
    /// If the histogram is not found, an error message is printed\n
    /// \n
    /// Modify a counter with one of the following methods\n
    /// \code
    ///     IncrementCounter(name)
    ///     IncrementCounter(name, delta)
    ///     DecrementCounter(name)
    ///     DecrementCounter(name, delta)
    ///     SetCounterValue(name, value)
    /// \endcode
    /// \n
    /// For use of fGeom, read DetectorAcceptance class.\n
    ///     WARNING: this class provides "exact" results, there is not tolerance. If the particle\n
    ///     passes in the sensitive volume of a detector it's considered as detected, wether it's close\n
    ///     to the edge or not. But as the class gives you the position of the passage point and the \n
    ///     estimated path length in the sensitive volume, you can apply further cuts from this\n
    ///     information at your convenience.\n
    /// \n
    /// To use the output of a different analyzer, use\n
    /// \code
    ///     outputType *var = GetOutput<outputType>("analyzerName.outputName", state);
    /// \endcode
    /// Where outputType is the variable type and state is of type outputState\n
    /// State is set with the state of the variable (kOUninit, kOInvalid ,kOValid). The value of the output should only be trusted if state == kOValid\n
    /// example :
    /// \code
    ///     TLorentzVector vertex = *(TLorentzVector*)GetOutput("simpleVertexAnalyzer.vertex", state);
    /// \endcode
    /// Before starting the processing of an event, the state flag of each output variable is reset to kOUninit\n
    /// When setting the value of an output variable, don't forget to set appropriately the state flag to either kOValid or kOInvalid\n
    /// to indicate if the value can/can't be used in other analyzer\n
    /// \code
    ///     SetOutputState("outputName", kOValid);
    /// \endcode
    /// If you want to append a candidate in one of your standard output Tree, use\n
    /// \code
    ///     KinePart *candidate = CreateStandardCandidate("treeName");
    /// \endcode
    /// and fill the properties of your candidate. It will be automatically written in the output tree.\n
    ///     \n
    /// If you want to save this event in your custom and standard TTrees (not the input tree replication), call\n
    /// \code
    ///     FillTrees();
    /// \endcode
    /// This will call the Fill method of every TTree created in this analyzer.\n
    ///     \n
    /// If you want to replicate this event in the output file, call\n
    /// \code
    ///     ExportEvent();
    /// \endcode
    /// The structure of all the trees that have been opened (by all Analyzer) will be copied in the output file\n
    /// and the events for which at least one analyzer called ExportEvent() will be replicated in the output trees.
    /// @see ROOT TParticlePDG for the particle properties
    /// @see ROOT TDatabasePDG for a list of PDG codes and particle naming convention
    /// \EndMemberDescr
    //if(fMCSimple.fStatus == MCSimple::kMissing){printIncompleteMCWarning(iEvent);return;}
    //      if(fMCSimple.fStatus == MCSimple::kEmpty){printNoMCWarning();return;}
    TRecoLKrEvent          *LKrEvent = (TRecoLKrEvent*)GetEvent("LKr");
    TRecoSpectrometerEvent *SpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");
    TRecoMUV1Event         *MUV1Event   = (TRecoMUV1Event*)GetEvent("MUV1");
    TRecoMUV2Event         *MUV2Event   = (TRecoMUV2Event*)GetEvent("MUV2","Reco");
    TRecoMUV3Event         *MUV3Event   = (TRecoMUV3Event*)GetEvent("MUV3");
    TRecoRICHEvent         *RICHEvent   = (TRecoRICHEvent*)GetEvent("RICH");
    TRecoCHODEvent         *CHODEvent   = (TRecoCHODEvent*)GetEvent("CHOD");
    TRecoCedarEvent        *CedarEvent = (TRecoCedarEvent*)GetEvent("Cedar");

    TRecoLKrCandidate*   LKrCluster;
    TRecoCHODCandidate*  CHODCandidate;
    TRecoCedarCandidate* CedarCandidate;
    TRecoRICHCandidate*  RingCandidate;
    TRecoMUV1Candidate*  MUV1Cluster;
    TRecoMUV2Candidate*  MUV2Cluster;
    TRecoMUV3Candidate*  MUV3Cluster;
    TRecoSpectrometerCandidate* Track;

    //Declaring the vectors for extrapolation of the the tracks
    //from the spectrometers to the corresponding detectors
    // All vectors with _extrap (example CHOD_extrap)!!!!!!
    TVector3 CHOD_extrap;
    TVector3 MUV1_extrap;
    TVector3 MUV2_extrap;
    TVector3 RICH_extrap;
    TVector3 LKr_extrap;
    TVector3 MUV3_extrap;

    //Slopes and Positions given by the spectrometer
    //Before the magnet (bdxdz and bdydz)
    TVector3 SlopesBefore;
    TVector3 PositionBefore;

    //Slopes given by the spectrometer after the magnet (dxdz and dydz)
    TVector3 SlopesAfter;
    TVector3 PositionAfter;

    //Beam position vector at trim5 and beam momentum
    //(at the moment no GTK assuming 75 GeV Definitions.h)
    TVector3 BeamTrim5Pos;
    TVector3 BeamP;

    //3Momentum of the charged track
    TVector3 TrackP;

    //3Momentum of the missing mass (BeamP - TrackP)
    TVector3 NuNubar;

    //Vertex using cda routine
    TVector3 Vertex;

    TVector2 MUV1Pos;
    TVector2 MUV2Pos;
    TVector2 CHODPos;
    TVector3 LkrPos;


    //Time Offset for all the detectors differences (ATM using only CHOD as reference)
    double LKrOffset   = +115; //Old 112.3 Use GetClusterTime new 115 Gia 121.2
    double CedarOffset =0.; //old 0 new 0 Gia 6.23
    double MUV1Offset  =0.;//new -5 old -8.5 Gia 6.08
    double MUV2Offset  =0.;//new -5 old -20.4 Gia 5.59
    double MUV3Offset  =0.;//old -11.35 Gia -2.75

    //Cuts used on the timedifferences for the Kmu2 selection
    double LKrOffsetCut  = 100.;//10;
    double RICHOffsetCut = 100.;//10;
    double MUV3OffsetCut = 100.;//5 ;
    double MUV1OffsetCut = 100.;//20;
    double MUV2OffsetCut = 100.;//20;
    double CedarOffsetCut= 100.;//2 ;

    //CUTComment:: Only one candidate in the STRAW
    if(SpectrometerEvent->GetNCandidates() != 1 ){return;}




    Track = ((TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0));

    //CUTComment:: Check if track is positive (K+ beam in NA62)
    if(Track->GetCharge() != 1){ return;}
    //cout << "MUV L0 trig == " << CHODEvent->GetL0TriggerType() << "STRAW L0 trig == " << MUV3Event->GetTriggerType() << endl;

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
    if(STRAW_chi2 > 20){ return;}
    if(STRAW_NC   < 3){ return;}

    for(int iCedarCand=0; iCedarCand < CedarEvent->GetNCandidates(); iCedarCand++){
        CedarCandidate = ((TRecoCedarCandidate*)CedarEvent->GetCandidate(iCedarCand));
        //CUTComment:: At least 5 sectors in CEDAR
        if(CedarCandidate->GetNSectors() < 5 ){return;}

    }


    //Getting the position vectors and the slopes of the tracks
    //before the magnet @ DCH1
    SlopesBefore.SetX(STRAW_bdxdz);
    SlopesBefore.SetY(STRAW_bdydz);
    SlopesBefore.SetZ(1.);
    PositionBefore = Track->GetPositionBeforeMagnet();

    //Building track momentum from STRAW slopes before magnet
    //and the position @ DCH1
    double norm = 1./sqrt(pow(SlopesBefore.X(),2) + pow(SlopesBefore.Y(),2) + 1.  );
    TrackP.SetXYZ(norm*SlopesBefore.X()*STRAW_P, norm*SlopesBefore.Y()*STRAW_P, norm*STRAW_P );


    //Beam momentum and position on Trim5 that will be used
    //for the making of the vertex
    double beam_norm = 1./sqrt(XAngle*XAngle + 1. );
    BeamP.SetXYZ(beam_norm*KEnergy*XAngle*1000,0,beam_norm*KEnergy*1000.);
    BeamTrim5Pos.SetXYZ(0., 0., Ztrim*1000.);

    //Calculating the intersection point between the kaon and
    //the track and getting the cda using the VertexCDA routine by Giuseppe
    double cda = 0;
    Vertex = VertexCDA(PositionBefore, TrackP, BeamTrim5Pos, BeamP, cda );


    //CUTComment:: Closest Approached Distance > 40.
    //Zvtx to be incide of the fiducial volume of NA62 detector (105 - 180 m)
    if(cda > 40.){return;}
    if( Vertex.Z() < 105000 || Vertex.Z() > 180000){return;}

    //Computing the NuNubar 3vector (Missing mass)
    NuNubar = BeamP - TrackP;

    //Making the 4Momentum of NuNuBar
    TLorentzVector PiP;
    TLorentzVector KP;
    TLorentzVector NuNubarP;
    double theta;

    //4Momentum of the track (assuming muon mass for Kmu2 selection)
    PiP.SetVect(TrackP);
    PiP.SetE( sqrt(pow(TrackP.Mag(),2) + pow(MuMass*1000,2)));

    //4Momentum of the Kaon (NO GTK yet)
    KP.SetVect(BeamP);
    KP.SetE(sqrt(pow(BeamP.Mag(),2) + pow(KMass*1000,2)));

    //4Momentum of the NuNubar system
    NuNubarP = KP - PiP;

    //Angle between the track and the Kaon
    theta = TrackP.Dot(BeamP)/(TrackP.Mag()*BeamP.Mag());


    //Getting the position vectors and the slopes of the tracks
    //after the magnet @ DCH4
    SlopesAfter.SetX(STRAW_dxdz);
    SlopesAfter.SetY(STRAW_dydz);
    SlopesAfter.SetZ(1.);
    PositionAfter = Track->GetPositionAfterMagnet();

    //Extrapolating the track to the other detectors
    RICH_extrap = PositionAfter + ( ZRICHStart*1000 - PositionAfter.Z() )*SlopesAfter;
    CHOD_extrap = PositionAfter + ( ZCHODStart*1000 - PositionAfter.Z() )*SlopesAfter;
    MUV1_extrap = PositionAfter + ( ZMUV1Start*1000 - PositionAfter.Z() )*SlopesAfter;
    MUV2_extrap = PositionAfter + ( ZMUV2Start*1000 - PositionAfter.Z() )*SlopesAfter;
    MUV3_extrap = PositionAfter + ( ZMUV3Start*1000 - PositionAfter.Z() )*SlopesAfter;
    LKr_extrap  = PositionAfter + ( ZLKrStart *1000 - PositionAfter.Z() )*SlopesAfter;

    //Energy scale correction and non-linearity correction for the LKr taken from Giuseppe
    Double_t fEScale = 1.03;
    for(int iLKrCand=0; iLKrCand<LKrEvent->GetNCandidates(); iLKrCand++){
        LKrCluster = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(iLKrCand));
        // ZS non linearity
        Double_t ue = LKrCluster->GetClusterEnergy();
        Double_t ce = ue;
        if (LKrCluster->GetNCells()>9) {
            if (ue<22) ce = ue/(0.7666+0.0573489*log(ue));
            if (ue>=22 && ue<65) ce = ue/(0.828962+0.0369797*log(ue));
            if (ue>=65) ce = ue/(0.828962+0.0369797*log(65));
        }

        LKrCluster->SetClusterEnergy(ce*fEScale);
        double LKrClusterEnergy = 1000*LKrCluster->GetClusterEnergy(); //  [MeV]

        //CUTComment:: MIP cluster requirement on number in LKr and the energy inside them
        if(LKrClusterEnergy > 800 ) {return ;} //  [MeV]
        if(LKrCluster->GetNCells()>5) {return;}
        if(LKrClusterEnergy < 300 ) {return;}
    }
    //CUTComment:: At least one track in MUV3 and in CEDAR
    if( MUV3Event->GetNCandidates() == 0 ){return;}
    if( CedarEvent->GetNCandidates() == 0 ){return;}


    //Passing the momentum of the particles together with the slopes after the magnet
    //to the function FindClosestCluster which gives you the index of the cluster and
    //the value of the closest distance between the extrapolated track and the cluster
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


    //CUTComment::At least one track associated with hit in the CHOD
    //Missing mass cut to select pure muons
    //Momentum measured by the STRAW to be between 10 and 65 GeV
    //Momentum given by the STRAW - Momentum given just by the pattern recognition
    // to be less than 20 GeV

    if(CHODClosestTrackIndex < 0. ){return;}
    if(MUV3TrackClusterIndex < 0. ){return;}
    if(fabs(NuNubarP.M2()*0.000001) > 0.01 ){return;}
    if(STRAW_P < 10000. || STRAW_P > 65000.){return;}
    if(fabs(STRAW_P - STRAW_Pbf) > 20000.){return;}
    double CHODR  = sqrt( pow(CHOD_extrap.X(),2) + pow(CHOD_extrap.Y(),2) );
    double STRAW4R  = sqrt( pow(PositionAfter.X(),2) + pow(PositionAfter.Y(),2) );
    double CD_CHODTime    = ((TRecoCHODCandidate*)CHODEvent->GetCandidate(CHODClosestTrackIndex))->GetTime();
    TVector2 CD_CHODPos   = ((TRecoCHODCandidate*)CHODEvent->GetCandidate(CHODClosestTrackIndex))->GetHitPosition();

    //CUTComment:: Distance to the extrapolated track bigger than 8 cm
    //Tracks to be in the CHOD and DCH 4 acceptance
    if(CHODdtrkcl_min > 80. ) {return;}
    if(CHODR < 100. || CHODR > 1200) {return;} //[mm]
    if(STRAW4R < 75. || STRAW4R > 1200) {return;} //[mm]
    //if(PositionAfter.Mag() < 120 || PositionAfter.Mag() > 1100) {return;}//[mm]
    //cout << PositionAfter.Mag() << endl;


    for(int iCHODCand=0; iCHODCand<CHODEvent->GetNCandidates(); iCHODCand++){

        CHODCandidate     = ((TRecoCHODCandidate*)CHODEvent->GetCandidate(iCHODCand));
        TRecoCHODCandidate* CHODNtCandidate = ((TRecoCHODCandidate*)CHODEvent->GetCandidate(CHODClosestTrackIndex));

        CHODPos           = CHODCandidate->GetHitPosition();
        TVector2 CHODNtPos= CHODNtCandidate->GetHitPosition();
        double CHODTime   = CHODCandidate->GetTime();
        double CHODntTime = CHODNtCandidate->GetTime();
        double CHOD_dtrk  = sqrt(pow(CHODPos.X()*10. - CHOD_extrap.X(), 2 ) + pow(CHODPos.Y()*10. - CHOD_extrap.Y(), 2 ) ) ;

        FillHisto("CHOD_trk_dist", CHOD_dtrk);
        FillHisto("CHOD_x_vs_y", CHODPos.X()*10.,CHODPos.Y()*10.);
        if(iCHODCand == CHODClosestTrackIndex){continue;}
        if(fabs(CHODntTime - CHODTime) < 5 && sqrt(pow(CHODPos.X()*10. - CHODNtPos.X()*10, 2 ) + pow(CHODPos.Y()*10. - CHODNtPos.Y()*10., 2 ) ) < 100){return;}
        FillHisto("CHOD_nt_timediff", CHODntTime - CHODTime);
        FillHisto("CHOD_nt_dtrk", sqrt(pow(CHODPos.X()*10. - CHODNtPos.X()*10, 2 ) + pow(CHODPos.Y()*10. - CHODNtPos.Y()*10., 2 ) ));

    }

//Do the same procedure for track cluster matching for CEDAR, but
    //matching only in time
    for(int iCedarCand=0; iCedarCand < CedarEvent->GetNCandidates(); iCedarCand++){
        CedarCandidate = ((TRecoCedarCandidate*)CedarEvent->GetCandidate(iCedarCand));
        double CedarCandidateTime = CedarCandidate->GetTime();
        double CedarTimeDiff      = CD_CHODTime - CedarCandidateTime + CedarOffset;
        CEDAR_STRAW_tdiff.push_back(CedarTimeDiff);
    }

    if(CEDAR_STRAW_tdiff.size() !=0){
        double CedarTime = *min_element(CEDAR_STRAW_tdiff.begin(), CEDAR_STRAW_tdiff.end());


        //CUTComment:: Cedar time difference cut
        if(fabs(CedarTime) > CedarOffsetCut){return;}
        FillHisto("CEDAR_timediff", CedarTime);
    }



    if(MUV1TrackClusterIndex > -1){
        TRecoMUV1Candidate* CD_MUV1Cluster = ((TRecoMUV1Candidate*)MUV1Event->GetCandidate(MUV1TrackClusterIndex));
        double CD_MUV1ClusterTime = ((TRecoMUV1Candidate*)MUV1Event->GetCandidate(MUV1TrackClusterIndex))->GetTime();
        double MUV1T0 = CD_CHODTime - CD_MUV1ClusterTime + MUV1Offset;
        //TClonesArray *MUV1Hits = MUV1Event->GetHits();
        //double MUV1Cluster_Charge=0;
        //double MUV1Cluster_SW=CD_MUV1Cluster->GetShowerWidth();



        // FillHisto("MUV1_Nhits", CD_MUV1Cluster->GetNHits() );
        // FillHisto("MUV1_timediff", MUV1T0);
        // FillHisto("MUV1_nearest_track_dtrkcl", MUV1dtrkcl_min);
        // FillHisto("MUV1_nearest_track_cluster_charge", MUV1Cluster_Charge);
        // FillHisto("MUV1_near_charge_vs_dtrkcl", MUV1Cluster_Charge, MUV1dtrkcl_min);
        // FillHisto("MUV1_nearest_track_x_vs_y", MUV1_extrap.X(), MUV1_extrap.Y() );
        // //Old Reco
        // //FillHisto("MUV1_cda_x_vs_y", CD_MUV1Cluster->GetPosition().X() - MUV1_extrap.X() + 60., CD_MUV1Cluster->GetPosition().Y() - MUV1_extrap.Y() + 60. );
        // //New Reco
        // FillHisto("MUV1_cda_x_vs_y", CD_MUV1Cluster->GetPosition().X() - MUV1_extrap.X(), CD_MUV1Cluster->GetPosition().Y() - MUV1_extrap.Y() );
        // FillHisto("MUV1_nt_SW", MUV1Cluster_SW);
        // FillHisto("MUV1_TrP_SW", STRAW_P ,MUV1Cluster_SW);


        if(fabs(MUV1T0) > MUV1OffsetCut){return;}
        if(MUV1dtrkcl_min > 100.) {return;}
        if(fabs(MUV1_extrap.X()) <= 130. && fabs(MUV1_extrap.Y()) <= 130.){return;}
        if(fabs(MUV1_extrap.X()) >= 1100. || fabs(MUV1_extrap.Y()) >= 1100.){return;}
        //Getting the closest cluster to the track, which is in a square with a side 40mm (+60 mm for old reconstruction)
        //Old Reco
        //if( fabs(CD_MUV1Cluster->GetPosition().X() - MUV1_extrap.X() + 60.) > 160. ||
        //    fabs(CD_MUV1Cluster->GetPosition().Y() - MUV1_extrap.Y() + 60.) > 160.){return;}
        //New Reco
        if( fabs(CD_MUV1Cluster->GetPosition().X() - MUV1_extrap.X()) > 160. ||
            fabs(CD_MUV1Cluster->GetPosition().Y() - MUV1_extrap.Y()) > 160.){return;}
        // if(MUV1_extrap.X() <= 90. || MUV1_extrap.Y() >= -24.){
        if(MUV1_extrap.X() <= 90. && MUV1_extrap.X() >= -24.){
            //cout << "GetPos.X == " << CD_MUV1Cluster->GetPosition().X() << "GetPos.Y == " << CD_MUV1Cluster->GetPosition().Y() << " WTFFF"<< endl;
            //cout << "GetVChannel == " << CD_MUV1Cluster->GetVerticalChannel() << "GetHChannel == " << CD_MUV1Cluster->GetHorizontalChannel() << " WTFFF" << endl;
        }
    }

    if(MUV2TrackClusterIndex > -1){

        TRecoMUV2Candidate* CD_MUV2Cluster = ((TRecoMUV2Candidate*)MUV2Event->GetCandidate(MUV2TrackClusterIndex));
        double CD_MUV2ClusterTime   = ((TRecoMUV2Candidate*)MUV2Event->GetCandidate(MUV2TrackClusterIndex))->GetTime();
        double MUV2T0 = CD_CHODTime - CD_MUV2ClusterTime + MUV2Offset;
        //double MUV2Cluster_Charge=0;
        //double MUV2Cluster_SW=CD_MUV2Cluster->GetShowerWidth();
        //TClonesArray *MUV2Hits = MUV2Event->GetHits();
        ////double MUV2R  = sqrt( pow(MUV2_extrap.X(),2) + pow(MUV2_extrap.Y(),2) );
        //for(int iMUV2Hit = 0; iMUV2Hit <CD_MUV2Cluster->GetNHits(); iMUV2Hit++ ){
        //    TRecoMUV2Hit* MUV2Hit = ((TRecoMUV2Hit*)MUV2Hits->At(iMUV2Hit));
        //    double MUV2HitCharge = MUV2Hit->GetCharge();
        //    MUV2Cluster_Charge  += MUV2HitCharge;
        //}

        // FillHisto("MUV2_Nhits", CD_MUV2Cluster->GetNHits() );
        // FillHisto("MUV2_timediff", MUV2T0);
        // FillHisto("MUV2_nearest_track_dtrkcl", MUV2dtrkcl_min);
        // FillHisto("MUV2_nearest_track_cluster_charge", MUV2Cluster_Charge);
        // FillHisto("MUV2_near_charge_vs_dtrkcl", MUV2Cluster_Charge, MUV2dtrkcl_min);
        // FillHisto("MUV2_nearest_track_x_vs_y", MUV2_extrap.X(), MUV2_extrap.Y() );
        // FillHisto("MUV2_cda_x_vs_y", CD_MUV2Cluster->GetPosition().X() - MUV2_extrap.X(), CD_MUV2Cluster->GetPosition().Y() - MUV2_extrap.Y());
        // FillHisto("MUV2_nt_SW", MUV2Cluster_SW);
        // FillHisto("MUV2_TrP_SW", STRAW_P ,MUV2Cluster_SW);

        if(fabs(MUV2_extrap.X()) <= 130.  && fabs(MUV2_extrap.Y()) <= 130.){return;}
        if(fabs(MUV2_extrap.X()) >= 1100. || fabs(MUV2_extrap.Y()) >= 1100.){return;}
        if(fabs(MUV2T0) > MUV2OffsetCut){return;}
        if(MUV2dtrkcl_min > 150.) {return;}
        if( fabs(CD_MUV2Cluster->GetPosition().X() - MUV2_extrap.X()) > 260. ||
            fabs(CD_MUV2Cluster->GetPosition().Y() - MUV2_extrap.Y()) > 260.){return;}
    }
    //Before
    //cout << "Before --" << endl;
    //cout << "MUV1 Cand == " <<  MUV1Event->GetNCandidates() <<  "MUV3 Cand == " <<  MUV3Event->GetNCandidates() << "MUV2 Cand == " <<  MUV2Event->GetNCandidates() << endl;
    //cout << "MUV1 == " << MUV1TrackClusterIndex << " MUV2 == " << MUV2TrackClusterIndex << " MUV3 == " << MUV3TrackClusterIndex << endl;
    //cout << "ENDOF" << endl;
    //if(MUV1Event->GetNHits() < 1 && ){return;}
    //if(MUV2Event->GetNHits() < 1){return;}
    if(LKrEvent->GetNHits() < 1){return;}
    if(MUV3TrackClusterIndex > -1){

        TRecoMUV3Candidate* CD_MUV3Cluster = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex));
        double CD_MUV3ClusterTime = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex))->GetTime();
        double MUV3T0 = CD_CHODTime - CD_MUV3ClusterTime + MUV3Offset;
        if(fabs(MUV3T0) > MUV3OffsetCut){return;}
        if(fabs(MUV3_extrap.X()) <= 130. && fabs(MUV3_extrap.Y()) <= 130.){return;}
        if(fabs(MUV3_extrap.X()) >= 1100. || fabs(MUV3_extrap.Y()) >= 1100.){return;}
        if( fabs(CD_MUV3Cluster->GetX() - MUV3_extrap.X()) > 200. ||
            fabs(CD_MUV3Cluster->GetY() - MUV3_extrap.Y()) > 200.){return;}
        //if(MUV2dtrkcl_min > 150.) {return;}
        //cout << "Channel one  == " << CD_MUV3Cluster->GetChannel1() << "Channel two == " << CD_MUV3Cluster->GetChannel2() << "Tile ID == " << CD_MUV3Cluster->GetTileID() << endl;
        //        cout << "ROChannel one  == " << CD_MUV3Cluster->GetROChannel1() << "ROChannel  two == " << CD_MUV3Cluster->GetROChannel2() << "Tile ID == " << CD_MUV3Cluster->GetTileID()  << endl;
        //FillHisto("MUV3_timediff", MUV3T0);
        //FillHisto("MUV3_nearest_track_dtrkcl", MUV3dtrkcl_min);
        //FillHisto("MUV3_nearest_track_x_vs_y", MUV3_extrap.X(), MUV3_extrap.Y() );
        //FillHisto("MUV3_cda_x_vs_y", CD_MUV3Cluster->GetX() - MUV3_extrap.X(), CD_MUV3Cluster->GetY() - MUV3_extrap.Y() );


    }

    if(LKrTrackClusterIndex > -1){
        TRecoLKrCandidate* CD_LKrCluster = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex));
        double CD_LKrClusterTime  = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterTime();
        double CD_LKrClusterDDead = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterDDeadCell();
        double LKrT0 = CD_CHODTime  - CD_LKrClusterTime + LKrOffset;
        double LKrR  = sqrt( pow(LKr_extrap.X(),2) + pow(LKr_extrap.Y(),2) );
        double Cluster_X = CD_LKrCluster->GetClusterX()*10;
        double Cluster_Y = CD_LKrCluster->GetClusterY()*10;


        //CUTComment:: LKr cluster Quality cuts
        //1. and 2. Detector acceptance ( 15cm < R < 110cm)
        //3. Distance between track and closest cluster
        //4. Time matching of the cluster
        //5. Distance to deadcell > 2cm
        if(LKrR < 150. ) {return;}
        if(LKrR > 1100. ) {return;}
        if(LKrdtrkcl_min  > 50. ) {return;}
        if(fabs(LKrT0) > LKrOffsetCut){return;}
        if(CD_LKrClusterDDead < 2.){return;}


        FillHisto("LKr_nearest_track_DDeadCell", CD_LKrClusterDDead );
        FillHisto("LKr_timediff" , LKrT0);
        FillHisto("LKr_nearest_track_dtrkcl", LKrdtrkcl_min );
        FillHisto("LKr_nearest_track_x_vs_y", Cluster_X, Cluster_Y );
        FillHisto("LKr_extrap_x_vs_y", LKr_extrap.X(), LKr_extrap.Y() );

    }
    //if(MUV3Event->GetBurstID() == 232 || MUV3Event->GetBurstID() == 389 || MUV3Event->GetBurstID() == 432 || MUV3Event->GetBurstID() == 855 ||
    //   MUV3Event->GetBurstID() == 772 || MUV3Event->GetBurstID() == 885 || MUV3Event->GetBurstID() == 1069|| MUV3Event->GetBurstID() == 1111||
    //   MUV3Event->GetBurstID() == 965){return;}
    //if(MUV3Event->GetBurstID() == 453 || MUV3Event->GetBurstID() == 901 || MUV3Event->GetBurstID() == 1038 || MUV3Event->GetBurstID() == 792){ return;}


    for(int iLKrCand=0; iLKrCand<LKrEvent->GetNCandidates(); iLKrCand++){
        LKrCluster = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(iLKrCand));
        TRecoLKrCandidate* LKrNtCluster = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex));

        LkrPos.SetX ( LKrCluster->GetClusterX() );
        LkrPos.SetY ( LKrCluster->GetClusterY() );
        TVector3 LkrNtPos;
        LkrNtPos.SetX ( LKrNtCluster->GetClusterX() );
        LkrNtPos.SetY ( LKrNtCluster->GetClusterY() );
        double LKrEcluster   = 1000*LKrCluster->GetClusterEnergy(); //  [MeV]
        double LKrEseed      = 1000*LKrCluster->GetClusterSeedEnergy(); //  [MeV]
        double LKrE77        = 1000*LKrCluster->GetCluster77Energy(); //  [MeV]
        int    LKrNcells     = LKrCluster->GetNCells();
        double ClusterTime   = LKrCluster->GetClusterTime();
        double ClusterNtTime = LKrNtCluster->GetClusterTime();

        FillHisto("LKr_x_vs_y", LkrNtPos.X()*10.,LkrNtPos.Y()*10.);
        FillHisto("LKr_cda_x_vs_y", LkrNtPos.X()*10. - LKr_extrap.X() , LkrNtPos.Y()*10. - LKr_extrap.Y());
        FillHisto("LKr_Ecl_vs_NCell", LKrEcluster, LKrNcells);
        FillHisto("LKr_EoP", LKrEcluster/STRAW_P );
        FillHisto("LKr_Ecl", LKrEcluster );
        FillHisto("LKr_Eseed_over_Ecl", LKrEseed/LKrEcluster );
        FillHisto("LKr_Es_Ecl_vs_E77_Ecl", LKrEseed/LKrEcluster , 1 - LKrE77/LKrEcluster );

        if(iLKrCand == LKrTrackClusterIndex){continue;}
        FillHisto("LKr_nt_timediff", ClusterNtTime - ClusterTime);
        FillHisto("LKr_nt_dtrk", sqrt(pow(LkrPos.X()*10. - LkrNtPos.X()*10, 2 ) + pow(LkrPos.Y()*10. - LkrNtPos.Y()*10., 2 ) ));
        if(fabs(ClusterNtTime - ClusterTime) < 5 && sqrt(pow(LkrPos.X()*10. - LkrNtPos.X()*10, 2 ) + pow(LkrPos.Y()*10. - LkrNtPos.Y()*10., 2 ) ) < 200){return;}

    }


    FillHisto("CHOD_cda_x_vs_y", CD_CHODPos.X()*10. - CHOD_extrap.X() , CD_CHODPos.Y()*10. - CHOD_extrap.Y());
    FillHisto("CHOD_nearest_track_dtrkcl", CHODdtrkcl_min);
    FillHisto("CHOD_nearest_track_x_vs_y", CD_CHODPos.X()*10.,CD_CHODPos.Y()*10. );
    FillHisto("CHOD_extrap_x_vs_y", CHOD_extrap.X(), CHOD_extrap.Y() );

    if(MUV1TrackClusterIndex > -1){
        TRecoMUV1Candidate* CD_MUV1Cluster = ((TRecoMUV1Candidate*)MUV1Event->GetCandidate(MUV1TrackClusterIndex));
        double CD_MUV1ClusterTime = ((TRecoMUV1Candidate*)MUV1Event->GetCandidate(MUV1TrackClusterIndex))->GetTime();
        double MUV1T0 = CD_CHODTime - CD_MUV1ClusterTime + MUV1Offset;
        TClonesArray *MUV1Hits = MUV1Event->GetHits();
        double MUV1Cluster_Charge=0;
        double MUV1Cluster_SW=CD_MUV1Cluster->GetShowerWidth();
        double MUV1Quality = CD_MUV1Cluster->GetQuality();
        TVector2 CD_MUV1Pos = CD_MUV1Cluster->GetPosition();
        for(int iMUV1Hit = 0; iMUV1Hit <CD_MUV1Cluster->GetNHits(); iMUV1Hit++ ){
            TRecoMUV1Hit* MUV1Hit = ((TRecoMUV1Hit*)MUV1Hits->At(iMUV1Hit));
            double MUV1HitCharge = MUV1Hit->GetCharge();
            MUV1Cluster_Charge  += MUV1HitCharge;
        }



        FillHisto("Quality", MUV1Quality);
        if(MUV1Quality==0){
            FillHisto("Q0_nearest_track_x_vs_y", CD_MUV1Pos.X(), CD_MUV1Pos.Y()  );
        }
        if(MUV1Quality==1){
            FillHisto("Q1_nearest_track_x_vs_y", CD_MUV1Pos.X(), CD_MUV1Pos.Y()  );
        }
        if(MUV1Quality==2){
            FillHisto("Q2_nearest_track_x_vs_y", CD_MUV1Pos.X(), CD_MUV1Pos.Y()  );
        }
        FillHisto("MUV1_Nhits", CD_MUV1Cluster->GetNHits() );
        FillHisto("MUV1_timediff", MUV1T0);
        FillHisto("MUV1_nearest_track_dtrkcl", MUV1dtrkcl_min);
        FillHisto("MUV1_nearest_track_cluster_charge", MUV1Cluster_Charge);
        FillHisto("MUV1_near_charge_vs_dtrkcl", MUV1Cluster_Charge, MUV1dtrkcl_min);
        FillHisto("MUV1_PvsQ", STRAW_P, MUV1Cluster_Charge);

        FillHisto("MUV1_nearest_track_x_vs_y", CD_MUV1Cluster->GetPosition().X(), CD_MUV1Cluster->GetPosition().Y() );
        FillHisto("MUV1_extrap_x_vs_y", MUV1_extrap.X(), MUV1_extrap.Y() );
        FillHisto("MUV1_nearest_track_VvsH", CD_MUV1Cluster->GetVerticalChannel(), CD_MUV1Cluster->GetHorizontalChannel() );
        //Old Reco
        //FillHisto("MUV1_cda_x_vs_y", CD_MUV1Cluster->GetPosition().X() - MUV1_extrap.X() + 60., CD_MUV1Cluster->GetPosition().Y() - MUV1_extrap.Y() + 60. );
        //New Reco
        FillHisto("MUV1_cda_x_vs_y", CD_MUV1Cluster->GetPosition().X() - MUV1_extrap.X(), CD_MUV1Cluster->GetPosition().Y() - MUV1_extrap.Y() );
        FillHisto("MUV1_nt_SW", MUV1Cluster_SW);
        FillHisto("MUV1_TrP_SW", STRAW_P ,MUV1Cluster_SW);

    }
    if(MUV2TrackClusterIndex > -1){

        TRecoMUV2Candidate* CD_MUV2Cluster = ((TRecoMUV2Candidate*)MUV2Event->GetCandidate(MUV2TrackClusterIndex));
        double CD_MUV2ClusterTime   = ((TRecoMUV2Candidate*)MUV2Event->GetCandidate(MUV2TrackClusterIndex))->GetTime();
        double MUV2T0 = CD_CHODTime - CD_MUV2ClusterTime + MUV2Offset;
        double MUV2Cluster_Charge=0;
        double MUV2Cluster_SW=CD_MUV2Cluster->GetShowerWidth();
        TClonesArray *MUV2Hits = MUV2Event->GetHits();
        //double MUV2R  = sqrt( pow(MUV2_extrap.X(),2) + pow(MUV2_extrap.Y(),2) );
        for(int iMUV2Hit = 0; iMUV2Hit <CD_MUV2Cluster->GetNHits(); iMUV2Hit++ ){
            TRecoMUV2Hit* MUV2Hit = ((TRecoMUV2Hit*)MUV2Hits->At(iMUV2Hit));
            double MUV2HitCharge = MUV2Hit->GetCharge();
            MUV2Cluster_Charge  += MUV2HitCharge;
        }


        FillHisto("MUV2_Nhits", CD_MUV2Cluster->GetNHits() );
        FillHisto("MUV2_timediff", MUV2T0);
        FillHisto("MUV2_nearest_track_dtrkcl", MUV2dtrkcl_min);
        FillHisto("MUV2_nearest_track_cluster_charge", MUV2Cluster_Charge);
        FillHisto("MUV2_near_charge_vs_dtrkcl", MUV2Cluster_Charge, MUV2dtrkcl_min);
        FillHisto("MUV2_PvsQ", STRAW_P, MUV2Cluster_Charge);
        FillHisto("MUV2_nearest_track_x_vs_y", CD_MUV2Cluster->GetPosition().X(), CD_MUV2Cluster->GetPosition().Y() );
        FillHisto("MUV2_nearest_track_VvsH", CD_MUV2Cluster->GetVerticalChannel(), CD_MUV2Cluster->GetHorizontalChannel() );
        FillHisto("MUV2_extrap_x_vs_y", MUV2_extrap.X(), MUV2_extrap.Y());
        FillHisto("MUV2_cda_x_vs_y", CD_MUV2Cluster->GetPosition().X() - MUV2_extrap.X(), CD_MUV2Cluster->GetPosition().Y() - MUV2_extrap.Y());
        FillHisto("MUV2_nt_SW", MUV2Cluster_SW);
        FillHisto("MUV2_TrP_SW", STRAW_P ,MUV2Cluster_SW);

    }

    if(MUV3TrackClusterIndex > -1){

        TRecoMUV3Candidate* CD_MUV3Cluster = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex));
        double CD_MUV3ClusterTime = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex))->GetTime();
        double MUV3T0 = CD_CHODTime - CD_MUV3ClusterTime + MUV3Offset;

        FillHisto("MUV3_timediff", MUV3T0);
        FillHisto("MUV3_nearest_track_dtrkcl", MUV3dtrkcl_min);
        FillHisto("MUV3_extrap_x_vs_y", MUV3_extrap.X(), MUV3_extrap.Y() );
        FillHisto("MUV3_nearest_track_x_vs_y", CD_MUV3Cluster->GetX(), CD_MUV3Cluster->GetY() );
        FillHisto("MUV3_cda_x_vs_y", CD_MUV3Cluster->GetX() - MUV3_extrap.X(), CD_MUV3Cluster->GetY() - MUV3_extrap.Y() );


    }

    double MUVs123=0;
    //Testung MUV candidates
    if(MUV3TrackClusterIndex > -1 &&
       MUV2TrackClusterIndex > -1 &&
       MUV1TrackClusterIndex > -1 ){
        TRecoMUV3Candidate* MUV3Cl_123 = (TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex);
        MUVs123 = 1;
        Int_t  MUV1_Hindex_MUV123 = MUV1Geometry::GetInstance()->GetScintillatorAt(MUV1_extrap.Y());
        Int_t  MUV1_Vindex_MUV123 = MUV1Geometry::GetInstance()->GetScintillatorAt(MUV1_extrap.X());
        Int_t  MUV2_Hindex_MUV123 = MUV2Geometry::GetInstance()->GetScintillatorAt(MUV2_extrap.Y());
        Int_t  MUV2_Vindex_MUV123 = MUV2Geometry::GetInstance()->GetScintillatorAt(MUV2_extrap.X());
        FillHisto("BadMUV1_HitMap_MUV123", MUV1_Vindex_MUV123, MUV1_Hindex_MUV123);
        FillHisto("BadMUV2_HitMap_MUV123", MUV2_Vindex_MUV123, MUV2_Hindex_MUV123);
        FillHisto("MUV3Hit_GoodEvent", MUV3Cl_123->GetTileID());
        FillHisto("BurstID_vs_MUV123", MUV3Event->GetBurstID(),MUVs123);
        FillHisto("MUV3_nearest_track_dtrkcl_MUV123", MUV3dtrkcl_min);
        FillHisto("MUV123", MUVs123);

        FillHisto("Nhits123_MUV1", MUV1Event->GetNHits());
        FillHisto("Nhits123_MUV2", MUV2Event->GetNHits());


        if(LKrTrackClusterIndex > -1){
            double CD_MUV3ClusterTime = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex))->GetTime();
            double CD_LKrClusterTime  = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterTime();
            FillHisto("MUV3_LKr_tdiff_MUV123", CD_MUV3ClusterTime - CD_LKrClusterTime + LKrOffset);
            FillHisto("Nhits123_LKr", LKrEvent->GetNHits());
        }
    }
    //return;
    double MUVs23=0;
    if(MUV3TrackClusterIndex > -1 &&
       MUV2TrackClusterIndex > -1 &&
       MUV1TrackClusterIndex < 0 ){
        MUVs23 = 1;
        TRecoMUV3Candidate* MUV3Cl_23 = (TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex);

        Int_t  MUV1_Hindex_MUV23 = MUV1Geometry::GetInstance()->GetScintillatorAt(MUV1_extrap.Y());
        Int_t  MUV1_Vindex_MUV23 = MUV1Geometry::GetInstance()->GetScintillatorAt(MUV1_extrap.X());
        Int_t  MUV2_Hindex_MUV23 = MUV2Geometry::GetInstance()->GetScintillatorAt(MUV2_extrap.Y());
        Int_t  MUV2_Vindex_MUV23 = MUV2Geometry::GetInstance()->GetScintillatorAt(MUV2_extrap.X());
        FillHisto("BadMUV1_HitMap_MUV23", MUV1_Vindex_MUV23, MUV1_Hindex_MUV23);
        FillHisto("BadMUV2_HitMap_MUV23", MUV2_Vindex_MUV23, MUV2_Hindex_MUV23);
        FillHisto("MUV1_Ncandidates_MUV23", MUV1Event->GetNCandidates());
        FillHisto("MUV2_Ncandidates_MUV23", MUV2Event->GetNCandidates());
        FillHisto("MUV3_Ncandidates_MUV23", MUV3Event->GetNCandidates());
        FillHisto("MUV3_nearest_track_dtrkcl_MUV23", MUV3dtrkcl_min);
        FillHisto("MUV3Hit_NoMUV1", MUV3Cl_23->GetTileID());

        FillHisto("BurstID_vs_MUV23", MUV3Event->GetBurstID(),MUVs23);
        FillHisto("TrackP_MUV23", STRAW_P);
        FillHisto("TrackPfit_TrackP_MUV23", STRAW_P - STRAW_Pbf);
        FillHisto("MUV23", MUVs23);
        FillHisto("Nhits0C23_MUV1", MUV1Event->GetNHits());
        FillHisto("Nhits0C23_MUV2", MUV2Event->GetNHits());
        //Checking MUV Nhits for the inefficient bursts
        if(MUV3Event->GetBurstID() == 232 || MUV3Event->GetBurstID() == 389 || MUV3Event->GetBurstID() == 432 || MUV3Event->GetBurstID() == 855 ||
           MUV3Event->GetBurstID() == 772 || MUV3Event->GetBurstID() == 885 || MUV3Event->GetBurstID() == 1069|| MUV3Event->GetBurstID() == 1111||
           MUV3Event->GetBurstID() == 965){
            FillHisto("Nhits0C23_MUV1_BB", MUV1Event->GetNHits());
        }
        TClonesArray *MUV1Hits = MUV1Event->GetHits();
        double MUV1_counter=0;
        double MUV1_Hcounter=0;
        double MUV1_Vcounter=0;

        for(int iMUV1Hit=0;  iMUV1Hit <  MUV1Event->GetNHits(); iMUV1Hit++){
            TRecoMUV1Hit* MUV1Hit = ((TRecoMUV1Hit*)MUV1Hits->At(iMUV1Hit));
            double Hit_CHOD_tdiff = CD_CHODTime -  MUV1Hit->GetTime() + MUV1Offset;

            if(MUV1Hit->GetChannelID()%100 < 50){
                //if(MUV1Hit->GetTime() -  < 50)
                FillHisto("0C_ChID_MUV1", MUV1Hit->GetChannelID());
                FillHisto("0C_VChID_diff_M1", (MUV1Hit->GetChannelID()%50) - MUV1_Vindex_MUV23);
                if(fabs((MUV1Hit->GetChannelID()%50) - MUV1_Vindex_MUV23) < 2){
                    FillHisto("0C_VM1_CHOD_t", Hit_CHOD_tdiff);
                    if(fabs(Hit_CHOD_tdiff) < 30){
                        MUV1_counter++;
                        MUV1_Vcounter++;
                        FillHisto("0C_VChID_MUV1", MUV1Hit->GetChannelID()%50);
                        FillHisto("MUV23_Vsaved_M1",(MUV1Hit->GetChannelID()%50 ));
                    }
                }

            }

            if(MUV1Hit->GetChannelID()%100 > 50){
                //if(MUV1Hit->GetChannelID()%50 - 1 < 50)
                FillHisto("0C_ChID_MUV1", MUV1Hit->GetChannelID());
                FillHisto("0C_HChID_diff_M1", (MUV1Hit->GetChannelID()%50) - MUV1_Hindex_MUV23);

                if(fabs(MUV1Hit->GetChannelID()%50 - MUV1_Hindex_MUV23) < 2){
                    FillHisto("0C_HM1_CHOD_t", Hit_CHOD_tdiff );

                    if( fabs(Hit_CHOD_tdiff) < 35 && fabs(Hit_CHOD_tdiff) > 15 ){

                        //TMUV1Digi* MUV1Digi = (TMUV1Digi*) MUV1Hit->GetDigi();
                        //Double_t *DigiSamples = MUV1Digi->GetAllSamples();
                        //for(int i = 0; i < MUV1Digi->GetNSamples(); i++){
                        //    cout << DigiSamples[i] << endl;
                        //}
                        FillHisto("ChannelID_25ns_away_M1",MUV1Hit->GetChannelID());

                    }

                    if(fabs(Hit_CHOD_tdiff) < 30){
                        MUV1_counter++;
                        MUV1_Hcounter++;
                        FillHisto("0C_HChID_MUV1", MUV1Hit->GetChannelID()%50);
                        FillHisto("MUV23_Hsaved_M1",(MUV1Hit->GetChannelID()%50));
                    }
                }
            }

            if(iMUV1Hit == (MUV1Event->GetNHits() - 1) ){

                FillHisto("RecVHits_M1", MUV1_Vcounter);
                FillHisto("RecHHits_M1", MUV1_Hcounter);
                //cout << "MUV1 ---------------------- " << endl;
                //cout << " VSaved == " << MUV1_Vcounter << " HSaved == " << MUV1_Hcounter << "CHOD_Tdiff == " << fabs(Hit_CHOD_tdiff) <<  "Test == " << CD_CHODTime -  MUV1Hit->GetTime() + MUV1Offset << endl;
                if(MUV1_Vcounter != 0  && MUV1_Hcounter != 0 // &&  fabs(Hit_CHOD_tdiff) < 30
                   ){
                    FillHisto("RecHits_M1", MUV1_counter);
                }

            }
        }

        if(LKrTrackClusterIndex > -1){
            double CD_MUV3ClusterTime = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex))->GetTime();
            double CD_LKrClusterTime  = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterTime();
            FillHisto("MUV3_LKr_tdiff_MUV23", CD_MUV3ClusterTime - CD_LKrClusterTime + LKrOffset);
            FillHisto("Nhits0C23_LKr", LKrEvent->GetNHits());
        }
    }

    double MUVs13=0;
    if(MUV3TrackClusterIndex > -1 &&
       MUV1TrackClusterIndex > -1 &&
       MUV2TrackClusterIndex < 0 ){
        MUVs13 = 1;
        TRecoMUV3Candidate* MUV3Cl_13 = (TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex);
        Int_t  MUV1_Hindex_MUV13 = MUV1Geometry::GetInstance()->GetScintillatorAt(MUV1_extrap.Y());
        Int_t  MUV1_Vindex_MUV13 = MUV1Geometry::GetInstance()->GetScintillatorAt(MUV1_extrap.X());
        Int_t  MUV2_Hindex_MUV13 = MUV2Geometry::GetInstance()->GetScintillatorAt(MUV2_extrap.Y());
        Int_t  MUV2_Vindex_MUV13 = MUV2Geometry::GetInstance()->GetScintillatorAt(MUV2_extrap.X());
        FillHisto("BadMUV1_HitMap_MUV13", MUV1_Vindex_MUV13, MUV1_Hindex_MUV13);
        FillHisto("BadMUV2_HitMap_MUV13", MUV2_Vindex_MUV13, MUV2_Hindex_MUV13);
        FillHisto("TrackP_MUV13", STRAW_P);
        FillHisto("TrackPfit_TrackP_MUV13", STRAW_P - STRAW_Pbf);
        FillHisto("MUV1_Ncandidates_MUV13", MUV1Event->GetNCandidates());
        FillHisto("MUV2_Ncandidates_MUV13", MUV2Event->GetNCandidates());
        FillHisto("MUV3_Ncandidates_MUV13", MUV3Event->GetNCandidates());
        FillHisto("MUV3_nearest_track_dtrkcl_MUV13", MUV3dtrkcl_min);
        FillHisto("MUV3Hit_NoMUV2", MUV3Cl_13->GetTileID());
        FillHisto("MUV13", MUVs13);
        FillHisto("BurstID_vs_MUV13", MUV3Event->GetBurstID(),MUVs13);
        FillHisto("Nhits0C13_MUV1", MUV1Event->GetNHits());
        FillHisto("Nhits0C13_MUV2", MUV2Event->GetNHits());
        //Checking MUV Nhits for the inefficient bursts
        if(MUV3Event->GetBurstID() == 453 || MUV3Event->GetBurstID() == 901 || MUV3Event->GetBurstID() == 1038 || MUV3Event->GetBurstID() == 792){
            FillHisto("Nhits0C13_MUV2_BB", MUV2Event->GetNHits());
        }

        TClonesArray *MUV2Hits = MUV2Event->GetHits();
        double MUV2_counter=0;
        double MUV2_Hcounter=0;
        double MUV2_Vcounter=0;
        for(int iMUV2Hit=0; iMUV2Hit < MUV2Event->GetNHits(); iMUV2Hit++){
            TRecoMUV2Hit* MUV2Hit = ((TRecoMUV2Hit*)MUV2Hits->At(iMUV2Hit));
            double Hit_M2CHOD_tdiff = CD_CHODTime -  MUV2Hit->GetTime() + MUV2Offset;
            FillHisto("0C_ChID_MUV2", MUV2Hit->GetChannelID());
            //FillHisto("0C_ChID_diff_M2", MUV2Hit->GetChannelID() - MUV1_);

            if(MUV2Hit->GetChannelID()%100 < 50){
                //if(MUV2Hit->GetChannelID()%50 - 1 < 50)
                FillHisto("0C_ChID_MUV2", MUV2Hit->GetChannelID());
                FillHisto("0C_VChID_diff_M2", (MUV2Hit->GetChannelID()%50) - MUV2_Vindex_MUV13);
                if( fabs((MUV2Hit->GetChannelID()%50) - MUV2_Vindex_MUV13) < 2){
                    FillHisto("0C_VM2_CHOD_t", Hit_M2CHOD_tdiff);

                    if( fabs(Hit_M2CHOD_tdiff) < 35 && fabs(Hit_M2CHOD_tdiff) > 15 ){


                           //Remove when over :D
        //TMUV2Event *Event = (TMUV2Event*)GetEvent("Digis","MUV2");
        FADCEvent *Event = (FADCEvent*)GetEvent("MUV2","Digis");
        //TClonesArray *MUV2Hits = MUV2Event->GetHits();
        //TMUV2Event* DetectorMUV2Event = (TMUV2Event*);
        // TClonesArray &Digis = (*(Event->GetHits()));
        TClonesArray* Digis = Event->GetHits();

        Int_t NDigis = Event->GetNHits();
        for(int iMUV2Hit=0;  iMUV2Hit <  NDigis; iMUV2Hit++){

            TMUV2Digi* MUV2Digi = (TMUV2Digi*) Digis->At(iMUV2Hit);
            Double_t* DigiSamples = MUV2Digi->GetAllSamples();
            Double_t Samples[12] = {1,2,3,4,5,6,7,8,9,10,11,12};

            Int_t DigiNSamples = MUV2Digi->GetNSamples();
            Int_t MUV2Quality = MUV2Digi->GetQuality();

            BookHisto(Form("Digi_%d",MUV2Event->GetID()),new TGraph(Form("Digi_%d",MUV2Event->GetID()),Form("Digi for event with ID = %d",MUV2Event->GetID())));

            for(int i = 0; i < MUV2Digi->GetNSamples(); i++){
                 // FillHisto(Form("Digi_%d",MUV2Event->GetID()), i, DigiSamples[i]);
                cout << "NS == " << DigiNSamples << "Sampe["<<i<<"]"<< DigiSamples[i]<< endl;
            }
        }


                        FillHisto("ChannelID_25ns_away_M2",MUV2Hit->GetChannelID());

                    }

                    if(fabs(Hit_M2CHOD_tdiff) < 30) {
                        MUV2_counter++;
                        MUV2_Vcounter++;
                        FillHisto("0C_VChID_MUV2", MUV2Hit->GetChannelID()%50);
                        FillHisto("MUV13_Vsaved_M2",(MUV2Hit->GetChannelID()%50));
                    }
                }

            }
            if(MUV2Hit->GetChannelID()%100 > 50){
                //if(MUV2Hit->GetChannelID()%50 - 1 < 50)
                FillHisto("0C_ChID_MUV2", MUV2Hit->GetChannelID());
                FillHisto("0C_HChID_diff_M2", (MUV2Hit->GetChannelID()%50) - MUV2_Hindex_MUV13);

                if( fabs((MUV2Hit->GetChannelID()%50) - MUV2_Hindex_MUV13) < 2){
                    FillHisto("0C_HM2_CHOD_t", Hit_M2CHOD_tdiff);

                    if(fabs(Hit_M2CHOD_tdiff) < 30) {
                        MUV2_counter++;
                        MUV2_Hcounter++;
                        FillHisto("0C_HChID_MUV2", MUV2Hit->GetChannelID()%50);
                        FillHisto("MUV13_Hsaved_M2",(MUV2Hit->GetChannelID()%50));
                    }

                }

            }

            if(iMUV2Hit == (MUV2Event->GetNHits() - 1) ){

                FillHisto("RecVHits_M2", MUV2_Vcounter);
                FillHisto("RecHHits_M2", MUV2_Hcounter);
                //cout << "MUV2 ---------------------- " << endl;
                //cout << " VSaved == " << MUV2_Vcounter << " HSaved == " << MUV2_Hcounter << "CHOD_Tdiff == " << fabs(Hit_M2CHOD_tdiff) << "Test == " << CD_CHODTime -  MUV2Hit->GetTime() + MUV2Offset << endl;
                if(MUV2_Vcounter != 0  && MUV2_Hcounter != 0 // &&  fabs(Hit_M2CHOD_tdiff) < 30
                   ){
                    FillHisto("RecHits_M2", MUV2_counter);
                }

            }

        }

        if(LKrTrackClusterIndex > -1){
            double CD_MUV3ClusterTime = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex))->GetTime();
            double CD_LKrClusterTime  = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterTime();
            FillHisto("Nhits0C13_LKr", LKrEvent->GetNHits());
            FillHisto("MUV3_LKr_tdiff_MUV13", CD_MUV3ClusterTime - CD_LKrClusterTime + LKrOffset);
        }

    }

    double MUVs3=0;
    if(MUV3TrackClusterIndex > -1 &&
       MUV2TrackClusterIndex < 0  &&
       MUV1TrackClusterIndex < 0 ){
        MUVs3 = 1;
        Int_t  MUV1_Hindex_MUV3only = MUV1Geometry::GetInstance()->GetScintillatorAt(MUV1_extrap.Y());
        Int_t  MUV1_Vindex_MUV3only = MUV1Geometry::GetInstance()->GetScintillatorAt(MUV1_extrap.X());
        Int_t  MUV2_Hindex_MUV3only = MUV2Geometry::GetInstance()->GetScintillatorAt(MUV2_extrap.Y());
        Int_t  MUV2_Vindex_MUV3only = MUV2Geometry::GetInstance()->GetScintillatorAt(MUV2_extrap.X());
        FillHisto("BadMUV1_HitMap_MUV3", MUV1_Vindex_MUV3only, MUV1_Hindex_MUV3only);
        FillHisto("BadMUV2_HitMap_MUV3", MUV2_Vindex_MUV3only, MUV2_Hindex_MUV3only);
        FillHisto("TrackP_MUV3", STRAW_P);
        FillHisto("TrackPfit_TrackP_MUV3", STRAW_P - STRAW_Pbf);
        FillHisto("MUV1_Ncandidates_MUV3", MUV1Event->GetNCandidates());
        FillHisto("MUV2_Ncandidates_MUV3", MUV2Event->GetNCandidates());
        FillHisto("MUV3_Ncandidates_MUV3", MUV3Event->GetNCandidates());
        FillHisto("MUV3_nearest_track_dtrkcl_MUV3", MUV3dtrkcl_min);
        FillHisto("MUV3Only", MUVs3);
        FillHisto("BurstID_vs_MUV3", MUV3Event->GetBurstID(),MUVs3);
        if(LKrTrackClusterIndex > -1){
            double CD_MUV3ClusterTime = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex))->GetTime();
            double CD_LKrClusterTime  = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterTime();
            FillHisto("MUV3_LKr_tdiff_MUV3", CD_MUV3ClusterTime - CD_LKrClusterTime + LKrOffset);
        }

        //Checking MUV Nhits for the inefficient bursts
        FillHisto("Nhits0C3_MUV2_BB", MUV2Event->GetNHits());
        FillHisto("Nhits0C3_MUV1_BB", MUV1Event->GetNHits());

        //for(int iMUV3Cand=0; iMUV3Cand < MUV3Event->GetNCandidates(); iMUV3Cand++){
        //    MUV3Cluster  = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(iMUV3Cand));
        //    double MUV3ClusterTime_MUV3  = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(iMUV3Cand))->GetTime();
        //
        //    if(LKrTrackClusterIndex > -1){
        //        //double CD_LKrClusterTime  = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterTime();
        //        FillHisto("MUV3_LKr_tdiff_MUV3", MUV3ClusterTime_MUV3);
        //    }
        //}
    }
    FillHisto("STRAW1_x_vs_y", PositionBefore.X() , PositionBefore.Y());
    FillHisto("STRAW4_x_vs_y", PositionAfter.X() , PositionAfter.Y());

    ///








    for(int iRICHCand=0; iRICHCand<RICHEvent->GetNRingCandidates(); iRICHCand++){ //loop su Ring Cand

        RingCandidate = RICHEvent->GetRingCandidate(iRICHCand);
        RingCandidate->SetEvent(RICHEvent);

        TVector2 RingCenter= RingCandidate->GetRingCenter();
        double RingCenterR = RingCandidate->GetRingRadius();
        double RingTime    = RingCandidate->GetRingTime();
        double NeonN       = 1.000067;
        double RICHAngle   = TMath::ATan( RingCenterR/17000. );
        double RICHMass    = STRAW_P*0.001*sqrt( pow(NeonN*TMath::Cos(RICHAngle), 2 ) - 1.);
        double slopex      = RingCenter.X()/17000.;
        double slopey      = RingCenter.Y()/17000.;
        //cout << "Beta == " << 1./(NeonN*TMath::Cos(RICHAngle)) << " 1/beta == " << NeonN*TMath::Cos(RICHAngle) << "Track_P == " << STRAW_P << endl;

        FillHisto("RICHRadius",RingCandidate->GetRingRadius());
        FillHisto("RICHAngle",RingCenterR/17000.);
        FillHisto("RICHMass",RICHMass);
        FillHisto("RICHMvsP", STRAW_P*0.001 , RICHMass);
        //FillHisto("RICHRvsp",RingCandidate->GetRingRadius(),STRAWCandidate->GetMomentum());

        FillHisto("RICH_STRAW_dxdzdiff", STRAW_dxdz - RingCenter.X()/17000. );
        FillHisto("RICH_STRAW_dydzdiff", STRAW_dydz - RingCenter.Y()/17000. );
        FillHisto("RICHRvsP", STRAW_P , RingCandidate->GetRingRadius());
        FillHisto("RICH_x_vs_y", RICH_extrap.X() , RICH_extrap.Y());
        FillHisto("RICH_cda_x_vs_y", RICH_extrap.X() - RingCenter.X() , RICH_extrap.Y() - RingCenter.Y());
        FillHisto("RICH_timediff", RingTime - CD_CHODTime + RICHOffsetCut);
        //if (RICH_extrap.X()>0) {
        //slopex+=0.00035;
        //slopey+=0.00057;
        FillHisto("RICH_dxdz_vs_dydz", slopex , slopey );
        // }
        //if (RICH_extrap.X()<0) {
        //slopex+=0.00049;
        //slopey+=0.00068;
        // FillHisto("RICH_dxdz_vs_dydz", slopex , slopey );
        //}
    }
    for(int iMUV1Cand=0; iMUV1Cand < MUV1Event->GetNCandidates(); iMUV1Cand++){

        MUV1Cluster =((TRecoMUV1Candidate*)MUV1Event->GetCandidate(iMUV1Cand));
        TRecoMUV1Candidate* MUV1NtCluster =((TRecoMUV1Candidate*)MUV1Event->GetCandidate(MUV1TrackClusterIndex));

        TVector2 MUV1PosOld       = MUV1Cluster->GetPosition();
        TVector2 MUV1NtPosOld     = MUV1NtCluster->GetPosition();
        int MUV1VerticalChannel   = MUV1Cluster->GetVerticalChannel();
        int MUV1NtVerticalChannel   = MUV1NtCluster->GetVerticalChannel();
        int MUV1HorizontalChannel = MUV1Cluster->GetHorizontalChannel();
        int MUV1NtHorizontalChannel = MUV1NtCluster->GetHorizontalChannel();
        double MUV1Time           = MUV1Cluster->GetTime();
        double MUV1ntTime         = MUV1NtCluster->GetTime();
        double MUV1SeedEnergy     = MUV1Cluster->GetSeedEnergy();
        double MUV1ClusterEnergy  = MUV1Cluster->GetEnergy();
        double MUV1SeedEnergyX    = MUV1Cluster->GetSeedEnergyHorizontal();
        double MUV1SeedEnergyY    = MUV1Cluster->GetSeedEnergyVertical();
        double MUV1_dtrk = sqrt(pow(MUV1Pos.X() - MUV1_extrap.X(),2 ) + pow(MUV1Pos.Y() - MUV1_extrap.Y(),2 ) ) ;

        // Old Reco!!!!!!!!!!!!!!!!!!!!!!
        //MUV1Pos   = MUV1PosOld + 60.;
        //New Reco
        MUV1Pos   = MUV1PosOld;
        TVector2 MUV1NtPos = MUV1NtPosOld + 60.;

        //cout << "Old positionX = " << MUV1PosOld.X() << "And new == " << MUV1Pos.X() << endl;
        FillHisto("MUV1_ClusterEnergy",MUV1ClusterEnergy);
        FillHisto("MUV1_Eseed_over_Ecl",MUV1SeedEnergy/MUV1ClusterEnergy);
        FillHisto("MUV1_SeedEnergy",MUV1SeedEnergy);
        FillHisto("MUV1_SeedEnergy_vs_ClusterEnergy",MUV1SeedEnergy, MUV1ClusterEnergy);
        FillHisto("MUV1_SeedEnergy_xvsy",MUV1SeedEnergyX,MUV1SeedEnergyY);
        FillHisto("MUV1xvsy", MUV1Pos.X(), MUV1Pos.Y());
        // FillHisto("MUV1_cda_x_vs_y", MUV1Pos.X() - MUV1_extrap.X(), MUV1Pos.Y() - MUV1_extrap.Y());
        FillHisto("MUV1_trk_dist",MUV1_dtrk);

        if(iMUV1Cand == MUV1TrackClusterIndex){continue;}

        FillHisto("MUV1_nt_timediff", MUV1ntTime - MUV1Time);
        FillHisto("MUV1_sc_HitMap", MUV1HorizontalChannel, MUV1VerticalChannel);
        FillHisto("MUV1_hz_nt_chdiff", MUV1NtHorizontalChannel - MUV1HorizontalChannel);
        FillHisto("MUV1_vt_nt_chdiff", MUV1NtVerticalChannel - MUV1VerticalChannel);
        FillHisto("MUV1_hz_vt_chdiff", MUV1NtHorizontalChannel - MUV1HorizontalChannel, MUV1NtVerticalChannel - MUV1VerticalChannel);
        if( fabs(MUV1NtHorizontalChannel - MUV1HorizontalChannel) < 4 || fabs(MUV1NtVerticalChannel - MUV1VerticalChannel) < 4 ){

            FillHisto("MUV1_zero_distance_timediff", MUV1ntTime - MUV1Time);
        }

    }

    for(int iMUV2Cand=0; iMUV2Cand < MUV2Event->GetNCandidates(); iMUV2Cand++){
        MUV2Cluster = ((TRecoMUV2Candidate*)MUV2Event->GetCandidate(iMUV2Cand));
        TRecoMUV2Candidate* MUV2NtCluster =((TRecoMUV2Candidate*)MUV2Event->GetCandidate(MUV2TrackClusterIndex));
        MUV2Pos           = MUV2Cluster->GetPosition();
        TVector2 MUV2ntPosOld     = MUV2NtCluster->GetPosition();
        int MUV2NtVerticalChannel   = MUV2NtCluster->GetVerticalChannel();
        int MUV2HorizontalChannel = MUV2Cluster->GetHorizontalChannel();
        int MUV2NtHorizontalChannel = MUV2NtCluster->GetHorizontalChannel();
        int MUV2VerticalChannel   = MUV2Cluster->GetVerticalChannel();
        double MUV2Time          = MUV2Cluster->GetTime();
        double MUV2ntTime        = MUV2NtCluster->GetTime();
        double MUV2SeedEnergy    = MUV2Cluster->GetSeedEnergy();
        double MUV2ClusterEnergy = MUV2Cluster->GetEnergy();
        double MUV2SeedEnergyX   = MUV2Cluster->GetSeedEnergyHorizontal();
        double MUV2SeedEnergyY   = MUV2Cluster->GetSeedEnergyVertical();
        double MUV2_dtrk = sqrt(pow(MUV2Pos.X() - MUV2_extrap.X(),2 ) + pow(MUV2Pos.Y() - MUV2_extrap.Y(),2 ) ) ;



        FillHisto("MUV2_ClusterEnergy",MUV2ClusterEnergy);
        FillHisto("MUV2_Eseed_over_Ecl",MUV2SeedEnergy/MUV2ClusterEnergy);
        FillHisto("MUV2_SeedEnergy",MUV2SeedEnergy);
        FillHisto("MUV2_SeedEnergy_xvsy",MUV2SeedEnergyX,MUV2SeedEnergyY);
        FillHisto("MUV2_SeedEnergy_vs_ClusterEnergy",MUV2SeedEnergy, MUV2ClusterEnergy);
        FillHisto("MUV2xvsy", MUV2Pos.X(), MUV2Pos.Y());
        // FillHisto("MUV2_cda_x_vs_y", MUV2Pos.X() - MUV2_extrap.X(), MUV2Pos.Y() - MUV2_extrap.Y());
        FillHisto("MUV2_trk_dist",MUV2_dtrk);
        if(iMUV2Cand == MUV2TrackClusterIndex){continue;}
        FillHisto("MUV2_nt_timediff", MUV2ntTime - MUV2Time);
        FillHisto("MUV2_sc_HitMap", MUV2HorizontalChannel, MUV2VerticalChannel);
        FillHisto("MUV2_hz_nt_chdiff", MUV2NtHorizontalChannel - MUV2HorizontalChannel);
        FillHisto("MUV2_vt_nt_chdiff", MUV2NtVerticalChannel - MUV2VerticalChannel);
        FillHisto("MUV2_hz_vt_chdiff", MUV2NtHorizontalChannel - MUV2HorizontalChannel, MUV2NtVerticalChannel - MUV2VerticalChannel);

        if( fabs(MUV2NtHorizontalChannel - MUV2HorizontalChannel) < 4 || fabs(MUV2NtVerticalChannel - MUV2VerticalChannel) < 4 ){
            FillHisto("MUV2_zero_distance_timediff", MUV2ntTime - MUV2Time);
        }
    }

    FillHisto("BurstID", MUV1Event->GetBurstID());

    FillHisto("BeamP",BeamP.Mag());
    FillHisto("Vertex_Z",Vertex.Z());
    FillHisto("Vertex_Y",Vertex.Y());
    FillHisto("Vertex_X",Vertex.X());
    FillHisto("Vertex_cda",cda);


    FillHisto("MUV1_Ncandidates",MUV1Event->GetNCandidates());
    FillHisto("MUV2_Ncandidates",MUV2Event->GetNCandidates());
    FillHisto("MUV3_Ncandidates",MUV3Event->GetNCandidates());
    FillHisto("CHOD_Ncandidates",CHODEvent->GetNCandidates());
    FillHisto("RICH_Ncandidates",RICHEvent->GetNCandidates());
    FillHisto("CEDAR_Ncandidates",CedarEvent->GetNCandidates());
    FillHisto("LKr_Ncandidates" ,LKrEvent ->GetNCandidates());
    FillHisto("STRAW_Nchambers",STRAW_NC);
    FillHisto("TrackChi2",STRAW_chi2);

    FillHisto("TrackPfit_TrackP", STRAW_P - STRAW_Pbf);
    FillHisto("TrackP",STRAW_P);
    //Missing mass squared
    FillHisto("MM2",NuNubarP.M2()*0.000001); //Converting MeV^2 to GeV^2
    FillHisto("Track_P_vs_MM2",TrackP.Mag()*0.001, NuNubarP.M2()*0.000001); //Converting MeV^2 to GeV^2
    FillHisto("Track_P_vs_Theta",TrackP.Mag()*0.001, TMath::ACos(theta) ); //Converting MeV^2 to GeV^2


}

void Kmu2::PostProcess(){
    /// \MemberDescr
    /// This function is called after an event has been processed by all analyzers. It could be used to free some memory allocated
    /// during the Process.
    /// \EndMemberDescr

}

void Kmu2::EndOfBurstUser(){
    /// \MemberDescr
    /// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the end of the last file\n
    /// Do here your start/end of burst processing if any
    /// \EndMemberDescr
}

void Kmu2::EndOfRunUser(){
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

void Kmu2::DrawPlot(){
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
TVector3 Kmu2::VertexCDA(TVector3 pos1, TVector3 p1, TVector3 pos2, TVector3 p2, Double_t &cda)
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

int Kmu2::FindClosestCluster( TRecoVEvent* Event, TVector3 Extrap_track, string detector_type, double& minimum){

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
