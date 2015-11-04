#include <stdlib.h>
#include <iostream>
#include <TChain.h>
#include "OneTrackSelection.hh"
#include "MCSimple.hh"
#include "functions.hh"
#include "Event.hh"
#include "Persistency.hh"
#include "TRecoVEvent.hh"
#include "Definition.h"

using namespace std;
using namespace NA62Analysis;
using namespace NA62Constants;


OneTrackSelection::OneTrackSelection(Core::BaseAnalysis *ba) : Analyzer(ba, "OneTrackSelection")
{

    RequestTree("LKr",new TRecoLKrEvent);
    RequestTree("Spectrometer",new TRecoSpectrometerEvent);
    RequestTree("MUV1",new TRecoMUV1Event);
    RequestTree("MUV2",new TRecoMUV2Event);
    RequestTree("MUV3",new TRecoMUV3Event);
    RequestTree("RICH",new TRecoRICHEvent);
    RequestTree("CHOD",new TRecoCHODEvent);
    RequestTree("Cedar",new TRecoCedarEvent);


}

void OneTrackSelection::InitOutput(){
}

void OneTrackSelection::InitHist(){
    BookHisto(new TH1F("CEDAR_timediff"," CEDAR_{time} - CHOD_{time} ; CEDAR_{time}- CHOD_{time} [ns]", 400, -100, 100.));
}

void OneTrackSelection::DefineMCSimple(){
}

void OneTrackSelection::StartOfRunUser(){
}

void OneTrackSelection::StartOfBurstUser(){
}

void OneTrackSelection::Process(int iEvent){
//    if(fMCSimple.fStatus == MCSimple::kMissing){printIncompleteMCWarning(iEvent);return;}
//    if(fMCSimple.fStatus == MCSimple::kEmpty){printNoMCWarning();return;}

    TRecoLKrEvent *LKrEvent = (TRecoLKrEvent*)GetEvent("LKr");
    TRecoSpectrometerEvent *SpectrometerEvent = (TRecoSpectrometerEvent*)GetEvent("Spectrometer");
    TRecoMUV1Event *MUV1Event = (TRecoMUV1Event*)GetEvent("MUV1");
    TRecoMUV2Event *MUV2Event = (TRecoMUV2Event*)GetEvent("MUV2");
    TRecoMUV3Event *MUV3Event = (TRecoMUV3Event*)GetEvent("MUV3");
    TRecoRICHEvent *RICHEvent = (TRecoRICHEvent*)GetEvent("RICH");
    TRecoCHODEvent *CHODEvent = (TRecoCHODEvent*)GetEvent("CHOD");
    TRecoCedarEvent *CedarEvent = (TRecoCedarEvent*)GetEvent("Cedar");

    //CUTComment:: Only one candidate in the STRAW
    if(SpectrometerEvent->GetNCandidates() != 1 ){return;}

    TRecoSpectrometerCandidate* Track = ((TRecoSpectrometerCandidate*)SpectrometerEvent->GetCandidate(0));

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

    //CUTComment:: At least one candidate in the Cedar
    if(CedarEvent->GetNCandidates() == 0){return;}

    for( int iCedarCand=0; iCedarCand < CedarEvent->GetNCandidates(); iCedarCand++){

        TRecoCedarCandidate* CedarCandidate = ((TRecoCedarCandidate*)CedarEvent->GetCandidate(iCedarCand));

        //CUTComment:: At least 5 sectors in the CEDAR
        if(CedarCandidate->GetNSectors() < 5 ){return;}

    }

    //Getting the position vectors and the slopes of the tracks
    //given by the spectrometer before the magnet @ DCH1 (dxdz and dydz)
    TVector3 SlopesBefore;
    TVector3 PositionBefore;

    SlopesBefore.SetX(STRAW_bdxdz);
    SlopesBefore.SetY(STRAW_bdydz);
    SlopesBefore.SetZ(1.);
    PositionBefore = Track->GetPositionBeforeMagnet();

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

    double CHODR  = sqrt( pow(CHOD_extrap.X(),2) + pow(CHOD_extrap.Y(),2) );

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
    if(CHODClosestTrackIndex < 0. ){return;}

    //Setting up the time of the track. STRAW time is not available, therefore
    //CHOD time will be used as the track time
    double STRAW4R  = sqrt( pow(PositionAfter.X(),2) + pow(PositionAfter.Y(),2) );
    double CD_TrackTime    = ((TRecoCHODCandidate*)CHODEvent->GetCandidate(CHODClosestTrackIndex))->GetTime();

    //CUTComment:: Distance to the extrapolated track in the CHOD bigger than 8 cm
    //Track to be in the CHOD geometrical  acceptance
    if(CHODdtrkcl_min > 80. ) {return;}
    if(CHODR < 100. || CHODR > 1200) {return;} //[mm]

    //Cedar event with the closest time to the track time selected

    for(int iCedarCand=0; iCedarCand < CedarEvent->GetNCandidates(); iCedarCand++){

        TRecoCedarCandidate* CedarCandidate = ((TRecoCedarCandidate*)CedarEvent->GetCandidate(iCedarCand));

        double CedarCandidateTime = CedarCandidate->GetTime();
        double CedarTimeDiff      = CD_TrackTime - CedarCandidateTime;

        CEDAR_STRAW_tdiff.push_back(CedarTimeDiff);

    }


    if(CEDAR_STRAW_tdiff.size() !=0){
        double CedarTime = *min_element(CEDAR_STRAW_tdiff.begin(), CEDAR_STRAW_tdiff.end());

        //CUTComment:: Cedar time difference cut
        if(fabs(CedarTime) > 3){return;}
        FillHisto("CEDAR_timediff", CedarTime);
    }

    //Quality of the LKr cluster, associated with the track (if any)
    if(LKrTrackClusterIndex > -1){
        TRecoLKrCandidate* CD_LKrCluster = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex));

        double CD_LKrClusterTime  = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterTime();
        double CD_LKrClusterDDead = ((TRecoLKrCandidate*)LKrEvent->GetCandidate(LKrTrackClusterIndex))->GetClusterDDeadCell();
        double LKrTrkTime = CD_TrackTime - CD_LKrClusterTime + 120;
        double LKrR  = sqrt( pow(LKr_extrap.X(),2) + pow(LKr_extrap.Y(),2) );
        double Cluster_X = CD_LKrCluster->GetClusterX()*10;
        double Cluster_Y = CD_LKrCluster->GetClusterY()*10;


        //CUTComment:: LKr cuts
        //1. Track to be inside of the LKr Geometrical acceptance
        //-----a) Inner diameter 15 cm
        //-----a) Outer diameter 110 cm
        //3. Distance between track and closest cluster 5cm
        //4. Distance to deadcell > 2cm
        //5. Cluster to be in time with the Track ( Current cut 10 ns)

        if(LKrR < 150. || LKrR > 1100. ) {return;}

        if(LKrdtrkcl_min  > 50. ) {return;}

        if(CD_LKrClusterDDead < 2.){return;}

        if(fabs(LKrTrkTime) > 10){return;}

    }

    //Quality of the MUV1 cluster, associated with the track (if any)
    if(MUV1TrackClusterIndex > -1){

        TRecoMUV1Candidate* CD_MUV1Cluster = ((TRecoMUV1Candidate*)MUV1Event->GetCandidate(MUV1TrackClusterIndex));

        double CD_MUV1ClusterTime = ((TRecoMUV1Candidate*)MUV1Event->GetCandidate(MUV1TrackClusterIndex))->GetTime();
        double MUV1TrkTime = CD_TrackTime - CD_MUV1ClusterTime;

        //CUTComment:: MUV1 cuts
        //1. Track to be inside of the MUV1 Geometrical acceptance
        //-----a) Inner Square 26 cm
        //-----b) Outer Square 220 cm
        //2. Track search square 12 cm ( 2 strips)
        //3. Distance between track and associated cluster 12 cm ??
        //4. Cluster to be in time with the Track ( Current cut 20 ns)

        if(fabs(MUV1_extrap.X()) <= 130. && fabs(MUV1_extrap.Y()) <= 130.){return;}
        if(fabs(MUV1_extrap.X()) >= 1100. || fabs(MUV1_extrap.Y()) >= 1100.){return;}

        if( fabs(CD_MUV1Cluster->GetPosition().X() - MUV1_extrap.X()) > 120. ||
            fabs(CD_MUV1Cluster->GetPosition().Y() - MUV1_extrap.Y()) > 120.){return;}

        if(fabs(MUV1TrkTime) > 20){return;}
        //if(MUV1dtrkcl_min > 120.) {return;}

    }

    //Quality of the MUV2 cluster, associated with the track (if any)
    if(MUV2TrackClusterIndex > -1){

        TRecoMUV2Candidate* CD_MUV2Cluster = ((TRecoMUV2Candidate*)MUV2Event->GetCandidate(MUV2TrackClusterIndex));

        double CD_MUV2ClusterTime   = ((TRecoMUV2Candidate*)MUV2Event->GetCandidate(MUV2TrackClusterIndex))->GetTime();
        double MUV2TrkTime = CD_TrackTime - CD_MUV2ClusterTime;

        //CUTComment:: MUV2 cuts
        //1. Track to be inside of the MUV2 Geometrical acceptance
        //-----a) Inner Square 26 cm
        //-----b) Outer Square 220 cm
        //2. Track search square 24 cm ( 2 strips)
        //3. Distance between track and associated cluster 24 cm ??
        //4. Cluster to be in time with the Track ( Current cut 20 ns)

        if(fabs(MUV2_extrap.X()) <= 130.  && fabs(MUV2_extrap.Y()) <= 130.){return;}
        if(fabs(MUV2_extrap.X()) >= 1100. || fabs(MUV2_extrap.Y()) >= 1100.){return;}

        if( fabs(CD_MUV2Cluster->GetPosition().X() - MUV2_extrap.X()) > 240. ||
            fabs(CD_MUV2Cluster->GetPosition().Y() - MUV2_extrap.Y()) > 240.){return;}

        //if(MUV2dtrkcl_min > 240.) {return;}
        if(fabs(MUV2TrkTime) > MUV2OffsetCut){return;}
    }

    //Quality of the MUV3 cluster, associated with the track (if any)
    if(MUV3TrackClusterIndex > -1){

        TRecoMUV3Candidate* CD_MUV3Cluster = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex));

        double CD_MUV3ClusterTime = ((TRecoMUV3Candidate*)MUV3Event->GetCandidate(MUV3TrackClusterIndex))->GetTime();
        double MUV3TrkTime = CD_TrackTime - CD_MUV3ClusterTime;

        //CUTComment:: MUV3 cuts
        //1. Track to be inside of the MUV3 Geometrical acceptance
        //-----a) Inner Square 26 cm
        //-----b) Outer Square 220 cm
        //2. Track search square 20 cm ( less than one block )
        //-----Reason: Interested in only straight tracks with no multiple scattering
        //3. Distance between track and associated cluster 24 cm ??
        //4. Cluster to be in time with the Track ( Current cut 5 ns)

        if(fabs(MUV3_extrap.X()) <= 130. && fabs(MUV3_extrap.Y()) <= 130.){return;}
        if(fabs(MUV3_extrap.X()) >= 1100. || fabs(MUV3_extrap.Y()) >= 1100.){return;}

        if( fabs(CD_MUV3Cluster->GetX() - MUV3_extrap.X()) > 200. ||
            fabs(CD_MUV3Cluster->GetY() - MUV3_extrap.Y()) > 200.){return;}

        if(fabs(MUV3TrkTime) > 5){return;}
        //if(MUV2dtrkcl_min > 150.) {return;}


    }



}

void OneTrackSelection::PostProcess(){

}

void OneTrackSelection::EndOfBurstUser(){
    /// \MemberDescr
    /// This method is called when a new file is opened in the ROOT TChain (corresponding to a start/end of burst in the normal NA62 data taking) + at the end of the last file\n
    /// Do here your start/end of burst processing if any
    /// \EndMemberDescr
}

void OneTrackSelection::EndOfRunUser(){
    SaveAllPlots();
}

void OneTrackSelection::DrawPlot(){
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
int OneTrackSelection::FindClosestCluster( TRecoVEvent* Event, TVector3 Extrap_track, string detector_type, double& minimum){

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
