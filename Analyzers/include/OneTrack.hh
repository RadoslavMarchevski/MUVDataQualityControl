#ifndef ONETRACK_HH
#define ONETRACK_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "TRecoVEvent.hh"
#include "DetectorAcceptance.hh"
#include <TCanvas.h>

class TH1I;
class TH2F;
class TGraph;
class TTree;


class OneTrack : public NA62Analysis::Analyzer
{
public:
    OneTrack(NA62Analysis::Core::BaseAnalysis *ba);
    void InitHist();
    void InitOutput();
    void DefineMCSimple();
    void Process(int iEvent);
    void StartOfBurstUser();
    void EndOfBurstUser();
    void StartOfRunUser();
    void EndOfRunUser();
    void PostProcess();
    void DrawPlot();
    TVector3 VertexCDA(TVector3 pos1, TVector3 p1, TVector3 pos2, TVector3 p2, Double_t &cda);
    int FindClosestCluster(TRecoVEvent* Event, TVector3 Extrap_track, string detector_type, double& dtrkcl_min);

protected:


};
#endif
