#ifndef ONETRACKSELECTION_HH
#define ONETRACKSELECTION_HH

#include <stdlib.h>
#include <vector>
#include "Analyzer.hh"
#include "MCSimple.hh"
#include "DetectorAcceptance.hh"
#include "TRecoVEvent.hh"
#include <TCanvas.h>

class TH1I;
class TH2F;
class TGraph;
class TTree;


class OneTrackSelection : public NA62Analysis::Analyzer
{
public:
    OneTrackSelection(NA62Analysis::Core::BaseAnalysis *ba);
    int FindClosestCluster( TRecoVEvent* Event, TVector3 Extrap_track, string detector_type, double& minimum);
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
protected:


};
#endif
