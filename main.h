#pragma once

#include <optional>

#include <TGraph.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>

#include "PaassRootStruct.hpp"

struct AnalysisResults
{
    Long64_t EventNumber = -1;
    Double_t PosX = -1;
    Double_t PosY = -1;

    // Add fit parameters for each channel
    struct ChannelFit
    {
        Double_t Amplitude = -1;
        Double_t PeakPosition = -1;
        Double_t DecayConstant = -1;
        Double_t RiseTimeConstant = -1;
        Double_t RisePower = -1;
        Double_t Baseline = -1;
    };

    std::map<std::string, ChannelFit> AnodeFits; // xa, xb, ya, yb

    struct DynodeFit
    {
        Double_t Amplitude = -1;
        Double_t PeakPosition = -1;
        Double_t FastDecay = -1;
        Double_t SlowDecay = -1;
        Double_t RiseTime = -1;
        Double_t UndershootAmp = -1;
        Double_t UndershootRecovery = -1;
        Double_t FastFraction = -1;
        Double_t Baseline = -1;
    } DynodeFitParams;
};

struct AnalysisHistograms
{
    std::map<std::string, std::vector<TH2D *> > ScatterPlots;
    std::map<std::string, std::vector<TProfile *> > Profiles;
    std::map<std::string, TProfile2D *> Profile2Ds;
    std::map<std::string, TH1D *> RisePowerHists;
    std::map<std::string, TH2D *> CountHist;
    std::map<std::string, TProfile2D *> FilteredRisePowerScatter;
};

struct ChannelFitData
{
    Double_t DecayTime;
    Double_t RiseTime;
    Double_t RisePower;
};

// RootInput
void LoadRequiredLibraries();

TFile *OpenRootFile(const char *FileName);

TTree *GetTree(TFile *InputFile, const char *TreeName);

std::pair<Int_t, Int_t> ExtractRunNumbers(const std::string &FileName);

std::string CreateTraceDirectory(const std::pair<Int_t, Int_t> &RunNumbers);

// EventSelection
std::vector<Long64_t> GetAllQualifyingEvents(TTree *TreeInput);

Bool_t MeetsSelectionCriteria(TTree *TreeInput, Long64_t Entry);

std::vector<Long64_t> ScanEvents(TTree *Tree, Long64_t MaxEventsToSave);

// TraceGraphs
void SaveTraceGraphs(TTree *TreeInput, Long64_t Entry, const char *ImagePath);

TGraph *CreateTraceGraph(const processor_struct::ROOTDEV &Device, const std::string &Title, Int_t DeviceIndex);

void SaveTraceGraphsWithFit(TTree *TreeInput, Long64_t Entry, const char *ImagePath);

void GraphFirstNEvents(TTree *TreeInput, const std::vector<Long64_t> &QualifyingEvents, Long64_t NumberOfEvents, const char *OutputPath);

std::optional<AnalysisResults::DynodeFit> ExtractDynodeFitParameters(const TF1 *FitFunc);

std::optional<AnalysisResults::ChannelFit> ExtractAnodeFitParameters(const TF1 *FitFunc);

std::optional<AnalysisResults> GetEventFitParameters(TTree *TreeInput, Long64_t Entry);

void SaveAnalysisResults(const std::vector<AnalysisResults> &Results, Int_t RunNumber, Int_t SubRunNumber);

// FitAnalysis
Double_t CalculateRisePower(const std::string &Channel, Double_t Position);

Double_t AnodePeakFunction(const Double_t *X, const Double_t *Parameters);

TF1 *FitPeakToTrace(TGraph *TraceGraph, Double_t FitRangeStart,
                    Double_t FitRangeEnd, const std::string &Channel,
                    Double_t PosX, Double_t PosY);

Double_t DynodePeakFunction(const Double_t *X, const Double_t *Parameters);

TF1 *FitDynodePeak(TGraph *TraceGraph, Double_t FitRangeStart, Double_t FitRangeEnd);

// AnalyseTraces
TH2D *CreateScatterPlot(const char *Name, const char *Title,
                        const char *XTitle, const char *YTitle,
                        Int_t NBinsX, Double_t XMin, Double_t XMax,
                        Int_t NBinsY, Double_t YMin, Double_t YMax);

AnalysisHistograms InitializeHistograms(const std::vector<std::string> &Channels);

void ProcessInputFiles(const std::vector<std::pair<Int_t, Int_t> > &RunsToAnalyze,
                       AnalysisHistograms &Histograms,
                       const std::vector<std::string> &Channels);

void CreateAndSaveChannelPlots(const AnalysisHistograms &Histograms,
                               const std::vector<std::string> &Channels,
                               const char *OutputDirectory);

void CleanupHistograms(AnalysisHistograms &Histograms);

void AnalyzePositionVsFitParameters(const std::vector<std::pair<Int_t, Int_t> > &RunsToAnalyze,
                                    const char *OutputDirectory);
