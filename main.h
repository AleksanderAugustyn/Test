#pragma once

#include <TGraph.h>
#include <TFitResult.h>
#include <TFile.h>

#include "PaassRootStruct.hpp"

// RootInput
void LoadRequiredLibraries();
TFile* OpenRootFile(const char* FileName);
TTree* GetTree(TFile* InputFile, const char* TreeName);

// EventSelection
std::vector<Long64_t> GetAllQualifyingEvents(TTree* TreeInput);
void GraphFirstNEvents(TTree* TreeInput, const std::vector<Long64_t>& QualifyingEvents,
    Long64_t NumberOfEvents, const char* OutputPath);
std::vector<Long64_t> ScanEvents(TTree* Tree, Long64_t MaxEventsToSave);
Bool_t MeetsSelectionCriteria(TTree* TreeInput, Long64_t Entry);

// TraceGraphs
void SaveTraceGraphs(TTree* TreeInput, Long64_t Entry, const char* ImagePath);
TGraph* CreateTraceGraph(const processor_struct::ROOTDEV& Device, const std::string& Title, Int_t DeviceIndex);
void SaveTraceGraphsWithFit(TTree* TreeInput, Long64_t Entry, const char* ImagePath);

// FitAnalysis
Double_t PeakFunction(const Double_t* X, const Double_t* Parameters);
TF1* FitPeakToTrace(TGraph* TraceGraph, Double_t FitRangeStart, Double_t FitRangeEnd);