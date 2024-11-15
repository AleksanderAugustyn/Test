#pragma once

#include <TGraph.h>
#include <TFitResult.h>
#include <TFile.h>

#include "PaassRootStruct.hpp"

// RootInput
void LoadRequiredLibraries();
TFile* OpenRootFile(const char* FileName);
TTree* GetTree(TFile* InputFile, const char* TreeName);
std::pair<Int_t, Int_t> ExtractRunNumbers(const std::string& FileName);
std::string CreateTraceDirectory(const std::pair<Int_t, Int_t>& RunNumbers);

// EventSelection
std::vector<Long64_t> GetAllQualifyingEvents(TTree* TreeInput);
Bool_t MeetsSelectionCriteria(TTree* TreeInput, Long64_t Entry);
std::vector<Long64_t> ScanEvents(TTree* Tree, Long64_t MaxEventsToSave);

// TraceGraphs
void SaveTraceGraphs(TTree* TreeInput, Long64_t Entry, const char* ImagePath);
TGraph* CreateTraceGraph(const processor_struct::ROOTDEV& Device, const std::string& Title, Int_t DeviceIndex);
void SaveTraceGraphsWithFit(TTree* TreeInput, Long64_t Entry, const char* ImagePath);
void GraphFirstNEvents(TTree* TreeInput, const std::vector<Long64_t>& QualifyingEvents, Long64_t NumberOfEvents, const char* OutputPath);

// FitAnalysis
Double_t AnodePeakFunction(const Double_t* X, const Double_t* Parameters);
TF1* FitPeakToTrace(TGraph* TraceGraph, Double_t FitRangeStart, Double_t FitRangeEnd);
Double_t DynodePeakFunction(const Double_t* X, const Double_t* Parameters);
TF1* FitDynodePeak(TGraph* TraceGraph, Double_t FitRangeStart, Double_t FitRangeEnd);