#include <vector>
#include <iostream>

#include <TFile.h>
#include <TTreeReaderValue.h>

#include "main.h"

int main()
{
    try
    {
        LoadRequiredLibraries();

        // 119 35
        // 055 24
        // const std::string InputFileName  = "pixie_bigrips_traces_119_35.root";
        const std::string InputFileName  = "pixie_bigrips_traces_055_24.root";

        TFile* InputFile = OpenRootFile(InputFileName.c_str());
        TTree* Tree = GetTree(InputFile, "pspmt");

        // Extract run numbers and create trace directory
        const auto RunNumbers = ExtractRunNumbers(InputFileName);
        const std::string OutputDirectory = CreateTraceDirectory(RunNumbers);

        // Get all qualifying events
        const std::vector<Long64_t> QualifyingEvents = GetAllQualifyingEvents(Tree);

        // Graph first N events
        constexpr Long64_t NumberOfEventsToProcess = 100;
        GraphFirstNEvents(Tree, QualifyingEvents, NumberOfEventsToProcess,  OutputDirectory.c_str());

        InputFile->Close();
        delete InputFile;

        return 0;
    }
    catch (const std::exception& Error)
    {
        std::cerr << "Error: " << Error.what() << std::endl;
        return 1;
    }
}