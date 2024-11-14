#include <vector>
#include <TFile.h>
#include <TTreeReaderValue.h>
#include <TSystem.h>

#include "main.h"

int main()
{
    try
    {
        LoadRequiredLibraries();

        TFile* InputFile = OpenRootFile("pixie_bigrips_traces_055_24.root");
        TTree* Tree = GetTree(InputFile, "pspmt");

        // Create output directory
        gSystem->Exec("mkdir -p traces");

        // Get all qualifying events
        const std::vector<Long64_t> QualifyingEvents = GetAllQualifyingEvents(Tree);

        // Graph first N events
        constexpr Long64_t NumberOfEventsToProcess = 10;
        GraphFirstNEvents(Tree, QualifyingEvents, NumberOfEventsToProcess, "traces");

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