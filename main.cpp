#include <iostream>
#include <string>
#include <stdexcept>

// ROOT headers
#include <TFile.h>
#include "TSystem.h"
#include "TTree.h"


int main()
{
    try
    {
        // rootlogon.C
        gSystem->Load("libPaassRootStruct.so");
        gSystem->Load("libyaml-cpp.so");
        gSystem->Load("libTraceAnalyzerLib.so");
        gSystem->Load("libmerger_data_dic.so");
        gSystem->Load("libMergerLib.so");

        // Open the ROOT file
        const auto Filename = "pixie_bigrips_traces_119_35.root";

        std::cout << "Opening the ROOT file: " << Filename << std::endl;

        TFile *File = TFile::Open(Filename);

        if ( File == nullptr)
        {
            throw std::runtime_error("Failed to open the ROOT file: " + std::string(Filename));
        }

        // Get the tree
        auto *Tree = dynamic_cast<TTree*>(File->Get("pspmt"));
        if (!Tree)
        {
            File->Close();
            throw std::runtime_error("Failed to get the tree from the ROOT file: " + std::string(Filename));
        }

        // Get number of entries
        const Long64_t NumberOfEntries = Tree->GetEntries();
        std::cout << "Number of events: " << NumberOfEntries << std::endl;

        // Setup variables to read the branches
        Int_t HighGainValid;
        Int_t LowGainValid;

        // Set branch addresses
        Tree->SetBranchAddress("high_gain_.valid_", &HighGainValid);
        Tree->SetBranchAddress("low_gain_.valid_", &LowGainValid);

        // Counters for statistics
        Long64_t OnlyHighValidCount = 0;    // Only high gain valid
        Long64_t OnlyLowValidCount = 0;     // Only low gain valid
        Long64_t EitherValidCount = 0;      // Either high or low valid
        Long64_t BothValidCount = 0;        // Both valid
        Long64_t NeitherValidCount = 0;     // Neither valid

        // Loop through all entries
        for (Long64_t i = 0; i < NumberOfEntries; i++)
        {
            Tree->GetEntry(i);

            const bool IsHighValid = (HighGainValid == 1);
            const bool IsLowValid = (LowGainValid == 1);

            // Count exclusive categories
            if (IsHighValid && !IsLowValid) OnlyHighValidCount++;
            if (!IsHighValid && IsLowValid) OnlyLowValidCount++;
            if (IsHighValid && IsLowValid) BothValidCount++;
            if (!IsHighValid && !IsLowValid) NeitherValidCount++;

            // Count inclusive category (either or both valid)
            if (IsHighValid || IsLowValid) EitherValidCount++;
        }

        // Print results with percentages
        std::cout << "\nValidity Statistics:\n" << std::endl;

        std::cout << "Events with only high gain valid: " << OnlyHighValidCount
                  << " (" << (100.0 * static_cast<double>(OnlyHighValidCount) / static_cast<double>(NumberOfEntries)) << "%)" << std::endl;

        std::cout << "Events with only low gain valid: " << OnlyLowValidCount
                  << " (" << (100.0 * static_cast<double>(OnlyLowValidCount) / static_cast<double>(NumberOfEntries)) << "%)" << std::endl;

        std::cout << "Events with either gain valid: " << EitherValidCount
                  << " (" << (100.0 * static_cast<double>(EitherValidCount) / static_cast<double>(NumberOfEntries)) << "%)" << std::endl;

        std::cout << "Events with both gains valid: " << BothValidCount
                  << " (" << (static_cast<double>(BothValidCount) * 100.0 / static_cast<double>(NumberOfEntries)) << "%)" << std::endl;

        std::cout << "Events with neither gain valid: " << NeitherValidCount
                  << " (" << (100.0 * static_cast<double>(NeitherValidCount) / static_cast<double>(NumberOfEntries)) << "%)" << std::endl;

        // Verification that numbers add up correctly
            if (const Long64_t TotalExclusiveEvents = OnlyHighValidCount + OnlyLowValidCount + BothValidCount + NeitherValidCount; TotalExclusiveEvents != NumberOfEntries)
        {
            std::cerr << "\nVerification:" << std::endl;
            std::cerr << "Sum of exclusive categories: " << TotalExclusiveEvents
                      << " (should equal total events: " << NumberOfEntries << ")" << std::endl;

            throw std::runtime_error("Error: Sum of exclusive categories does not equal total events!");
        }

        File->Close();
    }
    catch (const std::exception &E)
    {
        std::cerr << E.what() << std::endl;
        return -1;
    }


    return 0;
}
