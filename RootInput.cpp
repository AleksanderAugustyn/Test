#include <TROOT.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TSystem.h>

#include "main.h"

/**
 * Loads required ROOT libraries for the analysis
 * Throws runtime_error if any library fails to load
 */
void LoadRequiredLibraries()
{
    const char* RequiredLibraries[] = {
        "libPaassRootStruct.so",
        "libyaml-cpp.so",
        "libTraceAnalyzerLib.so",
        "libmerger_data_dic.so",
        "libMergerLib.so"
    };

    for (const auto& Library : RequiredLibraries)
    {
        if (gSystem->Load(Library) < 0)
        {
            throw std::runtime_error("Failed to load library: " + std::string(Library));
        }
    }
}

/**
 * Opens a ROOT file and returns the file pointer
 * @param FileName Path to the ROOT file
 * @return Pointer to the opened TFile
 * @throws runtime_error if file cannot be opened
 */
TFile* OpenRootFile(const char* FileName)
{
    // /home/aaugustyn/data/FileName
    auto* path = "/home/aaugustyn/data/";
    TFile* InputFile = TFile::Open((std::string(path) + FileName).c_str());

    if (!InputFile || InputFile->IsZombie())
    {
        throw std::runtime_error("Failed to open file: " + std::string(FileName));
    }

    return InputFile;
}

/**
 * Retrieves a tree from an open ROOT file
 * @param InputFile Pointer to the open ROOT file
 * @param TreeName Name of the tree to retrieve
 * @return Pointer to the retrieved TTree
 * @throws runtime_error if tree cannot be found
 */
TTree* GetTree(TFile* InputFile, const char* TreeName)
{
    if (!InputFile)
    {
        throw std::runtime_error("Invalid file pointer");
    }

    const auto Tree = dynamic_cast<TTree*>(InputFile->Get(TreeName));

    if (!Tree)
    {
        throw std::runtime_error("Failed to retrieve tree: " + std::string(TreeName));
    }

    return Tree;
}

