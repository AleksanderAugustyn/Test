#include "ROOTFileReader.h"
#include <iostream>
#include <filesystem>

const std::string ROOTFileReader::DATA_PATH = "/home/aaugustyn/data/";

ROOTFileReader::ROOTFileReader() : ROOTFile(nullptr), PSPMTTree(nullptr)
{
}

ROOTFileReader::~ROOTFileReader()
{
    CloseFile();
}

bool ROOTFileReader::OpenFile(const std::string &filename)
{
    // Close any existing file
    CloseFile();

    // Construct full path
    const std::string FullPath = DATA_PATH + filename;

    // Check if file exists
    if (!std::filesystem::exists(FullPath))
    {
        std::cerr << "Error: File " << FullPath << " does not exist!" << std::endl;
        return false;
    }

    // Try to open the ROOT file
    ROOTFile.reset(TFile::Open(FullPath.c_str(), "READ"));

    if (!ROOTFile || ROOTFile->IsZombie())
    {
        std::cerr << "Error: Could not open ROOT file " << FullPath << std::endl;
        ROOTFile.reset();
        return false;
    }

    // Get the pspmt tree
    PSPMTTree = dynamic_cast<TTree *>(ROOTFile->Get("pspmt"));

    if (!PSPMTTree)
    {
        std::cerr << "Error: Could not find TTree 'pspmt' in file " << FullPath << std::endl;
        CloseFile();
        return false;
    }

    return true;
}

void ROOTFileReader::CloseFile()
{
    if (ROOTFile)
    {
        ROOTFile->Close();
        ROOTFile.reset();
    }
    PSPMTTree = nullptr; // Tree pointer is owned by the file
}

TTree *ROOTFileReader::GetPSPMTTree() const
{
    return PSPMTTree;
}

bool ROOTFileReader::IsOpen() const
{
    return (ROOTFile != nullptr && !ROOTFile->IsZombie());
}
