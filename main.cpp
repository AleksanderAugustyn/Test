#include <iostream>
#include <string>
#include <stdexcept>

// Project headers
#include "main.h"

// ROOT headers
#include "ROOTFileReader.h"
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

        // Open the file
        ROOTFileReader FileReader;
        FileReader.OpenFile("pixie_bigrips_traces_119_35.root");
        // Get the tree
        TTree *Tree = FileReader.GetPSPMTTree();

    }
    catch (const std::exception &E)
    {
        std::cerr << E.what() << std::endl;
        return -1;
    }


    return 0;
}
