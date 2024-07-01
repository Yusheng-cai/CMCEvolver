#include "src/AtomGroupParsingStrategy.h"
#include "tools/Assert.h"
#include "tools/CommonTypes.h"
#include "tools/InputParser.h"
#include "tools/CommandLineArguments.h"
#include "tools/FileSystem.h"

#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <sstream>

int main(int argc, char** argv)
{
    ASSERT((argc > 1), "There must be at least 1 input argument which is the path to the input file.");
    CommandLineArguments cmd(argc, argv);
    std::string abspath;
    bool apathread = cmd.readString("apath", CommandLineArguments::Keys::Optional,abspath);

    std::string fname(argv[1]);
    InputParser ip;
    ParameterPack param;

    ip.ParseFile(fname, param);

    auto AGPacks = param.findParamPacks("atomgroup", ParameterPack::KeyType::Required);
    for (int i=0;i<AGPacks.size();i++)
    {
        auto& ag = AGPacks[i];
        std::vector<std::string> sel;
        std::string output;
        ag->ReadVectorString("selection", ParameterPack::KeyType::Required,sel);
        ag->ReadString("outputfile", ParameterPack::KeyType::Required,output);

        if (apathread)
        {
            sel[0] = FileSystem::joinPath(abspath, fname); 
        }

        GroFile g;
        AtomGroupParsingInput input = {g, sel};
        auto parsePtr = AtomGroupParsingRegistry::Factory::instance().create(sel[0],input);

        auto castPtr = dynamic_cast<IndexFileParsing*>(parsePtr);

        ASSERT((castPtr != nullptr), "The passed in argument is not for Index File parsing");

        std::ofstream ofs_(output);

        for (int i=0;i<castPtr->getNumFrames();i++)
        {
            ofs_ << castPtr -> getFrame()[i] << " ";
            std::vector<int> indices_;

            castPtr ->update(indices_);

            for (int j=0;j<indices_.size();j++)
            {
                ofs_ << indices_[j] << " ";
            }

            ofs_ << "\n";
        }

        ofs_.close();
    }
}