#pragma once
#include "CommonTypes.h"
#include "Assert.h"

#include <functional>
#include <vector>
#include <map>
#include <array>
#include <iostream>
#include <fstream>

class Output
{
    public:
        using outputFunc = std::function<void(std::string name)>;
        using perIteroutputFunc = std::function<void(std::ofstream& ofs)>;

        Output() = default;

        void registerOutputFunc(std::string name, outputFunc func);
        outputFunc& getOutputFuncByName(std::string name);

        void registerPerIterOutputFunc(std::string name, perIteroutputFunc func);
        perIteroutputFunc& getPerIterOutputFuncByName(std::string name);

        const std::vector<std::string>& getOutputNames() const { return OutputNames_;}

    private:
        std::map<std::string, outputFunc> MapNameToOutputFunc_;
        std::map<std::string, perIteroutputFunc> MapNameToPerIterOutputFunc_;
        std::vector<std::string> OutputNames_;
};