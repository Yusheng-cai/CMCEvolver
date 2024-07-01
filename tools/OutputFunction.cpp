#include "OutputFunction.h"

void Output::registerOutputFunc(std::string name, outputFunc func)
{
    auto it = MapNameToOutputFunc_.find(name);

    ASSERT((it == MapNameToOutputFunc_.end()), "The output with name " << name << " is already registered.");

    OutputNames_.push_back(name);

    MapNameToOutputFunc_.insert(std::make_pair(name, func));
}

Output::outputFunc& Output::getOutputFuncByName(std::string name)
{
    auto it = MapNameToOutputFunc_.find(name);

    ASSERT((it != MapNameToOutputFunc_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}

void Output::registerPerIterOutputFunc(std::string name, perIteroutputFunc func)
{
    auto it = MapNameToPerIterOutputFunc_.find(name);

    ASSERT((it == MapNameToPerIterOutputFunc_.end()), "The per iter output with name " << name << " is already registered.");

    MapNameToPerIterOutputFunc_.insert(std::make_pair(name, func));
}

Output::perIteroutputFunc& Output::getPerIterOutputFuncByName(std::string name)
{
    auto it = MapNameToPerIterOutputFunc_.find(name);

    ASSERT((it != MapNameToPerIterOutputFunc_.end()), "The output with name " << name << " is not registered.");

    return it -> second;
}