#include "Mesh.h"

MeshRefineStrategy::MeshRefineStrategy(MeshRefineStrategyInput& input)
:pack_(input.pack)
{
    pack_.ReadString("name", ParameterPack::KeyType::Required, name_);
}