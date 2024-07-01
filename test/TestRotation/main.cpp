#include "src/LinAlgTools.h"
#include "tools/InputParser.h"
#include "tools/Assert.h"
#include "tools/CommandLineArguments.h"
#include "tools/FileSystem.h"
#include "tools/CommonTypes.h"

#include <vector>
#include <array>
#include <string>

int main(int argc, char** argv)
{
    using Real = CommonTypes::Real;
    using Real3= CommonTypes::Real3;
    using Matrix = CommonTypes::Matrix;

    ASSERT(argc > 1, "The input argument must contain the input file.");
    CommandLineArguments cmd(argc, argv);

    Real tol = 1e-7;

    std::string fname(argv[1]); 

    InputParser ip;
    ParameterPack pack;

    ip.ParseFile(fname, pack);

    auto rotationPacks = pack.findParamPacks("rotation", ParameterPack::KeyType::Required);

    for (int i=0;i<rotationPacks.size();i++)
    {
        Real3 vector1;
        Real3 vector2;
        auto rotPack = rotationPacks[i];

        rotPack->ReadArrayNumber("vector1", ParameterPack::KeyType::Required, vector1);
        LinAlg3x3::normalize(vector1);

        rotPack->ReadArrayNumber("vector2", ParameterPack::KeyType::Required, vector2);
        LinAlg3x3::normalize(vector2);

        Matrix rotMat = LinAlg3x3::GetRotationMatrix(vector1, vector2);

        Real3 reconstructed = LinAlg3x3::MatrixDotVector(rotMat, vector1); 

        Real3 diff;

        for (int j=0;j<3;j++)
        {
            diff[j] = reconstructed[j] - vector2[j];
        }

        Real normdiff = LinAlg3x3::norm(diff);
        ASSERT((normdiff <= tol), "The difference is too large and it is " << normdiff);
    }
}