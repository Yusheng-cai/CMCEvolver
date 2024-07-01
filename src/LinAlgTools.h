#pragma once
#include "tools/CommonTypes.h"
#include "tools/Assert.h"
#include "tools/CommonOperations.h"
#include "tools/Constants.h"

#include <cmath>

namespace LinAlg3x3
{
    using Real = CommonTypes::Real;
    using Matrix = CommonTypes::Matrix;
    using Real3= CommonTypes::Real3;

    Real3 CrossProduct(const Real3& v1, const Real3& v2);
    Real DotProduct(const Real3& v1, const Real3& v2);
    Real norm(const Real3& v1);
    void normalize(Real3& v1);
    Matrix dyad(const Real3& v1, const Real3& v2);

    Real3 MatrixDotVector(const Matrix& mat1, const Real3& v1); 

    // This is the rotation matrix that rotates vector 1 onto vector 2
    Matrix GetRotationMatrix(const Real3& N1, const Real3& N2);

    // Rotate basis set
    void RotateBasisSet(Real3& N1, Real3& N2, const Real3& oldu1, const Real3& oldv1, Real3& u1, Real3& v1);

    // find cos(angle) between 2 vectors 
    Real findCosangle(Real3 vec1, Real3 vec2);

    // find sin(angle) between 2 vectors 
    Real findSinangle(Real3 vec1, Real3 vec2);

    // find angle
    Real findAngle(Real3 vec1, Real3 vec2);

    // find determinant
    Real MatrixDeterminant(Matrix& mat);
}