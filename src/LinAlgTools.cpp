#include "LinAlgTools.h"

LinAlg3x3::Real3 LinAlg3x3::CrossProduct(const Real3& v1, const Real3& v2){
    Real3 ret;
    ret.fill(0);

    ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
    ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
    ret[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return ret;
}

LinAlg3x3::Real LinAlg3x3::DotProduct(const Real3& v1, const Real3& v2){
    Real ret = 0.0;

    for (int i=0;i<3;i++){
        ret += v1[i]*v2[i];
    }

    return ret;
}

LinAlg3x3::Real LinAlg3x3::MatrixDeterminant(Matrix& mat){
    Real firstterm = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]);
    Real secondterm = mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]);
    Real thirdterm = mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

    return firstterm - secondterm + thirdterm;
}

LinAlg3x3::Real LinAlg3x3::norm(const Real3& v1){
    Real ret = DotProduct(v1,v1);
    ret = std::sqrt(ret);

    return ret;
}

LinAlg3x3::Matrix LinAlg3x3::dyad(const Real3& v1, const Real3& v2){
    Matrix ret;

    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            ret[i][j] = v1[i]*v2[j];
        }
    }

    return ret;
}

void LinAlg3x3::normalize(Real3& v1){
    Real sqr = norm(v1);
    v1 = v1 / sqr;
}

LinAlg3x3::Real3 LinAlg3x3::MatrixDotVector(const Matrix& A, const Real3& v1){
    Real3 ans;
    ans.fill(0);

    for (int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ans[i] += A[i][j]*v1[j];
        }
    }

    return ans;
}

LinAlg3x3::Matrix LinAlg3x3::GetRotationMatrix(const Real3& v1, const Real3& v2){
    Real3 crossProduct = LinAlg3x3::CrossProduct(v1, v2);
    Real cosine = LinAlg3x3::DotProduct(v1, v2);

    Matrix ret;
    Real denom = 1.0 + cosine;
    Real factor;

    if (std::abs(denom) < 1e-7)
    {
    }
    else{factor = 1.0/(1.0 + cosine);}


    ret[0][0] = 1 + factor*(-crossProduct[2]*crossProduct[2] - crossProduct[1]*crossProduct[1]);
    ret[0][1] = -crossProduct[2] + factor*crossProduct[0]*crossProduct[1];
    ret[0][2] = crossProduct[1] + factor*crossProduct[0]*crossProduct[2];
    ret[1][0] = crossProduct[2] + factor*crossProduct[0]*crossProduct[1];
    ret[1][1] = 1 + factor*(-crossProduct[2]*crossProduct[2] - crossProduct[0]*crossProduct[0]);
    ret[1][2] = -crossProduct[0] + factor*crossProduct[1]*crossProduct[2];
    ret[2][0] = -crossProduct[1] + factor*crossProduct[2]*crossProduct[0];
    ret[2][1] = crossProduct[0] + factor*crossProduct[2]*crossProduct[1];
    ret[2][2] = 1 + factor*(-crossProduct[1]*crossProduct[1] - crossProduct[0]*crossProduct[0]);

    return ret;
}

void LinAlg3x3::RotateBasisSet(Real3& N1, Real3& N2, const Real3& oldu1, const Real3& oldv1, Real3& u1, Real3& v1){
    LinAlg3x3::normalize(N1);
    LinAlg3x3::normalize(N2);

    Real3 crossProduct = LinAlg3x3::CrossProduct(N1, N2);
    Real cosine = LinAlg3x3::DotProduct(N1, N2);

    Matrix ret;
    Real denom = 1.0 + cosine;
    Real factor;

    if (std::abs(cosine) < 1e-5) { factor = 0.0;}
    else if (std::abs(denom) < 1e-5){
        for (int i=0;i<3;i++){
            u1[i] = -oldu1[i];
            v1[i] = -oldv1[i];
        }

        return;
    }
    else{
        factor = 1.0/(1.0+cosine);
    }

    ret[0][0] = 1 + factor*(-crossProduct[2]*crossProduct[2] - crossProduct[1]*crossProduct[1]);
    ret[0][1] = -crossProduct[2] + factor*crossProduct[0]*crossProduct[1];
    ret[0][2] = crossProduct[1] + factor*crossProduct[0]*crossProduct[2];
    ret[1][0] = crossProduct[2] + factor*crossProduct[0]*crossProduct[1];
    ret[1][1] = 1 + factor*(-crossProduct[2]*crossProduct[2] - crossProduct[0]*crossProduct[0]);
    ret[1][2] = -crossProduct[0] + factor*crossProduct[1]*crossProduct[2];
    ret[2][0] = -crossProduct[1] + factor*crossProduct[2]*crossProduct[0];
    ret[2][1] = crossProduct[0] + factor*crossProduct[2]*crossProduct[1];
    ret[2][2] = 1 + factor*(-crossProduct[1]*crossProduct[1] - crossProduct[0]*crossProduct[0]);

    u1 = MatrixDotVector(ret, oldu1);
    v1 = MatrixDotVector(ret, oldv1);

    LinAlg3x3::normalize(u1);
    LinAlg3x3::normalize(v1);

    return;
}

LinAlg3x3::Real LinAlg3x3::findCosangle(Real3 vec1, Real3 vec2){
    Real dot = DotProduct(vec1, vec2);
    Real normA = norm(vec1);
    Real normB = norm(vec2);

    return dot/(normA*normB);
}

LinAlg3x3::Real LinAlg3x3::findSinangle(Real3 vec1, Real3 vec2){
    Real3 cross = CrossProduct(vec1, vec2);
    Real normcross = norm(cross);
    Real normA = norm(vec1);
    Real normB = norm(vec2);

    return normcross/(normA * normB);
}

LinAlg3x3::Real LinAlg3x3::findAngle(Real3 vec1, Real3 vec2){
    Real sin_val = LinAlg3x3::norm(LinAlg3x3::CrossProduct(vec1, vec2));
    Real cos_val = LinAlg3x3::DotProduct(vec1, vec2);
    Real result = std::atan2(sin_val, cos_val);

    return std::fabs(result);
}