//
// Created by tianyan on 1/4/26.
//

#ifndef TINYVBD_MATH_HPP
#define TINYVBD_MATH_HPP
#include "Types.h"

namespace TY {

    [[nodiscard]] inline Mat3 Skew(const Vec3& v) {
    Mat3 m;
    m <<  0.0f,   -v.z(),  v.y(),
          v.z(),   0.0f,  -v.x(),
         -v.y(),   v.x(),  0.0f;
    return m;
    }

    [[nodiscard]] inline Mat3 ComputeNormalizedVectorDerivative(
        const float unnormalized_vec_length,
        const Vec3& normalized_vec,
        const Mat3& unnormalized_vec_derivative) {

        constexpr float eps = 1.0e-8f;
        if (unnormalized_vec_length < eps) {
            return Mat3::Zero();
        }

        const Mat3 I = Mat3::Identity();
        const Mat3 projection = I - (normalized_vec * normalized_vec.transpose()); // outer

        return (1.0f / unnormalized_vec_length) * projection * unnormalized_vec_derivative;
    }

    [[nodiscard]] inline Vec3 ComputeAngleDerivative(
    const Vec3& n1_hat,
    const Vec3& n2_hat,
    const Vec3& e_hat,
    const Mat3& dn1hat_dx,
    const Mat3& dn2hat_dx,
    const float sin_theta,
    const float cos_theta,
    const Mat3& skew_n1,
    const Mat3& skew_n2) {

        const Vec3 dsin_dx = (skew_n1 * dn2hat_dx - skew_n2 * dn1hat_dx).transpose() * e_hat;

        const Vec3 dcos_dx = dn1hat_dx.transpose() * n2_hat + dn2hat_dx.transpose() * n1_hat;

        return dsin_dx * cos_theta - dcos_dx * sin_theta;
    }

    inline Vec3 SolveSPDOrRegularize(Mat3& H, const Vec3& f) {
        using Scalar = Mat3::Scalar;
        H = static_cast<Scalar>(0.5) * (H + H.transpose());

        Eigen::LLT<Mat3> llt;

        Scalar tau = static_cast<Scalar>(1e-6) * H.diagonal().cwiseAbs().maxCoeff();
        if (!(tau > static_cast<Scalar>(0))) tau = static_cast<Scalar>(1e-6);

        for (int k = 0; k < 10; ++k) {
            Mat3 Hreg = H;
            Hreg.diagonal().array() += tau;

            llt.compute(Hreg);
            if (llt.info() == Eigen::Success)
                return llt.solve(f);

            tau *= static_cast<Scalar>(10);
        }

        return Vec3::Zero();
    }

}



#endif //TINYVBD_MATH_HPP