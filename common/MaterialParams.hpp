//
// Created by tianyan on 12/30/25.
//

#ifndef TINYVBD_MATERIALPARAMS_HPP
#define TINYVBD_MATERIALPARAMS_HPP

// avoid name conflict with raylib
struct MMaterial {
    float E() const {return _E;};
    float nu() const {return _nu;};
    float lambda() const {return _lambda;};
    float mu() const {return _mu;};

    MMaterial(const float E, const float nu) : _E{E}, _nu{nu} {
        update_lambda_mu();
    }

    void SetE(const float E) {
        _E = E;
        update_lambda_mu();
    }

    void SetNu(const float nu) {
        _nu = nu;
        update_lambda_mu();
    }

private:
    float _E;  // young's module
    float _nu;  // poisson's ratio
    float _lambda;  //first lame parameter
    float _mu;  // second lame parameter

    void update_lambda_mu() {
        _lambda = _E * _nu / (1.0f - _nu * _nu);
        _mu = _E / (2.0f * (1.0f + _nu));
    }

};

inline MMaterial default_cloth() {
    return {1e4f, 0.1f};
};


#endif //TINYVBD_MATERIALPARAMS_HPP