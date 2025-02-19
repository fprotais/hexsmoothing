#pragma once
#include "lbfgs.hpp"

#include <functional>

class LBFGS_wrapper {
public:
    using FuncEval = std::function<double(Eigen::VectorXd const &x, Eigen::VectorXd &g)>;
    // return true to stop the optimization
    using ProgressCallBack = std::function<bool(
        Eigen::VectorXd const &x, 
        Eigen::VectorXd const &g,
        double f,
        unsigned iter,
        unsigned nbEval
    )>;

    LBFGS_wrapper(const FuncEval f) : _f(f) {}

    FuncEval const _f;
    ProgressCallBack _callBack = [](Eigen::VectorXd const &, Eigen::VectorXd const &, double, unsigned, unsigned) { return false; };
    lbfgs::lbfgs_parameter_t _parameters;

    void setMaxIter(unsigned nbIter) {_parameters.max_iterations = nbIter; }
    
    int getStatus() const {return _status;}
    std::string getMessage() const {
        if (_status == -1) 
            return "Not run.";
        else
            return lbfgs::lbfgs_strerror(_status);
    }

    bool optimize(Eigen::VectorXd &x) {
        _status = lbfgs::lbfgs_optimize(x, _fval, staticFunctionCall, nullptr, staticProgressCall, static_cast<void*>(this), _parameters);
        return _status >= 0;
    }

    double minimum() const { return _fval; }

private: 
    int _status = -1;
    double _fval = 0.;
    static double staticFunctionCall(void *instance, Eigen::VectorXd const &x, Eigen::VectorXd &g) {
        LBFGS_wrapper const * wrapper = static_cast<LBFGS_wrapper const *>(instance);
        return wrapper->_f(x, g);
    }
    static int staticProgressCall(void *instance,
                                const Eigen::VectorXd &x,
                                const Eigen::VectorXd &g,
                                const double fx,
                                const double step,
                                const int k,
                                const int ls)
    {
        LBFGS_wrapper const * wrapper = static_cast<LBFGS_wrapper const *>(instance);
        return wrapper->_callBack(x,g,fx, k, ls);
    }   

};