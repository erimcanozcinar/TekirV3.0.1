#include <Eigen/Dense>
#include "mathoperations.hpp"

using namespace Eigen;
Differentiation::Differentiation(int size, double dt) : dt(dt) {
  X_prev = VectorXd::Zero(size);
  X_dprev = VectorXd::Zero(size);
}
Differentiation::Differentiation(int row, int column, double dt) : dt(dt) {
  X_prevM = MatrixXd::Zero(row, column);
  X_dprevM = MatrixXd::Zero(row, column);
}

VectorXd Differentiation::EulerDiff(const VectorXd &X) {
  VectorXd d_X = (X - X_prev) / dt;
  X_prev = X;

  return d_X;
}
MatrixXd Differentiation::EulerDiffMatrix(const MatrixXd &X) {
  MatrixXd d_X = (X - X_prevM) / dt;
  X_prevM = X;

  return d_X;
}
Integration::Integration(int size, double dt) : dt(dt) {
    X_prev = VectorXd::Zero(size);
    X_dprev = VectorXd::Zero(size);
    X_int = VectorXd::Zero(size);
}
VectorXd Integration::TrapezoidalInt(const VectorXd &X) {

    VectorXd X_sum = X_int + 0.5 * (X + X_prev) * dt;
    X_int = X_sum;
    X_prev = X;

    return X_sum;
}

VectorXd Integration::RungeKuttaInt(const VectorXd &X) {
    VectorXd k1 = X_dprev;
    VectorXd k2 = X;
    VectorXd k3 = X;
    VectorXd k4 = X;

    VectorXd X_sum = X_int + (k1 + 2 * k2 + 2 * k3 + k4) * (dt / 6.0);
    X_int = X_sum;
    X_dprev = X;

    return X_sum;
}
int signFuncCap(double input, double thresHold)
{
    int out;
    if (input > thresHold)
    {
        out = 1;
    }
    else if (input < -thresHold)
    {
        out = -1;
    }
    else
    {
        out = 0;
    }
    return out;
}



bool FuncEqual(double RealTime, double ta, double dt) {
    return (RealTime > (ta - dt * 0.25)) && (RealTime < (ta + dt * 0.25));
}

bool FuncGreater(double RealTime, double ta, double dt) {
    return RealTime > (ta + dt * 0.25);
}


bool FuncLess(double RealTime, double ta, double dt) {
    return RealTime < (ta - dt * 0.25);
}


bool FuncInterval(double RealTime, double ta, double ts, double dt) {
    return (RealTime > (ta - dt * 0.25)) && (RealTime < (ts - dt * 0.25));
}

bool FuncEqInterval(double RealTime, double ta, double ts, double dt) {
    return (RealTime > (ta - dt * 0.25)) && (RealTime < (ts + dt * 0.25));
}


double LPF(double Input, double prevOut, double freq, double dt)
{
return 0.0;
}




Matrix3d vec2SkewSym(const Vector3d& vec){

    Matrix3d SkewSym;
    SkewSym.setZero();

    SkewSym(0,1) = -vec(2);
    SkewSym(0,2) = vec(1);

    SkewSym(1,0) = vec(2);
    SkewSym(1,2) = -vec(0);

    SkewSym(2,0) = -vec(1);
    SkewSym(2,1) = vec(0);

    return SkewSym;
}
int signnum(int x) {
  if (x > 0) return 1;
  else return 0;
  return x;
}
double sigmoid(int x) {
    return 1.0 / (1.0 + std::exp(-x));
}
Eigen::VectorXd sigmoid(const Eigen::VectorXd& x) {
    // Apply the sigmoid function element-wise
    return 1.0 / (1.0 + (-x.array()).exp());
}
double calculateRMS(const std::vector<double>& data) {
    double sumOfSquares = 0.0;
    for (const double& x : data) {
        sumOfSquares += x * x;
    }
    double meanOfSquares = sumOfSquares / data.size();
    double rms = std::sqrt(meanOfSquares);
    return rms;
}

LPFilter::LPFilter(int size) : Out(Eigen::VectorXd::Zero(size)){};
Eigen::VectorXd LPFilter::LPFvec(const Eigen::VectorXd &Vecinput, double freq,
                                 double dt) {
  Out = (dt * freq * Vecinput + Out) / (1 + freq * dt);
  return Out;
}


  LPFilterSca::LPFilterSca(double initialOutput) : Out(initialOutput) {};

    double LPFilterSca::LPF(double input, double freq, double dt) {
        Out = (dt * freq * input + Out) / (1 + freq * dt);
        return Out;
    }


double detectContact::sigmoidLikeStep(double x, double threshold,
                                      double width) {
  // Calculate the distance from the threshold
  double distance = std::abs(x - threshold);

  // Apply a sigmoid-like transition using linear interpolation
  if (distance < width / 2.0) {
    // Inside the transition zone, interpolate
    return 0.5 + 0.5 * (distance / (width / 2.0));
  } else if (x < threshold) {
    // Below the threshold, output 0
    return 0.0;
  } else {
    // Above the threshold, output 1
    return 1.0;
  }
}
int detectContact::checkContact(double realTime, double Fconz) {
  double Eps = 5;

  if ((fabs(Fconz) < Eps) && ((realTime - lastTimeStamp) > 0.0005)) {
    statusContact = 0;
    lastTimeStamp = realTime;
  } else if ((fabs(Fconz) >= Eps) && ((realTime - lastTimeStamp) > 0.0005)) {
    statusContact = 1;
    lastTimeStamp = realTime;
  }
  return statusContact;
}
double zeroIfNegative(double num) {
    if (num < 0) {
        return 0;
    } else {
        return num;
    }
}
