#ifndef MATHOPREATIONS_HPP
#define MATHOPREATIONS_HPP
#include <Eigen/Dense>
#include <vector>


using namespace Eigen;


int signFuncCap(double input, double thresHold);
bool FuncEqual(double RealTime, double ta, double dt);
bool FuncGreater(double RealTime, double ta, double dt);
bool FuncLess(double RealTime, double ta, double dt);
bool FuncInterval(double RealTime, double ta, double ts, double dt);
bool FuncEqInterval(double RealTime, double ta, double ts, double dt);
double LPF(double Input, double prevOut, double freq, double dt);


class Integration {
private:
public:
  double dt;
  VectorXd X_prev, X_dprev, X_int;

  Integration(int size, double dt);
  VectorXd TrapezoidalInt(const VectorXd &X);
  VectorXd RungeKuttaInt(const VectorXd &X);
};

class Differentiation {
private:
public:
  double dt;
  VectorXd X_prev, X_dprev;
  MatrixXd X_prevM, X_dprevM;
  Differentiation(int size, double dt);
  Differentiation(int row, int column, double dt);

  VectorXd EulerDiff(const VectorXd &X);
  MatrixXd EulerDiffMatrix(const MatrixXd &X);
};
int signnum(int x) ;
double sigmoid(int x);
Eigen::VectorXd sigmoid(const Eigen::VectorXd& x);
double calculateRMS(const std::vector<double>& data) ;

Matrix3d vec2SkewSym(const Vector3d& vec);

class LPFilter {
private:
  Eigen::VectorXd Out;

public:
  LPFilter(int size);

  VectorXd LPFvec(const Eigen::VectorXd &Vecinput, double freq, double dt);
};

class detectContact{
private:
public:
    //double dt;
    double lastTimeStamp{};
    int statusContact{};

    double sigmoidLikeStep(double x, double threshold, double width);
    int checkContact(double realTime, double Fconz);
};
class LPFilterSca {
private:
    double Out;
public:
    LPFilterSca(double initialOutput );

    double LPF(double input, double freq, double dt);
    };
    double zeroIfNegative(double num);
#endif //MATHOPREATIONS_HPP
