#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "Eigen/Dense"
#include "unsupported/Eigen/MatrixFunctions"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cstdlib>

#define PI 3.14159265359
#define GRAVITY 9.81
#define MASS 36.100

bool AreDoubleSame(double dFirstVal, double dSecondVal);
double mapVal(double inVal, double inMax, double inMin, double outMax, double outMin);
int signFuncCap(double input, double thresHold);
double numIntegral(double In, double prevIn, double prevOut, double dt);
bool FuncEqual(double RealTime, double ta, double dt);
bool FuncGreater(double RealTime, double ta, double dt);
bool FuncInterval(double RealTime, double ta, double ts, double dt);
bool FuncLessOrEq(double RealTime, double ta, double dt);
double FuncAdmitFilt(double In, double PrevOut, double k, double b, double dt);
double LPF(double Input, double prevOut, double freq, double dt);
double HPF(double In, double prevIn, double prevOut, double g, double dt);
Eigen::VectorXd QuaternionRef(double Roll, double Pitch, double Yaw, double dRoll, double dPitch, double dYaw, double ddRoll, double ddPitch, double ddYaw);
Eigen::Matrix3d vec2SkewSym(Eigen::Vector3d V);
Eigen::Vector4d eulerToQuaternion(double yaw, double pitch, double roll);
Eigen::Matrix3d RotateRoll(double a);
Eigen::Matrix3d RotatePitch(double a);
Eigen::Matrix3d RotateYaw(double a);
Eigen::Vector3d quat2euler(double w, double x, double y, double z);
double Numdiff(double currX, double prevX, double dt);
Eigen::Vector2d FuncPoly3rd(double RealTime, double t_start, double t_end, double Z0, double dZ0, double Ze, double dZe);
Eigen::Vector2d FuncPoly4th(double RealTime, double t_start, double t_end, double Z0, double dZ0, double Ze, double dZe, double Fh);
Eigen::Vector3d FuncPoly5th(double RealTime, double t_start, double t_end, double Z0, double dZ0, double ddZ0, double Ze, double dZe, double ddZe);
Eigen::Vector3d FuncPoly6th(double RealTime, double t_start, double t_end, double Z0, double dZ0, double ddZ0, double Ze, double dZe, double ddZe, double Fh);
Eigen::Matrix3d JacTranspose_LF(double q0, double q1, double q2, double q3);
Eigen::Matrix3d JacTranspose_RF(double q0, double q1, double q2, double q3);
Eigen::Matrix3d JacTranspose_LB(double q0, double q1, double q2, double q3);
Eigen::Matrix3d JacTranspose_RB(double q0, double q1, double q2, double q3);
Eigen::Matrix<double, 4, 3> refForceCalc4(const Eigen::Vector3d& Rcf1, const Eigen::Vector3d& Rcf2, const Eigen::Vector3d& Rcf3, const Eigen::Vector3d& Rcf4, const Eigen::Vector3d& qlf, const Eigen::Vector3d& qrf, const Eigen::Vector3d& qlh, const Eigen::Vector3d& qrh, double dt);
Eigen::Matrix<double, 4, 3> VMC(const Eigen::Vector3d& Rcf1, const Eigen::Vector3d& Rcf2, const Eigen::Vector3d& Rcf3, const Eigen::Vector3d& Rcf4, const Eigen::Vector3d& Fc_LF, const Eigen::Vector3d& Fc_RF, const Eigen::Vector3d& Fc_LB, const Eigen::Vector3d& Fc_RB, Eigen::VectorXd C, double dt);
double SolveQuadCosSin(double ann, double bnn, double cnn, int mu);
Eigen::Vector3d fullBodyFK(Eigen::Vector3d torsoOrient, Eigen::Vector3d Rcom, Eigen::Vector3d Q, int u);
Eigen::Vector3d fullBodyFK2(Eigen::Matrix3d torsoOrient, Eigen::Vector3d Rcom, Eigen::Vector3d Q, int u);
Eigen::Vector3d fullBodyIKan(Eigen::Vector3d Rfoot, Eigen::Vector3d Rcom, Eigen::Vector3d torsoOrient, int u);
Eigen::Vector3d balanceControl(Eigen::Matrix3d Rd, Eigen::Matrix3d R);
/* Inverse Dynamics */
Eigen::Matrix3d quat2Rotmat(double qw, double qx, double qy, double qz);
Eigen::Vector3d funcNewtonEuler(Eigen::Vector3d rootAbsAcc, Eigen::Matrix3d rootOrient, Eigen::Vector3d rootAngVel, Eigen::Vector3d rootAngAcc, Eigen::Vector3d jPos, Eigen::Vector3d jVel, Eigen::Vector3d jAcc, Eigen::Vector3d grForce, Eigen::Vector3d Rcon, int legIndex);
Eigen::Vector3d funcNewtonEuler2(Eigen::Vector3d rootAbsAcc, Eigen::Matrix3d rootOrient, Eigen::Vector3d rootAngVel, Eigen::Vector3d rootAngAcc, Eigen::Vector3d jPos, Eigen::Vector3d jVel, Eigen::Vector3d jAcc, Eigen::Vector3d grForce, Eigen::Vector3d Rcon, int legIndex);
Eigen::Matrix<double, 12, 1> funNewtonEuler4Leg(Eigen::Vector3d rootAbsAcc, Eigen::Matrix3d rootOrient, Eigen::Vector3d rootAngVel, Eigen::Vector3d rootAngAcc, Eigen::VectorXd jPos, Eigen::VectorXd jVel, Eigen::VectorXd jAcc, Eigen::Vector3d grForce_LF, Eigen::Vector3d grForce_RF, Eigen::Vector3d grForce_LB, Eigen::Vector3d grForce_RB);

#endif