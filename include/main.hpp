#ifndef MAIN_HPP
#define MAIN_HPP

#include "raisim/World.hpp"
#include "raisim/RaisimServer.hpp"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include "Eigen/Dense"
#include <allegro5/allegro.h>
#include <qpOASES.hpp>
#include <allegro5/allegro_native_dialog.h>
#include "trajectory.hpp"

Eigen::Matrix3d Itorso;  

/* Vectors, matricesand variable declaration */
double t = 0.0, dt = 0.0;

Eigen::VectorXd initialConditions(19), genCoordinates(19);
Eigen::VectorXd F(18), genVelocity(18), genAcceleration(18), genAccelerationVec(18), prevgenVelocity(18);
Eigen::Vector3d Pf_LF, Pf_RF, Pf_LB, Pf_RB; 
Eigen::Vector3d Pcom, torsoRot, dtorsoRot, prev_torsoRot, ddtorsoRot, prev_dtorsoRot, imuRot, prev_imuRot, dimuRot;
Eigen::Vector3d torsoRot2, dtorsoRot2, prev_torsoRot2, ddtorsoRot2, prev_dtorsoRot2;

Eigen::Vector3d Rf_LF, Rf_RF, Rf_LB, Rf_RB;
Eigen::Vector3d Q_LF, Q_RF, Q_LB, Q_RB, prevQ_LF, prevQ_RF, prevQ_LB, prevQ_RB, dQ_LF, dQ_RF, dQ_LB, dQ_RB;
Eigen::Vector3d Fcon_LF, Fcon_RF, Fcon_LB, Fcon_RB, Fcon, q_LF, q_RF, q_LB, q_RB, dq_LF, dq_RF, dq_LB, dq_RB;
Eigen::Vector3d Pcon_LF, Pcon_RF, Pcon_LB, Pcon_RB;
Eigen::Vector3d Rmoment_LF, Rmoment_RF, Rmoment_LB, Rmoment_RB;
Eigen::VectorXd Pcap(8);

double Xi, Yi, m_L, m_R, b_L, b_R;

/* VMC parameters */
Eigen::Matrix <double, 3, 3> Rtorso;
Eigen::Matrix <double, 4, 3> Fmatrix;
Eigen::Vector3d Rc_h1, Rc_h2, Rc_h3, Rc_h4, F1cont, F2cont, F3cont, F4cont, Tau1_vmc, Tau2_vmc, Tau3_vmc, Tau4_vmc, Rcom;
Eigen::Vector3d Pbase, dPbase, des_Pbase, des_dPbase, prev_Pbase;
double Zcom_act = 0, prevZcom_act = 0, dZcom_act = 0;
Eigen::VectorXd Fvmc(6), Mvmc(3);
Eigen::Vector4d quat_ref;
Eigen::VectorXd orientControlRef(10);

/* Push Recovery parameters */
int capEnabled = 0;
double capTstart = 0.0;
double Xcapw = 0.0, Ycapw = 0.0;
Eigen::Vector2d accBody, prev_accBody, accBodyF, prev_accBodyF, velBody, prev_velBody;
double yaw = 0;

/* Inverse dynamics parameters */
Eigen::Matrix3d rootOrientation;
Eigen::Vector3d rootAbsposition, rootAbsvelocity, rootAbsacceleration, rootAngvelocity, rootAngacceleration, prevrootAngacceleration;
Eigen::VectorXd jffTorques(12), jffTorques2(12), jPositions(12), jVelocities(12), jAccelerations(12), Uctc(12);
Eigen::Vector3d ddQ_LF, ddQ_RF, ddQ_LB, ddQ_RB;
Eigen::Vector3d pre_dQ_LF, pre_dQ_RF, pre_dQ_LB, pre_dQ_RB;
Eigen::Vector3d Fref_LF, Fref_RF, Fref_LB, Fref_RB;

Eigen::Vector3d Tau1_JF, Tau2_JF, Tau3_JF, Tau4_JF;
Eigen::VectorXd JF(18), jointTorques(18);
/* Inverse dynamics parameters */

/* Joint position control gains */
Eigen::Vector3d Kp, Kd;

/* Initial Conditions */
// Note: GenCoord => (ComX, ComY, ComZ, W_quaternion, X_quaternion, Y_quaternion, Z_quaternion, Hip_AA_LF, Hip_FE_LF, Knee_FE_LF, Hip_AA_RF, Hip_FE_RF, Knee_FE_RF, Hip_AA_LB, Hip_FE_LB, Knee_FE_LB, Hip_AA_RB, Hip_FE_RB, Knee_FE_RB)
double initComX = 0, initComY = 0, initComZ = 0.53319176863337994221048177223565;
double initComRoll = 0*PI/180, initComPitch = 0*PI/180, initComYaw = 0*PI/180;
double initQuat_w = 1, initQuat_x = 0, initQuat_y = 0, initQuat_z = 0;
Eigen::Vector3d Tau_LF, Tau_RF, Tau_LB, Tau_RB;

/* Controller (Joystick) variables */
double cmdJoyF[22] = {0.0, 0.0, 0.0, 0.0, 0.0, initComZ, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double pre_cmdJoyF[22] = {0.0, 0.0, 0.0, 0.0, 0.0, initComZ, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double cmd_Vx = 0.0, cmd_Vy = 0.0;
double cmd_dyaw = 0.0;
double cmd_yaw = 0.0, cmd_pitch = 0.0, cmd_roll = 0.0;

/* Trajectory Output variables */
double ddZcom = 0.0;

/* Centrodial Momentum */
double Zfoot_offset = 0.0;
int Fcoef = 0;
Eigen::Vector3d Qcm_RF, Qcm_LB;
Eigen::VectorXd dQcm(6), prevdQcm(6), Qcm(6), prevQcm(6), Qcm_filtered(6), prevQcm_filtered(6);

/* Orientation Control Parameters */
double Zroll_left, Zroll_right, Zpitch_front, Zpitch_back;
double Zroll_LF, Zroll_RF, Zroll_LB, Zroll_RB;
double Zpitch_LF, Zpitch_RF, Zpitch_LB, Zpitch_RB;
double Kp_roll = 1.0, Kd_roll = 0.0;
double Kp_pitch = 1.0, Kd_pitch = 0.0;

#endif