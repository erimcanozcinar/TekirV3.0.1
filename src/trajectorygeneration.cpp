#include <eigen3/Eigen/Dense>
#include "trajectorygeneration.hpp"

using namespace Eigen;

double FuncPoly_1st(double RealTime,double  t_start,double t_end,double y0,double ye,double dt) {

    double y, dy,ddy;
    double v = (ye-y0)/(t_end-t_start); 


    if(FuncLess(RealTime, t_start, dt)==true){
        y = y0;
        dy = 0;}
    else if (FuncInterval(RealTime,t_start,t_end,dt)){
        y = y0 + v*(RealTime-t_start);
        dy = v;
        }

    else{
        y = ye;
        dy = 0;
    }

    ddy = 0;

    return y;
}

Vector3d FuncPoly4th(double RealTime, double t_start, double t_end, double Z0, double dZ0, double Fh, double Ze, double dZe, double dt){

    double tw = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel, Acc;
    double tw2 = tw*tw;
    double tw3 = tw2*tw;
    double tw4 = tw2*tw2;
    tw = t_end - t_start;

    Vector3d trajOut;

    double n0 = Z0;
    double n1 = dZ0;
    double n2 = -(11*Z0 - 16*Fh + 5*Ze + 4*dZ0*tw - dZe*tw)/tw2;
    double n3 = (18*Z0 - 32*Fh + 14*Ze + 5*dZ0*tw - 3*dZe*tw)/tw3;
    double n4 = -(2*(4*Z0 - 8*Fh + 4*Ze + dZ0*tw - dZe*tw))/tw4;
    
    if (FuncInterval(RealTime, t_start, t_end, dt) == true){
        Pos = n0 + Rt*(n1 + Rt*(n2 + Rt*(n3 + Rt*n4)));
        Vel = n1 + Rt*(2*n2 + Rt*(3*n3 + 4*Rt*n4));
        Acc = 2*n2 + Rt*(6*n3 + 12*Rt*n4);
    }
    else if (FuncLess(RealTime, t_start, dt) == true){
        Pos = Z0;
        Vel = dZ0;
        Acc = 0;
    }
    else{
        Pos = Ze;
        Vel = dZe;
        Acc = 0;
    }
    trajOut << Pos, Vel, Acc;

    return trajOut;
}

Vector3d FuncPoly5th(double RealTime, double t_start, double t_end, double z01, double v01, double a01, double z02, double v02, double a02, double dt){

    double tw = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel, Acc;
    double tw2 = tw*tw;
    double tw3 = tw2*tw;
    double tw4 = tw2*tw2;
    double tw5 = tw3*tw2;


    Vector3d trajOut;

    double p1 = -0.5*(12*z01-12*z02+6*v01*tw+6*v02*tw+a01*tw2-a02*tw2)/tw5;
    double p2 = 0.5*(30*z01-30*z02+16*v01*tw+14*v02*tw+3*a01*tw2-2*a02*tw2)/tw4;
    double p3 = -0.5*(20*z01-20*z02+12*v01*tw+8*v02*tw+3*a01*tw2-a02*tw2)/tw3;
    double p4 = 0.5*a01;
    double p5 = v01;
    double p6 = z01;

    if (FuncInterval(RealTime, t_start, t_end, dt) == true){
    
        Pos = p6 + Rt*(p5 + Rt*(p4 + Rt*(p3 + Rt*(p2 + Rt*p1))));
        Vel = p5 + Rt*(2*p4 + Rt*(3*p3 + Rt*(4*p2 + 5*Rt*p1)));
        Acc = (2*p4 + Rt*(6*p3 + Rt*(12*p2 + 20*Rt*p1)));

    }
    else if (FuncLess(RealTime, t_start, dt) == true){
        Pos = z01;
        Vel = v01;
        Acc = a01;
    }
    else{
        Pos = z02;
        Vel = v02;
        Acc = a02;
    }
    trajOut << Pos, Vel, Acc;
    return trajOut;
}

Vector3d FuncPoly6th(double RealTime, double t_start, double t_end, double z0, double dz0, double ddz0, double ze, double dze, double ddze, double zh, double dt){

    double Ts = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel, Acc;

    double Ts2 = Ts*Ts;
    double Ts3 = Ts2*Ts;
    double Ts4 = Ts2*Ts2;
    double Ts5 = Ts3*Ts2;
    double Ts6 = Ts3*Ts3;

    Vector3d trajOut;

    double n0 = z0;
    double n1 = dz0;
    double n2 =  ddz0/2;
    double n3 =  -(((5*ddz0)/2 + ddze/2)*Ts2 + (16*dz0 - 6*dze)*Ts + 42*z0 
    + 22*ze - 64*zh)/Ts3;
    double n4 = (((9*ddz0)/2 + 2*ddze)*Ts2 + (38*dz0 - 23*dze)*Ts + 111*z0 
    + 81*ze - 192*zh)/Ts4;
    double n5 = -(((7*ddz0)/2 + (5*ddze)/2)*Ts2 + (33*dz0 - 27*dze)*Ts 
    + 102*z0 + 90*ze - 192*zh)/Ts5;
    double n6 = ((ddz0 + ddze)*Ts2 + (10*dz0 - 10*dze)*Ts + 
    32*z0 + 32*ze - 64*zh)/Ts6;
    
    if (FuncInterval(RealTime, t_start, t_end, dt) == true){
    
        Pos = n0 + Rt*(n1 + Rt*(n2 + Rt*(n3 + Rt*(n4 + Rt*(n5 + Rt*n6)))));
        Vel = n1 + Rt*(2*n2 + Rt*(3*n3 + Rt*(4*n4 + Rt*(5*n5 + Rt*6*n6))));
        Acc = 2*n2 + Rt*(6*n3 + Rt*(12*n4 + Rt*(20*n5 + Rt*30*n6)));
    
    }

    else if (FuncLess(RealTime, t_start, dt) == true){
        
        Pos = z0;
        Vel = dz0;
        Acc = ddz0;
    }

    else{
        Pos = ze;
        Vel = dze;
        Acc = ddze;
    }
    trajOut << Pos, Vel, Acc;
    return trajOut;
}

//State Space Model for ZMP-LIPM with Triangular Discritizaion
void CoMPreviewTrajectory::StateSpace(double comPosZInit, double DT){

    // State Space Model
    A_sspace << 1, DT, pow(DT, 2) / 2, 0, 1, DT, 0, 0, 1;
    B_sspace << pow(DT, 3) / 6, pow(DT, 2) / 2, DT;
    C_sspace << 1, 0, -comPosZInit / GRAVITY;
    A_tilde << 1, C_sspace * A_sspace, MatrixXd::Zero(3, 1), A_sspace;
    B_tilde << C_sspace * B_sspace, B_sspace;
    C_tilde << 1, 0, 0, 0;
    I_tilde << 1, 0, 0, 0;
    //Augmented State Space Model for Preview Control

}

std::tuple < MatrixXd,MatrixXd>  dlqr_H(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R) {
  // Solve the discrete-time algebraic Riccati equation
  Eigen::MatrixXd P = Q;
  int iter = 0;
  while (iter < 200000) {
    Eigen::MatrixXd P_next = A.transpose() * P * A - A.transpose() * P * B *
      (R + B.transpose() * P * B).inverse() * B.transpose() * P * A + Q;
    if ((P_next - P).norm() < 1e-8) {
      break;
    }
    P = P_next;
    iter++;
  }

  // Compute the optimal feedback gain
  MatrixXd Kgain = (R + B.transpose() * P * B).inverse() * B.transpose() * P * A;
  return {Kgain,P};
}

void CoMPreviewTrajectory::CalculateZmpRef(double vx_mean,double vy_mean,double Ts,double Td, double comPosZInit){

  double  w = sqrt(GRAVITY / comPosZInit);
  double  Cxd = vx_mean * Ts / 2;
  double  Cx0 = (-Cxd);
  double  dCx0 = w * (-Cx0) * 1 / tanh(w * Ts / 2);
  double  Kx = dCx0 + w * (Cxd)*1 / tanh(w * Td / 2);

  double  Cyd = vy_mean * Ts / 2;
  double  Cy0 = (-Cyd);
  double  dCy0 = w * (-Cy0) * 1 / tanh(w * Ts / 2);
  double  Ky = dCy0 + w * (Cyd)*1 / tanh(w * Td / 2);



  Strx = 2 * Td * Kx;
  Stry = 2 * Td * Ky;

  XzmpRef=Strx/2;
  YzmpRef=Stry/2;
}

void CoMPreviewTrajectory::PreviewControlParameters( double Tpreview,double Rweight,double Qweight, double comPosZInit, double DT)
{
  StateSpace(comPosZInit, DT);
  // State Space Model
  Q = Eigen::DiagonalMatrix<double, 4>(Qweight, 0, 0, 0);
  VectorXd R(1);
  R <<Rweight;
  VectorXd Temp(1);
  Temp.setZero();
  MatrixXd Smatrix(4, 4);
  Smatrix.setZero();
  MatrixXd Kgain(1, 4);
  Kgain.setZero();
  std::tie(Kgain, Smatrix) = dlqr_H(A_tilde, B_tilde, Q, R);
  Gintegral = Kgain(0, 0);
  Gstate = Kgain.block<1, 3>(0, 1);
  int Nsteps = Tpreview / DT + 1;
  Gdelay.resize(Nsteps);
  Gdelay.setZero();
  Gdelay(0)=-Gintegral;
  MatrixXd Ac_tilde = A_tilde - B_tilde * Kgain;
  VectorXd X_tilde = -Ac_tilde.transpose() * Smatrix * I_tilde;
  for (int i = 1; i < Nsteps; i++) {
    Temp = ((R + B_tilde.transpose() * Smatrix * B_tilde).inverse()) *
          B_tilde.transpose() * X_tilde;
    Gdelay(i) = Temp(0);
    X_tilde = Ac_tilde.transpose() * X_tilde;
    Ac_tilde = A_tilde -
              B_tilde * (R + B_tilde.transpose() * Smatrix * B_tilde).inverse() *
                  B_tilde.transpose() * Smatrix * A_tilde;
  }
} 

void CoMPreviewTrajectory::comTrajPreview(double RealTime,double vx_mean,double vy_mean,int Nphase,double Ts,double Td,double Tpreview, double comPosZInit, double DT)
{

    double t0 = 0; // Beginning of the Universe(or Simulation)
    double t1 = t0 ; // Put robot on the ground
    double t2 = t1 + Ts; // First half step
    double tstart = t2 + Td; // Walking start
    double tend = tstart + (Ts + Td) * Nphase; // Walking end
    double t3 = tend+ Ts; // CoM stoping phase
    double Fh=0.15;
    double xzmpOffset=0.00;
    double yzmpOffset = 0.00;
    // PreviewControlParameters(TPreview,RWeight,QWeight);
    double rt4Fr = (RealTime-tstart) / (Ts + Td);
    double rt4Frcom = (RealTime-t0) / (Ts + Td);
    rt4Fr=zeroIfNegative(rt4Fr);
    rt4Frcom=zeroIfNegative(rt4Frcom);
    int k,kx,kcom;
    if (FuncInterval(RealTime, 0, t0, DT) == true) {
      k = 0;
      kx = 0;
      kcom = 0;
    } else if (FuncInterval(RealTime, t0, tstart, DT) == true) {
      k = 0;
      kcom = floor(round(rt4Frcom / DT) * DT);
      kx = 0;
    } else if (FuncGreater(RealTime, tend, DT) == true) {
      k = Nphase - 1;
      kcom = Nphase;
      kx = floor(k / (2 * DT)) * DT;
    } else {
      k = floor(round(rt4Fr / DT) * DT);
      kcom = floor(round(rt4Frcom / DT) * DT);
      kx = floor(round(rt4Fr / (2 * DT)) * DT);
    }

    double Rt = RealTime - (Ts + Td) * (kcom)-t0;
    CalculateZmpRef(vx_mean,vy_mean, Ts, Td, comPosZInit);
    // Left Foot Trajectory Generation

    if (FuncInterval(RealTime, 0, t1, DT) == true) {

      FootLx.setZero();
      FootLy.setZero();
      FootLz.setZero();
      FootRx.setZero();
      FootRz.setZero();
    } else if (FuncInterval(RealTime, t1, t2, DT) == true) {

      FootLx = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, XzmpRef, 0, 0, DT);
      FootLy = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, YzmpRef, 0, 0, DT);
      FootLz = FuncPoly6th(RealTime, t1, t2, 0, 0, 0, 0, 0, 0, Fh, DT);
      FootRx.setZero();
      FootRz.setZero();

    }
    else if (FuncInterval(RealTime, t2, tstart, DT) == true) {

      FootLx << XzmpRef, 0, 0;
      FootLy << YzmpRef, 0, 0;
      FootLz.setZero();
      FootRx.setZero();
      FootRz.setZero();
    } else if (FuncInterval(RealTime, tstart, tend, DT) == true) {

      if (((k + 1) % 2) == 0) {
        double ts_l = tstart + (Ts + Td) * k;
        double ts_r = tstart + (Ts + Td) * k;
        FootLx = FuncPoly5th(RealTime, ts_l, ts_l + Ts, XzmpRef + Strx * kx, 0,
                             0, XzmpRef + Strx * (kx + 1), 0, 0, DT);
        FootLy = FuncPoly5th(RealTime, ts_l, ts_l + Ts, YzmpRef + Stry * kx, 0,
                             0, YzmpRef + Stry * (kx + 1), 0, 0, DT);
        FootLz =
            FuncPoly6th(RealTime, ts_l, ts_l + Ts, 0, 0, 0, 0, 0, 0, Fh, DT);

        FootRx = FuncPoly5th(RealTime, ts_l, ts_r + Ts, Strx * (kx + 1), 0, 0,
                             Strx * (kx + 1), 0, 0, DT);
        FootRy = FuncPoly5th(RealTime, ts_l, ts_r + Ts, Stry * (kx + 1), 0, 0,
                             Stry * (kx + 1), 0, 0, DT);
        FootRz =
            FuncPoly6th(RealTime, ts_r, ts_r + Ts, 0, 0, 0, 0, 0, 0, 0, DT);
      }
      else {
        double ts_l = tstart + (Ts + Td) * k;
        double ts_r = tstart + (Ts + Td) * k;
        FootLx = FuncPoly5th(RealTime, ts_l, ts_l + Ts, XzmpRef + Strx * kx, 0,
                             0, XzmpRef + Strx * kx, 0, 0, DT);
        FootLy = FuncPoly5th(RealTime, ts_l, ts_l + Ts, YzmpRef + Stry * kx, 0,
                             0, YzmpRef + Stry * kx, 0, 0, DT);
        FootLz =
            FuncPoly6th(RealTime, ts_l, ts_l + Ts, 0, 0, 0, 0, 0, 0, 0, DT);

        FootRx = FuncPoly5th(RealTime, ts_r, ts_r + Ts, Strx * kx, 0, 0,
                             Strx * (kx + 1), 0, 0, DT);
        FootRy = FuncPoly5th(RealTime, ts_r, ts_r + Ts, Stry * kx, 0, 0,
                             Stry * (kx + 1), 0, 0, DT);
        FootRz =
            FuncPoly6th(RealTime, ts_r, ts_r + Ts, 0, 0, 0, 0, 0, 0, Fh, DT);

      }
    } else if (FuncInterval(RealTime, tend, t3, DT) == true) {
      if ((Nphase % 2) == 0) {
        FootLx << Strx / 2 + Strx * (kx + 1), 0, 0;
        FootLy << Stry / 2 + Stry * (kx + 1), 0, 0;
        FootLz.setZero();
        FootRx = FuncPoly5th(RealTime, tend, t3, Strx * (kx + 1), 0, 0,
                             Strx / 2 + Strx * (kx + 1), 0, 0, DT);
        FootRy = FuncPoly5th(RealTime, tend, t3, Stry * (kx + 1), 0, 0,
                             Stry / 2 + Stry * (kx + 1), 0, 0, DT);                             
        FootRz =
            FuncPoly6th(RealTime, tend, t3, 0, 0, 0, 0, 0, 0, Fh, DT);
      }
      else {

        FootLx = FuncPoly5th(RealTime, tend, t3, Strx / 2 + Strx * (kx), 0, 0,
                             Strx + Strx * (kx), 0, 0, DT);
        FootLy = FuncPoly5th(RealTime, tend, t3, Stry / 2 + Stry * (kx), 0, 0,
                             Stry + Stry * (kx), 0, 0, DT);
        FootLz = FuncPoly6th(RealTime, tend, t3, 0, 0, 0, 0, 0, 0, Fh, DT);

        FootRx << Strx * (kx + 1), 0, 0;
        FootRy << Stry * (kx + 1), 0, 0;
        FootRz.setZero();
      }
    } 
    else{

      if ((Nphase % 2) == 0) {

        FootLx << Strx / 2 + Strx * (kx + 1), 0, 0;
        FootLy << Stry / 2 + Stry * (kx + 1), 0, 0;
        FootLz =
            FuncPoly6th(RealTime, tend, t3, 0, 0, 0, 0, 0, 0, Fh, DT);

        FootRx << Strx / 2 + Strx * (kx + 1), 0, 0;
        FootRy << Stry / 2 + Stry * (kx + 1), 0, 0;
        FootRz.setZero();

      }

      else {
        FootLx << Strx * (kx + 1), 0, 0;
        FootLy << Stry * (kx + 1), 0, 0;
        FootLz.setZero();

        FootRx << Strx * (kx + 1), 0, 0;
        FootRy << Stry * (kx + 1), 0, 0;
        FootRz.setZero();
      }
    }
    if (FuncInterval(Rt, 0, Ts, DT) == true) {
      px_ref = XzmpRef*kcom;
      py_ref = YzmpRef*kcom;


    }

    else if (FuncInterval(Rt, Ts, Ts + Td, DT) == true) {
      px_ref = FuncPoly_1st(Rt, Ts, Ts + Td, XzmpRef*kcom, xzmpOffset+XzmpRef * (kcom + 1), DT);
      py_ref = FuncPoly_1st(Rt, Ts, Ts + Td, YzmpRef*kcom, yzmpOffset+YzmpRef * (kcom + 1), DT);

    } 
    // else if (FuncEqual(Rt, Ts + Td - DT, DT) == true) {

    //   px = px_ref;
    //   py = py_ref;

    // }

    PreviewactionX = 0;
    PreviewactionY=0;
    px_previewref = px_ref;
    px_preview = px;
    py_previewref = py_ref;
    py_preview = py;
    for (int i = 0; i < Tpreview / DT ; i++) {



      double Rtpreview = DT * (i) + RealTime;
      int Rtfloor= round(Rtpreview/DT);
    //   if (Rtpreview >= 0 && Rtpreview <= 30) {
    //     if (Rtpreview >= 5 && Rtpreview <= 10) {
    //    vx_mean=FuncStep(Rtpreview, 5, 0, 0.10);
    //     } else if (Rtpreview >= 10 && Rtpreview <= 15) {
    //    vx_mean=FuncStep(Rtpreview, 10, 0, 0.11);
    //     } else if (Rtpreview >= 15 && Rtpreview <= 20) {
    //        vx_mean=FuncStep(Rtpreview, 4, 0, 0.12);
    //     } el
    //       vxmean = 0.12;
    //     }
    //   }
    //   vx_mean = FuncPoly_1st(Rtpreview, 3, 10, 0, 0.6, DT);

    //   vx_mean = FuncPoly_1st(Rtpreview, 3, 3.001, 0, 0.4, DT);
    // vx_mean = FuncPoly_1st(Rtpreview, 8, 8.001, 0.4, 0.5, DT);
//  vx_mean = FuncPoly_1st(Rtpreview, 3, 10, 0, 0.4, DT);

      CalculateZmpRef(vx_mean,vy_mean, Ts, Td, comPosZInit);
      double Rt4Frpreview = (Rtpreview - t0) / (Ts + Td);
      Rt4Frpreview = zeroIfNegative(Rt4Frpreview);
      double kpreview = floor(round(Rt4Frpreview / DT) * DT);

    
        if (FuncInterval(Rtpreview, t0, tstart, DT) == true){
            kpreview = floor(round(Rt4Frpreview / DT) *DT);
        }
        else if (FuncGreater(Rtpreview, tend, DT) == true){

            kpreview = Nphase;
        }
        else{

            kpreview = floor(round(Rt4Frpreview / DT) *DT);
        }
      Rt = Rtpreview - (Ts + Td) * (kpreview)-t0;

      if (FuncInterval(Rt, 0, Ts, DT) == true){
        px_previewref = XzmpRef*kpreview;
        py_previewref = YzmpRef*kpreview;

      }
      else if(FuncInterval(Rt, Ts, Ts + Td, DT) == true) {
      px_previewref =FuncPoly_1st(Rt, Ts, Ts + Td, XzmpRef*kpreview, xzmpOffset+(XzmpRef) * (kpreview + 1),DT);
      py_previewref = FuncPoly_1st(Rt, Ts, Ts + Td, YzmpRef*kpreview,
                                   yzmpOffset+YzmpRef * (kpreview + 1), DT);

      }

    //   else if(FuncEqual(Rt, Ts + Td - DT, DT) == true){

    //     px_previewref = px_previewref;
    //     py_previewref = py_previewref;
    // }

          PreviewactionX = PreviewactionX + Gdelay(i) * px_previewref;
          PreviewactionY = PreviewactionY + Gdelay(i) * py_previewref;

    }


    double Ux =  +Gintegral*Xzmperror- Gstate * xstatesX - PreviewactionX;
    double Uy = -Gintegral * Yzmperror - Gstate * xstatesY - PreviewactionY;

    Yx = C_sspace * xstatesX;
    Yy = C_sspace * xstatesY;
    Xzmperror += px_ref - Yx;

    Yzmperror = py_ref - Yy;
    xstatesX = A_sspace * xstatesX + B_sspace * Ux;
    xstatesY = A_sspace * xstatesY + B_sspace * Uy;
}

void CoMPreviewTrajectory::trajGeneration_preview(double RealTime, bool walkEnable, double command_Vx, double command_Vy, double prev_height, double height, double dt)
{
    if(RealTime < dt) { PreviewControlParameters(2, 10e-6, 1, height, dt); }
    else if(!AreDoubleSame(prev_height, height)){ PreviewControlParameters(2, 10e-6, 1, height, dt); }

    // std::cout << traj.Gdelay << std::endl;

    Eigen::VectorXd outVals(13);
    double footPx_R = 0.0, footPy_R = 0, footPz_R = 0.0;
    double footPx_L = 0.0, footPy_L = 0, footPz_L = 0.0;
    

    if(!walk_enabled && walkEnable) w_start_time = RealTime;
    walk_enabled = walkEnable;

    if(walk_enabled){
        can_switch = (AreDoubleSame(FootRz(0),Fc)) || (AreDoubleSame(FootLz(0),Fc));
    } else {
        can_switch = false;
    }

    if(abs(Vx_mean) > 0 || abs(Vy_mean) > 0){
        can_stop = false;
    } else {
        can_stop = (AreDoubleSame(FootRz(0),0)) || (AreDoubleSame(FootLz(0),0));
    }
     

    if(can_switch && walk_enabled){
        Vx_mean = command_Vx; 
        Vy_mean = command_Vy; 
    }

    if(prev_vx_mean != Vx_mean) Comx += xstatesX(0);
    if(prev_vy_mean != Vy_mean) Comy += xstatesY(0);

    if((prev_vx_mean != Vx_mean) || (prev_vy_mean != Vy_mean)){
        if(AreDoubleSame(FootRz(0),Fc)){
            w_start_time = RealTime - Ts - Td - Ts/2.0;
            comTrajPreview(RealTime-w_start_time, Vx_mean, Vy_mean, Nphase, Ts, Td, 2, height, dt);
            if(prev_vx_mean != Vx_mean) Comx -= Strx/2;
            if(prev_vy_mean != Vy_mean) Comy -= Stry/2;
        }else if(AreDoubleSame(FootLz(0),Fc)){
            w_start_time = RealTime - Ts - Td - Ts - Td - Ts/2.0 ;
            comTrajPreview(RealTime-w_start_time, Vx_mean, Vy_mean, Nphase, Ts, Td, 2, height, dt);
            if(prev_vx_mean != Vx_mean) Comx -= Strx;
            if(prev_vy_mean != Vy_mean) Comy -= Stry;
        }

        if(prev_vx_mean != Vx_mean) prev_vx_mean = Vx_mean;
        if(prev_vy_mean != Vy_mean) prev_vy_mean = Vy_mean;
    }
    
    if(walk_enabled){
        Nphase += 1;
    }
    else if(!walk_enabled && can_stop){
        Nphase = 0;
    }
    
    comTrajPreview(RealTime-w_start_time, Vx_mean, Vy_mean, Nphase, Ts, Td, 2, height, dt);
    Xc = xstatesX(0)+Comx; Yc = xstatesY(0);
    dXc = xstatesX(1); dYc = xstatesY(1); 
    ddXc = xstatesX(2); ddYc = xstatesY(2);
    Pfoot_R << FootRx(0)+Comx, FootRy(0)+Comy, FootRz(0);
    Pfoot_L << FootLx(0)+Comx, FootLy(0)+Comy, FootLz(0);
}