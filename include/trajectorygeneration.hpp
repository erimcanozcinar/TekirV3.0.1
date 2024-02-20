#ifndef TRAJECTORYGENERATION_HPP
#define TRAJECTORYGENERATION_HPP

#include <Eigen/Dense>
#include <allegro5/allegro.h>
#include <allegro5/allegro_native_dialog.h>
#include "functions.hpp"
using namespace Eigen;

double FuncPoly_1st(double RealTime, double t_start, double t_end, double y0, double ye, double dt);
Vector3d FuncPoly4th(double RealTime, double t_start, double t_end, double Z0, double dZ0, double Fh, double Ze, double dZe, double dt);
Vector3d FuncPoly5th(double RealTime, double t_start, double t_end, double z01, double v01, double a01, double z02, double v02, double a02, double dt);
Vector3d FuncPoly6th(double RealTime, double t_start, double t_end, double z0, double dz0, double ddz0, double ze, double dze, double ddze, double zh, double dt);
double SolveQuadCosSin(double ann, double bnn, double cnn, int mu);

class controller {
    protected:
        double Zc = 0.53319176863337994221048177223565;
        double Kv = 0.1;
        double MIN_BODY_HEIGHT = Zc*0.6;
        double MAX_BODY_HEIGHT = 0.58;
        double Vx_mean = 0.0, Vy_mean = 0.0;
        double cmdZc = Zc;
        double roll = 0.0, pitch = 0.0, yaw = 0.0;
        double decreaseHeight = 0.0, increaseHeight = 0.0;;

    public:
        double joyCmd[22] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.53319176863337994221048177223565, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};;
        bool close = false;
        bool walkEnable = false;

        void xboxController(ALLEGRO_EVENT_QUEUE *event_queue, ALLEGRO_EVENT event, ALLEGRO_JOYSTICK* joyStick);
        void dualShockController(ALLEGRO_EVENT_QUEUE *event_queue, ALLEGRO_EVENT ev);
};

class CoMPreviewTrajectory{
    protected:
        double prev_vx_mean = 0.0, prev_vy_mean = 0.0;
        bool walk_enabled = false;
        bool can_switch = false;
        bool can_stop = false;
        double w_start_time = 0.0;

        double Zc = 0.53319176863337994221048177223565;
        double Vx_mean = 0.0;					            // Mean velocity of CoM in X-axis
        double Vy_mean = 0.0;					            // Mean velocity of CoM in X-axis
        int Nphase = 0;						                // Total number of phase(Each Ts + Td is one phase)
        // double px = 0.0;                                    // X ZMP
        // double py = 0.0;                                    // Y ZMP
        double Ts = 0.28;						            // Single support phase period
        double Td = 0.12;						            // Double support phase period
        double Fc = 0.15;
        double Xzmp = 0, Yzmp = 0, Kphase;

        double Cx, dCx, ddCx;
        double Cy, dCy, ddCy;

        double Comx, Comy;
        // double Strx = 0.0, Stry = 0.0;

    private:
        //Initialized State Space Model Variables
        Matrix3d A_sspace{};
        Vector3d B_sspace{};
        RowVector3d C_sspace{};
        Matrix4d A_tilde{};
        Vector4d B_tilde{};
        RowVector4d C_tilde{};
        Vector4d I_tilde{};
        Matrix4d Q{};




        double XzmpRef{};

        double px{};
        double px_previewref{},px_preview{};


        double YzmpRef{};

        double py{};
        double py_previewref{},py_preview{};



    public:
        double Gintegral{};
        VectorXd Gdelay{};
        RowVector3d Gstate{};
        Vector3d xstatesX{};
        Vector3d xstatesY{};
        Vector3d FootLx{},FootLy{}, FootLz{};
        Vector3d FootRx{},FootRy{},  FootRz{};
        double px_ref{};
        double py_ref{};
        double Yx{};
        double Yy{};
        double PreviewactionX{};
        double PreviewactionY{};
        double Strx{};
        double Stry{};
        double Xzmperror{};
        double Yzmperror{};
        double Xcop, Ycop;
        Eigen::Vector3d Footz_L, Footz_R, Footx_L, Footx_R, Footy_L, Footy_R;
        double Xc = 0.0, Yc = 0.0;
        double dXc = 0.0, dYc = 0.0;
        double ddXc = 0.0, ddYc = 0.0;
        double height = Zc;
        Eigen::Vector3d Pfoot_R, Pfoot_L;

        void StateSpace(double comPosZInit, double DT);
        void PreviewControlParameters(double Tpreview,double Rweight,double Qweight, double comPosZInit, double DT);
        void comTrajPreview(double RealTime,double vx_mean,double vy_mean,int Nphase,double Ts,double Td,double Tpreview, double comPosZInit, double DT);
        void CalculateZmpRef(double vx_mean,double vy_mean,double Ts,double Td, double comPosZInit);
        void trajGeneration_preview(double RealTime, bool walkEnable, double command_Vx, double command_Vy, double prev_height, double height, double dt);
};

#endif // TRAJECTORYGENERATION_HPP