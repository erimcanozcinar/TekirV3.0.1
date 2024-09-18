#include "functions.hpp"

bool AreDoubleSame(double dFirstVal, double dSecondVal)
{
    return std::fabs(dFirstVal - dSecondVal) <  1E-6;
}

double mapVal(double inVal, double inMax, double inMin, double outMax, double outMin)
{
    double mappedVal = (inVal - inMin)*((outMax-outMin)/(inMax - inMin)) + outMin;
    return mappedVal;
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

double numIntegral(double In, double prevIn, double prevOut, double dt)
{
    double Out = prevOut + (In + prevIn) * dt / 2;
    return Out;
}

bool FuncEqual(double RealTime, double ta, double dt)
{
    double Lim1 = ta - dt * 0.05;
    double Lim2 = ta + dt * 0.05;
    if ((RealTime > Lim1) && (RealTime < Lim2))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool FuncGreater(double RealTime, double ta, double dt)
{
    double Lim1 = ta + dt * 0.25;

    if (RealTime > Lim1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool FuncInterval(double RealTime, double ta, double ts, double dt)
{
    double Lim1 = ta - dt * 0.25;
    double Lim2 = ts + dt * 0.25;

    if ((RealTime > Lim1) && (RealTime < Lim2)) // t = [Lim1 Lim2]
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool FuncLessOrEq(double RealTime, double ta, double dt)
{
    double Lim1 = ta + dt * 0.25;
    if (RealTime < Lim1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

double FuncAdmitFilt(double In, double PrevOut, double k, double b, double dt)
{
    double Bk = b / dt;
    double Out = (In + PrevOut * Bk) / (k + Bk);
    // Out(i,1) = (u(i,1) + Out((i-1),1)*Ba/dt)/den;
    // den = Ka + Ba/dt;
    return Out;
}

double LPF(double Input, double prevOut, double freq, double dt)
{
    double out = (dt * freq * Input + prevOut) / (1 + freq * dt);
    return out;
}

double HPF(double In, double prevIn, double prevOut, double g, double dt)
{
    double Out = (In - prevIn + prevOut) / (1 + g * dt);
    return Out;
}

Eigen::VectorXd QuaternionRef(double Roll, double Pitch, double Yaw, double dRoll, double dPitch, double dYaw, double ddRoll, double ddPitch, double ddYaw) 
{   
    Eigen::VectorXd out(10);
    Eigen::Matrix3d R;
    double Wx = -dPitch*sin(Yaw)+dRoll*cos(Pitch)*cos(Yaw);
    double Wy = dPitch*cos(Yaw)+dRoll*cos(Pitch)*sin(Yaw);
    double Wz = dYaw-dRoll*sin(Pitch);

    double dWx = -dYaw*(dPitch*cos(Yaw)+dRoll*cos(Pitch)*sin(Yaw))-ddPitch*sin(Yaw)+ddRoll*cos(Pitch)*cos(Yaw)-dRoll*dPitch*cos(Yaw)*sin(Pitch);
    double dWy = -dYaw*(dPitch*sin(Yaw)-dRoll*cos(Pitch)*cos(Yaw))+ddPitch*cos(Yaw)+ddRoll*cos(Pitch)*sin(Yaw)-dRoll*dPitch*sin(Pitch)*sin(Yaw);
    double dWz = ddYaw-ddRoll*sin(Pitch)-dRoll*dPitch*cos(Pitch);

    /* Rotation Matrix to Quaternion */
    R(0,0) = cos(Pitch)*cos(Yaw);
    R(0,1) = -cos(Roll)*sin(Yaw)+cos(Yaw)*sin(Roll)*sin(Pitch);
    R(0,2) = sin(Roll)*sin(Yaw)+cos(Roll)*cos(Yaw)*sin(Pitch);
    R(1,0) = cos(Pitch)*sin(Yaw);
    R(1,1) = cos(Roll)*cos(Yaw)+sin(Roll)*sin(Pitch)*sin(Yaw);
    R(1,2) = -cos(Yaw)*sin(Roll)+cos(Roll)*sin(Pitch)*sin(Yaw);
    R(2,0) = -sin(Pitch);
    R(2,1) = cos(Pitch)*sin(Roll);
    R(2,2) = cos(Roll)*cos(Pitch);
    double qw = 0.5 * sqrt(1 + R(0,0) + R(1,1) + R(2,2));
    double qx = (R(2,1) - R(1,2)) / (4 * qw);
    double qy = (R(0,2) - R(2,0)) / (4 * qw);
    double qz = (R(1,0) - R(0,1)) / (4 * qw);

    out << qw,qx,qy,qz,Wx,Wy,Wz,dWx,dWy,dWz;
    return out;
}

Eigen::Matrix3d vec2SkewSym(Eigen::Vector3d V)
{
    Eigen::Matrix3d M;
    M << 0, -V(2), V(1),
         V(2), 0, -V(0),
         -V(1), V(0), 0;
    return M;
}

Eigen::Vector4d eulerToQuaternion(double yaw, double pitch, double roll) {
    double w, x, y, z;
    Eigen::Vector4d quat;

    double cy = std::cos(yaw * 0.5);
    double sy = std::sin(yaw * 0.5);
    double cp = std::cos(pitch * 0.5);
    double sp = std::sin(pitch * 0.5);
    double cr = std::cos(roll * 0.5);
    double sr = std::sin(roll * 0.5);

    w = cy * cp * cr + sy * sp * sr;
    x = cy * cp * sr - sy * sp * cr;
    y = cy * sp * cr + sy * cp * sr;
    z = sy * cp * cr - cy * sp * sr;

    quat << w, x, y, z;

    return quat;
}

Eigen::Matrix3d RotateRoll(double a)
{
    double RotX[3][3];
    Eigen::Matrix3d rollMatrix;

    RotX[0][0] = 1;
    RotX[0][1] = 0;
    RotX[0][2] = 0;

    RotX[1][0] = 0;
    RotX[1][1] = cos(a);
    RotX[1][2] = -sin(a);

    RotX[2][0] = 0;
    RotX[2][1] = sin(a);
    RotX[2][2] = cos(a);

    rollMatrix << RotX[0][0], RotX[0][1], RotX[0][2], RotX[1][0], RotX[1][1], RotX[1][2], RotX[2][0], RotX[2][1], RotX[2][2];
    return rollMatrix;
}

Eigen::Matrix3d RotatePitch(double a)
{
    double RotY[3][3];
    Eigen::Matrix3d pitchMatrix;

    RotY[0][0] = cos(a);
    RotY[0][1] = 0;
    RotY[0][2] = sin(a);

    RotY[1][0] = 0;
    RotY[1][1] = 1;
    RotY[1][2] = 0;

    RotY[2][0] = -sin(a);
    RotY[2][1] = 0;
    RotY[2][2] = cos(a);

    pitchMatrix << RotY[0][0], RotY[0][1], RotY[0][2], RotY[1][0], RotY[1][1], RotY[1][2], RotY[2][0], RotY[2][1], RotY[2][2];
    return pitchMatrix;
}

Eigen::Matrix3d RotateYaw(double a)
{
    double RotZ[3][3];
    Eigen::Matrix3d yawMatrix;

    RotZ[0][0] = cos(a);
    RotZ[0][1] = -sin(a);
    RotZ[0][2] = 0;

    RotZ[1][0] = sin(a);
    RotZ[1][1] = cos(a);
    RotZ[1][2] = 0;

    RotZ[2][0] = 0;
    RotZ[2][1] = 0;
    RotZ[2][2] = 1;

    yawMatrix << RotZ[0][0], RotZ[0][1], RotZ[0][2], RotZ[1][0], RotZ[1][1], RotZ[1][2], RotZ[2][0], RotZ[2][1], RotZ[2][2];
    return yawMatrix;
}

Eigen::Vector3d quat2euler(double w, double x, double y, double z)
{
    double roll, pitch, yaw;
    Eigen::Vector3d orientation;
    roll = atan2(2 * (w * x + y * z), 1 - 2 * (x * x + y * y));
    // pitch = asin(2 * (w * y - z * x));
    pitch = 2*atan2(sqrt(1+2*(w*y-x*z)),sqrt(1-2*(w*y-x*z))) - M_PI_2;
    yaw = atan2(2 * (w * z + x * y), 1 - 2 * (y * y + z * z));
    orientation << roll, pitch, yaw;
    return orientation;
}

double Numdiff(double currX, double prevX, double dt)
{
    return (currX - prevX) / dt;
}

Eigen::Vector2d FuncPoly3rd(double RealTime, double t_start, double t_end, double Z0, double dZ0, double Ze, double dZe)
{
    double tw = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel;
    Eigen::Vector2d trajOut;

    double n0 = Z0;
    double n1 = dZ0;
    double n2 = (3 * Ze) / pow(tw, 2) - (3 * Z0) / pow(tw, 2) - (2 * dZ0) / tw - dZe / tw;
    double n3 = (2 * Z0) / pow(tw, 3) - (2 * Ze) / pow(tw, 3) + dZ0 / pow(tw, 2) + dZe / pow(tw, 2);

    Pos = n0 + n1 * Rt + n2 * pow(Rt, 2) + n3 * pow(Rt, 3);
    Vel = n1 + 2 * n2 * Rt + 3 * n3 * pow(Rt, 2);

    if (RealTime >= t_start && RealTime <= t_end)
    {
        Pos = n0 + n1 * Rt + n2 * pow(Rt, 2) + n3 * pow(Rt, 3);
        Vel = n1 + 2 * n2 * Rt + 3 * n3 * pow(Rt, 2);
    }
    else if (RealTime < t_start)
    {
        Pos = Z0;
        Vel = dZ0;
    }
    else if (RealTime > t_end)
    {
        Pos = Ze;
        Vel = dZe;
    }
    trajOut << Pos, Vel;
    return trajOut;
}

Eigen::Vector2d FuncPoly4th(double RealTime, double t_start, double t_end, double Z0, double dZ0, double Ze, double dZe, double Fh)
{
    double tw = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel;
    Eigen::Vector2d trajOut;

    double n0 = Z0;
    double n1 = dZ0;
    double n2 = -(11 * Z0 - 16 * Fh + 5 * Ze + 4 * dZ0 * tw - dZe * tw) / pow(tw, 2);
    double n3 = (18 * Z0 - 32 * Fh + 14 * Ze + 5 * dZ0 * tw - 3 * dZe * tw) / pow(tw, 3);
    double n4 = -(2 * (4 * Z0 - 8 * Fh + 4 * Ze + dZ0 * tw - dZe * tw)) / pow(tw, 4);

    Pos = n0 + n1 * Rt + n2 * pow(Rt, 2) + n3 * pow(Rt, 3) + n4 * pow(Rt, 4);
    Vel = n1 + 2 * n2 * Rt + 3 * n3 * pow(Rt, 2) + 4 * n4 * pow(Rt, 3);

    if (RealTime >= t_start && RealTime <= t_end)
    {
        Pos = n0 + n1 * Rt + n2 * pow(Rt, 2) + n3 * pow(Rt, 3) + n4 * pow(Rt, 4);
        Vel = n1 + 2 * n2 * Rt + 3 * n3 * pow(Rt, 2) + 4 * n4 * pow(Rt, 3);
    }
    else if (RealTime < t_start)
    {
        Pos = Z0;
        Vel = dZ0;
    }
    else if (RealTime > t_end)
    {
        Pos = Ze;
        Vel = dZe;
    }
    trajOut << Pos, Vel;
    return trajOut;
}

Eigen::Vector3d FuncPoly5th(double RealTime, double t_start, double t_end, double Z0, double dZ0, double ddZ0, double Ze, double dZe, double ddZe)
{
    double tw = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel, Acc;
    Eigen::Vector3d trajOut;

    double n0 = Z0;
    double n1 = dZ0;
    double n2 = ddZ0 / 2;
    double n3 = -(20 * Z0 - 20 * Ze + 12 * dZ0 * tw + 8 * dZe * tw + 3 * ddZ0 * pow(tw, 2) - ddZe * pow(tw, 2)) / (2 * pow(tw, 3));
    double n4 = (30 * Z0 - 30 * Ze + 16 * dZ0 * tw + 14 * dZe * tw + 3 * ddZ0 * pow(tw, 2) - 2 * ddZe * pow(tw, 2)) / (2 * pow(tw, 4));
    double n5 = -(12 * Z0 - 12 * Ze + 6 * dZ0 * tw + 6 * dZe * tw + ddZ0 * pow(tw, 2) - ddZe * pow(tw, 2)) / (2 * pow(tw, 5));

    if (RealTime >= t_start && RealTime <= t_end)
    {
        Pos = n0 + n1 * Rt + n2 * pow(Rt, 2) + n3 * pow(Rt, 3) + n4 * pow(Rt, 4) + n5 * pow(Rt, 5);
        Vel = n1 + 2 * n2 * Rt + 3 * n3 * pow(Rt, 2) + 4 * n4 * pow(Rt, 3) + 5 * n5 * pow(Rt, 4);
        Acc = 2 * n2 + 6 * n3 * Rt + 12 * n4 * pow(Rt, 2) + 20 * n5 * pow(Rt, 3);
    }
    else if (RealTime < t_start)
    {
        Pos = Z0;
        Vel = dZ0;
        Acc = ddZ0;
    }
    else if (RealTime > t_end)
    {
        Pos = Ze;
        Vel = dZe;
        Acc = ddZe;
    }
    trajOut << Pos, Vel, Acc;
    return trajOut;
}

Eigen::Vector3d FuncPoly6th(double RealTime, double t_start, double t_end, double Z0, double dZ0, double ddZ0, double Ze, double dZe, double ddZe, double Fh)
{
    double tw = t_end - t_start;
    double Rt = RealTime - t_start;
    double Pos, Vel, Acc;
    Eigen::Vector3d trajOut;

    double n0 = Z0;
    double n1 = dZ0;
    double n2 = ddZ0 / 2;
    double n3 = -(84 * Z0 - 128 * Fh + 44 * Ze + 32 * dZ0 * tw - 12 * dZe * tw + 5 * ddZ0 * pow(tw, 2) + ddZe * pow(tw, 2)) / (2 * pow(tw, 3));
    double n4 = (222 * Z0 - 384 * Fh + 162 * Ze + 76 * dZ0 * tw - 46 * dZe * tw + 9 * ddZ0 * pow(tw, 2) + 4 * ddZe * pow(tw, 2)) / (2 * pow(tw, 4));
    double n5 = -(204 * Z0 - 384 * Fh + 180 * Ze + 66 * dZ0 * tw - 54 * dZe * tw + 7 * ddZ0 * pow(tw, 2) + 5 * ddZe * pow(tw, 2)) / (2 * pow(tw, 5));
    double n6 = (32 * Z0 - 64 * Fh + 32 * Ze + 10 * dZ0 * tw - 10 * dZe * tw + ddZ0 * pow(tw, 2) + ddZe * pow(tw, 2)) / pow(tw, 6);

    if (RealTime >= t_start && RealTime <= t_end)
    {
        Pos = n0 + n1 * Rt + n2 * pow(Rt, 2) + n3 * pow(Rt, 3) + n4 * pow(Rt, 4) + n5 * pow(Rt, 5) + n6 * pow(Rt, 6);
        Vel = n1 + 2 * n2 * Rt + 3 * n3 * pow(Rt, 2) + 4 * n4 * pow(Rt, 3) + 5 * n5 * pow(Rt, 4) + 6 * n6 * pow(Rt, 5);
        Acc = 2 * n2 + 6 * n3 * Rt + 12 * n4 * pow(Rt, 2) + 20 * n5 * pow(Rt, 3) + 30 * n6 * pow(Rt, 4);
    }
    else if (RealTime < t_start)
    {
        Pos = Z0;
        Vel = dZ0;
        Acc = ddZ0;
    }
    else if (RealTime > t_end)
    {
        Pos = Ze;
        Vel = dZe;
        Acc = ddZe;
    }
    trajOut << Pos, Vel, Acc;
    return trajOut;
}

Eigen::Matrix3d JacTranspose_LF(double q0, double q1, double q2, double q3)
{
    /* Left Front 3x3 Jacobian Matrix Transpose */

    double A0[3][3];
    Eigen::Matrix3d Out;

    // Link Lengths
    double d0 = 0.079;
    double d1 = 0.0262;
    double d2 = 0.1102;
    double d3 = -0.27;
    double d4 = -0.224282;
    double d5 = -0.107131;

    double t2 = cos(q0);
    double t3 = cos(q1);
    double t4 = sin(q0);
    double t5 = sin(q1);
    double t6 = q1 + q2;
    double t7 = cos(t6);
    double t8 = q3 + t6;
    double t9 = sin(t6);
    double t10 = cos(t8);
    double t11 = sin(t8);
    double t12 = d4 * t7;
    double t14 = d4 * t2 * t9;
    double t15 = d4 * t4 * t9;
    double t13 = d5 * t10;
    double t16 = d5 * t2 * t11;
    double t17 = d5 * t4 * t11;
    double t18 = -t14;
    double t19 = -t16;
    A0[0][0] = 0.0;
    A0[0][1] = -d1 * t4 - d2 * t4 - t2 * t12 - t2 * t13 - d3 * t2 * t3;
    A0[0][2] = d1 * t2 + d2 * t2 - t4 * t12 - t4 * t13 - d3 * t3 * t4;
    A0[1][0] = t12 + t13 + d3 * t3;
    A0[1][1] = t15 + t17 + d3 * t4 * t5;
    A0[1][2] = t18 + t19 - d3 * t2 * t5;
    A0[2][0] = t12 + t13;
    A0[2][1] = t15 + t17;
    A0[2][2] = t18 + t19;

    Out << A0[0][0], A0[0][1], A0[0][2],
        A0[1][0], A0[1][1], A0[1][2],
        A0[2][0], A0[2][1], A0[2][2];
    return Out;
}

Eigen::Matrix3d JacTranspose_RF(double q0, double q1, double q2, double q3)
{
    /* Right Front 3x3 Jacobian Matrix Transpose */

    double A0[3][3];
    Eigen::Matrix3d Out;

    // Link Lengths
    double d0 = 0.079;
    double d1 = 0.0262;
    double d2 = 0.1102;
    double d3 = -0.27;
    double d4 = -0.224282;
    double d5 = -0.107131;

    double t2 = cos(q0);
    double t3 = cos(q1);
    double t4 = sin(q0);
    double t5 = sin(q1);
    double t6 = q1 + q2;
    double t7 = cos(t6);
    double t8 = q3 + t6;
    double t9 = sin(t6);
    double t10 = cos(t8);
    double t11 = sin(t8);
    double t12 = d4 * t7;
    double t14 = d4 * t2 * t9;
    double t15 = d4 * t4 * t9;
    double t13 = d5 * t10;
    double t16 = -t12;
    double t17 = d5 * t2 * t11;
    double t18 = d5 * t4 * t11;
    double t20 = -t14;
    double t19 = -t13;
    double t21 = -t17;
    A0[0][0] = 0.0;
    A0[0][1] = d1 * t4 + d2 * t4 + t2 * t16 + t2 * t19 - d3 * t2 * t3;
    A0[0][2] = -d1 * t2 - d2 * t2 + t4 * t16 + t4 * t19 - d3 * t3 * t4;
    A0[1][0] = t16 + t19 - d3 * t3;
    A0[1][1] = t15 + t18 + d3 * t4 * t5;
    A0[1][2] = t20 + t21 - d3 * t2 * t5;
    A0[2][0] = t16 + t19;
    A0[2][1] = t15 + t18;
    A0[2][2] = t20 + t21;

    Out << A0[0][0], A0[0][1], A0[0][2],
        A0[1][0], A0[1][1], A0[1][2],
        A0[2][0], A0[2][1], A0[2][2];
    return Out;
}

Eigen::Matrix3d JacTranspose_LB(double q0, double q1, double q2, double q3)
{
    /* Left Back 3x3 Jacobian Matrix Transpose */

    double A0[3][3];
    Eigen::Matrix3d Out;

    // Link Lengths
    double d0 = 0.079;
    double d1 = 0.0262;
    double d2 = 0.1102;
    double d3 = -0.27;
    double d4 = -0.224282;
    double d5 = -0.107131;

    double t2 = cos(q0);
    double t3 = cos(q1);
    double t4 = sin(q0);
    double t5 = sin(q1);
    double t6 = q1 + q2;
    double t7 = cos(t6);
    double t8 = q3 + t6;
    double t9 = sin(t6);
    double t10 = cos(t8);
    double t11 = sin(t8);
    double t12 = d4 * t7;
    double t14 = d4 * t2 * t9;
    double t15 = d4 * t4 * t9;
    double t13 = d5 * t10;
    double t16 = d5 * t2 * t11;
    double t17 = d5 * t4 * t11;
    double t18 = -t14;
    double t19 = -t15;
    double t20 = -t16;
    double t21 = -t17;
    A0[0][0] = 0.0;
    A0[0][1] = -d1 * t4 - d2 * t4 + t2 * t12 + t2 * t13 + d3 * t2 * t3;
    A0[0][2] = -d1 * t2 - d2 * t2 - t4 * t12 - t4 * t13 - d3 * t3 * t4;
    A0[1][0] = t12 + t13 + d3 * t3;
    A0[1][1] = t19 + t21 - d3 * t4 * t5;
    A0[1][2] = t18 + t20 - d3 * t2 * t5;
    A0[2][0] = t12 + t13;
    A0[2][1] = t19 + t21;
    A0[2][2] = t18 + t20;

    Out << A0[0][0], A0[0][1], A0[0][2],
        A0[1][0], A0[1][1], A0[1][2],
        A0[2][0], A0[2][1], A0[2][2];
    return Out;
}

Eigen::Matrix3d JacTranspose_RB(double q0, double q1, double q2, double q3)
{
    /* Right Back 3x3 Jacobian Matrix Transpose */

    double A0[3][3];
    Eigen::Matrix3d Out;

    // Link Lengths
    double d0 = 0.079;
    double d1 = 0.0262;
    double d2 = 0.1102;
    double d3 = -0.27;
    double d4 = -0.224282;
    double d5 = -0.107131;

    double t2 = cos(q0);
    double t3 = cos(q1);
    double t4 = sin(q0);
    double t5 = sin(q1);
    double t6 = q1 + q2;
    double t7 = cos(t6);
    double t8 = q3 + t6;
    double t9 = sin(t6);
    double t10 = cos(t8);
    double t11 = sin(t8);
    double t12 = d4 * t7;
    double t14 = d4 * t2 * t9;
    double t15 = d4 * t4 * t9;
    double t13 = d5 * t10;
    double t16 = -t12;
    double t17 = d5 * t2 * t11;
    double t18 = d5 * t4 * t11;
    double t20 = -t14;
    double t21 = -t15;
    double t19 = -t13;
    double t22 = -t17;
    double t23 = -t18;
    A0[0][0] = 0.0;
    A0[0][1] = d1 * t4 + d2 * t4 + t2 * t12 + t2 * t13 + d3 * t2 * t3;
    A0[0][2] = d1 * t2 + d2 * t2 + t4 * t16 + t4 * t19 - d3 * t3 * t4;
    A0[1][0] = t16 + t19 - d3 * t3;
    A0[1][1] = t21 + t23 - d3 * t4 * t5;
    A0[1][2] = t20 + t22 - d3 * t2 * t5;
    A0[2][0] = t16 + t19;
    A0[2][1] = t21 + t23;
    A0[2][2] = t20 + t22;

    Out << A0[0][0], A0[0][1], A0[0][2],
        A0[1][0], A0[1][1], A0[1][2],
        A0[2][0], A0[2][1], A0[2][2];
    return Out;
}

Eigen::Matrix<double, 4, 3> refForceCalc4(const Eigen::Vector3d& Rcf1, const Eigen::Vector3d& Rcf2, const Eigen::Vector3d& Rcf3, const Eigen::Vector3d& Rcf4, const Eigen::Vector3d& qlf, const Eigen::Vector3d& qrf, const Eigen::Vector3d& qlh, const Eigen::Vector3d& qrh, double dt)
{
    double F1x, F2x, F3x, F4x;
    double F1y, F2y, F3y, F4y;
    double F1z, F2z, F3z, F4z;
    Eigen::Matrix<double, 4, 3> Out;

    double q01 = qlf(0);
    double q11 = qlf(1);
    double q21 = qlf(2);
    double q02 = qrf(0);
    double q12 = qrf(1);
    double q22 = qrf(2);
    double q03 = qlh(0);
    double q13 = qlh(1);
    double q23 = qlh(2);
    double q04 = qrh(0);
    double q14 = qrh(1);
    double q24 = qrh(2);
    double qf_f = 30.0 * PI / 180;
    double qf_h = 30.0 * PI / 180;


    if ((Rcf1(2)) > 0 || (Rcf4(2)) > 0)
    {
        F1z = 0; F2z = MASS * GRAVITY / 2; F3z = MASS * GRAVITY / 2; F4z = 0;
    }
    else if ((Rcf2(2)) > 0 || (Rcf3(2)) > 0)
    {
        F1z = MASS * GRAVITY / 2; F2z = 0; F3z = 0; F4z = MASS * GRAVITY / 2;
    }
    else
    {
        F1z = MASS * GRAVITY / 4; F2z = MASS * GRAVITY / 4; F3z = MASS * GRAVITY / 4; F4z = MASS * GRAVITY / 4;
    }

    
    F1x = F1z * tan(q11 + q21 + qf_f) * (1 / cos(q01));
    F2x = F2z * tan(q12 + q22 - qf_f) * (1 / cos(q02));
    F3x = F3z * tan(q13 + q23 + qf_h) * (1 / cos(q03));
    F4x = F4z * tan(q14 + q24 - qf_h) * (1 / cos(q04));

    F1y = -F1z * tan(q01);
    F2y = -F2z * tan(q02);
    F3y = -F3z * tan(q03);
    F4y = -F4z * tan(q04);

    Out << F1x, F1y, F1z,
        F2x, F2y, F2z,
        F3x, F3y, F3z,
        F4x, F4y, F4z;
    return Out;
}

Eigen::Matrix<double, 4, 3> VMC(const Eigen::Vector3d& Rcf1, const Eigen::Vector3d& Rcf2, const Eigen::Vector3d& Rcf3, const Eigen::Vector3d& Rcf4, const Eigen::Vector3d& Fc_LF, const Eigen::Vector3d& Fc_RF, const Eigen::Vector3d& Fc_LB, const Eigen::Vector3d& Fc_RB, Eigen::VectorXd C, double dt)
{
    double F1x, F2x, F3x, F4x;
    double F1y, F2y, F3y, F4y;
    double F1z, F2z, F3z, F4z;
    Eigen::Matrix<double, 4, 3> Out;

    int contactState_L, contactState_R, contactState_LB, contactState_RB;

    // Foot position wrt CoM
    double r1xc = Rcf1(0); double r1yc = Rcf1(1); double r1zc = Rcf1(2); // LF
    double r2xc = Rcf2(0); double r2yc = Rcf2(1); double r2zc = Rcf2(2); // RF
    double r3xc = Rcf3(0); double r3yc = Rcf3(1); double r3zc = Rcf3(2); // LB
    double r4xc = Rcf4(0); double r4yc = Rcf4(1); double r4zc = Rcf4(2); // RB

    
    Eigen::Matrix <double, 6, 12> P;
    Eigen::VectorXd F(12);
    F.setZero();

    if(Fc_LF(2) > 0 || Fc_RB(2) > 20) { contactState_L = 1; }
    else { contactState_L = 0; }
    if(Fc_RF(2) > 0 || Fc_LB(2) > 20) { contactState_R = 1; }
    else { contactState_R = 0; }

    P << Eigen::MatrixXd::Identity(3,3)*contactState_L, Eigen::MatrixXd::Identity(3,3)*contactState_R, Eigen::MatrixXd::Identity(3,3)*contactState_R, Eigen::MatrixXd::Identity(3,3)*contactState_L,
         vec2SkewSym(Rcf1)*contactState_L, vec2SkewSym(Rcf2)*contactState_R, vec2SkewSym(Rcf3)*contactState_R, vec2SkewSym(Rcf4)*contactState_L;
    F = P.completeOrthogonalDecomposition().solve(C);
    
    F1x = F(0); F2x = F(3); F3x = F(6); F4x = F(9);
    F1y = F(1); F2y = F(4); F3y = F(7); F4y = F(10);
    F1z = F(2); F2z = F(5); F3z = F(8); F4z = F(11);

    Out << F1x, F1y, F1z,
        F2x, F2y, F2z,
        F3x, F3y, F3z,
        F4x, F4y, F4z;
    return Out;
}

double SolveQuadCosSin(double ann, double bnn, double cnn, int mu)
{
    // ann = bnn * cos(x) + cnn * sin(x)
    // asinx + bcosx = c WEB
    double nx = sqrt(pow(bnn, 2) + pow(cnn, 2));
    double an = ann / nx;
    double bn = bnn / nx;
    double cn = cnn / nx;
    double dn = sqrt(pow(bn, 2) + pow(cn, 2) - pow(an, 2));
    double Out1 = atan2((mu * dn * bn + an * cn), (bn * an - mu * dn * cn));
    ///// Out = mu * atan2(dn, an) + atan2(cn, bn);% Onceki calisan kısım
    double Out = -2 * PI * 0 + Out1;
    return Out;
    // hnn = ann ^ 2 - cnn ^ 2;
    //
    // Out = atan2(((cnn * bnn + mu * ann * sqrt(bnn ^ 2 - hnn))), hnn);
}

Eigen::Vector3d fullBodyFK(Eigen::Vector3d torsoOrient, Eigen::Vector3d Rcom, Eigen::Vector3d Q, int u)
{
    int m, n, f;
    double L, W, H;

    double R = torsoOrient(0);
    double P = torsoOrient(1);
    double Y = torsoOrient(2);

    double Xc = Rcom(0);
    double Yc = Rcom(1);
    double Zc = Rcom(2);
    
    double q0 = Q(0);
    double q1 = Q(1);
    double q2 = Q(2);
    double q3 = 30*PI/180;

    Eigen::Vector3d P7;

    switch (u)
    {
    case 1:
        m = 1;
        n = 1;
        f = 1;
        L = 0.3102;
        W = 0.105;
        H = 0.002;
        break;
    case 2:
        m = 1;
        n = -1;
        f = 1;
        L = 0.3102;
        W = -0.105;
        H = 0.002;
        break;
    case 3:
        m = -1;
        n = 1;
        f = -1;
        L = -0.3102;
        W = 0.105;
        H = 0.002;
        break;
    case 4:
        m = -1;
        n = -1;
        f = -1;
        L = -0.3102;
        W = -0.105;
        H = 0.002;
        break;
    }

    // Link Lengths
    double d0 = m * 0.079;
    double d1 = n * 0.0262;
    double d2 = n * 0.1102;
    double d3 = -0.27;
    double d4 = -0.224282;
    double d5 = -0.107131; //-0.107131;

    P7[0] = Xc+d3*(cos(q1)*(sin(P)*cos(q0)+cos(P)*sin(Y)*sin(q0))+cos(P)*cos(Y)*sin(q1))+d5*(cos(f*q3)*(cos(q2)*(cos(q1)*(sin(P)*cos(q0)+cos(P)*sin(Y)*sin(q0))+cos(P)*cos(Y)*sin(q1))-sin(q2)*(sin(q1)*(sin(P)*cos(q0)+cos(P)*sin(Y)*sin(q0))-cos(P)*cos(Y)*cos(q1)))-sin(f*q3)*(cos(q2)*(sin(q1)*(sin(P)*cos(q0)+cos(P)*sin(Y)*sin(q0))-cos(P)*cos(Y)*cos(q1))+sin(q2)*(cos(q1)*(sin(P)*cos(q0)+cos(P)*sin(Y)*sin(q0))+cos(P)*cos(Y)*sin(q1))))+d4*(cos(q2)*(cos(q1)*(sin(P)*cos(q0)+cos(P)*sin(Y)*sin(q0))+cos(P)*cos(Y)*sin(q1))-sin(q2)*(sin(q1)*(sin(P)*cos(q0)+cos(P)*sin(Y)*sin(q0))-cos(P)*cos(Y)*cos(q1)))+H*sin(P)+d1*(sin(P)*sin(q0)-cos(P)*sin(Y)*cos(q0))+d2*(sin(P)*sin(q0)-cos(P)*sin(Y)*cos(q0))+L*cos(P)*cos(Y)+d0*cos(P)*cos(Y)-W*cos(P)*sin(Y);
    P7[1] = Yc+d4*(cos(q2)*(sin(q1)*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))-cos(q1)*(sin(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))+cos(P)*sin(R)*cos(q0)))+sin(q2)*(cos(q1)*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))+sin(q1)*(sin(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))+cos(P)*sin(R)*cos(q0))))+L*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))+d5*(cos(f*q3)*(cos(q2)*(sin(q1)*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))-cos(q1)*(sin(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))+cos(P)*sin(R)*cos(q0)))+sin(q2)*(cos(q1)*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))+sin(q1)*(sin(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))+cos(P)*sin(R)*cos(q0))))+sin(f*q3)*(cos(q2)*(cos(q1)*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))+sin(q1)*(sin(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))+cos(P)*sin(R)*cos(q0)))-sin(q2)*(sin(q1)*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))-cos(q1)*(sin(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))+cos(P)*sin(R)*cos(q0)))))+W*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))+d0*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))+d3*(sin(q1)*(cos(R)*sin(Y)+cos(Y)*sin(P)*sin(R))-cos(q1)*(sin(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))+cos(P)*sin(R)*cos(q0)))+d1*(cos(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))-cos(P)*sin(R)*sin(q0))+d2*(cos(q0)*(cos(R)*cos(Y)-sin(P)*sin(R)*sin(Y))-cos(P)*sin(R)*sin(q0))-H*cos(P)*sin(R);
    P7[2] = Zc+d5*(cos(f*q3)*(cos(q2)*(sin(q1)*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))-cos(q1)*(sin(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))-cos(P)*cos(R)*cos(q0)))+sin(q2)*(cos(q1)*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))+sin(q1)*(sin(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))-cos(P)*cos(R)*cos(q0))))+sin(f*q3)*(cos(q2)*(cos(q1)*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))+sin(q1)*(sin(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))-cos(P)*cos(R)*cos(q0)))-sin(q2)*(sin(q1)*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))-cos(q1)*(sin(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))-cos(P)*cos(R)*cos(q0)))))+d4*(cos(q2)*(sin(q1)*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))-cos(q1)*(sin(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))-cos(P)*cos(R)*cos(q0)))+sin(q2)*(cos(q1)*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))+sin(q1)*(sin(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))-cos(P)*cos(R)*cos(q0))))+L*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))+W*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))+d3*(sin(q1)*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))-cos(q1)*(sin(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))-cos(P)*cos(R)*cos(q0)))+d0*(sin(R)*sin(Y)-cos(R)*cos(Y)*sin(P))+d1*(cos(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))+cos(P)*cos(R)*sin(q0))+d2*(cos(q0)*(cos(Y)*sin(R)+cos(R)*sin(P)*sin(Y))+cos(P)*cos(R)*sin(q0))+H*cos(P)*cos(R);

    return P7;
}

Eigen::Vector3d fullBodyIKan(Eigen::Vector3d Rfoot, Eigen::Vector3d Rcom, Eigen::Vector3d torsoOrient, int u)
{
    int m, n, f;
    double L, W, H;

    switch (u)
    {
    case 1:
        m = 1;
        n = 1;
        f = 1;
        L = 0.3102;
        W = 0.105;
        H = 0.002;
        break;
    case 2:
        m = 1;
        n = -1;
        f = 1;
        L = 0.3102;
        W = -0.105;
        H = 0.002;
        break;
    case 3:
        m = -1;
        n = 1;
        f = -1;
        L = -0.3102;
        W = 0.105;
        H = 0.002;
        break;
    case 4:
        m = -1;
        n = -1;
        f = -1;
        L = -0.3102;
        W = -0.105;
        H = 0.002;
        break;
    }

    // Link Lengths
    double d0 = m * 0.079;
    double d1 = n * 0.0262;
    double d2 = n * 0.1102;
    double d3 = -0.27;
    double d4 = -0.224282;
    double d5 = -0.107131; //-0.107131;

    // Hip AA position wrt CoM
    Eigen::Vector3d Rcom2hip;
    //Rcom2hip << L / 2 + dx_hip, W / 2 - dy_hip, 0;
    Rcom2hip << L, W, H;

    // Euler Angles (torsoOrient: torso orinetatiın vector) (torsoOrient = [torsoRoll, torsoPitch, torsoYaw])
    Eigen::Matrix3d torsoRotation;
    // Eigen::Vector4d quat;
    // quat = eulerToQuaternion(torsoOrient(2),torsoOrient(1),torsoOrient(0));
    // torsoRotation = quat2Rotmat(quat(0),quat(1),quat(2),quat(3));
    torsoRotation = RotateYaw(torsoOrient(2))*RotatePitch(torsoOrient(1))*RotateRoll(torsoOrient(0));

    // Rcom: CoM position vector
    // Rfoot: Foot position vector

    // Foot position wrt Hip AA
    Eigen::Vector3d Rhip2foot;
    Rhip2foot = torsoRotation.transpose() * (Rfoot - Rcom) - Rcom2hip;
    double Pfx = Rhip2foot(0);
    double Pfy = Rhip2foot(1);
    double Pfz = Rhip2foot(2);

    Eigen::Vector3d Q;

    double q3 = f * 30.0 * PI / 180;
    double ann2 = (pow(Pfx - d0, 2) + pow(Pfy, 2) + pow(Pfz, 2) - pow(d1 + d2, 2) - pow(d3, 2) - pow(d4, 2) - pow(d5, 2) - 2 * cos(q3) * d4 * d5) / d3;
    double d6 = 2 * d4 + 2 * cos(q3) * d5;
    double q2 = SolveQuadCosSin(ann2, d6, -2 * sin(q3) * d5, f);
    double s2 = sin(q2);
    double c2 = cos(q2);

    double ann1 = Pfx - d0;
    double bnn1 = d4 * s2 + d5 * c2 * sin(q3) + d5 * cos(q3) * s2; // d6* s2 / 2 + d5 * c2 / 2
    double cnn1 = d3 + d4 * c2 - d5 * s2 * sin(q3) + d5 * cos(q3) * c2; // d3 + d6 * c2 / 2 - d5 * s2 / 2
    double q1 = SolveQuadCosSin(ann1, bnn1, cnn1, 1);


    double An = d1 + d2;
    double Bn = d4 * cos(q1 + q2) + d3 * cos(q1) + d5 * cos(q1 + q2 + q3);
    // Cosq0 = (Pfy * A + Pfz * B) / (A ^ 2 + B ^ 2);
    // Sinq0 = sqrt(1 - Cosq0 ^ 2);
    // q0 = atan2(Sinq0, Cosq0);
    // q0 = SolveQuadCosSin(Pfz, Bn, An, -1);
    double q0 = SolveQuadCosSin(Pfy, An, -Bn, -1);

    Q << m*q0, n*q1, n*q2;
    return Q;
}

Eigen::Vector3d balanceControl(Eigen::Matrix3d Rd, Eigen::Matrix3d R)
{
    Eigen::Matrix3d Re, angleSkw;
    Eigen::Vector3d angles;
    double thetha;
    Re = Rd*R.transpose();

    // thetha = acos((Re.trace() - 1)*0.5);
    // angleSkw = (Re - Re.transpose())*thetha/(2*sin(thetha));

    angleSkw = Re.log();

    angles << angleSkw(2,1), angleSkw(0,2), angleSkw(1,0);
    return angles;
}

/* Inverse Dynamics */
Eigen::Matrix3d quat2Rotmat(double qw, double qx, double qy, double qz)
{
    double r31, r32, r33, r11, r12, r13, r21, r22, r23;
    Eigen::Matrix3d rMat;
    r11 = 1 - 2 * qy * qy - 2 * qz * qz;
    r12 = 2 * (qx * qy - qz * qw);
    r13 = 2 * (qx * qz + qy * qw);
    r21 = 2 * (qx * qy + qz * qw);
    r22 = 1 - 2 * qx * qx - 2 * qz * qz;
    r23 = 2 * (qy * qz - qx * qw);
    r31 = 2 * (qx * qz - qy * qw);
    r32 = 2 * (qy * qz + qx * qw);
    r33 = 1 - 2 * qx * qx - 2 * qy * qy;
    rMat << r11, r12, r13, r21, r22, r23, r31, r32, r33;
    return rMat;
}

Eigen::Vector3d funcNewtonEuler(Eigen::Vector3d rootAbsAcc, Eigen::Matrix3d rootOrient, Eigen::Vector3d rootAngVel, Eigen::Vector3d rootAngAcc, Eigen::Vector3d jPos, Eigen::Vector3d jVel, Eigen::Vector3d jAcc, Eigen::Vector3d grForce, Eigen::Vector3d Rcon, int legIndex)
{
    Eigen::Vector3d JTorques;
    Eigen::Matrix3d RG0, R01, R12, R23, R34, R45, RG5;
    Eigen::Vector3d P01, P12, P23, P34, P45;
    Eigen::Vector3d Pc0, Pc1, Pc2, Pc3, Pc4;
    Eigen::Vector3d Jvel0, Jvel1, Jvel2, Jvel3;
    Eigen::Vector3d Jacc0, Jacc1, Jacc2, Jacc3;
    Eigen::Vector3d w0, w1, w2, w3, w4;
    Eigen::Vector3d dw0, dw1, dw2, dw3, dw4;
    Eigen::Vector3d dv0, dv1, dv2, dv3, dv4;
    Eigen::Vector3d dvc0, dvc1, dvc2, dvc3, dvc4;
    Eigen::Vector3d Temp;
    Eigen::Vector3d F0, F1, F2, F3, F4, N0, N1, N2, N3, N4;
    Eigen::Vector3d f0, f1, f2, f3, f4, f5, n0, n1, n2, n3, n4, n5;
    Eigen::Matrix3d  Inertia0, Inertia1, Inertia2, Inertia3, Inertia4;
    double q0, q1, q2, q3;
    double dq0, dq1, dq2, dq3;
    double ddq0, ddq1, ddq2, ddq3;
    double m0 = 20.516;
    double m1 = 1.681; // Hip AA 2 Hip FE
    double m2 = 2.005; // Hip FE 2 Knee FE
    double m3 = 0.21; // Knee FE 2 Ankle
    double m4 = 0;
    Eigen::Vector3d gravityVec;
    gravityVec << 0, 0, GRAVITY;

    double m, n, k;
    switch (legIndex) {
    case 1: // LF
        m = 1; // For roll and x axis
        n = 1; // For pitch and y axis
        q3 = 30*PI/180;
        break;
    case 2: // RF
        m = 1; // For roll and x axis
        n = -1; // For pitch and y axis
        q3 = -30*PI/180;
        break;
    case 3: // LB
        m = -1; // For roll and x axis
        n = 1; // For pitch and y axis
        q3 = -30*PI/180;
        break;
    case 4: // RB
        m = -1; // For roll and x axis
        n = -1; // For pitch and y axis
        q3 = 30*PI/180;
        break;
    }

    // Joint angle
    q0 = m*jPos(0);
    q1 = n*jPos(1);
    q2 = n*jPos(2);
    q3 = n*q3;
    // Joint angular velocity
    dq0 = m*jVel(0);
    dq1 = n*jVel(1);
    dq2 = n*jVel(2);
    dq3 = 0;
    Jvel0 << dq0, 0, 0;
    Jvel1 << 0, dq1, 0;
    Jvel2 << 0, dq2, 0;
    Jvel3 << 0, 0, 0;
    // Joint angular acceleration
    ddq0 = m*jAcc(0);
    ddq1 = n*jAcc(1);
    ddq2 = n*jAcc(2);
    ddq3 = 0;
    Jacc0 << ddq0, 0, 0;
    Jacc1 << 0, ddq1, 0;
    Jacc2 << 0, ddq2, 0;
    Jacc3 << 0, 0, 0;

    // Joint position wrt previous joint
    P01 << m*0.3102, n*0.105, 0.002;
    P12 << m*0.079, n*0.0262, 0.0;
    P23 << m*0.0, n*0.1102, -0.27;
    P34 << m*0.0, n*0.0, -0.224282;
    P45 << m*0, n*0, -0.107131;
    //P45 << m * sin(q3)*(-0.066), n * 0, cos(q3)*(-0.066);

    // Link center of mass positions wrt previous joint
    Pc0 << 0.00145, 0.003271, 0.01792;
    Pc1 << m*(-0.003717), n*(-0.015862), -0.000015;
    Pc2 << m*(-0.003917), n*(0.076506), -0.013793;
    Pc3 << m*(-0.015346), n*(-0.000102), -0.186842;
    Pc4 = 0.5 * P45;

    // Link inertia
    Inertia0 << 0.410, 0, 0,
                0, 0.908, 0,
                0, 0, 1.1923;
    Inertia1 << 0.003, 0, 0,
                0, 0.003, 0,
                0, 0, 0.003;
    Inertia2 << 0.008, 0, 0,
                0, 0.007, 0,
                0, 0, 0.004;
    Inertia3 << 0.002, 0, 0,
                0, 0.002, 0,
                0, 0, 1.086e-4;
    Inertia4 << 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0;
    // Inertia1 << 0.003, 9.851e-5, -1.213e-6,
    //             9.851e-5, 0.003, -1.327e-6,
    //             -1.213e-6, -1.327e-6, 0.003;
    // Inertia2 << 0.008, -1.818e-4, 5.307e-4,
    //             -1.818e-4, 0.007, 0.001,
    //             5.307e-4, 0.001, 0.004;
    // Inertia3 << 0.002, -7.048e-7, -2.682e-4,
    //             -7.048e-7, 0.002, -1.088e-5,
    //             -2.682e-4, -1.088e-5, 1.086e-4;
    // Inertia4 << 0.0, 0.0, 0.0,
    //             0.0, 0.0, 0.0,
    //             0.0, 0.0, 0.0;

    RG0 = rootOrient;
    R01 = RotateRoll(q0);
    R12 = RotatePitch(q1);
    R23 = RotatePitch(q2);
    R34 = RotatePitch(q3);
    R45.setIdentity();


    /*********** Outward iteration **********/

    // root (Joint 0)
    w0 = RG0.transpose() * rootAngVel;
    dw0 = RG0.transpose() * rootAngAcc;
    // Hip AA (Joint 1)
    Temp = R01.transpose() * w0;
    w1 = Temp + Jvel0;
    dw1 = R01.transpose() * dw0 + Temp.cross(Jvel0) + Jacc0;
    // Hip FE (Joint 2)
    Temp = R12.transpose() * w1;
    w2 = Temp + Jvel1;
    dw2 = R12.transpose() * dw1 + Temp.cross(Jvel1) + Jacc1;
    // Knee FE (Joint 3)
    Temp = R23.transpose() * w2;
    w3 = Temp + Jvel2;
    dw3 = R23.transpose() * dw2 + Temp.cross(Jvel2) + Jacc2;
    // Ankle (Joint 4)
    Temp = R34.transpose() * w3;
    w4 = Temp + Jvel3;
    dw4 = R34.transpose() * dw3 + Temp.cross(Jvel3) + Jacc3;

    // root (Joint 0)
    dv0 = RG0.transpose() * (rootAbsAcc - gravityVec);
    // Hip AA (Joint 1)
    Temp = w0.cross(P01);
    Temp = w0.cross(Temp);
    Temp = dw0.cross(P01) + Temp + dv0;
    dv1 = R01.transpose() * Temp;
    // Hip FE (Joint 2)
    Temp = w1.cross(P12);
    Temp = w1.cross(Temp);
    Temp = dw1.cross(P12) + Temp + dv1;
    dv2 = R12.transpose() * Temp;
    // Knee FE (Joint 3)
    Temp = w2.cross(P23);
    Temp = w2.cross(Temp);
    Temp = dw2.cross(P23) + Temp + dv2;
    dv3 = R23.transpose() * Temp;
    // Ankle (Joint 4)
    Temp = w3.cross(P34);
    Temp = w3.cross(Temp);
    Temp = dw3.cross(P34) + Temp + dv3;
    dv4 = R34.transpose() * Temp;

    // Root (Torso)
    Temp = rootAngVel.cross(Pc0);
    Temp = rootAngVel.cross(Temp);
    dvc0 = rootAngAcc.cross(Pc0) + Temp + (rootAbsAcc - gravityVec);
    // Link 1 (Between HAA and HFE joints)
    Temp = w1.cross(Pc1);
    Temp = w1.cross(Temp);
    dvc1 = dw1.cross(Pc1) + Temp + dv1;
    // Link 2 (Between HFE and KFE joints)
    Temp = w2.cross(Pc2);
    Temp = w2.cross(Temp);
    dvc2 = dw2.cross(Pc2) + Temp + dv2;
    // Link 3 (Between KFE and ankle joints)
    Temp = w3.cross(Pc3);
    Temp = w3.cross(Temp);
    dvc3 = dw3.cross(Pc3) + Temp + dv3;
    // Link 4 (Between Ankle and foot tip (frame 5) joints)
    Temp = w4.cross(Pc4);
    Temp = w4.cross(Temp);
    dvc4 = dw4.cross(Pc4) + Temp + dv4;

    F0 = m0 * dvc0;
    F1 = m1 * dvc1;
    F2 = m2 * dvc2;
    F3 = m3 * dvc3;
    F4 = m4 * dvc4;
    N0 = Inertia0 * rootAngAcc + rootAngVel.cross(Inertia0 * rootAngVel);
    N1 = Inertia1 * dw1 + w1.cross(Inertia1 * w1);
    N2 = Inertia2 * dw2 + w2.cross(Inertia2 * w2);
    N3 = Inertia3 * dw3 + w3.cross(Inertia3 * w3);
    N4 = Inertia4 * dw4 + w4.cross(Inertia4 * w4);

    /*********** Inward iteration **********/
    RG5 = RG0 * R01 * R12 * R23 * R34 * R45;
    f5 = RG5.transpose() * (-grForce);
    // f5 << 0, 0, 0;
    f4 = R45 * f5 + F4;
    f3 = R34 * f4 + F3;
    f2 = R23 * f3 + F2;
    f1 = R12 * f2 + F1;
    f0 = RG0 * f1 + F0;

    // Eigen::Vector3d contOffset;
    // contOffset = compContPoint(RG5);
    // Temp = contOffset.cross(grForce);
    // n5 = RG5.transpose() * Temp;

    // Temp = minPFunRotm(RG5);
    // n5 = Temp.cross(grForce);
    // n5 = -RG5.transpose()*n5;

    n5 << 0, 0, 0;
    /*n5 = Rcon.cross(grForce);
    n5 = -RG5.transpose()*n5;*/
    n4 = N4 + R45 * n5 + Pc4.cross(F4) + P45.cross(R45 * f5);
    n3 = N3 + R34 * n4 + Pc3.cross(F3) + P34.cross(R34 * f4);
    n2 = N2 + R23 * n3 + Pc2.cross(F2) + P23.cross(R23 * f3);
    n1 = N1 + R12 * n2 + Pc1.cross(F1) + P12.cross(R12 * f2);
    n0 = N0 + R01 * n1 + Pc0.cross(F0) + P01.cross(R01 * f1);

    JTorques << m*n1(0), n*n2(1), n*n3(1);
    return JTorques;
}

Eigen::Vector3d funcNewtonEuler2(Eigen::Vector3d rootAbsAcc, Eigen::Matrix3d rootOrient, Eigen::Vector3d rootAngVel, Eigen::Vector3d rootAngAcc, Eigen::Vector3d jPos, Eigen::Vector3d jVel, Eigen::Vector3d jAcc, Eigen::Vector3d grForce, Eigen::Vector3d Rcon, int legIndex)
{
    Eigen::Vector3d JTorques;
    double q0, q1, q2, q3, q4, q5, q6, qf;
    double dq0, dq1, dq2, dq3, dq4, dq5, dq6;
    double ddq0, ddq1, ddq2, ddq3, ddq4, ddq5, ddq6;
    Eigen::Matrix3d R10, R21, R32, R43, R54, R65, R60;
    Eigen::Vector3d Jvel0, Jvel1, Jvel2, Jvel3, Jvel4, Jvel5, Jvel6;
    Eigen::Vector3d Jacc0, Jacc1, Jacc2, Jacc3, Jacc4, Jacc5, Jacc6;
    Eigen::Vector3d P01, P12, P23, P34, P45, P56;
    Eigen::Vector3d Pc1, Pc2, Pc3, Pc4, Pc5;
    double m1, m2, m3, m4, m5;
    Eigen::Matrix3d Inertia1, Inertia2, Inertia3, Inertia4, Inertia5;

    Eigen::Vector3d w0, w1, w2, w3, w4, w5;
    Eigen::Vector3d dw0, dw1, dw2, dw3, dw4, dw5;
    Eigen::Vector3d dv0, dv1, dv2, dv3, dv4, dv5;
    Eigen::Vector3d dvc0, dvc1, dvc2, dvc3, dvc4, dvc5;
    Eigen::Vector3d F1, F2, F3, F4, F5; 
    Eigen::Vector3d N1, N2, N3, N4, N5;
    Eigen::Vector3d f1, f2, f3, f4, f5, f6; 
    Eigen::Vector3d n1, n2, n3, n4, n5, n6;

    Eigen::Vector3d gravityVec;
    gravityVec << 0, 0, GRAVITY;

    int m, n;
    switch (legIndex) {
    case 1: // LF
        m = 1; // For roll and x axis
        n = 1; // For pitch and y axis
        qf = 30*PI/180;
        break;
    case 2: // RF
        m = 1; // For roll and x axis
        n = -1; // For pitch and y axis
        qf = -30*PI/180;
        break;
    case 3: // LB
        m = -1; // For roll and x axis
        n = 1; // For pitch and y axis
        qf = -30*PI/180;
        break;
    case 4: // RB
        m = -1; // For roll and x axis
        n = -1; // For pitch and y axis
        qf = 30*PI/180;
        break;
    }

    /* #region: Joint position & orientation matrix */
    q0 = 0.0; // world joint (fixed)
    q1 = 0.0; // Torso orientation (q1 is not used rootOrient matrix is used instead) 
    q2 = m*jPos(0); // Hip AA joint
    q3 = n*jPos(1); // Hip FE joint
    q4 = n*jPos(2); //Knee FE joint
    q5 = n*qf; // Foot Frame joint
    q6 = 0.0; // Tip of foot

    // R10: Rotates points in frame {1} to frame {0}
    R10 = rootOrient; // Rotation from Torso frame to World frame
    R21 = RotateRoll(q2); // Rotation from Hip AA joint frame to Torso frame
    R32 = RotatePitch(q3); // Rotation from Hip FE joint frame to Hip AA joint frame
    R43 = RotatePitch(q4); // Rotation from Knee FE joint frame to Hip FE joint frame
    R54 = RotatePitch(q5); // Rotation from Foot frame to Knee FE joint frame
    R65.setIdentity(); // Rotation from Tip of foot frame to Foot frame
    /* #endregion */

    /* #region: Joint angular velocity */
    dq0 = 0.0; // world joint (fixed)
    dq1 = 0.0; // Torso angular velocity (dq1 is not used rootAngVel vector is used instead) 
    dq2 = m*jVel(0); // Hip AA joint
    dq3 = n*jVel(1); // Hip FE joint
    dq4 = n*jVel(2); //Knee FE joint
    dq5 = 0.0; // Foot Frame joint
    dq6 = 0.0; // Tip of foot

    Jvel0 << 0, 0, 0;
    Jvel1 = rootAngVel;
    Jvel2 << dq2, 0, 0;
    Jvel3 << 0, dq3, 0;
    Jvel4 << 0, dq4, 0;
    Jvel5 << 0, dq5, 0;
    Jvel6 << 0, dq6, 0;
    /* #endregion */

    /* #region: Joint angular acceleration */
    ddq0 = 0.0; // world joint (fixed)
    ddq1 = 0.0; // Torso angular acceleration (ddq1 is not used rootAngAcc vector is used instead) 
    ddq2 = m*jAcc(0); // Hip AA joint
    ddq3 = n*jAcc(1); // Hip FE joint
    ddq4 = n*jAcc(2); //Knee FE joint
    ddq5 = 0.0; // Foot Frame joint
    ddq6 = 0.0; // Tip of foot

    Jacc0 << 0, 0, 0;
    Jacc1 = rootAngAcc;
    Jacc2 << ddq2, 0, 0;
    Jacc3 << 0, ddq3, 0;
    Jacc4 << 0, ddq4, 0;
    Jacc5 << 0, ddq5, 0;
    Jacc6 << 0, ddq6, 0;
    /* #endregion */

    /* #region: Joint position wrt prev. joint */
    P01 << 0.0000, 0.0000, 0.0000;
    P12 << m*0.3102, n*0.1050, 0.0020;
    P23 << m*0.0790, n*0.0262, 0.0000;
    P34 << m*0.0000, n*0.1102, -0.2700;
    P45 << m*0.0000, n*0.0000, -0.224282;
    P56 << m*0.0000, n*0.0000, -0.107131;
    /* #endregion */

    /* #region: Link center of mass positions wrt prev. joint */
    Pc1 << 0.00145, 0.003271, 0.01792; // Link 1 (Torso COM)
    Pc2 << m*(-0.003717), n*(-0.015862), -0.000015; // Link 2
    Pc3 << m*(-0.003917), n*(0.076506), -0.013793; // Link 3
    Pc4 << m*(-0.015346), n*(-0.000102), -0.186842; // Link 4
    Pc5 << m*(0.0), n*(0.0), -0.003546; // Link 5
    /* #endregion */

    /* #region: Link mass */
    m1 = 20.516; // Torso (Joint 1-2)
    m2 = 1.681; // Hip FE&AA motors (Joint 2-3)
    m3 = 2.005; // Thigh (Joint 3-4)
    m4 = 0.21; // Shank (Joint 4-5) 
    m5 = 0.0; // Foot (Joint 5-6)
    /* #endregion */

    /* #region: Link inertia */
    Inertia1 << 0.410, 0.003, -0.012,  
                0.003, 0.908, 0.007,
                -0.012, 0.007, 1.192; // Torso (Joint 1-2)
    Inertia2 << 0.003, 9.851e-5, -1.213e-6,
                9.851e-5, 0.003, -1.327e-6,
                -1.213e-6, -1.327e-6, 0.003; // Hip FE&AA motors (Joint 2-3)
    Inertia3 << 0.008, -1.818e-4, 5.307e-4,
                -1.818e-4, 0.007, 0.001,
                5.307e-4, 0.001, 0.004; // Thigh (Joint 3-4)
    Inertia4 << 0.002, -7.048e-7, -2.682e-4,
                -7.048e-7, 0.002, -1.088e-5,
                -2.682e-4, -1.088e-5, 1.086e-4; // Shank (Joint 4-5) 
    Inertia5 << 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0; // Foot (Joint 5-6)
    /* #endregion */
    

    /*********** Outward iteration **********/
    // World (Joint 0)
    w0.setZero();
    dw0.setZero();
    dv0 = gravityVec;
    dvc0.setZero();
    // Torso (Joint 1)
    w1 = R10.transpose()*w0 + Jvel1;
    dw1 = R10.transpose()*dw0 + (R10.transpose()*w0).cross(Jvel1) + Jacc1;
    dv1 = R10.transpose()*(dw0.cross(P01) + w0.cross(w0.cross(P01)) + dv0 + rootAbsAcc);
    dvc1 = dw1.cross(Pc1) + w1.cross(w1.cross(Pc1)) + dv1;
    F1 = m1*dvc1;
    N1 = Inertia1*dw1 + w1.cross(Inertia1*w1);
    // Hip AA (Joint 2)
    w2 = R21.transpose()*w1 + Jvel2;
    dw2 = R21.transpose()*dw1 + (R21.transpose()*w1).cross(Jvel2) + Jacc2;
    dv2 = R21.transpose()*(dw1.cross(P12) + w1.cross(w1.cross(P12)) + dv1);
    dvc2 = dw2.cross(Pc2) + w2.cross(w2.cross(Pc2)) + dv2;
    F2 = m2*dvc2;
    N2 = Inertia2*dw2 + w2.cross(Inertia2*w2);
    // Hip FE (Joint 3)
    w3 = R32.transpose()*w2 + Jvel3;
    dw3 = R32.transpose()*dw2 + (R32.transpose()*w2).cross(Jvel3) + Jacc3;
    dv3 = R32.transpose()*(dw2.cross(P23) + w2.cross(w2.cross(P23)) + dv2);
    dvc3 = dw3.cross(Pc3) + w3.cross(w3.cross(Pc3)) + dv3;
    F3 = m3*dvc3;
    N3 = Inertia3*dw3 + w3.cross(Inertia3*w3);
    // Knee FE (Joint 4)
    w4 = R43.transpose()*w3 + Jvel4;
    dw4 = R43.transpose()*dw3 + (R43.transpose()*w3).cross(Jvel4) + Jacc4;
    dv4 = R43.transpose()*(dw3.cross(P34) + w3.cross(w3.cross(P34)) + dv3);
    dvc4 = dw4.cross(Pc4) + w4.cross(w4.cross(Pc4)) + dv4;
    F4 = m4*dvc4;
    N4 = Inertia4*dw4 + w4.cross(Inertia4*w4);
    // Footframe (Joint 5)
    w5 = R54.transpose()*w4 + Jvel5;
    dw5 = R54.transpose()*dw4 + (R54.transpose()*w4).cross(Jvel5) + Jacc5;
    dv5 = R54.transpose()*(dw4.cross(P45) + w4.cross(w4.cross(P45)) + dv4);
    dvc5 = dw5.cross(Pc5) + w5.cross(w5.cross(Pc5)) + dv5;
    F5 = m5*dvc5;
    N5 = Inertia5*dw5 + w5.cross(Inertia5*w5);

    /*********** Inward iteration **********/
    R60 = R10*R21*R32*R43*R54*R65;
    f6 = R60.transpose()*(-grForce); // Force acting on tip of foot (Joint 6)    
    f5 = R65*f6 + F5; // Force acting on Footframe (Joint 5) 
    f4 = R54*f5 + F4; // Force acting on Knee FE (Joint 4)
    f3 = R43*f4 + F3; // Force acting on Hip FE (Joint 3)
    f2 = R32*f3 + F2; // Force acting on Hip AA (Joint 2)
    f1 = R21*f2 + F1; // Force acting on Torso (Joint 1)
    
    n6 << 0, 0, 0; // Moment acting on tip of foot (Joint 6)    
    n5 = N5 + R65*n6 + Pc5.cross(F5) + P56.cross(R65*f6); // Moment acting on Footframe (Joint 5)
    n4 = N4 + R54*n5 + Pc4.cross(F4) + P45.cross(R54*f5); // Moment acting on Knee FE (Joint 4)
    n3 = N3 + R43*n4 + Pc3.cross(F3) + P34.cross(R43*f4); // Moment acting on Hip FE (Joint 3)
    n2 = N2 + R32*n3 + Pc2.cross(F2) + P23.cross(R32*f3); // Moment acting on Hip AA (Joint 2)
    n1 = N1 + R21*n2 + Pc1.cross(F1) + P12.cross(R21*f2); // Moment acting on Torso (Joint 1)

    JTorques << m*n2(0), n*n3(1), n*n4(1);
    return JTorques;
}

Eigen::Matrix<double, 12, 1> funNewtonEuler4Leg(Eigen::Vector3d rootAbsAcc, Eigen::Matrix3d rootOrient, Eigen::Vector3d rootAngVel, Eigen::Vector3d rootAngAcc, Eigen::VectorXd jPos, Eigen::VectorXd jVel, Eigen::VectorXd jAcc, Eigen::Vector3d grForce_LF, Eigen::Vector3d grForce_RF, Eigen::Vector3d grForce_LB, Eigen::Vector3d grForce_RB)
{
    Eigen::Vector3d jPos_LF, jPos_RF, jPos_LB, jPos_RB;
    Eigen::Vector3d jVel_LF, jVel_RF, jVel_LB, jVel_RB;
    Eigen::Vector3d jAcc_LF, jAcc_RF, jAcc_LB, jAcc_RB;
    Eigen::Vector3d Tau_LF, Tau_RF, Tau_LB, Tau_RB;
    Eigen::VectorXd jTorques(12);
    Eigen::Matrix3d Rfoot_LF, Rfoot_RF, Rfoot_LB, Rfoot_RB;
    Eigen::Vector3d Rcon_LF, Rcon_RF, Rcon_LB, Rcon_RB;

    // Joint Positions
    jPos_LF << jPos(0), jPos(1), jPos(2);
    jPos_RF << jPos(3), jPos(4), jPos(5);
    jPos_LB << jPos(6), jPos(7), jPos(8);
    jPos_RB << jPos(9), jPos(10), jPos(11);
    // Joint Velocities
    jVel_LF << jVel(0), jVel(1), jVel(2);
    jVel_RF << jVel(3), jVel(4), jVel(5);
    jVel_LB << jVel(6), jVel(7), jVel(8);
    jVel_RB << jVel(9), jVel(10), jVel(11);
    // Joint Accelerations
    jAcc_LF << jAcc(0), jAcc(1), jAcc(2);
    jAcc_RF << jAcc(3), jAcc(4), jAcc(5);
    jAcc_LB << jAcc(6), jAcc(7), jAcc(8);
    jAcc_RB << jAcc(9), jAcc(10), jAcc(11);

    Rfoot_LF = rootOrient * RotateRoll(jPos_LF(0)) * RotatePitch(jPos_LF(1)) * RotatePitch(jPos_LF(2)) * RotatePitch(30 * PI / 180);
    Rfoot_RF = rootOrient * RotateRoll(jPos_RF(0)) * RotatePitch(jPos_RF(1)) * RotatePitch(jPos_RF(2)) * RotatePitch(-30 * PI / 180);
    Rfoot_LB = rootOrient * RotateRoll(jPos_LB(0)) * RotatePitch(jPos_LB(1)) * RotatePitch(jPos_LB(2)) * RotatePitch(30 * PI / 180);
    Rfoot_RB = rootOrient * RotateRoll(jPos_RB(0)) * RotatePitch(jPos_RB(1)) * RotatePitch(jPos_RB(2)) * RotatePitch(-30 * PI / 180);

    Tau_LF = funcNewtonEuler2(rootAbsAcc, rootOrient, rootAngVel, rootAngAcc, jPos_LF, jVel_LF, jAcc_LF, grForce_LF, Rcon_LF, 1);
    Tau_RF = funcNewtonEuler2(rootAbsAcc, rootOrient, rootAngVel, rootAngAcc, jPos_RF, jVel_RF, jAcc_RF, grForce_RF, Rcon_RF, 2);
    Tau_LB = funcNewtonEuler2(rootAbsAcc, rootOrient, rootAngVel, rootAngAcc, jPos_LB, jVel_LB, jAcc_LB, grForce_LB, Rcon_LB, 3);
    Tau_RB = funcNewtonEuler2(rootAbsAcc, rootOrient, rootAngVel, rootAngAcc, jPos_RB, jVel_RB, jAcc_RB, grForce_RB, Rcon_RB, 4);

    jTorques << Tau_LF(0), Tau_LF(1), Tau_LF(2), Tau_RF(0), Tau_RF(1), Tau_RF(2), Tau_LB(0), Tau_LB(1), Tau_LB(2), Tau_RB(0), Tau_RB(1), Tau_RB(2);

    return jTorques;
}
