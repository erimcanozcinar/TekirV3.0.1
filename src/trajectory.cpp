#include "trajectory.hpp"


/* Functions */
void controller::xboxController(ALLEGRO_EVENT_QUEUE *event_queue, ALLEGRO_EVENT event, ALLEGRO_JOYSTICK* joyStick)
{
    
    /* Gamepad routines */
    al_wait_for_event_timed(event_queue, &event, 0.0001);
    switch (event.type) {
        case ALLEGRO_EVENT_JOYSTICK_CONFIGURATION:
            printf("Joystick configuration changed.\n");
            break;
        case ALLEGRO_EVENT_JOYSTICK_AXIS:
            // Handle gamepad axis events
            switch (event.joystick.stick) {
                case 0: /* Left Stick */
                    switch (event.joystick.axis) {
                        case 0: /* y axis */
                            Vy_mean = -0.15*event.joystick.pos;
                            Vy_mean = floor(Vy_mean*10000)/10000; /* Truncate to 3 decimal */
                            if(abs(Vy_mean)<0.03) { Vy_mean = 0.0; }
                            break;
                        case 1: /* x axis */
                            Vx_mean = -Kv*event.joystick.pos;                                
                            Vx_mean = floor(Vx_mean*10000)/10000; /* Truncate to 3 decimal */
                            if(abs(Vx_mean)<0.03) { Vx_mean = 0.0; }
                            break;
                    }
                    break;                    
                case 1: /* Right Stick */
                    switch (event.joystick.axis) {
                        case 0: /* y axis */
                            std::cout << event.joystick.axis << ": " << event.joystick.pos << std::endl;
                            break;
                        case 1: /* x axis */
                            std::cout << event.joystick.axis << ": " << event.joystick.pos << std::endl;
                            break;
                    }
                    break; 
                case 2: /* D-pad */
                    switch (event.joystick.axis) {
                        case 0: /* D-pad x axis */                    
                            break;                    
                        case 1: /* D-pad y axis */
                            if(event.joystick.pos == -1) /* D-pad Up */
                            {
                                if(Kv > 0.6){Kv = 0.6;}                                
                                else { Kv = Kv + 0.2; }
                                std::cout << "Max. Vel. Set to:" << Kv << std::endl;
                            }
                            else if(event.joystick.pos == 1) /* D-pad Down */
                            {
                                if(Kv <= 0.2){Kv = 0.2;}                                
                                else { Kv = Kv - 0.2; }
                                std::cout << "Max. Vel. Set to:" << Kv << std::endl;
                            }
                            else
                            {
                                Kv = Kv;
                            }
                            break;
                    }
                    break; 
            }
            break;
        case ALLEGRO_EVENT_JOYSTICK_BUTTON_DOWN:
            switch (event.joystick.button){
                case 0:
                    /* Y Button */
                    std::cout << "Y: " << event.joystick.button << std::endl;
                    break;
                case 1:
                    /* B Button */
                    std::cout << "B: " << event.joystick.button << std::endl;
                    break;
                case 2:
                    /* A Button */
                    std::cout << "A: " << event.joystick.button << std::endl;
                    break;
                case 3:
                    /* X Button */
                    std::cout << "X: " << event.joystick.button << std::endl;
                    break;
                case 4:
                    /* LB Button */
                    std::cout << "LB: " << event.joystick.button << std::endl;
                    break;
                case 5:
                    /* RB Button */
                    std::cout << "RB: " << event.joystick.button << std::endl;
                    break;
                case 6:
                    /* LT Button */
                    std::cout << "LT: " << event.joystick.button << std::endl;
                    break;
                case 7:
                    /* RT Button */
                    std::cout << "RT: " << event.joystick.button << std::endl;
                    break;
                case 8:
                    /* Back Button */
                    std::cout << "Back: " << event.joystick.button << std::endl;
                    close = true;
                    break;
                case 9:
                    /* Start Button */
                    if(walkEnable){ 
                        walkEnable = false; 
                    }
                    else{ 
                        walkEnable = true; 
                    }            
                    break;
                case 10:
                    /* Back Button */
                    std::cout << "LSB: " << event.joystick.button << std::endl;
                    break;
                case 11:
                    /* Start Button */
                    std::cout << "RSB: " << event.joystick.button << std::endl;
                    break;
                
                default:
                    break;
            }
            break;
        case ALLEGRO_EVENT_JOYSTICK_BUTTON_UP:
            // Handle gamepad button up events
            // event.joystick.button contains the button ID
            // Add your game logic here
            break;
        // Add more cases for other event types if needed
        default:
            // Handle other event types or do nothing
            break;
    }
}

void controller::dualShockController(ALLEGRO_EVENT_QUEUE *event_queue, ALLEGRO_EVENT ev)
{
    ALLEGRO_JOYSTICK *joyStick = al_get_joystick(0);
    ALLEGRO_JOYSTICK_STATE joyState;
    if(joyStick) {al_get_joystick_state(joyStick, &joyState);}   
    

    al_wait_for_event_timed(event_queue, &ev, 0.0001);
    switch (ev.type) {
        case ALLEGRO_EVENT_JOYSTICK_CONFIGURATION:
            al_reconfigure_joysticks();
        case ALLEGRO_EVENT_JOYSTICK_AXIS:
            switch (ev.joystick.stick) {
                case 0: /* Left-Stick */                        
                    switch (ev.joystick.axis) {
                        case 0: /* y-axis */
                            // std::cout << "Left stick y-axis: " << ev.joystick.pos << std::endl;
                            Vy_mean = -0.15*ev.joystick.pos;
                            Vy_mean = floor(Vy_mean*10000)/10000; /* Truncate to 3 decimal */
                            if(abs(Vy_mean)<0.03) { Vy_mean = 0.0; }
                            joyCmd[1] = Vy_mean;
                            break;
                        case 1: /* x-axis */
                            // std::cout << "Left stick x-axis: " << ev.joystick.pos << std::endl;
                            Vx_mean = -Kv*ev.joystick.pos;                                
                            Vx_mean = floor(Vx_mean*10000)/10000; /* Truncate to 3 decimal */
                            if(abs(Vx_mean)<0.03) { Vx_mean = 0.0; }
                            joyCmd[0] = Vx_mean;
                            break;
                    }
                    break;                    
                case 1: /* Right stick y-axis and L2 axis */
                    switch (ev.joystick.axis) {
                        case 0: /* L2 axis */
                            cmdZc = mapVal(ev.joystick.pos, 1.0, -1.0, MIN_BODY_HEIGHT, initZc);
                            joyCmd[5] = cmdZc;
                            break;
                        case 1:
                            // std::cout << "Right stick y-axis: " << ev.joystick.pos << std::endl;
                            if(walkEnable) { 
                                yaw = -0.0873*ev.joystick.pos;
                                roll = 0.0; 
                            }
                            else {
                                if(joyState.button[4]) { roll = 0.274*ev.joystick.pos; }
                                else { yaw = -0.35*ev.joystick.pos; }                                 
                            }
                            if(abs(yaw) < 0.03) { yaw = 0.0;}
                            if(abs(roll) < 0.03) { roll = 0.0;}
                            joyCmd[2] = yaw;
                            joyCmd[4] = roll;
                            break;
                    }
                    break;
                case 2: /* Right stick x-axis and R2 axis */
                    switch (ev.joystick.axis) {
                        case 0: /* Right stick x-axis */
                            // std::cout << "Right stick x-axis: " << ev.joystick.pos << std::endl;
                            if(!walkEnable) { 
                                pitch = 0.1*ev.joystick.pos;
                            }                               
                            else { 
                                pitch = 0.0;
                            }
                            if(abs(pitch) < 0.03) { pitch = 0.0;}
                            joyCmd[3] = pitch;
                            break;
                        case 1: /* R2 axis */
                            cmdZc = mapVal(ev.joystick.pos, 1.0, -1.0, MAX_BODY_HEIGHT, initZc);
                            joyCmd[5] = cmdZc;
                            break;
                    }
                    break;
                case 3: /* D-pad */
                    switch (ev.joystick.axis) {
                        case 0: /* D-pad R/L */
                            // std::cout << "D-pad R/L: " << ev.joystick.pos << std::endl;
                            break;
                        case 1: /* D-pad U/D */
                            Kv = Kv - 0.1*ev.joystick.pos;
                            if(Kv >= 1.0) { Kv = 1.0;}
                            else if(Kv <= 0.1) { Kv = 0.1; }
                            if(abs(ev.joystick.pos)) { std::cout << "Max. Vel. set to: " << Kv << std::endl; }                                              
                            break;
                    }
                    break;
            }
            break;
        case ALLEGRO_EVENT_JOYSTICK_BUTTON_DOWN:
            switch (ev.joystick.button) {
                // case 0: /* X Button */
                //     std::cout << "X: " << ev.joystick.button << std::endl;
                //     break;
                // case 1: /* Circle Button */
                //     std::cout << "Circle: " << ev.joystick.button << std::endl;
                //     break;
                // case 2: /* Triangle Button */
                //     std::cout << "Triangle: " << ev.joystick.button << std::endl;
                //     break;
                // case 3: /* Square Button */
                //     std::cout << "Square: " << ev.joystick.button << std::endl;
                //     break;
                // case 4: /* L1 Button */
                //     std::cout << "L1: " << ev.joystick.button << std::endl;                                        
                //     break;
                // case 5: /* R1 Button */
                //     std::cout << "R1: " << ev.joystick.button << std::endl;
                //     break;
                case 6: /* L2 Button */
                    // std::cout << "L2: " << ev.joystick.button << std::endl;
                    break;
                case 7: /* R2 Button */
                    // std::cout << "R2: " << ev.joystick.button << std::endl;
                    break;
                case 8: /* Back Button */
                    // std::cout << "Share: " << ev.joystick.button << std::endl;
                    close = true;
                    al_release_joystick(joyStick);
                    break;
                case 9: /* Options */
                    if(walkEnable){ 
                        walkEnable = false;
                        std::cout << "Walking mode disabled!" << std::endl; 
                    }
                    else{ 
                        walkEnable = true;
                        std::cout << "Walking mode enabled!" << std::endl;
                    } 
                    break;
                // case 10: /* PS Button */
                //     std::cout << "PS: " << ev.joystick.button << std::endl;
                //     break;
                // case 11: /* LS Button */
                //     std::cout << "LS: " << ev.joystick.button << std::endl;
                //     break;
                // case 12: /* RS Button */
                //     std::cout << "RS: " << ev.joystick.button << std::endl;
                //     break;
                
                default:
                    break;
            }
            break;
        case ALLEGRO_EVENT_JOYSTICK_BUTTON_UP:
            // Handle gamepad button up events
            // ev.joystick.button contains the button ID
            // std::cout << ev.joystick.button << std::endl;
            // Add your game logic here
            break;
        // Add more cases for other event types if needed
        default:
            // Handle other event types or do nothing
            break;
    }
}

void trajectory::comTrajectory(double RealTime, double Ts, double Td, int Nphase, double px, double py, double vx_mean, double vy_mean, double yaw, double Cz, double Fh, double dt)
{    
    double t0 = 0.0; // Beginning of the Universe(or Simulation)
    double t1 = t0 + 0; // Put robot on the ground
    double t2 = t1 + Ts; // First half step
    double tstart_w = t2 + Td; // Walking start
    double tend_w = tstart_w + (Ts + Td) * Nphase; // Walking end
    double t3 = tend_w + Ts; // CoM stoping phase

    double g = 9.81; // Gravitational acceleration
    double w = sqrt(g / Cz); // Natural freq of equivalent pendulum

    Eigen::Vector3d trajComX, trajComY, trajComYaw;
    Eigen::VectorXd returnVals(16);

    comVel << vx_mean, vy_mean, 0;
    comVel = RotateYaw(yaw)*comVel;   

    /* #region: Phase Calculation Start(Each Ts + Td is one phase) */
    int k, kx;
    double Rt4Fr = (RealTime - tstart_w) / (Ts + Td);
    if (FuncInterval(RealTime, t0, tstart_w, dt) == true)
    {
        k = 0;
        kx = 0;
    }
    else if (FuncGreater(RealTime, tend_w, dt) == true)
    {
        k = Nphase - 1;
        kx = floor(k / 2);
    }
    else
    {
        k = floor(Rt4Fr);
        kx = floor(Rt4Fr / 2);
    }
    /* #endregion: Phase Calculation Start(Each Ts + Td is one phase) */
        
    /* #region: Xcom parameters */
    double Cxd = (px + comVel(0) * Ts / 2);
    double Cx0 = (2 * px - Cxd);
    double dCx0 = w * (px - Cx0) * 1 / tanh(w * Ts / 2);
    double ddCx0 = w * w * (Cx0 - px);
    double d_Cx0 = Cxd;
    double d_dCx0 = dCx0;

    double Kx = d_dCx0 + w * (Cxd - px) * 1 / tanh(w * Td / 2);
    Strx = 2 * Td * Kx;
    double Cxd0 = Cx0 + Strx / 2;
    /* #endregion: Xcom parameters */

    /* #region: Ycom parameters */
    double Cyd = (py + comVel(1) * Ts / 2);
    double Cy0 = (2 * py - Cyd);
    double dCy0 = w * (py - Cy0) * 1 / tanh(w * Ts / 2);
    double ddCy0 = w * w * (Cy0 - py);
    double d_Cy0 = Cyd;
    double d_dCy0 = dCy0;

    double Ky = d_dCy0 + w * (Cyd - py) * 1 / tanh(w * Td / 2);
    Stry = 2 * Td * Ky;
    double Cyd0 = Cy0 + Stry / 2;
    /* #endregion: Ycom parameters */
    
    /* #region: Yaw */
    offsetPf_LF << + Pfx_offset, + Pfy_offset + LatOut, Pfz_offset;
    offsetPf_RF << + Pfx_offset, - Pfy_offset - LatOut, Pfz_offset;
    offsetPf_LB << - Pfx_offset, + Pfy_offset + LatOut, Pfz_offset;
    offsetPf_RB << - Pfx_offset, - Pfy_offset - LatOut, Pfz_offset;
    yawStr_LF = (RotateYaw(yaw*(kx+1))*offsetPf_LF - offsetPf_LF);
    yawStr_RF = (RotateYaw(yaw*(kx+1))*offsetPf_RF - offsetPf_RF);
    yawStr_LB = (RotateYaw(yaw*(kx+1))*offsetPf_LB - offsetPf_LB);
    yawStr_RB = (RotateYaw(yaw*(kx+1))*offsetPf_RB - offsetPf_RB);
    pre_yawStr_LF = (RotateYaw(yaw*(kx))*offsetPf_LF - offsetPf_LF);
    pre_yawStr_RF = (RotateYaw(yaw*(kx))*offsetPf_RF - offsetPf_RF);
    pre_yawStr_LB = (RotateYaw(yaw*(kx))*offsetPf_LB - offsetPf_LB);
    pre_yawStr_RB = (RotateYaw(yaw*(kx))*offsetPf_RB - offsetPf_RB);
    /*  #endregion */
    
    /* #region: CoM Trajectory */
    if (FuncInterval(RealTime, t0, t1, dt) == true) // Put robot on the ground
    {
        Cx = 0; dCx = 0; ddCx = 0;
        Cy = 0; dCy = 0; ddCy = 0;
        trajComYaw.setZero();
    }
    else if (FuncInterval(RealTime, t1, t2, dt) == true) // Initialize CoM position (First half step)
    {
        trajComX = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, Cxd, dCx0, -ddCx0);
        trajComY = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, Cyd, dCy0, -ddCy0);
        Cx = trajComX(0); dCx = trajComX(1); ddCx = trajComX(2);
        Cy = trajComY(0); dCy = trajComY(1); ddCy = trajComY(2);
        trajComYaw.setZero();
    }
    else if (FuncInterval(RealTime, t2, tstart_w, dt) == true) // Initialite CoM trajectory (Just before walking starts)
    {
        trajComX = FuncPoly5th(RealTime, t2, tstart_w, Cxd, dCx0, -ddCx0, Cxd0, dCx0, ddCx0);
        trajComY = FuncPoly5th(RealTime, t2, tstart_w, Cyd, dCy0, -ddCy0, Cyd0, dCy0, ddCy0);
        Cx = trajComX(0); dCx = trajComX(1); ddCx = trajComX(2);
        Cy = trajComY(0); dCy = trajComY(1); ddCy = trajComY(2);
        trajComYaw = FuncPoly5th(RealTime, t2, tstart_w, 0, 0, 0, 0, 0, 0);
    }
    else if (FuncInterval(RealTime, tstart_w, tend_w, dt) == true) // Start walking
    {
        double Rt = RealTime - tstart_w - (Ts + Td) * (k);
        if (FuncInterval(Rt, 0, Ts, dt) == true) // Single Support Phase
        {
            double s_T = Rt;
            Cx = (Cx0 - px) * cosh(w * s_T) + (dCx0 / w) * sinh(w * s_T) + px;
            dCx = w * (Cx0 - px) * sinh(w * s_T) + dCx0 * cosh(w * s_T);
            ddCx = w * w * (Cx0 - px) * cosh(w * s_T) + w * dCx0 * sinh(w * s_T);
            Cy = (Cy0 - py) * cosh(w * s_T) + (dCy0 / w) * sinh(w * s_T) + py;
            dCy = w * (Cy0 - py) * sinh(w * s_T) + dCy0 * cosh(w * s_T);
            ddCy = w * w * (Cy0 - py) * cosh(w * s_T) + w * dCy0 * sinh(w * s_T);
            Cx = Cx + (k + 1) * Strx / 2;
            Cy = Cy + (k + 1) * Stry / 2;
            trajComYaw << yaw*k; 0; 0;
            Cyaw = trajComYaw(0);
        }
        else if (FuncInterval(Rt, Ts, Ts + Td, dt) == true) // Double Support Phase
        {
            double d_T = Rt - Ts;
            Cx = (d_Cx0 - px) * cosh(w * d_T) + ((d_dCx0 - Kx) / w) * sinh(w * d_T) + Kx * d_T + px;
            dCx = w * (d_Cx0 - px) * sinh(w * d_T) + (d_dCx0 - Kx) * cosh(w * d_T) + Kx;
            ddCx = w * w * (d_Cx0 - px) * cosh(w * d_T) + w * (d_dCx0 - Kx) * sinh(w * d_T);
            Cy = (d_Cy0 - py) * cosh(w * d_T) + ((d_dCy0 - Ky) / w) * sinh(w * d_T) + Ky * d_T + py;
            dCy = w * (d_Cy0 - py) * sinh(w * d_T) + (d_dCy0 - Ky) * cosh(w * d_T) + Ky;
            ddCy = w * w * (d_Cy0 - py) * cosh(w * d_T) + w * (d_dCy0 - Ky) * sinh(w * d_T);
            Cx = Cx + (k + 1) * Strx / 2;
            Cy = Cy + (k + 1) * Stry / 2;
            
            trajComYaw = FuncPoly5th(Rt, Ts, Ts+Td, yaw*k, 0, 0, yaw*(k+1), 0, 0);
            Cyaw = trajComYaw(0);
        }
    }
    else if (FuncInterval(RealTime, tend_w, t3, dt) == true) // Stop CoM trajectory(After walking ends)
    {
        double PosX_start = Cxd0 + (k + 1) * Strx / 2;
        double PosX_end = Strx / 2 + (k + 1) * Strx / 2;
        double PosY_start = Cyd0 + (k + 1) * Stry / 2;
        double PosY_end = Stry / 2 + (k + 1) * Stry / 2;
        trajComX = FuncPoly5th(RealTime,  tend_w,  t3,  PosX_start,  dCx0,  ddCx0, PosX_end, 0, 0);
        trajComY = FuncPoly5th(RealTime,  tend_w,  t3,  PosY_start,  dCy0,  ddCy0, PosY_end, 0, 0);
        Cx = trajComX(0); dCx = trajComX(1); ddCx = trajComX(2);
        Cy = trajComY(0); dCy = trajComY(1); ddCy = trajComY(2);

        trajComYaw = FuncPoly5th(RealTime, tend_w, t3, yaw*(k+1), 0, 0, yaw*(k+1), 0, 0);
        Cyaw = trajComYaw(0);
    }
    else
    {
        Cx = Strx / 2 + (k + 1) * Strx / 2; dCx = 0; ddCx = 0;
        Cy = Stry / 2 + (k + 1) * Stry / 2; dCy = 0; ddCy = 0;
        trajComYaw << yaw*(k+1); 0; 0;
        Cyaw = trajComYaw(0);
    }    

    Xcop = Cx - Cz * ddCx / g;
    Ycop = Cy - Cz * ddCy / g;
    /* #endregion: CoM Trajectory End */

    /* #region: Feet Trajectory */
    if (FuncInterval(RealTime, t0, t1, dt) == true) // Put robot on the ground
    {
        Footx_LF.setZero(); Footy_LF.setZero(); Footz_LF.setZero();       
        Footx_RF.setZero(); Footy_RF.setZero(); Footz_RF.setZero();
        Footx_LB.setZero(); Footy_LB.setZero(); Footz_LB.setZero();       
        Footx_RB.setZero(); Footy_RB.setZero(); Footz_RB.setZero();       
    }
    else if (FuncInterval(RealTime, t1, t2, dt) == true) // Feet trajectory initialization(First half step)
    {
        // Footx_LF = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, Strx*0.5, 0, 0);
        // Footy_LF = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, Stry*0.5, 0, 0);
        // Footz_LF = FuncPoly6th(RealTime, t1, t2, 0, 0, 0, 0, 0, 0, Fh);
        Footx_LF = FuncPoly3rd(RealTime, t1, t2, 0, 0, Strx*0.5, 0);
        Footy_LF = FuncPoly3rd(RealTime, t1, t2, 0, 0, Stry*0.5, 0);
        Footz_LF = FuncPoly4th(RealTime, t1, t2, 0, 0, 0, 0, Fh);

        Footx_RF.setZero(); Footy_RF.setZero(); Footz_RF.setZero();
        Footx_LB.setZero(); Footy_LB.setZero(); Footz_LB.setZero();

        // Footx_RB = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, Strx*0.5, 0, 0);
        // Footy_RB = FuncPoly5th(RealTime, t1, t2, 0, 0, 0, Stry*0.5, 0, 0);
        // Footz_RB = FuncPoly6th(RealTime, t1, t2, 0, 0, 0, 0, 0, 0, Fh); 
        Footx_RB = FuncPoly3rd(RealTime, t1, t2, 0, 0, Strx*0.5, 0);
        Footy_RB = FuncPoly3rd(RealTime, t1, t2, 0, 0, Stry*0.5, 0);
        Footz_RB = FuncPoly4th(RealTime, t1, t2, 0, 0, 0, 0, Fh);                
    }
    else if (FuncInterval(RealTime, t2, tstart_w, dt) == true)
    {
        Footx_LF << Strx*0.5, 0, 0;
        Footy_LF << Stry*0.5, 0, 0;
        Footz_LF << 0, 0, 0;

        Footx_RF.setZero(); Footy_RF.setZero(); Footz_RF.setZero();
        Footx_LB.setZero(); Footy_LB.setZero(); Footz_LB.setZero();

        Footx_RB << Strx*0.5, 0, 0;
        Footy_RB << Stry*0.5, 0, 0;
        Footz_RB << 0, 0, 0;        
    }
    else if (FuncInterval(RealTime, tstart_w, tend_w, dt) == true) // Feet trajectory after initialization
    {
        double ts = tstart_w + (Ts + Td) * k;
        if (((k + 1) % 2) == 0) // Left foot swing, right foot stand
        {
            Pstart_LF << Strx*0.5 + Strx*kx + pre_yawStr_LF(0), Stry*0.5 + Stry*kx + pre_yawStr_LF(1), 0;
            Pend_LF << Strx*0.5 + Strx*(kx + 1) + yawStr_LF(0), Stry*0.5 + Stry*(kx + 1) + yawStr_LF(1), 0;
            // Footx_LF = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_LF(0), 0, 0, Pend_LF(0), 0, 0);
            // Footy_LF = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_LF(1), 0, 0, Pend_LF(1), 0, 0);
            // Footz_LF = FuncPoly6th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0, 0, Fh);
            Footx_LF = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_LF(0), 0, Pend_LF(0), 0);
            Footy_LF = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_LF(1), 0, Pend_LF(1), 0);
            Footz_LF = FuncPoly4th(RealTime, ts, ts + Ts, 0, 0, 0, 0, Fh);

            Pstart_RF << Strx*(kx + 1) + yawStr_RF(0), Stry*(kx + 1) + yawStr_RF(1), 0;
            Pend_RF << Strx*(kx + 1) + yawStr_RF(0), Stry*(kx + 1) + yawStr_RF(1), 0;
            // Footx_RF = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_RF(0), 0, 0, Pend_RF(0), 0, 0);
            // Footy_RF = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_RF(1), 0, 0, Pend_RF(1), 0, 0);
            // Footz_RF = FuncPoly6th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0, 0, 0);
            Footx_RF = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_RF(0), 0, Pend_RF(0), 0);
            Footy_RF = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_RF(1), 0, Pend_RF(1), 0);
            Footz_RF = FuncPoly4th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0);

            Pstart_LB << Strx*(kx + 1) + yawStr_LB(0), Stry*(kx + 1) + yawStr_LB(1), 0;
            Pend_LB << Strx*(kx + 1) + yawStr_LB(0), Stry*(kx + 1) + yawStr_LB(1), 0;
            // Footx_LB = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_LB(0), 0, 0, Pend_LB(0), 0, 0);
            // Footy_LB = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_LB(1), 0, 0, Pend_LB(1), 0, 0);
            // Footz_LB = FuncPoly6th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0, 0, 0);
            Footx_LB = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_LB(0), 0, Pend_LB(0), 0);
            Footy_LB = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_LB(1), 0, Pend_LB(1), 0);
            Footz_LB = FuncPoly4th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0);

            Pstart_RB << Strx*0.5 + Strx*kx + pre_yawStr_RB(0), Stry*0.5 + Stry*kx + pre_yawStr_RB(1), 0;
            Pend_RB << Strx*0.5 + Strx*(kx + 1) + yawStr_RB(0), Stry*0.5 + Stry*(kx + 1) + yawStr_RB(1), 0;
            // Footx_RB = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_RB(0), 0, 0, Pend_RB(0), 0, 0);
            // Footy_RB = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_RB(1), 0, 0, Pend_RB(1), 0, 0);
            // Footz_RB = FuncPoly6th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0, 0, Fh);
            Footx_RB = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_RB(0), 0, Pend_RB(0), 0);
            Footy_RB = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_RB(1), 0, Pend_RB(1), 0);
            Footz_RB = FuncPoly4th(RealTime, ts, ts + Ts, 0, 0, 0, 0, Fh);
        }
        else // Right foot swing, left foot stand
        {
            Pstart_LF << Strx*0.5 + Strx*kx + pre_yawStr_LF(0), Stry*0.5 + Stry*kx + pre_yawStr_LF(1), 0;
            Pend_LF << Strx*0.5 + Strx*kx + pre_yawStr_LF(0), Stry*0.5 + Stry*kx + pre_yawStr_LF(1), 0;
            // Footx_LF = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_LF(0), 0, 0, Pend_LF(0), 0, 0);
            // Footy_LF = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_LF(1), 0, 0, Pend_LF(1), 0, 0);
            // Footz_LF = FuncPoly6th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0, 0, 0);
            Footx_LF = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_LF(0), 0, Pend_LF(0), 0);
            Footy_LF = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_LF(1), 0, Pend_LF(1), 0);
            Footz_LF = FuncPoly4th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0);

            Pstart_RF << Strx*kx + pre_yawStr_RF(0), Stry*kx + pre_yawStr_RF(1), 0;
            Pend_RF << Strx*(kx + 1) + yawStr_RF(0), Stry*(kx + 1) + yawStr_RF(1), 0;
            // Footx_RF = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_RF(0), 0, 0, Pend_RF(0), 0, 0);
            // Footy_RF = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_RF(1), 0, 0, Pend_RF(1), 0, 0);
            // Footz_RF = FuncPoly6th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0, 0, Fh);
            Footx_RF = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_RF(0), 0, Pend_RF(0), 0);
            Footy_RF = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_RF(1), 0, Pend_RF(1), 0);
            Footz_RF = FuncPoly4th(RealTime, ts, ts + Ts, 0, 0, 0, 0, Fh);

            Pstart_LB << Strx*kx + pre_yawStr_LB(0), Stry*kx + pre_yawStr_LB(1), 0;
            Pend_LB << Strx*(kx + 1) + yawStr_LB(0), Stry*(kx + 1) + yawStr_LB(1), 0;
            // Footx_LB = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_LB(0), 0, 0, Pend_LB(0), 0, 0);
            // Footy_LB = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_LB(1), 0, 0, Pend_LB(1), 0, 0);
            // Footz_LB = FuncPoly6th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0, 0, Fh);
            Footx_LB = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_LB(0), 0, Pend_LB(0), 0);
            Footy_LB = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_LB(1), 0, Pend_LB(1), 0);
            Footz_LB = FuncPoly4th(RealTime, ts, ts + Ts, 0, 0, 0, 0, Fh);

            Pstart_RB << Strx*0.5 + Strx*kx + pre_yawStr_RB(0), Stry*0.5 + Stry*kx + pre_yawStr_RB(1), 0;
            Pend_RB << Strx*0.5 + Strx*kx + pre_yawStr_RB(0), Stry*0.5 + Stry*kx + pre_yawStr_RB(1), 0;
            // Footx_RB = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_RB(0), 0, 0, Pend_RB(0), 0, 0);
            // Footy_RB = FuncPoly5th(RealTime, ts, ts + Ts, Pstart_RB(1), 0, 0, Pend_RB(1), 0, 0);
            // Footz_RB = FuncPoly6th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0, 0, 0);
            Footx_RB = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_RB(0), 0, Pend_RB(0), 0);
            Footy_RB = FuncPoly3rd(RealTime, ts, ts + Ts, Pstart_RB(1), 0, Pend_RB(1), 0);
            Footz_RB = FuncPoly4th(RealTime, ts, ts + Ts, 0, 0, 0, 0, 0);
        }

    }
    else if (FuncInterval(RealTime, tend_w, t3, dt) == true)
    {
        if (((k + 1) % 2) == 0)
        {
            Pend_LF << Strx*0.5 + Strx*(kx + 1) + yawStr_LF(0), Stry*0.5 + Stry*(kx + 1) + yawStr_LF(1), 0;
            Footx_LF << Pend_LF(0), 0, 0;
            Footy_LF << Pend_LF(1), 0, 0;
            Footy_LF << 0, 0, 0;

            Pstart_RF << Strx*(kx + 1) + yawStr_RF(0), Stry*(kx + 1) + yawStr_RF(1), 0;
            Pend_RF << Strx*0.5 + Strx*(kx + 1) + yawStr_RF(0), Stry*0.5 + Stry*(kx + 1) + yawStr_RF(1), 0;
            // Footx_RF = FuncPoly5th(RealTime, tend_w, t3, Pstart_RF(0), 0, 0, Pend_RF(0), 0, 0);
            // Footy_RF = FuncPoly5th(RealTime, tend_w, t3, Pstart_RF(1), 0, 0, Pend_RF(1), 0, 0);
            // Footz_RF = FuncPoly6th(RealTime, tend_w, t3, 0, 0, 0, 0, 0, 0, Fh);
            Footx_RF = FuncPoly3rd(RealTime, tend_w, t3, Pstart_RF(0), 0, Pend_RF(0), 0);
            Footy_RF = FuncPoly3rd(RealTime, tend_w, t3, Pstart_RF(1), 0, Pend_RF(1), 0);
            Footz_RF = FuncPoly4th(RealTime, tend_w, t3, 0, 0, 0, 0, Fh);

            Pstart_LB << Strx*(kx + 1) + yawStr_LB(0), Stry*(kx + 1) + yawStr_LB(1), 0;
            Pend_LB << Strx*0.5 + Strx*(kx + 1) + yawStr_LB(0), Stry*0.5 + Stry*(kx + 1) + yawStr_LB(1), 0;
            // Footx_LB = FuncPoly5th(RealTime, tend_w, t3, Pstart_LB(0), 0, 0, Pend_LB(0), 0, 0);
            // Footy_LB = FuncPoly5th(RealTime, tend_w, t3, Pstart_LB(1), 0, 0, Pend_LB(1), 0, 0);
            // Footz_LB = FuncPoly6th(RealTime, tend_w, t3, 0, 0, 0, 0, 0, 0, Fh);
            Footx_LB = FuncPoly3rd(RealTime, tend_w, t3, Pstart_LB(0), 0, Pend_LB(0), 0);
            Footy_LB = FuncPoly3rd(RealTime, tend_w, t3, Pstart_LB(1), 0, Pend_LB(1), 0);
            Footz_LB = FuncPoly4th(RealTime, tend_w, t3, 0, 0, 0, 0, Fh);

            Pend_RB << Strx*0.5 + Strx*(kx + 1) + yawStr_RB(0), Stry*0.5 + Stry*(kx + 1) + yawStr_RB(1), 0;
            Footx_RB << Pend_RB(1), 0, 0;
            Footy_RB << Pend_RB(1), 0, 0;
            Footy_RB << 0, 0, 0;
        }
        else
        {
            Pstart_LF << Strx*0.5 + Strx*kx + pre_yawStr_LF(0), Stry*0.5 + Stry*kx + pre_yawStr_LF(1), 0;
            Pend_LF << Strx + Strx*kx + yawStr_LF(0), Stry + Stry*kx + yawStr_LF(1), 0;
            // Footx_LF = FuncPoly5th(RealTime, tend_w, t3, Pstart_LF(0), 0, 0, Pend_LF(0), 0, 0);
            // Footy_LF = FuncPoly5th(RealTime, tend_w, t3, Pstart_LF(1), 0, 0, Pend_LF(1), 0, 0);
            // Footz_LF = FuncPoly6th(RealTime, tend_w, t3, 0, 0, 0, 0, 0, 0, Fh);
            Footx_LF = FuncPoly3rd(RealTime, tend_w, t3, Pstart_LF(0), 0, Pend_LF(0), 0);
            Footy_LF = FuncPoly3rd(RealTime, tend_w, t3, Pstart_LF(1), 0, Pend_LF(1), 0);
            Footz_LF = FuncPoly4th(RealTime, tend_w, t3, 0, 0, 0, 0, Fh);

            Pend_RF << Strx*(kx + 1) + yawStr_RF(0), Stry*(kx + 1) + yawStr_RF(1), 0;
            Footx_RF << Pend_RF(0), 0, 0;
            Footy_RF << Pend_RF(1), 0, 0;
            Footz_RF << 0, 0, 0;

            Pend_LB << Strx*(kx + 1) + yawStr_LB(0), Stry*(kx + 1) + yawStr_LB(1), 0;
            Footx_LB << Pend_LB(0), 0, 0;
            Footy_LB << Pend_LB(1), 0, 0;
            Footz_LB << 0, 0, 0;

            Pstart_RB << Strx*0.5 + Strx*kx + pre_yawStr_RB(0), Stry*0.5 + Stry*kx + pre_yawStr_RB(1), 0;
            Pend_RB << Strx + Strx*kx + yawStr_RB(0), Stry + Stry*kx + yawStr_RB(1), 0;
            // Footx_RB = FuncPoly5th(RealTime, tend_w, t3, Pstart_RB(0), 0, 0, Pend_RB(0), 0, 0);
            // Footy_RB = FuncPoly5th(RealTime, tend_w, t3, Pstart_RB(1), 0, 0, Pend_RB(1), 0, 0);
            // Footz_RB = FuncPoly6th(RealTime, tend_w, t3, 0, 0, 0, 0, 0, 0, Fh);
            Footx_RB = FuncPoly3rd(RealTime, tend_w, t3, Pstart_RB(0), 0, Pend_RB(0), 0);
            Footy_RB = FuncPoly3rd(RealTime, tend_w, t3, Pstart_RB(1), 0, Pend_RB(1), 0);
            Footz_RB = FuncPoly4th(RealTime, tend_w, t3, 0, 0, 0, 0, Fh);
        }       
    }
    else
    {
        if (((k + 1) % 2) == 0)
        {
            Footx_LF << Strx*0.5 + Strx*(kx + 1) + yawStr_LF(0), 0, 0;
            Footy_LF << Stry*0.5 + Stry*(kx + 1) + yawStr_LF(1), 0, 0;
            Footz_LF << 0, 0, 0;

            Footx_RF << Strx*0.5 + Strx*(kx + 1) + yawStr_RF(0), 0, 0;
            Footy_RF << Stry*0.5 + Stry*(kx + 1) + yawStr_RF(1), 0, 0;
            Footz_RF << 0, 0, 0;

            Footx_LB << Strx*0.5 + Strx*(kx + 1) + yawStr_LB(0), 0, 0;
            Footy_LB << Stry*0.5 + Stry*(kx + 1) + yawStr_LB(1), 0, 0;
            Footz_LB << 0, 0, 0;

            Footx_RB << Strx*0.5 + Strx*(kx + 1) + yawStr_RB(0), 0, 0;
            Footy_RB << Stry*0.5 + Stry*(kx + 1) + yawStr_RB(1), 0, 0;
            Footz_RB << 0, 0, 0;
        }
        else
        {
            Footx_LF << Strx*(kx + 1) + yawStr_LF(0), 0, 0;
            Footy_LF << Stry*(kx + 1) + yawStr_LF(1), 0, 0;
            Footz_LF << 0, 0, 0;

            Footx_RF << Strx*(kx + 1) + yawStr_RF(0), 0, 0;
            Footy_RF << Stry*(kx + 1) + yawStr_RF(1), 0, 0;
            Footz_RF << 0, 0, 0;

            Footx_LB << Strx*(kx + 1) + yawStr_LB(0), 0, 0;
            Footy_LB << Stry*(kx + 1) + yawStr_LB(1), 0, 0;
            Footz_LB << 0, 0, 0;

            Footx_RB << Strx*(kx + 1) + yawStr_RB(0), 0, 0;
            Footy_RB << Stry*(kx + 1) + yawStr_RB(1), 0, 0;
            Footz_RB << 0, 0, 0;
        }
    }
    /* #endregion: Feet Trajectory End */
   
}

void trajectory::trajGeneration(double RealTime, bool walkEnable, double command_Vx, double command_Vy, double command_Yaw, double height, double dt)
{   

    if(!walk_enabled && walkEnable) w_start_time = RealTime;
    walk_enabled = walkEnable;

    if(walk_enabled){
        can_switch = (AreDoubleSame(Footz_RF(0),Fc)) || (AreDoubleSame(Footz_LF(0),Fc));
    } else {
        can_switch = false;
    }


    if(abs(Vx_mean) > 0 || abs(Vy_mean) > 0){
        can_stop = false;
    } else {
        can_stop = (AreDoubleSame(Footz_RF(0),0)) && (AreDoubleSame(Footz_LF(0),0));
    }

    if(can_switch && walk_enabled){
        Vx_mean = command_Vx; 
        Vy_mean = command_Vy;
        Yaw_mean = command_Yaw;
    }


    if(prev_vx_mean != Vx_mean) Comx += Cx;
    if(prev_vy_mean != Vy_mean) Comy += Cy;
    if(prev_Yaw_mean != Yaw_mean) ComYaw += Cyaw; 

    if((prev_vx_mean != Vx_mean) || (prev_vy_mean != Vy_mean) || (prev_Yaw_mean != Yaw_mean)){
        if((Vx_mean > prev_vx_mean) && AreDoubleSame(prev_vx_mean,0.0)) w_start_time = RealTime;
        if(AreDoubleSame(Footz_RF(0),Fc)){
            w_start_time = RealTime - Ts - Td - Ts*0.5;
            comTrajectory(RealTime-w_start_time, Ts, Td, Nphase, px, py, Vx_mean, Vy_mean, Yaw_mean, height, Fc, dt);
            if(prev_vx_mean != Vx_mean) Comx -= Strx*0.5;
            if(prev_vy_mean != Vy_mean) Comy -= Stry*0.5;
            if(prev_Yaw_mean != Yaw_mean) ComYaw -= Cyaw;
        }else if(AreDoubleSame(Footz_LF(0),Fc)){
            w_start_time = RealTime - Ts - Td - Ts - Td - Ts*0.5;
            comTrajectory(RealTime-w_start_time, Ts, Td, Nphase, px, py, Vx_mean, Vy_mean, Yaw_mean, height, Fc, dt);
            if(prev_vx_mean != Vx_mean) Comx -= Strx;
            if(prev_vy_mean != Vy_mean) Comy -= Stry;
            if(prev_Yaw_mean != Yaw_mean) ComYaw -= Cyaw;
        }

        if(prev_vx_mean != Vx_mean) prev_vx_mean = Vx_mean;
        if(prev_vy_mean != Vy_mean) prev_vy_mean = Vy_mean;
        if(prev_Yaw_mean != Yaw_mean) prev_Yaw_mean = Yaw_mean;
    }
    
    if(walk_enabled){
        Nphase += 1;
    }
    else if(!walk_enabled && can_stop){
        Nphase = 0;
    }

    
    comTrajectory(RealTime-w_start_time, Ts, Td, Nphase, px, py, Vx_mean, Vy_mean, Yaw_mean, height, Fc, dt);
    Zc = height;
    Xc = Cx+Comx; Yc = Cy+Comy; 
    dXc = dCx; dYc = dCy; 
    ddXc = ddCx; ddYc = ddCy;
    if(walk_enabled) { Yawc = Cyaw/2; }
    else { Yawc = command_Yaw; }  

    Pfoot_LF << Footx_LF(0) + offsetPf_LF(0) + Comx, Footy_LF(0) + offsetPf_LF(1) + Comy, Footz_LF(0) + offsetPf_LF(2);  
    Pfoot_RF << Footx_RF(0) + offsetPf_RF(0) + Comx, Footy_RF(0) + offsetPf_RF(1) + Comy, Footz_RF(0) + offsetPf_RF(2);  
    Pfoot_LB << Footx_LB(0) + offsetPf_LB(0) + Comx, Footy_LB(0) + offsetPf_LB(1) + Comy, Footz_LB(0) + offsetPf_LB(2);  
    Pfoot_RB << Footx_RB(0) + offsetPf_RB(0) + Comx, Footy_RB(0) + offsetPf_RB(1) + Comy, Footz_RB(0) + offsetPf_RB(2);  
}
