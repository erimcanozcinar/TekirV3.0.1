#include "main.hpp"

int main(int argc, char** argv) {   
    /* #region: Raisim */
    auto binaryPath = raisim::Path::setFromArgv(argv[0]);
    raisim::World::setActivationKey(binaryPath.getDirectory() + "\\activation.raisim");
    raisim::World world;

    world.setTimeStep(0.001);
    world.setMaterialPairProp("steel", "steel", 0.95, 0.95, 0.001, 0.95, 0.001);
    world.setMaterialPairProp("steel", "rubber", 0.95, 0.15, 0.001, 0.95, 0.001);
    auto ground = world.addGround(0, "steel");
    // ground->setAppearance("hidden");

    auto quadruped = world.addArticulatedSystem("/home/erim/RaiSim_Simulations/TekirV3.0.1/rsc/urdf/tekir3mesh_new.urdf");
    // auto quadruped = world.addArticulatedSystem("/home/erim/raisim_ws/rsc/Tekir/urdf/tekir3.urdf");
    quadruped->getCollisionBody("Foot_lf/0").setMaterial("rubber");
    quadruped->getCollisionBody("Foot_rf/0").setMaterial("rubber");
    quadruped->getCollisionBody("Foot_lb/0").setMaterial("rubber");
    quadruped->getCollisionBody("Foot_rb/0").setMaterial("rubber");
    /* #endregion */

    /* #region: Create Log file */
    FILE* fp0;
    fp0 = fopen("/home/erim/RaiSim_Simulations/TekirV3.0.1/Log/dataLog.txt", "w");
    /* #endregion */

    /* #region: Allegro */
    al_init();
    al_install_joystick();
    ALLEGRO_EVENT_QUEUE *event_queue = al_create_event_queue();
    al_register_event_source(event_queue, al_get_joystick_event_source());
    std::cout << "Gamepad routines enabled." << std::endl;
    ALLEGRO_EVENT event;
    /* #endregion */

    /* #region: Init Eigen variables */
    Itorso << 0.410, 0, 0,
              0, 0.908, 0,
              0, 0, 1.192;

    prevgenVelocity.setZero();
    dQ_LF.setZero(); dQ_RF.setZero(); dQ_LB.setZero(); dQ_RB.setZero();
    Fcon_LF.setZero(); Fcon_RF.setZero(); Fcon_LB.setZero(); Fcon_RB.setZero();

    prev_Pbase.setZero();

    accBody.setZero(); velBody.setZero();
    jPositions.setZero(); jVelocities.setZero(); jAccelerations.setZero(); jffTorques.setZero();

    Kp << 120, 70, 70;
    Kd << 12, 7, 7;

    Tau_LF.setZero(); Tau_RF.setZero(); Tau_LB.setZero(); Tau_RB.setZero();
    
    prevdQcm.setZero(); prevQcm.setZero(); prevQcm_filtered.setZero();
    /* #endregion */

    /* #region: Init classes*/
    controller jStick;
    CoMPreviewTrajectory traj;
    /* #endregion */

    /* #region: Initialize Robot */
    Pf_LF << Pfx_f, Pfy + LatOut, Pfz;
    Pf_RF << Pfx_f, -Pfy - LatOut, Pfz;
    Pf_LB << -Pfx_b, Pfy + LatOut, Pfz;
    Pf_RB << -Pfx_b, -Pfy - LatOut, Pfz;
    Pcom << ComX, ComY, ComZ;
    torsoRot << ComRoll, ComPitch, ComYaw;
    Q_LF = fullBodyIKan(Pf_LF, Pcom, torsoRot, 1);
    Q_RF = fullBodyIKan(Pf_RF, Pcom, torsoRot, 2);
    Q_LB = fullBodyIKan(Pf_LB, Pcom, torsoRot, 3);
    Q_RB = fullBodyIKan(Pf_RB, Pcom, torsoRot, 4);
    initialConditions << ComX, ComY, ComZ, quat_w, quat_x, quat_y, quat_z, Q_LF(0), Q_LF(1), Q_LF(2), Q_RF(0), Q_RF(1), Q_RF(2), Q_LB(0), Q_LB(1), Q_LB(2), Q_RB(0), Q_RB(1), Q_RB(2);
    quadruped->setGeneralizedCoordinate(initialConditions);
    /* #endregion */

    /* #region: Launch raisim server for visualization.Can be visualized on raisimUnity */
    raisim::RaisimServer server(&world);
    server.setMap("default");
    server.focusOn(quadruped);
    server.launchServer();
    /* #endregion */

    while (!jStick.close) {
        auto start_time = std::chrono::high_resolution_clock::now();
        RS_TIMED_LOOP(int(world.getTimeStep()*1e6));
        t = world.getWorldTime();
        dt = world.getTimeStep();

        /* #region: TRAJECTORY GENERATION WITH CONTROLLER */
        jStick.dualShockController(event_queue, event);
        for(int i=0; i<22; i++)
        {
            cmdJoyF[i] = LPF(jStick.joyCmd[i], pre_cmdJoyF[i], 2*PI*0.2, dt);
            pre_cmdJoyF[i] = cmdJoyF[i];
        }
        cmd_Vx = jStick.joyCmd[0]; cmd_Vy = jStick.joyCmd[1];
        cmd_yaw = cmdJoyF[2];  cmd_pitch = cmdJoyF[3]; cmd_roll = cmdJoyF[4];

        traj.trajGeneration_preview(t, jStick.walkEnable, cmd_Vx, cmd_Vy, pre_cmdJoyF[5], cmdJoyF[5], dt);
        /* #endregion */
        
        /* #region: CONTACT DEFINITION START */
        for (auto& contact : quadruped->getContacts()) // LF:3, RF:2, LB:1, RB:0
        {
            if (contact.skip()) continue;
            if (quadruped->getBodyIdx("shank_lf") == contact.getlocalBodyIndex())
            {
                Fcon_LF = -contact.getContactFrame().e().transpose() * contact.getImpulse().e() / dt;
                Pcon_LF = contact.getPosition().e().transpose();
            }
            else if (quadruped->getBodyIdx("shank_rf") == contact.getlocalBodyIndex())
            {
                Fcon_RF = -contact.getContactFrame().e().transpose() * contact.getImpulse().e() / dt;
                Pcon_RF = contact.getPosition().e().transpose();
            }
            else if (quadruped->getBodyIdx("shank_lb") == contact.getlocalBodyIndex())
            {
                Fcon_LB = -contact.getContactFrame().e().transpose() * contact.getImpulse().e() / dt;
                Pcon_LB = contact.getPosition().e().transpose();
            }
            else if (quadruped->getBodyIdx("shank_rb") == contact.getlocalBodyIndex())
            {
                Fcon_RB = -contact.getContactFrame().e().transpose() * contact.getImpulse().e() / dt;
                Pcon_RB = contact.getPosition().e().transpose();
            }
        }
        /* #endregion */

        /* #region: READ ACTUAL DATA START */
        genCoordinates = quadruped->getGeneralizedCoordinate().e();
        genVelocity = quadruped->getGeneralizedVelocity().e();
        for (int i = 0; i < 18; i++)
        {
            genAcceleration(i) = Numdiff(genVelocity(i), prevgenVelocity(i), dt);
            prevgenVelocity(i) = genVelocity(i);
        }

        imuRot = quat2euler(genCoordinates(3), genCoordinates(4), genCoordinates(5), genCoordinates(6));
        dimuRot << Numdiff(imuRot(0), prev_imuRot(0), dt), Numdiff(imuRot(1), prev_imuRot(1), dt), Numdiff(imuRot(2), prev_imuRot(2), dt);
        prev_imuRot = imuRot;

        q_LF << genCoordinates(7), genCoordinates(8), genCoordinates(9);
        q_RF << genCoordinates(10), genCoordinates(11), genCoordinates(12);
        q_LB << genCoordinates(13), genCoordinates(14), genCoordinates(15);
        q_RB << genCoordinates(16), genCoordinates(17), genCoordinates(18);
        dq_LF << genVelocity(6), genVelocity(7), genVelocity(8);
        dq_RF << genVelocity(9), genVelocity(10), genVelocity(11);
        dq_LB << genVelocity(12), genVelocity(13), genVelocity(14);
        dq_RB << genVelocity(15), genVelocity(16), genVelocity(17);

        rootOrientation = quat2Rotmat(genCoordinates(3), genCoordinates(4), genCoordinates(5), genCoordinates(6));
        rootAngvelocity << genVelocity(3), genVelocity(4), genVelocity(5);
        rootAngacceleration << genAcceleration(3), genAcceleration(4), genAcceleration(5);

        rootAbsposition << genCoordinates(0), genCoordinates(1), genCoordinates(2);
        rootAbsvelocity << genVelocity(0), genVelocity(1), genVelocity(2);
        rootAbsacceleration << genAcceleration(0), genAcceleration(1), genAcceleration(2);

        jPositions << q_LF(0), q_LF(1), q_LF(2), q_RF(0), q_RF(1), q_RF(2), q_LB(0), q_LB(1), q_LB(2), q_RB(0), q_RB(1), q_RB(2);
        jVelocities << dq_LF(0), dq_LF(1), dq_LF(2), dq_RF(0), dq_RF(1), dq_RF(2), dq_LB(0), dq_LB(1), dq_LB(2), dq_RB(0), dq_RB(1), dq_RB(2);
        jAccelerations << ddQ_LF(0), ddQ_LF(1), ddQ_LF(2), ddQ_RF(0), ddQ_RF(1), ddQ_RF(2), ddQ_LB(0), ddQ_LB(1), ddQ_LB(2), ddQ_RB(0), ddQ_RB(1), ddQ_RB(2);    
        /* #endregion */
         
        /* #region: ORIENTATION CONTROL */
        // Zpitch_front = Kp_pitch*(0 - imuRot(1)) + Kd_pitch*(0 - dimuRot(1));
        // Zpitch_back = -Kp_pitch*(0 - imuRot(1)) - Kd_pitch*(0 - dimuRot(1));
        // Zroll_left = -Kp_roll*(0 - imuRot(0)) - Kd_roll*(0 - dimuRot(0));
        // Zroll_right = Kp_roll * (0 - imuRot(0)) + Kd_roll * (0 - dimuRot(0));
        // Zpitch_LF = Zpitch_front; Zpitch_RF = Zpitch_front; Zpitch_LB = Zpitch_back; Zpitch_RB = Zpitch_back;
        // Zroll_LF = Zroll_left; Zroll_RF = Zroll_right; Zroll_LB = Zroll_left; Zroll_RB = Zroll_right;
        // if (traj.Pfoot_L(2) > 0.0)
        // {
        //     Zroll_LF = 0.0; Zpitch_LF = 0.0;
        //     Zroll_RB = 0.0; Zpitch_RB = 0.0;
        // }
        // if (traj.Pfoot_R(2) > 0.0)
        // {
        //     Zroll_RF = 0.0; Zpitch_RF = 0.0;
        //     Zroll_LB = 0.0; Zpitch_LB = 0.0;
        // }
        /* #endregion */

        /* #region: CENTRODIAL MOMENTUM */
        // if (t>5)
        // {
        //     Zfoot_offset = 0.1;
        //     Fcoef = 0;
        // }else{
        //     Zfoot_offset = 0.0;
        //     Fcoef = 1;
        // }

        // dQcm = centrodialMomentum(quadruped);
        // for(int i = 0; i < 6; i++)
        // {
        //     Qcm(i) = numIntegral(dQcm(i), prevdQcm(i), prevQcm(i), dt);
        //     Qcm_filtered(i) = HPF(Qcm(i), prevQcm(i), prevQcm_filtered(i), 2*PI*0.1, dt);
        //     prevdQcm(i) = dQcm(i);
        //     prevQcm(i) = Qcm(i);
        //     prevQcm_filtered(i) = Qcm_filtered(i);
        // }
        // Qcm_RF << Qcm_filtered(0), Qcm_filtered(1), Qcm_filtered(2);
        // Qcm_LB << Qcm_filtered(3), Qcm_filtered(4), Qcm_filtered(5);
        /* #endregion */

        /* #region: DESIRED FOOT POSITION */
        Pcom << traj.Xc, traj.Yc, traj.height;
        Rcom << genCoordinates(0), genCoordinates(1), genCoordinates(2);
        torsoRot << cmd_roll, cmd_pitch, cmd_yaw; // Torso Orientation       

        Pf_LF << traj.Pfoot_L(0) + Pfx_f, traj.Pfoot_L(1) + Pfy + LatOut, Pfz + traj.Pfoot_L(2);
        Pf_RF << traj.Pfoot_R(0) + Pfx_f, traj.Pfoot_R(1) - Pfy - LatOut, Pfz + traj.Pfoot_R(2);
        Pf_LB << traj.Pfoot_R(0) - Pfx_b, traj.Pfoot_R(1) + Pfy + LatOut, Pfz + traj.Pfoot_R(2);
        Pf_RB << traj.Pfoot_L(0) - Pfx_b, traj.Pfoot_L(1) - Pfy - LatOut, Pfz + traj.Pfoot_L(2);
        /* #endregion */

        /* #region: ANALYTIC INVERSE KINEMATIC */
        pre_dQ_LF = dQ_LF;
        pre_dQ_RF = dQ_RF;
        pre_dQ_LB = dQ_LB;
        pre_dQ_RB = dQ_RB;
        prevQ_LF = Q_LF;
        prevQ_RF = Q_RF;
        prevQ_LB = Q_LB;
        prevQ_RB = Q_RB;
        Q_LF = fullBodyIKan(Pf_LF, Pcom, torsoRot, 1);
        Q_RF = fullBodyIKan(Pf_RF, Pcom, torsoRot, 2);
        Q_LB = fullBodyIKan(Pf_LB, Pcom, torsoRot, 3);
        Q_RB = fullBodyIKan(Pf_RB, Pcom, torsoRot, 4);

        dQ_LF << Numdiff(Q_LF(0), prevQ_LF(0), dt), Numdiff(Q_LF(1), prevQ_LF(1), dt), Numdiff(Q_LF(2), prevQ_LF(2), dt);
        dQ_RF << Numdiff(Q_RF(0), prevQ_RF(0), dt), Numdiff(Q_RF(1), prevQ_RF(1), dt), Numdiff(Q_RF(2), prevQ_RF(2), dt);
        dQ_LB << Numdiff(Q_LB(0), prevQ_LB(0), dt), Numdiff(Q_LB(1), prevQ_LB(1), dt), Numdiff(Q_LB(2), prevQ_LB(2), dt);
        dQ_RB << Numdiff(Q_RB(0), prevQ_RB(0), dt), Numdiff(Q_RB(1), prevQ_RB(1), dt), Numdiff(Q_RB(2), prevQ_RB(2), dt);

        ddQ_LF << Numdiff(dQ_LF(0), pre_dQ_LF(0), dt), Numdiff(dQ_LF(1), pre_dQ_LF(1), dt), Numdiff(dQ_LF(2), pre_dQ_LF(2), dt);
        ddQ_RF << Numdiff(dQ_RF(0), pre_dQ_RF(0), dt), Numdiff(dQ_RF(1), pre_dQ_RF(1), dt), Numdiff(dQ_RF(2), pre_dQ_RF(2), dt);
        ddQ_LB << Numdiff(dQ_LB(0), pre_dQ_LB(0), dt), Numdiff(dQ_LB(1), pre_dQ_LB(1), dt), Numdiff(dQ_LB(2), pre_dQ_LB(2), dt);
        ddQ_RB << Numdiff(dQ_RB(0), pre_dQ_RB(0), dt), Numdiff(dQ_RB(1), pre_dQ_RB(1), dt), Numdiff(dQ_RB(2), pre_dQ_RB(2), dt);
        /* #endregion */

        /* #region: VMC CONTROLLER FOR TORSO */
        Rf_LF = fullBodyFK(imuRot, {0,0,0}, q_LF, 1);
        Rf_RF = fullBodyFK(imuRot, {0,0,0}, q_RF, 2);
        Rf_LB = fullBodyFK(imuRot, {0,0,0}, q_LB, 3);
        Rf_RB = fullBodyFK(imuRot, {0,0,0}, q_RB, 4);

        dtorsoRot << Numdiff(torsoRot(0), prev_torsoRot(0), dt), Numdiff(torsoRot(1), prev_torsoRot(1), dt), Numdiff(torsoRot(2), prev_torsoRot(2), dt);
        ddtorsoRot << Numdiff(dtorsoRot(0), prev_dtorsoRot(0), dt), Numdiff(dtorsoRot(1), prev_dtorsoRot(1), dt), Numdiff(dtorsoRot(2), prev_dtorsoRot(2), dt);
        prev_torsoRot = torsoRot;
        prev_dtorsoRot = dtorsoRot;

        orientControlRef = QuaternionRef(torsoRot(0), torsoRot(1), torsoRot(2), dtorsoRot(0), dtorsoRot(1), dtorsoRot(2), ddtorsoRot(0), ddtorsoRot(1), ddtorsoRot(2));
        Eigen::Vector3d Wref, dWref, quatVecRef, quatVec;
        Wref << orientControlRef(4),orientControlRef(5),orientControlRef(6);
        dWref << orientControlRef(7),orientControlRef(8),orientControlRef(9);
        quatVecRef << orientControlRef(1), orientControlRef(2), orientControlRef(3);
        quatVec << genCoordinates(4), genCoordinates(5), genCoordinates(6);

        Mvmc = Itorso*(dWref + sqrt(200)*(Wref - rootAngvelocity) - 200*(orientControlRef(0)*quatVec - genCoordinates(3)*quatVecRef + 0*vec2SkewSym(quatVecRef)*quatVec));

        Fvmc << MASS*(traj.ddXc + (traj.ddXc-genAcceleration(0))*0.1), MASS*(traj.ddYc + (traj.ddYc-genAcceleration(1))*0.1), MASS*(GRAVITY + (ddZcom-genAcceleration(2))*0.1), Mvmc(0), Mvmc(1), Mvmc(2);
        Fmatrix = VMC(Rf_LF, Rf_RF, Rf_LB, Rf_RB, Fcon_LF, Fcon_RF, Fcon_LB, Fcon_RB, Fvmc, dt);
        F1cont << Fmatrix(0, 0), Fmatrix(0, 1), Fmatrix(0, 2);
        F2cont << Fmatrix(1, 0), Fmatrix(1, 1), Fmatrix(1, 2);
        F3cont << Fmatrix(2, 0), Fmatrix(2, 1), Fmatrix(2, 2);
        F4cont << Fmatrix(3, 0), Fmatrix(3, 1), Fmatrix(3, 2);
        /* #endregion */
        
        /* #region: INVERSE DYNAMICS */
        jffTorques = funNewtonEuler4Leg(rootAbsacceleration, rootOrientation, rootAngvelocity, rootAngacceleration, jPositions, jVelocities, jAccelerations, rootOrientation.transpose()*F1cont, rootOrientation.transpose()*F2cont, rootOrientation.transpose()*F3cont, rootOrientation.transpose()*F4cont);
        jffTorques2 = funNewtonEuler4Leg(rootAbsacceleration, rootOrientation, rootAngvelocity*0, rootAngacceleration, jPositions, jVelocities*0, jAccelerations, F1cont*0, F2cont*0, F3cont*0, F4cont*0);
        /* #endregion */

        /* #region: INVERSE DYNAMICS WITH RAISIM*/
        Tau1_JF = JacTranspose_LF(q_LF(0), q_LF(1), q_LF(2), 30*PI/180) * (rootOrientation.transpose() * F1cont);
        Tau2_JF = JacTranspose_RF(q_RF(0), q_RF(1), q_RF(2), -30*PI/180) * (rootOrientation.transpose() * F2cont);
        Tau3_JF = JacTranspose_LB(q_LB(0), q_LB(1), q_LB(2), -30*PI/180) * (rootOrientation.transpose() * F3cont);
        Tau4_JF = JacTranspose_RB(q_RB(0), q_RB(1), q_RB(2), 30*PI/180) * (rootOrientation.transpose() * F4cont);

        JF << 0, 0, 0, 0, 0, 0, Tau1_JF(0), Tau1_JF(1), Tau1_JF(2), Tau2_JF(0), Tau2_JF(1), Tau2_JF(2), Tau3_JF(0), Tau3_JF(1), Tau3_JF(2), Tau4_JF(0), Tau4_JF(1), Tau4_JF(2);
        genAccelerationVec << genAcceleration(0), genAcceleration(1), genAcceleration(2), rootAngacceleration(0), rootAngacceleration(1), rootAngacceleration(2), ddQ_LF(0), ddQ_LF(1), ddQ_LF(2), ddQ_RF(0), ddQ_RF(1), ddQ_RF(2), ddQ_LB(0), ddQ_LB(1), ddQ_LB(2), ddQ_RB(0), ddQ_RB(1), ddQ_RB(2);
        jointTorques = quadruped->getMassMatrix().e() * genAccelerationVec + quadruped->getNonlinearities(world.getGravity()).e() - 0*JF;
        /* #endregion */

        /* #region: PD CONTROLLER */
        for (int i = 0; i < 3; i++)
        {
            // i = 0: Hip AA, i = 1: Hip FE, i = 2: Knee FE
            Tau_LF(i) = Kp(i)*(Q_LF(i) - q_LF(i)) + Kd(i)*(dQ_LF(i) - dq_LF(i)) + jffTorques(i);
            Tau_RF(i) = Kp(i)*(Q_RF(i) - q_RF(i)) + Kd(i)*(dQ_RF(i) - dq_RF(i)) + jffTorques(i+3);
            Tau_LB(i) = Kp(i)*(Q_LB(i) - q_LB(i)) + Kd(i)*(dQ_LB(i) - dq_LB(i)) + jffTorques(i+6);
            Tau_RB(i) = Kp(i)*(Q_RB(i) - q_RB(i)) + Kd(i)*(dQ_RB(i) - dq_RB(i)) + jffTorques(i+9);
        }
        /* #endregion */

        /* #region: SEND COMMEND TO THE ROBOT */
        F << 0, 0, 0, 0, 0, 0, Tau_LF(0), Tau_LF(1), Tau_LF(2), Tau_RF(0), Tau_RF(1), Tau_RF(2), Tau_LB(0), Tau_LB(1), Tau_LB(2), Tau_RB(0), Tau_RB(1), Tau_RB(2);
        quadruped->setGeneralizedForce(F);
        /* #endregion */

        /* #region: PUSH ROBOT */
        Eigen::Vector3d F_ex;
        if (t > 6.13 && t < 6.23)
        {
            F_ex << 0, 100, 0;
        }
        else 
        {
            F_ex << 0, 0, 0;
        }
        //quadruped->setExternalForce(quadruped->getBodyIdx("torso"), F_ex);
        /* #endregion */
                
        // server.setCameraPositionAndLookAt({genCoordinates(0)-1.3,genCoordinates(1),genCoordinates(2)+0.6}, {genCoordinates(0),genCoordinates(1),genCoordinates(2)});
               
                 
        server.integrateWorldThreadSafe(); 

        auto end_time = std::chrono::high_resolution_clock::now();
        auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        // std::cout << traj.ddYc << " microseconds" << std::endl;  

        /* Log data (fp0:NewtonEulerTorques, fp1:PD_torques, fp2: InvDynTorques ) */
        fprintf(fp0, "%f %f %f %f %f %f %f %f\n", t, traj.Xc, traj.Yc, traj.dXc, traj.dYc, traj.Pfoot_R(0), traj.Pfoot_L(0), cmdJoyF[5]);

    }
    server.killServer();
    fclose(fp0);
    al_destroy_event_queue(event_queue);
    std::cout << "end of simulation" << std::endl;
}