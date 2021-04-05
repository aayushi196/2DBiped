#include <MatrixMath.h>
#include <math.h>
// out means left
#define left_hip_dir 52
#define left_hip_pwm 9
#define right_hip_dir 50
#define right_hip_pwm 6

#define left_knee_dir 53
#define left_knee_pwm 10
#define right_knee_dir 51
#define right_knee_pwm 7

#define left_hip_encA 20
#define left_hip_encB 44
#define right_hip_encA 18
#define right_hip_encB 42

#define left_knee_encA 21
#define left_knee_encB 45
#define right_knee_encA 19
#define right_knee_encB 43

#define limit_switch_lh 34
#define limit_switch_rh 32
#define limit_switch_lk 35
#define limit_switch_rk 33

//double theta_hip[42]= [1.20061, 1.20747, 1.21564, 1.22496, 1.23530, 1.24661, 1.25878, 1.271805, 1.28562, 1.30022, 1.31559, 1.33172, 1.34863, 1.36633,
//                       1.38484, 1.40420, 1.42447, 1.44569, 1.46796, 1.49138, 1.51609, 1.530914, 1.55239, 1.54884, 1.50927, 1.43563, 1.33883, 1.23476,
//                       1.14029, 1.06929, 1.02862, 1.01416, 1.00964, 1.01121, 1.02322, 1.045936, 1.07661, 1.11073, 1.14320, 1.16959, 1.18735, 1.19703];
//
//double theta_knee[42]= [0.46714, 0.48853, 0.50684, 0.52244, 0.53559, 0.54651, 0.55534, 0.56221, 0.56720, 0.57037, 0.57176, 0.57138, 0.56924, 0.56530,
//                        0.55951, 0.55182, 0.54211, 0.53026, 0.51609, 0.49937, 0.47981, 0.46714, 0.45061, 0.47025, 0.54234, 0.66351, 0.81686, 0.97808,
//                        1.12156, 1.22645, 1.28282, 1.29774, 1.29510, 1.26465, 1.18718, 1.06366, 0.90946, 0.74829, 0.60599, 0.50449, 0.45560, 0.45497];
mtx_type dT_dt_dthetaDot[4][4];
mtx_type dT_dtheta[4][4];
mtx_type B_q[4][4];
mtx_type C_q_qdot[4][4];
mtx_type D_x[4][1];
mtx_type E_x[4][1];
mtx_type G_q[4][1];
mtx_type F_ext[4][1];
mtx_type Torque_dash[4][1];
mtx_type Torque[4][1];
mtx_type D_E_xdot[4][1];
mtx_type theta_dot_column[4][1];
mtx_type alpha_torque_dash[4][1];
mtx_type beta1[4][1];
mtx_type beta[4][1];

double pulley_ratio = 0.6667;
volatile long enc_lh = 0, enc_rh = 0, enc_lk = 0, enc_rk = 0;

double l[4] = {0.32, 0.32, 0.32, 0.32};
double lc[4] = {0.041, 0.142, 0.041, 0.142};
double m[4] = {1.094, 0.83, 1.094, 0.83};
double Ic[4] = {0.02, 0.001, 0.02, 0.001};

double theta[4] = { 1.200 , 0.467, 1.481 , 0.467}, theta_dot[4] = {0, 0, 0, 0}, theta_ddot[4] = {0, 0, 0, 0};
double prev_theta[4] = { 0, 0, 0, 0}, prev_theta_dot[4] = { 0, 0, 0, 0};
double err[4] = {0, 0, 0, 0}, err_dot[4] = {0, 0, 0, 0}, err_sum[4] = {0, 0, 0, 0};

double theta_lk = 0, theta_rk = 0, theta_lh = 0, theta_rh = 0;
double err_lk = 0.0, err_rk = 0.0, err_lh = 0.0, err_rh = 0.0;

double thetaD[4] = { 0 , 0, 0 , 0}, thetaD_dot[4] = {0, 0, 0, 0}, thetaD_ddot[4] = {0, 0, 0, 0};
double thetaD_stance[2] = {0, 0}, thetaD_stance_dot[2] = {0, 0}, thetaD_stance_ddot[2] = {0, 0};
double thetaD_swing[2] = {0, 0}, thetaD_swing_dot[2] = {0, 0}, thetaD_swing_ddot[2] = {0, 0};
double thetaD_stance_vel_acc[4] = {0, 0, 0, 0};

double Kp = 25.0, Kv = 10.0, Ki = 0.1;

//COM parameters
double z = 0.62;
double g = 9.81;
double Tc = sqrt(z / g);
double origin[2] = {0, 0};
double stride_length = 0.15;
double x_dot0 = 0.30;
double x_0 = -stride_length / 2;
double h = 0.05;
double stride_time = Tc * log((-stride_length / 2 - Tc*x_dot0) / (stride_length / 2 - Tc*x_dot0));
double P_com[2] = {0 , 0}, P_com_dot[2] = {0 , 0}, P_com_ddot[2] = {0 , 0};
double com_dot[2] = {0, 0};

const double Po[2] = { -x_0, z}, Pm[2] = {0, z - h}, Pf[2] = {x_0, z};
const double Po_dot[2] = { -x_dot0, 0}, Pm_dot[2] = {0.2, 0}, Pf_dot[2] = { -x_dot0, 0};
const double Po_ddot[2] = { -g*x_0 / z, 0}, Pm_ddot[2] = {0.2, 0}, Pf_ddot[2] = { -g*x_0 / z, 0};

double theta_o[2] = {0.0, 0.0}, theta_m[2] = {0, 0}, theta_f[2] = {0, 0};
double theta_o_vel_acc[4] = {0, 0, 0, 0}, theta_m_vel_acc[4] = {0, 0, 0, 0}, theta_f_vel_acc[4] = {0, 0, 0, 0};
double theta_1o[3] = {0, 0, 0}, theta_1m[3] = {0, 0, 0}, theta_1f[3] = {0, 0, 0}, theta_2o[3] = {0, 0, 0},  theta_2m[3] = {0, 0, 0},  theta_2f[3] = {0, 0, 0};
double a1_im[6], a1_mf[6], a2_im[6], a2_mf[6];

unsigned long   t_init = 0.0;
double t = 0, t_prev = 0.0, t1 = 0, dt = 0 ;
int flag = 1;

float Current[4] = { 0, 0, 0, 0}, Voltage[4] = { 0, 0, 0, 0} ;
int pwm[4] = {0, 0, 0, 0};
float Resistance = 2.2, K_omega = 0.4, K_torque = 0.4;
//float Resistance = 2.2, K_omega = 0.52, K_torque = 0.3;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
//  Serial.print(" Stride_time= ");
//  Serial.println(stride_time);
  motor_initialize();
  encoder_initialize();
  //  encoder_reset();
  inverse_kinematics(  Po, l, theta_o);
  inverse_jacobian(Po_dot, Po_ddot, theta_o, l, theta_o_vel_acc );
  inverse_kinematics(  Pm, l, theta_m);
  inverse_jacobian(Pm_dot, Pm_ddot, theta_m, l, theta_m_vel_acc );
  inverse_kinematics(  Pf, l, theta_f);
  inverse_jacobian(Pf_dot, Pf_ddot, theta_f, l, theta_f_vel_acc );

  theta_1o[0] = theta_o[0];
  theta_1o[1] = theta_o_vel_acc[0];
  theta_1o[2] = theta_o_vel_acc[2];
  theta_1m[0] = theta_m[0];
  theta_1m[1] = theta_m_vel_acc[0];
  theta_1m[2] = theta_m_vel_acc[2];
  theta_1f[0] = theta_f[0];
  theta_1f[1] = theta_f_vel_acc[0];
  theta_1f[2] = theta_f_vel_acc[2];

  theta_2o[0] = theta_o[1];
  theta_2o[1] = theta_o_vel_acc[1];
  theta_2o[2] = theta_o_vel_acc[3];
  theta_2m[0] = theta_m[1];
  theta_2m[1] = theta_m_vel_acc[1];
  theta_2m[2] = theta_m_vel_acc[3];
  theta_2f[0] = theta_f[1];
  theta_2f[1] = theta_f_vel_acc[1];
  theta_2f[2] = theta_f_vel_acc[3];

  trajectory_generator_im( theta_1o, theta_1m, theta_2o, theta_2m, stride_time / 2);
  trajectory_generator_mf( theta_1m, theta_1f, theta_2m, theta_2f, stride_time / 2);

  //  t_init = millis();
  //  Serial.println(t_init);
}

void loop() {
  // put your main code here, to run repeatedly:
  if (t < stride_time) {
    //  t = (millis()- t_init)%636;
    //  t = t/1000.00;
    //  if (t > stride_time) t = 0.0;
    t = t + 0.01;
    //  Serial.print(" Time= ");
    //  Serial.println(t);

    P_com[0] = x_0 * cosh(t / Tc) + Tc * x_dot0 * sinh(t / Tc);
    P_com[1] = z;
    P_com_dot[0] = -(x_0 * sinh(t / Tc) / Tc + x_dot0 * cosh(t / Tc));
    P_com_dot[1] = 0;
    P_com_ddot[0] =-(g / z * P_com[0]);
    P_com_ddot[1] = 0;

    //swing leg trajectory
    inverse_kinematics(  P_com, l, thetaD_stance);
    inverse_jacobian( P_com_dot, P_com_ddot, thetaD_stance, l, thetaD_stance_vel_acc );

    thetaD_stance_dot[0] = thetaD_stance_vel_acc[0];
    thetaD_stance_dot[1] = thetaD_stance_vel_acc[1];
    thetaD_stance_ddot[0] = thetaD_stance_vel_acc[2];
    thetaD_stance_ddot[1] = thetaD_stance_vel_acc[3];

    //swing leg trajectory
    if (t < stride_time / 2) {
      thetaD_swing[0] = a1_im[0] + a1_im[1] * t + a1_im[2] * pow(t, 2) + a1_im[3] * pow(t, 3) + a1_im[4] * pow(t, 4) + a1_im[5] * pow(t, 5);
      thetaD_swing[1] = a2_im[0] + a2_im[1] * t + a2_im[2] * pow(t, 2) + a2_im[3] * pow(t, 3) + a2_im[4] * pow(t, 4) + a2_im[5] * pow(t, 5);

      thetaD_swing_dot[0] = a1_im[1] + 2 * a1_im[2] * t + 3 * a1_im[3] * pow(t, 2) + 4 * a1_im[4] * pow(t, 3) + 5 * a1_im[5] * pow(t, 4);
      thetaD_swing_dot[1] = a2_im[1] + 2 * a2_im[2] * t + 3 * a2_im[3] * pow(t, 2) + 4 * a2_im[4] * pow(t, 3) + 5 * a2_im[5] * pow(t, 4);

      thetaD_swing_ddot[0] = 2 * a1_im[2] + 6 * a1_im[3] * t + 12 * a1_im[4] * pow(t, 2) + 20 * a1_im[5] * pow(t, 3);
      thetaD_swing_ddot[1] = 2 * a2_im[2] + 6 * a2_im[3] * t + 12 * a2_im[4] * pow(t, 2) + 20 * a2_im[5] * pow(t, 3);
    }
    else {
      t1 = t - stride_time / 2;
      thetaD_swing[0] = a1_mf[0] + a1_mf[1] * t1 + a1_mf[2] * pow(t1, 2) + a1_mf[3] * pow(t1, 3) + a1_mf[4] * pow(t1, 4) + a1_mf[5] * pow(t1, 5);
      thetaD_swing[1] = a2_mf[0] + a2_mf[1] * t1 + a2_mf[2] * pow(t1, 2) + a2_mf[3] * pow(t1, 3) + a2_mf[4] * pow(t1, 4) + a2_mf[5] * pow(t1, 5);

      thetaD_swing_dot[0] = a1_mf[1] + 2 * a1_mf[2] * t1 + 3 * a1_mf[3] * pow(t1, 2) + 4 * a1_mf[4] * pow(t1, 3) + 5 * a1_mf[5] * pow(t1, 4);
      thetaD_swing_dot[1] = a2_mf[1] + 2 * a2_mf[2] * t1 + 3 * a2_mf[3] * pow(t1, 2) + 4 * a2_mf[4] * pow(t1, 3) + 5 * a2_mf[5] * pow(t1, 4);

      thetaD_swing_ddot[0] = 2 * a1_mf[2] + 6 * a1_mf[3] * t1 + 12 * a1_mf[4] * pow(t1, 2) + 20 * a1_mf[5] * pow(t1, 3);
      thetaD_swing_ddot[1] = 2 * a2_mf[2] + 6 * a2_mf[3] * t1 + 12 * a2_mf[4] * pow(t1, 2) + 20 * a2_mf[5] * pow(t1, 3);
    }

    if (flag == 1) {
      thetaD[0] = thetaD_stance[0];
      thetaD[1] = thetaD_stance[1];
      thetaD[2] = thetaD_swing[0];
      thetaD[3] = thetaD_swing[1];
      thetaD_dot[0] = thetaD_stance_dot[0];
      thetaD_dot[1] = thetaD_stance_dot[1];
      thetaD_dot[2] = thetaD_swing_dot[0];
      thetaD_dot[3] = thetaD_swing_dot[1];
      thetaD_ddot[0] = thetaD_stance_ddot[0];
      thetaD_ddot[1] = thetaD_stance_ddot[1];
      thetaD_ddot[2] = thetaD_swing_ddot[0];
      thetaD_ddot[3] = thetaD_swing_ddot[1];
    }
    else
    {
      thetaD[0] = thetaD_swing[0];
      thetaD[1] = thetaD_swing[1];
      thetaD[2] = thetaD_stance[0];
      thetaD[3] = thetaD_stance[1];
      thetaD_dot[0] = thetaD_swing_dot[0];
      thetaD_dot[1] = thetaD_swing_dot[1];
      thetaD_dot[2] = thetaD_stance_dot[0];
      thetaD_dot[3] = thetaD_stance_dot[1];
      thetaD_ddot[0] = thetaD_swing_ddot[0];
      thetaD_ddot[1] = thetaD_swing_ddot[1];
      thetaD_ddot[2] = thetaD_stance_ddot[0];
      thetaD_ddot[3] = thetaD_stance_ddot[1];
    }
    //Inverse Dynamics and and controller
    com_dot[0] = -P_com_dot[0]; com_dot[1] = -P_com_ddot[0];
    inverse_dynamics( theta, theta_dot, l, lc, com_dot, m, Ic );

    for (int i = 0; i < 4; i++) {
      err[i] = thetaD[i] - theta[i];
      err_dot[i] = thetaD_dot[i] - theta_dot[i];
      err_sum[i] = err_sum[i] + err[i];
      //torque_dash using PID
      Torque_dash[i][0] = thetaD_ddot[i] + Kp * err[i] + Kv * err_dot[i] + Ki * err_sum[i];
      theta_dot_column[i][0] = theta_dot[i];
    }
    //Torque = B_q*Torque_dash' + C_q_qdot*theta_dot' + F_ext;
    Matrix.Multiply((mtx_type*)B_q, (mtx_type*)Torque_dash, 4, 4, 1, (mtx_type*)alpha_torque_dash);
    Matrix.Multiply((mtx_type*)C_q_qdot, (mtx_type*)theta_dot_column, 4, 4, 1, (mtx_type*)beta1);
    Matrix.Add((mtx_type*) beta1, (mtx_type*)F_ext , 4, 1, (mtx_type*)beta );
    Matrix.Add((mtx_type*) alpha_torque_dash, (mtx_type*)beta , 4, 1, (mtx_type*)Torque );

    //pwm using motor dynamics
    if (flag == 1) {
      for (int i = 0; i < 4; i++) {
        Current[i] = Torque[i][0] * 1 / K_torque;
        Voltage[i] = Current[i] * Resistance + K_omega * theta_dot[i];
        pwm[i] = round(Voltage[i] * 255 / 12);
        //        if (pwm[i] > 150) pwm[i] = 150;
        //        if (pwm[i] < -150) pwm[i] = -150;

      }
      //    if (P_com[0] > stride_length / 2) {
      if (t >= stride_time) {
        flag = 0;
      }
    }

    else {
      //pwm using motor dynamics
      for (int i = 0; i < 4; i++) {
        Current[i] = Torque[i][0] * 1 / K_torque;
        Voltage[i] = Current[i] * Resistance + K_omega * theta_dot[i];
        pwm[i] = round(Voltage[i] * 255 / 12);
        //        if (pwm[i] > 150) pwm[i] = 150;
        //        if (pwm[i] < -150) pwm[i] = -150;

      }
      //  if (P_com[0] > stride_length / 2) {
      if (t >= stride_time) {
        flag = 1;
      }
    }
//            Serial.print(" thetaD[0]= ");
//            Serial.print(thetaD[0]);
//            Serial.print(" thetaD[1]= ");
//            Serial.print(thetaD[1]);
//            Serial.print(" thetaD[2]= ");
//            Serial.print(thetaD[2]);
//            Serial.print(" thetaD[3]= ");
//            Serial.println(thetaD[3]);

    //        Serial.print(" theta[0]= ");
    //        Serial.print(theta[0]);
    //        Serial.print(" theta[1]= ");
    //        Serial.print(theta[1]);
    //        Serial.print(" theta[2]= ");
    //        Serial.print(theta[2]);
    //        Serial.print(" theta[3]= ");
    //        Serial.println(theta[3]);
    //
    //
//        Serial.print(" err[0]= ");
//        Serial.print(err[0]);
//        Serial.print(" err[1]= ");
//        Serial.print(err[1]);
//        Serial.print(" err[2]= ");
//        Serial.print(err[2]);
//        Serial.print(" err[3]= ");
//        Serial.println(err[3]);
    
    Serial.print(" Torque[0]= ");
    Serial.print(Torque[0][0]);
    Serial.print(" Torque[1]= ");
    Serial.print(Torque[1][0]);
    Serial.print(" Torque[2]= ");
    Serial.print(Torque[2][0]);
    Serial.print(" Torque[3]= ");
    Serial.println(Torque[3][0]);

    motor_run(left_hip_dir, left_hip_pwm, pwm[0]);
    motor_run(left_knee_dir, left_knee_pwm, pwm[1]);
    motor_run(right_hip_dir, right_hip_pwm, pwm[2]);
    motor_run(right_knee_dir, right_knee_pwm, pwm[3]);
    //    update_thetas_velocities();
    simulate_theta();

  }
  else {

    motor_run(left_hip_dir, left_hip_pwm, 0);
    motor_run(left_knee_dir, left_knee_pwm, 0);
    motor_run(right_hip_dir, right_hip_pwm, 0);
    motor_run(right_knee_dir, right_knee_pwm, 0);
  }

}
