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
mtx_type thetax[4][2];

double pulley_ratio = 0.6667;
volatile long enc_lh = 0, enc_rh = 0, enc_lk = 0, enc_rk = 0;

double l[4] = {0.32, 0.32, 0.32, 0.32};
double lc[4] = {0.041, 0.142, 0.041, 0.142};
double m[4] = {1.094, 0.83, 1.094, 0.83};
double Ic[4] = {0.02, 0.001, 0.02, 0.001};
double P0[2] = {0 , 0}, P1[2] = {0 , 0}, P2[2] = {0 , 0}, P4[2] = {0 , 0}, P5[2] = { 0, 0};
double theta[4] = { 1.200 , 0.467, 1.481 , 0.467}, theta_dot[4] = {0, 0, 0, 0}, theta_ddot[4] = {0, 0, 0, 0};
double prev_theta[4] = { 0, 0, 0, 0}, prev_theta_dot[4] = { 0, 0, 0, 0};
double theta_lk = 0, theta_rk = 0, theta_lh = 0, theta_rh = 0;
double err[4] = {0, 0, 0, 0}, err_dot[4] = {0, 0, 0, 0}, err_sum[4] = {0, 0, 0, 0};
double prev_err[4] = {0, 0, 0, 0};
double thetaD[4] = { 0 , 0, 0 , 0}, thetaD_dot[4] = {0, 0, 0, 0}, thetaD_ddot[4] = {0, 0, 0, 0};
double thetaD_stance[2] = {0, 0}, thetaD_stance_dot[2] = {0, 0}, thetaD_stance_ddot[2] = {0, 0};
double thetaD_swing[2] = {0, 0}, thetaD_swing_dot[2] = {0, 0}, thetaD_swing_ddot[2] = {0, 0};
double thetaD_stance_vel_acc[4] = {0, 0, 0, 0};

double pwm_lk = 0.0, pwm_rk = 0.0, pwm_lh = 0.0, pwm_rh = 0.0;
double err_lk = 0.0, err_rk = 0.0, err_lh = 0.0, err_rh = 0.0;
double corr_lk = 0.0, corr_rk = 0.0, corr_lh = 0.0, corr_rh = 0.0;
double prev_err_lk = 0.0, prev_err_rk = 0.0, prev_err_lh = 0.0, prev_err_rh = 0.0;
double diff_err_lk = 0.0, diff_err_rk = 0.0, diff_err_lh = 0.0, diff_err_rh = 0.0;
double err_sum_lk = 0.0, err_sum_rk = 0.0, err_sum_lh = 0.0, err_sum_rh = 0.0;
double Kp = 50, Kv = 10, Ki = 0.1;


//COM parameters
double z = 0.62;
double g = 9.81;
double Tc = sqrt(z / g);
double origin[2] = {0, 0};
double stride_length = 0.15;
double x_dot0 = 0.35;
double x_0 = -stride_length / 2;
double h = 0.06;
double stride_time = Tc * log((-stride_length / 2 - Tc*x_dot0) / (stride_length / 2 - Tc*x_dot0));
double P_com[2] = {0 , 0}, P_com_dot[2] = {0 , 0}, P_com_ddot[2] = {0 , 0};
double com_dot[2] = {0, 0};



int j = 0;
//int flag=1;
//double dt= 0.01;

const double Po[2] = { -x_0, z}, Pm[2] = {0, z - h}, Pf[2] = {x_0, z};
const double Po_dot[2] = { -x_dot0, 0}, Pm_dot[2] = {0.3, 0}, Pf_dot[2] = { -x_dot0, 0};
const double Po_ddot[2] = { -g*x_0 / z, 0}, Pm_ddot[2] = {0.3, 0}, Pf_ddot[2] = { -g*x_0 / z, 0};

double theta_o[2] = {0.0, 0.0}, theta_m[2] = {0, 0}, theta_f[2] = {0, 0};
double theta_o_vel_acc[4] = {0, 0, 0, 0}, theta_m_vel_acc[4] = {0, 0, 0, 0}, theta_f_vel_acc[4] = {0, 0, 0, 0};
double theta_1o[3] = {0, 0, 0}, theta_1m[3] = {0, 0, 0}, theta_1f[3] = {0, 0, 0}, theta_2o[3] = {0, 0, 0},  theta_2m[3] = {0, 0, 0},  theta_2f[3] = {0, 0, 0};
double a1_im[6], a1_mf[6], a2_im[6], a2_mf[6];
unsigned long   t_init = 0.0;
double t = 0, t_prev = 0.0, t1 = 0, dt = 0 ;
int flag = 1;
float Current[4] = { 0, 0, 0, 0}, Voltage[4] = { 0, 0, 0, 0} ;
int pwm[4] = {0, 0, 0, 0};
float Resistance = 2.2, K_omega = 0.2, K_torque = 0.2;

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  Serial.print(" Stride_time= ");
  Serial.println(stride_time);
  motor_initialize();
  encoder_initialize();
  encoder_reset();
  inverse_kinematics(  Po, l, theta_o);
  inverse_jacobian(Po_dot, Po_ddot, theta_o, l, theta_o_vel_acc );
  inverse_kinematics(  Pm, l, theta_m);
  inverse_jacobian(Pm_dot, Pm_ddot, theta_m, l, theta_m_vel_acc );
  inverse_kinematics(  Pf, l, theta_f);
  inverse_jacobian(Pf_dot, Pf_ddot, theta_f, l, theta_f_vel_acc );

  thetax[0][0] = 1.23;
  thetax[0][1] = 0.44;
  thetax[1][0] = 1.35;
  thetax[1][1] = 1.01;
  thetax[2][0] = 1.48;
  thetax[2][1] = 0.44;
  thetax[3][0] = 1.23;
  thetax[3][1] = 0.66;

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

  //  Serial.println(theta_1o[0]);
  //  Serial.println(theta_2o[0]);
  //  Serial.println(theta_1m[0]);
  //  Serial.println(theta_2m[0]);
  //  Serial.println(theta_1f[0]);
  //  Serial.println(theta_2f[0]);



  trajectory_generator_im( theta_1o, theta_1m, theta_2o, theta_2m, stride_time / 2);
  trajectory_generator_mf( theta_1m, theta_1f, theta_2m, theta_2f, stride_time / 2);

  t_init = millis();
  Serial.println(t_init);
}

void loop() {
  // put your main code here, to run repeatedly:
  //
  //     t = (millis()- t_init)%636;
  //     t = t/1000.00;
  //
  //     Serial.print(" Time= ");
  //     Serial.println(t);

  for (int i = 0; i < 4; i++) {

    if (j < 2) {
      j = i + 2;
    }
    else
    { j = 3 - i;
    }

    theta_lk = (enc_lk * 0.3 * PI) / 180;
    theta_rk = (enc_rk * 0.3 * PI) / 180;
    theta_lh = (enc_lh * 0.2 * PI) / 180;
    theta_rh = (enc_lh * 0.2 * PI) / 180;

    err_lk = theta_lk - thetax[i][1];
    err_rk = theta_rk - thetax[j][1];
    err_lh = thetax[i][0] - theta_lh;
    err_rh = thetax[j][0] - theta_rh;

    while ((abs(err_lk) > 0.005) && (abs(err_rk) > 0.005) && (abs(err_lh) > 0.005) && (abs(err_rh) > 0.005)) {
      theta_lk = (enc_lk * 0.3 * PI) / 180;
      theta_rk = (enc_rk * 0.3 * PI) / 180;
      theta_lh = (enc_lh * 0.2 * PI) / 180;
      theta_rh = (enc_lh * 0.2 * PI) / 180;

      prev_err_lk = err_lk;
      prev_err_rk = err_rk;
      prev_err_lh = err_lh;
      prev_err_rh = err_rh;

      err_lk = theta_lk - thetax[i][1];
      err_rk = theta_rk - thetax[j][1];
      err_lh = thetax[i][0] - theta_lh;
      err_rh = thetax[j][0] - theta_rh;

      diff_err_lk = err_lk - prev_err_lk;
      diff_err_rk = err_rk - prev_err_rk;
      diff_err_lh = err_lh - prev_err_lh;
      diff_err_rh = err_rh - prev_err_rh;

      err_sum_lk = err_lk + err_sum_lk;
      err_sum_rk = err_rk + err_sum_rk;
      err_sum_lh = err_lh + err_sum_lh;
      err_sum_rh = err_rh + err_sum_rh;

      corr_lk =  Kp * err_lk + Kv * diff_err_lk + Ki * err_sum_lk;
      corr_rk =  Kp * err_rk + Kv * diff_err_rk + Ki * err_sum_rk;
      corr_lh =  Kp * err_lh + Kv * diff_err_lh + Ki * err_sum_lh;
      corr_rh =  Kp * err_rh + Kv * diff_err_rh + Ki * err_sum_rh;

      if (corr_lh > 40) corr_lh = 50;
      if (corr_lh < -40) corr_lh = -50;
      if (corr_rh > 40) corr_rh = 50;
      if (corr_rh < -40) corr_rh = -50;

      if (corr_lk > 40) corr_lk = 50;
      if (corr_lk < -40) corr_lk = -50;
      if (corr_rk > 40) corr_rk = 50;
      if (corr_rk < -40) corr_rk = -50;

      //KNEE.............................................................................................
      //     if (abs(corr_lk) > 2) {
      if (err_lk > 0.01) {
        pwm_lk = 60 + corr_lk;
        digitalWrite(left_knee_dir, LOW);
        analogWrite(left_knee_pwm, 60);
      }
      if ( err_lk < -0.01) {
        pwm_lk = -60 + corr_lk;
        digitalWrite(left_knee_dir, HIGH);
        analogWrite(left_knee_pwm, 60);
      }
      //      }
      //      else {
      //        digitalWrite(left_knee_dir, HIGH);
      //        analogWrite(left_knee_pwm, 0);
      //        //        Serial.print(" ///");
      //      }
      //      if (abs(corr_rk) > 2) {
      if (err_rk > 0.01)
      {
        pwm_rk = 60 + corr_rk;
        digitalWrite(right_knee_dir, LOW);
        analogWrite(right_knee_pwm, 60);
      }
      if (err_rk < -0.01)
      {
        pwm_rk = -60 + corr_rk;
        digitalWrite(right_knee_dir, HIGH);
        analogWrite(right_knee_pwm, 60);
      }
      //      }
      //      else {
      //        digitalWrite(right_knee_dir, HIGH);
      //        analogWrite(right_knee_pwm, 0);
      //        //        Serial.print(" /// ");
      //      }
      //HIP.............................................................................................


      //      if ( abs(corr_lh) > 2) {
      if (err_lh > 0.01)
      {
        pwm_lh = 50 + corr_lh;
        digitalWrite(left_hip_dir, HIGH);
        analogWrite(left_hip_pwm, 65);
      }

      if (err_lh < -0.01)
      {
        pwm_lh = -50 + corr_lh;
        digitalWrite(left_hip_dir, LOW);
        analogWrite(left_hip_pwm, 65);
      }
      //      }
      //      else {
      //        digitalWrite(left_hip_dir, LOW);
      //        analogWrite(left_hip_pwm, 0);
      //        //        Serial.print(" ///");
      //      }

      //      if ( abs(corr_rh) > 2) {
      if (err_rh > 0.01)
      {
        pwm_rh = 50 + corr_rh;
        digitalWrite(right_hip_dir, LOW);
        analogWrite(right_hip_pwm, 65);
      }
      if (err_rh < -0.01)
      {
        pwm_rh = -50 + corr_rh;
        digitalWrite(right_hip_dir, HIGH);
        analogWrite(right_hip_pwm, 65);
      }
      //      }
      //      else
      //      {
      //        digitalWrite(right_hip_dir, HIGH);
      //        analogWrite(right_hip_pwm, 0);
      //        //        Serial.print(" /// ");
      //      }

      //      Serial.print(" ");
      //      Serial.print(corr_lk);
      //      Serial.print(" ");
      //      Serial.print(corr_rk);
      //      Serial.print(" ");
      //      Serial.print(corr_lh);
      //      Serial.print(" ");
      //      Serial.println(corr_rh);
      delay(500);
    }
    Serial.println(i);
    digitalWrite(left_knee_dir, HIGH);
    analogWrite(left_knee_pwm, 0);
    digitalWrite(right_knee_dir, HIGH);
    analogWrite(right_knee_pwm, 0);
    digitalWrite(left_hip_dir, LOW);
    analogWrite(left_hip_pwm, 0);
    digitalWrite(right_hip_dir, HIGH);
    analogWrite(right_hip_pwm, 0);
    
  }
}
