void motor_initialize(){
  pinMode(right_hip_pwm, OUTPUT);
  pinMode(right_hip_dir, OUTPUT);
  pinMode(left_hip_pwm, OUTPUT);
  pinMode(left_hip_dir, OUTPUT);

  pinMode(right_knee_pwm, OUTPUT);
  pinMode(right_knee_dir, OUTPUT);
  pinMode(left_knee_pwm, OUTPUT);
  pinMode(left_knee_dir, OUTPUT);
 }

  void update_thetas_velocities(){
    theta[0] = (enc_lh*0.3*pulley_ratio*PI)/180; 
    theta[1] = (enc_lk*0.3*PI)/180;
    theta[2] = (enc_rh*0.3*pulley_ratio*PI)/180;
    theta[3] = (enc_rk*0.3*PI)/180;

    theta_dot[0]= theta[0]-prev_theta[0];
    theta_dot[1]= theta[1]-prev_theta[1];
    theta_dot[2]= theta[2]-prev_theta[2];
    theta_dot[3]= theta[3]-prev_theta[3];
    
    theta_ddot[0]= theta_dot[0]-prev_theta_dot[0];
    theta_ddot[1]= theta_dot[1]-prev_theta_dot[1];
    theta_ddot[2]= theta_dot[2]-prev_theta_dot[2];
    theta_ddot[3]= theta_dot[3]-prev_theta_dot[3];

    prev_theta[0]= theta[0];
    prev_theta[1]= theta[1];
    prev_theta[2]= theta[2];
    prev_theta[3]= theta[3];

    prev_theta_dot[0]= theta_dot[0];
    prev_theta_dot[1]= theta_dot[1];
    prev_theta_dot[2]= theta_dot[2];
    prev_theta_dot[3]= theta_dot[3];
    }
 
void motor_run(int motor_dir_pin,int motor_pwm_pin,double pwm){
  if(pwm>=0)
 {
  digitalWrite(motor_dir_pin,HIGH);
  analogWrite(motor_pwm_pin,abs(pwm));
 }
 else
 {
   digitalWrite(motor_dir_pin,LOW);
   analogWrite(motor_pwm_pin,abs(pwm)); 
 }
 }
