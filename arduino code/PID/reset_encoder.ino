void encoder_reset() {
  // Initial Resetting of encoders
  // for knee
  while (digitalRead(limit_switch_lk) == LOW)
  { digitalWrite(left_knee_dir, HIGH); // The direction where LIMIT SWITCH is pressed is LOW
    analogWrite(left_knee_pwm, 60);
  }
  enc_lk = 108 / 0.3;
  theta_lk = (enc_lk * 0.3 * PI) / 180;
  err_lk = theta_lk - theta[1];

  while (err_lk > 0.01)
  {
    theta_lk = (enc_lk * 0.3 * PI) / 180;
    err_lk = theta_lk - theta[1];
    digitalWrite(left_knee_dir, LOW);
    analogWrite(left_knee_pwm, 50);
  }
  digitalWrite(left_knee_dir, LOW);
  analogWrite(left_knee_pwm, 0);

  while (digitalRead(limit_switch_rk) == LOW)
  { digitalWrite(right_knee_dir, HIGH); // The direction where LIMIT SWITCH is pressed is LOW
    analogWrite(right_knee_pwm, 60);
  }

  enc_rk = 104 / 0.3;
  theta_rk = (enc_rk * 0.3 * PI) / 180;
  err_rk = theta_rk - theta[3];

  while (err_rk > 0.01)
  {
    theta_rk = (enc_rk * 0.3 * PI) / 180;
    err_rk = theta_rk - theta[3];
    digitalWrite(right_knee_dir, LOW);
    analogWrite(right_knee_pwm, 50);
  }

  digitalWrite(right_knee_dir, LOW);
  analogWrite(right_knee_pwm, 0);

    // for hip
    while (digitalRead(limit_switch_lh) == LOW)
    { digitalWrite(left_hip_dir, LOW); // The direction where LIMIT SWITCH is pressed is LOW
      analogWrite(left_hip_pwm, 30);
    }
    enc_lh = 145;
    theta_lh = (enc_lh * 0.2 * PI) / 180;
    err_lh = theta[0] - theta_lh;
  
    while (err_lh > 0.01)
    {
      theta_lh = (enc_lh * 0.2 * PI) / 180;
      err_lh = theta[0] - theta_lh;
      digitalWrite(left_hip_dir, HIGH);
      analogWrite(left_hip_pwm, 60);
    }
    digitalWrite(left_hip_dir, HIGH);
    analogWrite(left_hip_pwm, 0);

  while (digitalRead(limit_switch_rh) == LOW)
  { digitalWrite(right_hip_dir, LOW); // The direction where LIMIT SWITCH is pressed is LOW
    analogWrite(right_hip_pwm, 80);
  }
  enc_rh = 160;
  theta_rh = (enc_rh * 0.2 * PI) / 180;
  err_rh = theta[2] - theta_rh;
  
  while (err_rh > 0.01)
  {
    theta_rh = (enc_rh * 0.2 * PI) / 180;
    err_rh = theta[2] - theta_rh;
    digitalWrite(right_hip_dir, HIGH);
    analogWrite(right_hip_pwm, 60);
  }

  digitalWrite(right_hip_dir, HIGH);
  analogWrite(right_hip_pwm, 0);

}
