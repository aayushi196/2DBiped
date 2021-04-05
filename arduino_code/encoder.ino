void encoder_initialize(){
 pinMode(right_hip_encA, INPUT);
  pinMode(right_hip_encB, INPUT);
  pinMode(left_hip_encA, INPUT);
  pinMode(left_hip_encB, INPUT);

  pinMode(right_knee_encA, INPUT);
  pinMode(right_knee_encB, INPUT);
  pinMode(left_knee_encA, INPUT);
  pinMode(left_knee_encB, INPUT);

  digitalWrite(right_hip_encA, HIGH);
  digitalWrite(right_hip_encB, HIGH);
  digitalWrite(left_hip_encA, HIGH);
  digitalWrite(left_hip_encB, HIGH);

  digitalWrite(right_knee_encA, HIGH);
  digitalWrite(right_knee_encB, HIGH);
  digitalWrite(left_knee_encA, HIGH);
  digitalWrite(left_knee_encB, HIGH);
  
  interrupts(); 
 
  attachInterrupt(digitalPinToInterrupt(left_hip_encA),int_left_hip,CHANGE);
  attachInterrupt(digitalPinToInterrupt(right_hip_encA),int_right_hip,CHANGE);
  attachInterrupt(digitalPinToInterrupt(left_knee_encA),int_left_knee,CHANGE);
  attachInterrupt(digitalPinToInterrupt(right_knee_encA),int_right_knee,CHANGE);
  
}

void int_left_hip()
{  
  if(digitalRead(left_hip_encA) == HIGH) 
    enc_lh= (digitalRead(left_hip_encB) == LOW)?enc_lh-1:enc_lh+1;
  if(digitalRead(left_hip_encA) == LOW)                                        
    enc_lh= (digitalRead(left_hip_encB) == HIGH)?enc_lh-1:enc_lh+1;
  
}

void int_right_hip()
{
  if(digitalRead(right_hip_encA) == HIGH) 
    enc_rh= (digitalRead(right_hip_encB) == LOW)?enc_rh+1:enc_rh-1;
  if(digitalRead(right_hip_encA) == LOW)                                        
    enc_rh= (digitalRead(right_hip_encB) == HIGH)?enc_rh+1:enc_rh-1;
  
}

void int_left_knee()
{
  if(digitalRead(left_knee_encA) == HIGH) 
    enc_lk= (digitalRead(left_knee_encB) == LOW)?enc_lk+1:enc_lk-1;
  if(digitalRead(left_knee_encA) == LOW)                                        
    enc_lk= (digitalRead(left_knee_encB) == HIGH)?enc_lk+1:enc_lk-1;
  
}

void int_right_knee()
{
  if(digitalRead(right_knee_encA) == HIGH) 
    enc_rk= (digitalRead(right_knee_encB) == LOW)?enc_rk-1:enc_rk+1;
  if(digitalRead(right_knee_encA) == LOW)                                        
    enc_rk= (digitalRead(right_knee_encB) == HIGH)?enc_rk-1:enc_rk+1;
  
}
