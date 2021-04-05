%Thigh and leg length
clc
clear all
l= [0.32 0.32 0.32 0.32];
lc= [0.041 0.142 0.041 0.142];
m= [1.094 0.83 1.094 0.83]; 
Ic= [0.02 0.001 0.02 0.001];
P0= [0 , 0];
P1= [0 , 0]; P2= [0 , 0]; P4= [0 , 0]; P5= [ 0, 0];
theta = [ 1.200 0.467 1.481 0.467];
% theta = [ 1.45 0.40 1.45 0.40];
theta_dot = [0 0 0 0];
theta_ddot = [0 0 0 0];

%controller 
err= [0 0 0 0];
err_dot= [0 0 0 0];
err_sum= [0 0 0 0];
prev_err= [0 0 0 0];
x_real=0.0;
K= [ 25 10 0.1];
% K= [ 0 0 0];


%motor parameters
K_torque=0.35; %K_torque=0.65*K_omega
K_omega= 0.35;
Resistance= 2.2;
Current= [0 0 0 0];
Voltage= [0 0 0 0];
pwm= [0 0 0 0];

axis(gca, 'equal');
axis([-1 3 -1 1]);
grid on;

%COM parameters
z = 0.62;
g = 9.81;
Tc = sqrt(z/g);
origin =[0 0];
stride_length =0.15;
x_dot0 = 0.35;
x_0 = -stride_length/2;
h = 0.06;
stride_time = Tc*log((-stride_length/2 -Tc*x_dot0)/(stride_length/2 -Tc*x_dot0));
i=0;
flag=1;
dt= 0.01;

%swing leg tajectories
Po = [-x_0, z]; Pm= [0, z-h]; Pf=[x_0, z];
Po_dot = [-x_dot0, 0]; Pm_dot = [0.2, 0]; Pf_dot = [-x_dot0, 0];
Po_ddot = [-g*x_0/z, 0]; Pm_ddot = [0.2, 0]; Pf_ddot = [-g*x_0/z, 0];

[theta_1o, theta_2o] = inverse_kinematics_stance( Po(1), Po(2), l );
[theta_1odot, theta_2odot, theta_1oddot, theta_2oddot] = inverse_jacobian( Po_dot(1), Po_dot(2), Po_ddot(1), Po_ddot(2), theta_1o, theta_2o, l);
[theta_1m, theta_2m] = inverse_kinematics_stance( Pm(1), Pm(2), l);
[theta_1mdot, theta_2mdot, theta_1mddot, theta_2mddot] = inverse_jacobian( Pm_dot(1), Pm_dot(2), Pm_ddot(1), Pm_ddot(2), theta_1m, theta_2m, l);
[theta_1f, theta_2f] = inverse_kinematics_stance( Pf(1), Pf(2), l);
[theta_1fdot, theta_2fdot, theta_1fddot, theta_2fddot] = inverse_jacobian( Pf_dot(1), Pf_dot(2), Pf_ddot(1), Pf_ddot(2), theta_1f, theta_2f, l);

theta_1o=[theta_1o, theta_1odot, theta_1oddot];
theta_1m=[theta_1m, theta_1mdot, theta_1mddot];
theta_1f=[theta_1f, theta_1fdot, theta_1fddot];
theta_2o=[theta_2o, theta_2odot, theta_2oddot];
theta_2m=[theta_2m, theta_2mdot, theta_2mddot];
theta_2f=[theta_2f, theta_2fdot, theta_2fddot];

for step= 1:1:4
  i=0;
for t=0:dt:10
i =i+1;
%stance leg trajectory
x(i)= x_0*cosh(t/Tc) + Tc*x_dot0*sinh(t/Tc);
x_dot(i) = x_0*sinh(t/Tc)/Tc +x_dot0*cosh(t/Tc);
x_ddot(i) = g/z*x(i);
% x_real(i)= x(i);
% x_dot_real(i)= x_real(i);
% x_ddot_real(i)= x_ddot(i);
[thetaD_stance(1), thetaD_stance(2)] = inverse_kinematics_stance(x(i), z, l);
[thetaD_dot_stance(1), thetaD_dot_stance(2), thetaD_ddot_stance(1), thetaD_ddot_stance(2)] = inverse_jacobian(-x_dot(i), 0,-x_ddot(i), 0, thetaD_stance(1), thetaD_stance(2),l);

%swing leg trajectory
if (t< stride_time/2)
[a3, a4] = trajectory_generator(theta_1o, theta_1m, theta_2o, theta_2m, stride_time/2);
thetaD_swing(1) = a3(1) + a3(2)*t +a3(3)*t^2 + a3(4)*t^3 + a3(5)*t^4 +a3(6)*t^5;
thetaD_swing(2) = a4(1) + a4(2)*t +a4(3)*t^2 + a4(4)*t^3 + a4(5)*t^4 +a4(6)*t^5;

thetaD_dot_swing(1) = a3(2) +2*a3(3)*t + 3*a3(4)*t^2 + 4*a3(5)*t^3 +5*a3(6)*t^4;
thetaD_dot_swing(2) = a4(2) +2*a4(3)*t + 3*a4(4)*t^2 + 4*a4(5)*t^3 +5*a4(6)*t^4;

thetaD_ddot_swing(1) = 2*a3(3) + 6*a3(4)*t + 12*a3(5)*t^2 +20*a3(6)*t^3;
thetaD_ddot_swing(2) = 2*a4(3) + 6*a4(4)*t + 12*a4(5)*t^2 +20*a4(6)*t^3;

else
[a3, a4] = trajectory_generator(theta_1m, theta_1f, theta_2m, theta_2f, stride_time/2);
t1= t- stride_time/2;
thetaD_swing(1) = a3(1) + a3(2)*t1 +a3(3)*t1^2 + a3(4)*t1^3 + a3(5)*t1^4 +a3(6)*t1^5;
thetaD_swing(2) = a4(1) + a4(2)*t1 +a4(3)*t1^2 + a4(4)*t1^3 + a4(5)*t1^4 +a4(6)*t1^5;

thetaD_dot_swing(1) = a3(2) +2*a3(3)*t1 + 3*a3(4)*t1^2 + 4*a3(5)*t1^3 +5*a3(6)*t1^4;
thetaD_dot_swing(2) = a4(2) +2*a4(3)*t1 + 3*a4(4)*t1^2 + 4*a4(5)*t1^3 +5*a4(6)*t1^4;

thetaD_ddot_swing(1) = 2*a3(3) + 6*a3(4)*t1 + 12*a3(5)*t1^2 +20*a3(6)*t1^3;
thetaD_ddot_swing(2) = 2*a4(3) + 6*a4(4)*t1 + 12*a4(5)*t1^2 +20*a4(6)*t1^3;
end

if (flag ==1)
%Inverse Dynamics and and controller
thetaD(:,1:2)= thetaD_stance;
thetaD(:,3:4)= thetaD_swing;
thetaD_dot(:,1:2)= thetaD_dot_stance;
thetaD_dot(:,3:4)= thetaD_dot_swing;
thetaD_ddot(:,1:2)= thetaD_ddot_stance;
thetaD_ddot(:,3:4)= thetaD_ddot_swing;

% prev_x_real(i)=x_real(i);
% prev_x_dot_real(i)= x_dot_real(i);
% x_real(i) = -l(1)*cos(theta(1))-l(2)*cos(theta(1)+ theta(2));
% x_dot_real(i)= x_real(i)-prev_x_real(i);
% x_ddot_real(i)= x_dot_real(i)-prev_x_dot_real(i);

com_dot(1) = x_dot(i); com_dot(2) = x_ddot(i);
[B_q, C_q_qdot, G_q, D_x, E_x] = inverse_dynamics( theta, theta_dot, l, lc, com_dot, m, Ic );
F_ext= G_q + D_x*com_dot(2) + E_x*com_dot(1);
[ err, err_dot, err_sum, prev_err] = calculate_error( thetaD, theta, thetaD_dot, theta_dot, err_sum, prev_err);

% [Torque_dash, Torque, correction]= fuzzy_controller( B_q, C_q_qdot, F_ext, err, err_dot, thetaD_ddot, theta_dot);
[Torque_dash, Torque, correction]= controller( B_q, C_q_qdot, F_ext, err, err_dot, err_sum, thetaD_ddot, theta_dot, K);

%Calculating Thetas
theta_ddot= Torque_dash;
% theta_ddot = (inv(B_q)*(Torque - C_q_qdot*theta_dot' - F_ext))';
theta_dot = theta_dot + theta_ddot*dt;
theta= theta + theta_dot*dt;
 
%current and voltage
Current = Torque*1/K_torque;
Voltage = Current*Resistance + K_omega*theta_dot;
pwm= floor(Voltage*255/12);

%forward kinematics
P0 = [origin(1)+x(i) z];
P1 = P0 + [l(1)*cos(theta(1)) , -l(1)*sin(thetaD(1))];
P2 = P0 + [l(1)*cos(theta(1)) + l(2)*cos(theta(1)+theta(2)) , -l(1)*sin(theta(1)) - l(2)*sin(theta(1) + theta(2))];
P4 = P0 + [l(1)*cos(theta(3)) , -l(1)*sin(theta(3))];
P5 = P0 + [l(1)*cos(theta(3)) + l(2)*cos(theta(3)+theta(4)) , -l(1)*sin(theta(3))-l(2)*sin(theta(3) + theta(4))];
simulation( P0, P1, P2, P4, P5, origin, x(i), z );

if(x(i)>stride_length/2)
    flag =0;
    
     origin = origin + [stride_length 0];
     i=0;
     break;
end 
else
%Inverse Dynamics and and controller
thetaD(:,1:2)= thetaD_swing;
thetaD(:,3:4)= thetaD_stance;
thetaD_dot(:,1:2)= thetaD_dot_swing;
thetaD_dot(:,3:4)= thetaD_dot_stance;
thetaD_ddot(:,1:2)= thetaD_ddot_swing;
thetaD_ddot(:,3:4)= thetaD_ddot_stance;

% prev_x_real(i)=x_real(i);
% prev_x_dot_real(i)= x_dot_real(i);
% x_real(i) = -l(1)*cos(theta(1))-l(2)*cos(theta(1)+ theta(2));
% x_dot_real(i)= x_real(i)-prev_x_real(i);
% x_ddot_real(i)= x_dot_real(i)-prev_x_dot_real(i);

com_dot(1) =x_dot(i); com_dot(2) = x_ddot(i);
[B_q, C_q_qdot, G_q, D_x, E_x] = inverse_dynamics( theta, theta_dot, l, lc, com_dot, m, Ic );
F_ext= G_q + D_x*com_dot(2) + E_x*com_dot(1);
[ err, err_dot, err_sum, prev_err] = calculate_error( thetaD, theta, thetaD_dot, theta_dot, err_sum, prev_err);

%  [Torque_dash, Torque, correction]= fuzzy_controller( B_q, C_q_qdot, F_ext, err, err_dot, thetaD_ddot, theta_dot);
[Torque_dash, Torque, correction]= controller( B_q, C_q_qdot, F_ext, err, err_dot, err_sum, thetaD_ddot, theta_dot, K);

%Calculating Thetas
theta_ddot= Torque_dash;
% theta_ddot = (inv(B_q)*(Torque - C_q_qdot*theta_dot' - F_ext))';
theta_dot = theta_dot + theta_ddot*dt;
theta= theta + theta_dot*dt;

%current and voltage
Current = Torque*1/K_torque;
Voltage = Current*Resistance+K_omega*theta_dot;
pwm= round(Voltage*255/12);

%forward kinematics
P0 = [origin(1)+x(i) z];
P1 =  P0 + [l(1)*cos(theta(1)) , -l(1)*sin(theta(1))];
P2 = P0 +[l(1)*cos(theta(1)) + l(2)*cos(theta(1)+theta(2)) , -l(1)*sin(theta(1)) - l(2)*sin(theta(1) + theta(2))];
P4 = P0 + [l(1)*cos(theta(3)) , -l(1)*sin(theta(3))];
P5 = P0 + [l(1)*cos(theta(3)) + l(2)*cos(theta(3)+theta(4)) , -l(1)*sin(theta(3))-l(2)*sin(theta(3) + theta(4))];
simulation( P0, P1, P2, P4, P5, origin, x(i), z );
if(x(i)>stride_length/2)
     origin = origin + [stride_length 0]; 
     i=0;
    flag = 1;
    break; 
end
end
%plot velocities and acceleration
time(i)=t;
% time(64+i)=0.63 + t;

% origin_1(j,1) =P0(1);
theta_3(i)= thetaD(3);
theta_4(i)= thetaD(4);
theta_1(i)= thetaD(1);
theta_2(i)= thetaD(2);

% thet_dot3(i) = thetaD_dot(3);
% thet_dot4(i) = thetaD_dot(4);
% thet_dot1(i) = thetaD_dot(1);
% thet_dot2(i) = thetaD_dot(2);
% 
% thet_ddot3(i) = thetaD_ddot(3);
% thet_ddot4(i) = thetaD_ddot(4);
% thet_ddot1(i) = thetaD_ddot(1);
% thet_ddot2(i) = thetaD_ddot(2);

torq_1(i)=Torque(1);
torq_2(i)=Torque(2);
torq_3(i)=Torque(3);
torq_4(i)=Torque(4);

% torq_dash_1(i)=Torque_dash(1);
% torq_dash_2(i)=Torque_dash(2);
% torq_dash_3(i)=Torque_dash(3);
% torq_dash_4(i)=Torque_dash(4);

curr_1(i)=Current(1);
curr_2(i)=Current(2);
curr_3(i)=Current(3);
curr_4(i)=Current(4);

volt_1(i)=pwm(1);
volt_2(i)=pwm(2);
volt_3(i)=pwm(3);
volt_4(i)=pwm(4);

err_1(i) = err(1);
err_2(i) = err(2);
err_3(i) = err(3);
err_4(i) = err(4);

err_dot_1(i) = err_dot(1);
err_dot_2(i) = err_dot(2);
err_dot_3(i) = err_dot(3);
err_dot_4(i) = err_dot(4);

corr_1(i) = correction(1);
corr_2(i) = correction(2);
corr_3(i) = correction(3);
corr_4(i) = correction(4);

end

end
%  plot (time, torq_1, 'Color', 'blue');
%  hold on
%  plot (time, torq_2, 'Color', 'red');
%  hold on
%  plot (time, torq_3, 'Color', 'green');
%  hold on
%  plot (time, torq_4, 'Color', 'yellow');
%  hold on

  subplot(2,2,1);  
  plot(time, err_1, 'LineWidth',3 , 'Color', 'blue' );
  title('Subplot 1: stance-hip');
  ylabel('error in radians'); 
  xlabel('time in seconds') ;
  subplot(2,2,2);
  plot(time, err_2, 'LineWidth',3, 'Color', 'red' );
  title('Subplot 2:  stance-knee');
  ylabel('error in radians'); 
  xlabel('time in seconds') ;
  subplot(2,2,3);
  plot(time, err_3, 'LineWidth',3, 'Color', 'green' );
  title('Subplot 3:  swing-hip');
  ylabel('error in radians'); 
  xlabel('time in seconds') ;
  subplot(2,2,4);
  plot(time, err_4, 'LineWidth',3 , 'Color', 'yellow');
  title('Subplot 4:  swing-knee');
  ylabel('error in radians'); 
  xlabel('time in seconds') ;
%   subplot(2,2,1);  
%   plot(time, torq_1, 'LineWidth',3 , 'Color', 'blue' );
%   title('Subplot 1: stance-hip');
%   ylabel('torque in N-m'); 
%   xlabel('time in seconds') ;
%   subplot(2,2,2);
%   plot(time, torq_2, 'LineWidth',3, 'Color', 'red' );
%   title('Subplot 2:  stance-knee');
%   ylabel('torque in N-m'); 
%   xlabel('time in seconds') ;
%   subplot(2,2,3);
%   plot(time, torq_3, 'LineWidth',3, 'Color', 'green' );
%   title('Subplot 3:  swing-hip');
%   ylabel('torque in N-m'); 
%   xlabel('time in seconds') ;
%   subplot(2,2,4);
%   plot(time, torq_4, 'LineWidth',3 , 'Color', 'yellow');
%   title('Subplot 4:  swing-knee');
%   ylabel('torque in N-m'); 
%   xlabel('time in seconds') ;
