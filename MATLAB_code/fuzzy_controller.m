function [ Torque_dash, Torque, correction] = fuzzy_controller( B_q, C_q_qdot, F_ext, err, err_dot,  thetaD_ddot, theta_dot)

% fis = readfis('Triangular.fis');
fis = readfis('Gaussian.fis');
for j= 1:1:4
 correction(j)= evalfis(fis,[err(j) err_dot(j)]);
% correction(j) = 0.0;
end

Torque_dash = thetaD_ddot + correction;
Torque = B_q*Torque_dash' + C_q_qdot*theta_dot' + F_ext;

% for i= 1:1:4
%  if(Torque_dash(i)>3) 
%     Torque_dash(i)=3 ;
%  end
%  if(Torque_dash(i)<-3) 
%     Torque_dash(i)=-3 ;
%  end
% end
end
