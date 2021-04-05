function [ Torque_dash, Torque,correction] = controller( B_q, C_q_qdot, F_ext, err, err_dot, err_sum,  thetaD_ddot, theta_dot, K)
Kp= K(1);
Kv= K(2);
Ki= K(3);
correction = Kp*err + Kv*err_dot + Ki*err_sum;
Torque_dash = thetaD_ddot + Kp*err + Kv*err_dot + Ki*err_sum;
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

