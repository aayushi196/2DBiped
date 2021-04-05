void inverse_kinematics(  double P[2], double l[4], double theta[2] ) {
  double num, den, sin_1, cos_1, thet_1, thet_2, thet1, x, y;
  x = P[0];
  y = P[1];
  num = (sq(x) + sq(y) - sq(l[0]) - sq(l[1]));
  den = 2 * l[0] * l[1];
  thet_2 = acos(num/ den);

  sin_1 = ((l[1] + l[0] * cos(thet_2)) * y - l[0] * sin(thet_2) * x) / ( sq(x) + sq(y));
  cos_1 = ((l[1] + l[0] * cos(thet_2)) * x + l[0] * sin(thet_2) * y) / ( sq(x) + sq(y));
  thet_1 = atan2(sin_1, cos_1);
  thet1 = PI - thet_1 - thet_2;

  theta[0] = thet1;
  theta[1] = thet_2;

}

void inverse_jacobian(double P_dot[2], double P_ddot[2], double thet[2], double l[4], double theta_vel_acc[4] ) {

  mtx_type J[2][2];
  mtx_type J_dot[2][2];
  mtx_type velocity[2][1];
  mtx_type acc[2][1];
  mtx_type theta_dot[2][1];
  mtx_type theta_ddot[2][1];
  mtx_type J_dot_theta_dot[2][1];
  mtx_type acc_J_dot_theta_dot[2][1];
  double theta_1, theta_2, theta_1dot, theta_2dot;
  double x_dot, y_dot, x_ddot, y_ddot;

  theta_1 = thet[0];
  theta_2 = thet[1];
  x_dot = P_dot[0];
  y_dot = P_dot[1];
  x_ddot = P_ddot[0];
  y_ddot = P_ddot[1];

  J[0][0] = -(l[0] * sin(theta_1)) - (l[1] * sin(theta_1 + theta_2));
  J[0][1] = -(l[1] * sin(theta_1 + theta_2));
  J[1][0] = -(l[0] * cos(theta_1)) - (l[1] * cos(theta_1 + theta_2));
  J[1][1] = -(l[1] * cos(theta_1 + theta_2));

  Matrix.Invert((mtx_type*)J, 2); //inverts the matrix and saves in the same variable  Matrix.Invert((mtx_type*)A, N)
  velocity[0][0] = x_dot;
  velocity[1][0] = y_dot;
  //theta_dot = inv(J)*velocity;
  Matrix.Multiply((mtx_type*)J, (mtx_type*)velocity, 2, 2, 1, (mtx_type*)theta_dot);
  
  theta_1dot = theta_dot[0][0];
  theta_2dot = theta_dot[1][0];
  J_dot[0][0] = -(l[0] * cos(theta_1) * theta_1dot) - (l[0] * cos(theta_1 + theta_2) * (theta_1dot + theta_2dot));
  J_dot[0][1] = -(l[1] * cos(theta_1 + theta_2) * (theta_1dot + theta_2dot));
  J_dot[1][0] = (l[0] * sin(theta_1) * theta_1dot) + (l[1] * sin(theta_1 + theta_2) * (theta_1dot + theta_2dot));
  J_dot[1][1] = (l[1] * sin(theta_1 + theta_2) * (theta_1dot + theta_2dot));

  acc[0][0] = x_ddot;
  acc[1][0] = y_ddot;
  //theta_ddot= inv(J)*(acc - J_dot*theta_dot);
  Matrix.Multiply((mtx_type*)J_dot, (mtx_type*)theta_dot, 2, 2, 1, (mtx_type*)J_dot_theta_dot);
  Matrix.Subtract((mtx_type*) acc, (mtx_type*)J_dot_theta_dot , 2, 1, (mtx_type*)acc_J_dot_theta_dot );
  Matrix.Multiply((mtx_type*)J, (mtx_type*)acc_J_dot_theta_dot, 2, 2, 1, (mtx_type*)theta_ddot);
//  Matrix.Print((mtx_type*)theta_ddot, 2, 1, "theta_ddot");
  
  theta_vel_acc[0] = theta_1dot;
  theta_vel_acc[1] = theta_2dot;
  theta_vel_acc[2] = theta_ddot[0][0];
  theta_vel_acc[3] = theta_ddot[1][0];

}
