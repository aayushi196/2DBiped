void inverse_dynamics( double theta[4], double theta_dot[4], double l[4], double lc[4], double com_dot[2], double m[4], double Ic[4] ) {
  //short forms for cosine and sine
  
  double c1 = cos(theta[0]);
  double c2 = cos(theta[1]);
  double c12 = cos(theta[0] + theta[1]);
  double c3 = cos(theta[2]);
  double c4 = cos(theta[3]);
  double c34 = cos(theta[2] + theta[3]);

  double s1 = sin(theta[0]);
  double s2 = sin(theta[1]);
  double s12 = sin(theta[0] + theta[1]);
  double s3 = sin(theta[2]);
  double s4 = sin(theta[3]);
  double s34 = sin(theta[2] + theta[3]);

  //parameters for B matrix
  double A1 = (m[0] + m[1]) / 2;
  double A2 = Ic[0] + Ic[1] + m[0] * sq(lc[0]) + m[1] * sq(l[1]) + m[1] * sq(lc[1]) / 2 + m[1] * l[0] * lc[1] * c2;
  double A3 = (Ic[1] + m[1] * sq(lc[2])) / 2;

  double b1 = -(m[0] * lc[0] * s1 + m[1] * l[0] * s1 + m[1] * lc[1] * s12);
  double b2 = (Ic[1] + m[1] * sq(lc[1]) + m[1] * l[0] * lc[1] * c2);
  double b3 = -m[1] * lc[1] * s12;

  double C1 = (m[2] + m[3]) / 2;
  double C2 = (Ic[2] + Ic[3] + m[2] * sq(lc[2]) + m[3] * sq(l[3]) + m[3] * sq(lc[3])) / 2 + m[3] * l[2] * lc[3] * c4;
  double C3 = (Ic[3] + m[3] * sq(lc[3])) / 2;

  double D1 = -(m[2] * lc[2] * s3 + m[3] * l[2] * s3 + m[3] * lc[3] * s34);
  double D2 = (Ic[3] + m[3] * sq(lc[3]) + m[3] * l[2] * lc[3] * c4);
  double D3 = -m[3] * lc[3] * s34;

  //parameters for C matrix
  double dB1_dtheta1 = -m[0] * lc[0] * c1 - m[1] * l[0] * c1 - m[1] * lc[1] * c12;
  double dB3_dtheta1 = -m[1] * l[1] * c12;

  double dA2_dtheta2 = -m[1] * l[0] * lc[1] * s2;
  double dB1_dtheta2 = -m[1] * lc[1] * c12;
  double dB2_dtheta2 = -m[1] * l[0] * lc[1] * s2;
  double dB3_dtheta2 = -m[1] * l[1] * c12;

  double dD1_dtheta3 = -(m[2] * lc[2] * c3 + m[3] * l[2] * c3 + m[3] * lc[3] * c34);
  double dD3_dtheta3 = -m[3] * l[3] * c34;

  double dC2_dtheta4 = -m[3] * l[2] * lc[3] * s4;
  double dD1_dtheta4 = -m[3] * lc[3] * c34;
  double dD2_dtheta4 = -m[3] * l[2] * lc[3] * s4;
  double dD3_dtheta4 = -m[3] * l[3] * c34;


  double dA1_dt = 0 ;
  double dA2_dt = -m[1] * l[0] * lc[1] * s2 * theta_dot[1];
  double dA3_dt = 0;

  double dB1_dt = -(theta_dot[0] * (m[0] * lc[0] * c1 + m[1] * l[0] * c1 + m[1] * lc[1] * c12) + theta_dot[1] * (m[1] * lc[1] * c12));
  double dB2_dt = -m[1] * l[0] * lc[1] * s2 * theta_dot[1];
  double dB3_dt = -m[1] * l[1] * c12 * (theta_dot[0] + theta_dot[1]);

  double dC1_dt = 0;
  double dC2_dt = -m[3] * l[2] * lc[3] * s4 * theta_dot[3];
  double dC3_dt = 0;

  double dD1_dt = -(theta_dot[2] * (m[2] * lc[2] * c3 + m[3] * l[2] * c3 + m[3] * lc[3] * c34) + theta_dot[3] * (m[3] * lc[3] * c34));
  double dD2_dt = -m[3] * l[2] * lc[3] * s4 * theta_dot[3];
  double dD3_dt = -m[3] * l[3] * c34 * (theta_dot[2] + theta_dot[3]);

  //computation of the matrices
  dT_dt_dthetaDot[0][0] = 2 * dA2_dt;
  dT_dt_dthetaDot[0][1] = dB2_dt;
  dT_dt_dthetaDot[0][2] = 0;
  dT_dt_dthetaDot[0][3] = 0;
  dT_dt_dthetaDot[1][0] = dB2_dt;
  dT_dt_dthetaDot[1][1] = 2 * dA3_dt,
  dT_dt_dthetaDot[1][2] = 0;
  dT_dt_dthetaDot[1][3] = 0;
  dT_dt_dthetaDot[2][0] = 0;
  dT_dt_dthetaDot[2][1] = 0;
  dT_dt_dthetaDot[2][2] = 2 * dC2_dt;
  dT_dt_dthetaDot[2][3] = dD2_dt;
  dT_dt_dthetaDot[3][0] = 0;
  dT_dt_dthetaDot[3][1] = 0;
  dT_dt_dthetaDot[3][2] = dD2_dt;
  dT_dt_dthetaDot[3][3] = 2 * dC3_dt;

  dT_dtheta[0][0] = com_dot[0] * dB1_dtheta1;
  dT_dtheta[0][1] = com_dot[0] * dB3_dtheta1;
  dT_dtheta[0][2] = 0;
  dT_dtheta[0][3] = 0;
  dT_dtheta[1][0] = theta_dot[0] * dA2_dtheta2 + com_dot[0] * dB1_dtheta2;
  dT_dtheta[1][1] = theta_dot[0] * dB2_dtheta2 + com_dot[0] * dB3_dtheta2;
  dT_dtheta[1][2] = 0;
  dT_dtheta[1][3] = 0;
  dT_dtheta[2][0] = 0;
  dT_dtheta[2][1] = 0;
  dT_dtheta[2][2] = com_dot[0] * dD1_dtheta3;
  dT_dtheta[2][3] = com_dot[0] * dD3_dtheta3;
  dT_dtheta[3][0] = 0;
  dT_dtheta[3][1] = 0;
  dT_dtheta[3][2] = theta_dot[2] * dC2_dtheta4 + com_dot[0] * dD1_dtheta4;
  dT_dtheta[3][3] = theta_dot[2] * dD2_dtheta4 + com_dot[0] * dD3_dtheta4;

  B_q[0][0] = 2 * A2;
  B_q[0][1] = b2;
  B_q[0][2] = 0;
  B_q[0][3] = 0;
  B_q[1][0] = b2;
  B_q[1][1] = 2 * A3;
  B_q[1][2] = 0;
  B_q[1][3] = 0;
  B_q[2][0] = 0;
  B_q[2][1] = 0;
  B_q[2][2] = 2 * C2;
  B_q[2][3] = D2;
  B_q[3][0] = 0;
  B_q[3][1] = 0;
  B_q[3][2] = D2;
  B_q[3][3] = 2 * C3;

  D_x[0][0] = b1 * com_dot[1];
  D_x[1][0] = b3 * com_dot[1];
  D_x[2][0] = D1 * com_dot[1];
  D_x[3][0] = D3 * com_dot[1];

  //C_q_qdot = dT_dt_dthetaDot -dT_dtheta;
  Matrix.Subtract((mtx_type*)dT_dt_dthetaDot , (mtx_type*)dT_dtheta , 4, 4, (mtx_type*)C_q_qdot );

  E_x[0][0] = dB1_dt * com_dot[0];
  E_x[1][0] = dB3_dt * com_dot[0];
  E_x[2][0] = dD1_dt * com_dot[0];
  E_x[3][0] = dD3_dt * com_dot[0];

  G_q[0][0] = -g * (m[0] * lc[0] * c1 + m[1] * l[0] * c1 + m[1] * lc[1] * c12);
  G_q[1][0] = -g * m[1] * lc[1] * c12;
  G_q[2][0] = -g * (m[0] * lc[2] * c3 + m[1] * l[2] * c3 + m[1] * lc[3] * c34);
  G_q[3][0] = -g * m[1] * lc[3] * c34;

  //F_ext= G_q + D_x*com_dot[1] + E_x*com_dot[0];
  Matrix.Add((mtx_type*) D_x, (mtx_type*)E_x , 4, 1, (mtx_type*)D_E_xdot );
  Matrix.Add((mtx_type*) G_q, (mtx_type*)D_E_xdot , 4, 1, (mtx_type*)F_ext );
  //    F_ext[0][0] = G_q[0][0];
  //    F_ext[1][0] = G_q[1][0];
  //    F_ext[2][0] = G_q[2][0];
  //    F_ext[3][0] = G_q[3][0];
}
