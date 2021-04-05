void trajectory_generator_im(double theta1_init[3], double theta1_middle[3], double theta2_init[3], double theta2_middle[3], double tf){
a1_im[0] = theta1_init[0];
a1_im[1] = theta1_init[1];
a1_im[2] = theta1_init[2]/2;
a1_im[3] = (20*(theta1_middle[0] -theta1_init[0]) -(8*theta1_middle[1] + 12*theta1_init[1])*tf -(3*theta1_init[2]-theta1_middle[2])*sq(tf))/(2*pow(tf,3));
a1_im[4] = (30*(theta1_init[0] -theta1_middle[0]) +(14*theta1_middle[1] + 16*theta1_init[1])*tf +(3*theta1_init[2]-2*theta1_middle[2])*sq(tf))/(2*pow(tf,4));
a1_im[5] = (12*(theta1_middle[0] -theta1_init[0]) -6*(theta1_middle[1] + theta1_init[1])*tf -(theta1_init[2] -theta1_middle[2])*sq(tf))/(2*pow(tf,5));

a2_im[0] = theta2_init[0];
a2_im[1] = theta2_init[1];
a2_im[2] = theta2_init[2]/2;
a2_im[3] = (20*(theta2_middle[0] -theta2_init[0]) -(8*theta2_middle[1] + 12*theta2_init[1])*tf -(3*theta2_init[2]-theta2_middle[2])*sq(tf))/(2*pow(tf,3));
a2_im[4] = (30*(theta2_init[0] -theta2_middle[0]) +(14*theta2_middle[1] + 16*theta2_init[1])*tf +(3*theta2_init[2]-2*theta2_middle[2])*sq(tf))/(2*pow(tf,4));
a2_im[5] = (12*(theta2_middle[0] -theta2_init[0]) -6*(theta2_middle[1] + theta2_init[1])*tf -(theta2_init[2] -theta2_middle[2])*sq(tf))/(2*pow(tf,5));

}

void trajectory_generator_mf(double  theta1_middle[3], double theta1_final[3], double theta2_middle[3], double theta2_final[3], double tf){
a1_mf[0] = theta1_middle[0];
a1_mf[1] = theta1_middle[1];
a1_mf[2] = theta1_middle[2]/2;
a1_mf[3] = (20*(theta1_final[0] -theta1_middle[0]) -(8*theta1_final[1] + 12*theta1_middle[1])*tf -(3*theta1_middle[2]-theta1_final[2])*sq(tf))/(2*pow(tf,3));
a1_mf[4] = (30*(theta1_middle[0] -theta1_final[0]) +(14*theta1_final[1] + 16*theta1_middle[1])*tf +(3*theta1_middle[2]-2*theta1_final[2])*sq(tf))/(2*pow(tf,4));
a1_mf[5] = (12*(theta1_final[0] -theta1_middle[0]) -6*(theta1_final[1] + theta1_middle[1])*tf -(theta1_middle[2] -theta1_final[2])*sq(tf))/(2*pow(tf,5));

a2_mf[0] = theta2_middle[0];
a2_mf[1] = theta2_middle[1];
a2_mf[2] = theta2_middle[2]/2;
a2_mf[3] = (20*(theta2_final[0] -theta2_middle[0]) -(8*theta2_final[1] + 12*theta2_middle[1])*tf -(3*theta2_middle[2]-theta2_final[2])*sq(tf))/(2*pow(tf,3));
a2_mf[4] = (30*(theta2_middle[0] -theta2_final[0]) +(14*theta2_final[1] + 16*theta2_middle[1])*tf +(3*theta2_middle[2]-2*theta2_final[2])*sq(tf))/(2*pow(tf,4));
a2_mf[5] = (12*(theta2_final[0] -theta2_middle[0]) -6*(theta2_final[1] + theta2_middle[1])*tf -(theta2_middle[2] -theta2_final[2])*sq(tf))/(2*pow(tf,5));

}
