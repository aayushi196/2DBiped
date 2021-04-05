function [ err, err_dot, err_sum, prev_err] = calculate_error( thetaD, theta, thetaD_dot, theta_dot, err_sum, prev_err )
err= thetaD - theta;
err_dot= thetaD_dot - theta_dot;
derr= err- prev_err;
err_sum= err_sum + err;
prev_err= err;
end

