function [q_est, P] = ekf_attitude(q_prev, P, gyro, accel, dt, Q, R, use_multiplicative)
    % EKF Attitude Estimation using Quaternions
    % Inputs:
    % q_prev - Previous quaternion estimate [qw, qx, qy, qz]'
    % P - Previous covariance matrix
    % gyro - Gyroscope measurements [wx, wy, wz] (rad/s)
    % accel - Accelerometer measurements [ax, ay, az] (m/s^2)
    % dt - Time step (s)
    % Q - Process noise covariance matrix
    % R - Measurement noise covariance matrix
    % use_multiplicative - Boolean flag to toggle between additive and multiplicative EKF
    % Outputs:
    % q_est - Estimated quaternion
    % P - Updated covariance matrix
    
    % Normalize accelerometer
    accel = accel / norm(accel);
    
    % Predict step
    omega = [0, -gyro(1), -gyro(2), -gyro(3);
             gyro(1), 0, gyro(3), -gyro(2);
             gyro(2), -gyro(3), 0, gyro(1);
             gyro(3), gyro(2), -gyro(1), 0];
    
    q_pred = q_prev + 0.5 * dt * omega * q_prev;
    q_pred = q_pred / norm(q_pred);
    
    % Jacobian of the process model
    F = eye(4) + 0.5 * dt * omega;
    P = F * P * F' + Q;
    
    % Measurement update
    g_ref = [0; 0; 0; 1];  % Reference gravity in quaternion form
    h = quatrotate(q_pred', g_ref(2:4));  % Rotate reference gravity
    
    H = 2 * [q_pred(3), -q_pred(4), q_pred(1), -q_pred(2);
             -q_pred(2), -q_pred(1), -q_pred(4), -q_pred(3);
             -q_pred(1), q_pred(2), q_pred(3), q_pred(4)];
    
    K = P * H' / (H * P * H' + R);
    q_update = K * (accel - h);
    
    if use_multiplicative
        % Convert correction to a small quaternion
        delta_q = [1; 0.5 * q_update];
        delta_q = delta_q / norm(delta_q);
        q_est = quatmultiply(delta_q', q_pred');
        q_est = q_est';
    else
        q_est = q_pred + q_update;
        q_est = q_est / norm(q_est);
    end
    
    P = (eye(4) - K * H) * P;
end
