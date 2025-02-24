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
% use_multiplicative - Boolean flag for additive vs. multiplicative EKF
% Outputs:
% q_est - Estimated quaternion
% P - Updated covariance matrix

% Predict Step: Quaternion Propagation using Gyroscope
omega_mat = 0.5 * [ 0,    -gyro(1), -gyro(2), -gyro(3);
    gyro(1),  0,     gyro(3), -gyro(2);
    gyro(2), -gyro(3),  0,    gyro(1);
    gyro(3),  gyro(2), -gyro(1),  0];

q_pred = q_prev + dt * omega_mat * q_prev;
q_pred = q_pred / norm(q_pred); % Normalize quaternion
% Define gyroscope noise covariance matrix
Q_gyro = diag([0.015^2, 0.015^2, 0.015^2])/0.1; % Gyro noise

% Compute process noise for quaternion update
G = 0.5 * [ -q_pred(2), -q_pred(3), -q_pred(4);
    q_pred(1), -q_pred(4),  q_pred(3);
    q_pred(4),  q_pred(1), -q_pred(2);
    -q_pred(3),  q_pred(2),  q_pred(1)];

% Compute Q in quaternion space
Q = G * Q_gyro * G' * dt^2;

% Covariance update
% P = F * P * F' + Q;
% Jacobian of the Process Model (State Transition Matrix)
F = eye(4) + dt * omega_mat;
P = F * P * F' + Q; % Covariance update

% Measurement Update: Accelerometer
g_ref = [0; 0; -9.81];  % Gravity in world frame (real-world units)
h = quatrotate(q_pred', g_ref')';  % Rotate gravity into body frame

% Jacobian of Measurement Model
% H = 2 * [-q_pred(3),  q_pred(4), -q_pred(1), q_pred(2);
%           q_pred(2),  q_pred(1),  q_pred(4), q_pred(3);
%           q_pred(1), -q_pred(2), -q_pred(3), q_pred(4)];
H = -2 * [ q_pred(3), -q_pred(4),  q_pred(1), -q_pred(2);
    -q_pred(2), -q_pred(1), -q_pred(4), -q_pred(3);
    -q_pred(1),  q_pred(2),  q_pred(3),  q_pred(4)];


% Kalman Gain Computation
S = H * P * H' + R;
K = P * H' / S;

% Compute Innovation (Residual)
y = accel - h;  % Now accel is in m/s², consistent with g_ref

% Update Quaternion Estimate
q_update = K * y;

if use_multiplicative
    % Convert correction to small quaternion delta
    delta_q = [1; 0.5 * q_update(1:3)]; % Small rotation quaternion
    delta_q = delta_q / norm(delta_q);  % Normalize correction

    % Apply quaternion multiplication for correction
    q_est = quatmultiply(delta_q', q_pred')';
else
    % Additive Update
    q_est = q_pred + q_update;
end

% Normalize Quaternion
q_est = q_est / norm(q_est);

% Covariance Update
P = (eye(4) - K * H) * P;
end
