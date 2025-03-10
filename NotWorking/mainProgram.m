clc; clear; close all;

% Load IMU data (Assuming CSV with time, gyro_x, gyro_y, gyro_z, accel_x, accel_y, accel_z)
imu_data = readmatrix('G:\REPOS\ESEKF_IMU\data\imu_noise.txt');
% load ground truth
gt_data = readmatrix('G:\REPOS\ESEKF_IMU\data\traj_gt.txt');;
time = imu_data(:,1);
gyro = imu_data(:,2:4);
accel = imu_data(:,5:7);
dt = mean(diff(time)); % Assuming uniform sampling

% Initialize EKF
% q_est = [1; 0; 0; 0]; % Initial quaternion (no rotation)
q_est = [gt_data(1,8); gt_data(1,5:7)']; %set initial to GT
P = eye(4) * 1; % Initial covariance matrix
% Q = diag([0.019,0.019,0.019,0.019])*1e-12; % Process noise covariance
Q=diag([0.015^2,0.015^2,0.015^2,0.015^2]);
R = eye(3) * 0.015^2; % Measurement noise covariance
use_multiplicative = true; % Toggle between additive and multiplicative EKF

% Storage for results
q_history = zeros(length(time), 4);

% Run EKF
for i = 1:length(time)
    [q_est, P] = ekf_attitude(q_est, P, gyro(i,:)', accel(i,:)', dt, Q, R, use_multiplicative);
    q_history(i, :) = q_est';
end

% Convert quaternions to Euler angles (roll, pitch, yaw)
% euler_angles = quat2eul(q_history, 'ZYX');
euler_angles = quat2eul([q_history(:,4) q_history(:,1:3)], 'ZYX');
euler_anglesGT = rad2deg(quat2eul([gt_data(:,8) gt_data(:,5:7)],'ZYX'));
% Plot results
figure;
subplot(3,1,1);
plot(time, euler_angles(:,3) * 180/pi, 'r', time,euler_anglesGT(:,3), 'k' ); 
ylabel('Roll (deg)'); grid on;
subplot(3,1,2);
plot(time, euler_angles(:,2) * 180/pi, 'g', time,euler_anglesGT(:,2), 'k' ); 
ylabel('Pitch (deg)'); grid on;
subplot(3,1,3);
plot(time, euler_angles(:,1) * 180/pi, 'b', time,euler_anglesGT(:,1), 'k' ); 
ylabel('Yaw (deg)'); xlabel('Time (s)'); grid on;

title('Estimated Attitude using EKF');