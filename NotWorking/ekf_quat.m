% Extended Kalman Filter (EKF) for Attitude Estimation using Quaternions

clc;
clear;
close all;

% Constants
GYRO_NOISE = 0.015; % rad/sec
ACC_NOISE = 1.0; % m/s^2
g = 9.81; % gravity

% Load IMU data
imu_data = readmatrix('G:\REPOS\ESEKF_IMU\data\imu_noise.txt');
time = imu_data(:,1);
gyro = imu_data(:,2:4);
accel = imu_data(:,5:7);

% Load ground truth
gt_data = readmatrix('G:\REPOS\ESEKF_IMU\data\traj_gt.txt');
q = cell(size(time));
q{1} = gt_data(1, 2:5)'; % Use ground truth as the first state

% Time configuration
t0 = time(1);
tf = time(end);
dt = mean(diff(time));
N = length(time);

% Generate sensor data
w = gyro;
a = accel;

% EKF estimation
[estimated_q, estimated_b, b_err] = estimate(time, w, a);

% Convert to Euler angles
rpy = quaternion_to_euler_degrees(cell2mat(q)');
estimated_rpy = quaternion_to_euler_degrees(estimated_q);

% Plot results
plot_results(time, w, a, rpy, estimated_rpy);

%% Functions

function [estimated_q, estimated_b, b_err] = estimate(time, w, a)
    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);
    roll = atan2(2 * (qw * qx + qy * qz), 1 - 2 * (qx^2 + qy^2));
    pitch = asin(2 * (qw * qy - qz * qx));
    yaw = atan2(2 * (qw * qz + qx * qy), 1 - 2 * (qy^2 + qz^2));
end
