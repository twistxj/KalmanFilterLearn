
    % Load IMU parameters
    addpath(genpath('G:\REPOS\readyaml'));
    imu_parameters = load_imu_parameters();
    
    % Load IMU and ground truth data
    imu_data = load('./data/imu_noise.txt');
    gt_data = load('./data/traj_gt.txt');
    
    % Initialize state
    init_nominal_state = zeros(19,1);
    init_nominal_state(1:10) = gt_data(1, 2:11)'; % Init p, q, v
    init_nominal_state(11:13) = 0;  % Init accelerometer bias
    init_nominal_state(14:16) = 0;  % Init gyroscope bias
    init_nominal_state(17:19) = [0; 0; -9.81]; % Init gravity
    
    estimator = ESEKF(init_nominal_state, imu_parameters);
    
    % Filter data
    test_duration = [0, 61];
    start_time = imu_data(1,1);
    mask_imu = imu_data(:,1) >= start_time + test_duration(1) & imu_data(:,1) <= start_time + test_duration(2);
    mask_gt = gt_data(:,1) >= start_time + test_duration(1) & gt_data(:,1) <= start_time + test_duration(2);
    imu_data = imu_data(mask_imu, :);
    gt_data = gt_data(mask_gt, :);
    
    traj_est = zeros(size(gt_data,1), 8);
    traj_est(1,:) = gt_data(1,1:8);
    update_ratio = 10;
    sigma_measurement_p = 0.02;
    sigma_measurement_q = 0.015;
    sigma_measurement = diag([sigma_measurement_p^2 * ones(1,3), sigma_measurement_q^2 * ones(1,3)]);
    
    for i = 2:size(imu_data,1)
        timestamp = imu_data(i,1);
        estimator = estimator.predict(imu_data(i, :));
        
        if mod(i, update_ratio) == 0
            gt_pose = gt_data(i, 2:8)';
            gt_pose(1:3) = gt_pose(1:3) + randn(3,1) * sigma_measurement_p;
            noise_axis = randn(3,1);
            noise_angle = norm(noise_axis) * sigma_measurement_q;
            noise_axis = noise_axis / norm(noise_axis);
            qn = axang2quat([noise_axis' noise_angle]);
            gt_pose(4:7) = quatmultiply(qn, gt_pose(4:7)');
            estimator = estimator.update(gt_pose, sigma_measurement);
        end
        
        traj_est(i,:) = [timestamp, estimator.nominal_state(1:7)'];
    end
    
    % Save estimated trajectory
    save('./data/traj_esekf_out.txt', 'traj_est', '-ascii');

    %% Plotting routines
time = traj_est(:,1);
euler_anglesGT = quat2eul([gt_data(:,8) gt_data(:,5:7)]); %GT
euler_anglesEst = quat2eul([traj_est(:,8) traj_est(:,5:7)]); % Est
euler_anglesGT=flip(euler_anglesGT,2); %reverse order since it comes out as YPR
euler_anglesEst=flip(euler_anglesEst,2); %reverse order since it comes out as YPR
colororder({'r', 'b'}) %mapped colors
figure;

subplot(3,1,1);
plot(time,rad2deg(euler_anglesGT(:,1)),'k','LineWidth',1);
hold on; plot(time,rad2deg(euler_anglesEst(:,1)),'b-.','LineWidth',2);
ylabel('Roll (Degrees)'); xlabel('time (s)');

subplot(3,1,2);
plot(time,rad2deg(euler_anglesGT(:,2)),'k','LineWidth',1);
hold on; plot(time,rad2deg(euler_anglesEst(:,2)),'b-.','LineWidth',2);
ylabel('Pitch (Degrees)'); xlabel('time (s)');

subplot(3,1,3);
plot(time,rad2deg(euler_anglesGT(:,3)),'k','LineWidth',1);
hold on; plot(time,rad2deg(euler_anglesEst(:,3)),'b-.','LineWidth',2);
ylabel('Yaw (Degrees)'); xlabel('time (s)');
legend('Ground Truth','Estimate');

%errors
figure;
subplot(3,1,1);
plot(time,rad2deg(euler_anglesGT(:,1))-rad2deg(euler_anglesEst(:,1)),'b','LineWidth',2);
ylabel('Roll Error (Degrees)'); xlabel('time (s)');

subplot(3,1,2);
plot(time,rad2deg(euler_anglesGT(:,2))-rad2deg(euler_anglesEst(:,2)),'b','LineWidth',2);
ylabel('Pitch Error (Degrees)'); xlabel('time (s)');

subplot(3,1,3);
plot(time,rad2deg(euler_anglesGT(:,3))-rad2deg(euler_anglesEst(:,3)),'b','LineWidth',2);
ylabel('Yaw Error (Degrees)'); xlabel('time (s)');
title("Errors'")
%% Load IMU params
function imu_params = load_imu_parameters()
    % Load IMU parameters from a YAML file
    yamlFile = './data/testparams.yaml';
    
    if exist(yamlFile, 'file')
        params = readyaml(yamlFile); % Requires MATLAB YAML parser
        imu_params = ImuParameters();
        imu_params.frequency = params.IMU.frequency;
        imu_params.sigma_a_n = params.IMU.acc_noise_sigma;
        imu_params.sigma_w_n = params.IMU.gyro_noise_sigma;
        imu_params.sigma_a_b = params.IMU.acc_bias_sigma;
        imu_params.sigma_w_b = params.IMU.gyro_bias_sigma;
    else
        error('YAML file not found: %s', yamlFile);
    end
end
