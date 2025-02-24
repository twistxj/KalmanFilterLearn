classdef ESEKF
    properties
        nominal_state;
        error_covar;
        noise_covar;
        last_predict_time = 0;
    end

    methods
        function obj = ESEKF(init_nominal_state, imu_params)
            obj.nominal_state = init_nominal_state;
            if obj.nominal_state(4) < 0
                obj.nominal_state(4:7) = obj.nominal_state(4:7) * -1;
            end

            % Initialize noise covariance
            noise_covar = zeros(12, 12);
            noise_covar(1:3, 1:3) = (imu_params.sigma_a_n^2) * eye(3);
            noise_covar(4:6, 4:6) = (imu_params.sigma_w_n^2) * eye(3);
            noise_covar(7:9, 7:9) = (imu_params.sigma_a_b^2) * eye(3);
            noise_covar(10:12, 10:12) = (imu_params.sigma_w_b^2) * eye(3);

            G = zeros(18, 12);
            G(4:6, 4:6) = -eye(3);
            G(7:9, 1:3) = -eye(3);
            G(10:12, 7:9) = eye(3);
            G(13:15, 10:12) = eye(3);

            obj.noise_covar = G * noise_covar * G';
            obj.error_covar = 0.01 * obj.noise_covar;

        end
        function obj = predict(obj, imu_measurement)
            if obj.last_predict_time == imu_measurement(1)
                return;
            end
            obj = obj.predictErrorCovar(imu_measurement);
            obj = obj.predictNominalState(imu_measurement);
            obj.last_predict_time = imu_measurement(1);
        end
        function obj = update(obj, gt_measurement, measurement_covar)
            % Update step
            H = zeros(6, 18);
            H(1:3, 1:3) = eye(3);
            H(4:6, 4:6) = eye(3);
            PHt = obj.error_covar * H'; % 18x6
            K = PHt / (H * PHt + measurement_covar); % 18x6

            % Update error covariance matrix
            obj.error_covar = (eye(18) - K * H) * obj.error_covar;
            obj.error_covar = 0.5 * (obj.error_covar + obj.error_covar'); % Force symmetry

            % Compute measurements based on nominal state and ground truth
            if gt_measurement(4) < 0
                gt_measurement(4:7) = -gt_measurement(4:7);
            end
            gt_p = gt_measurement(1:3);
            gt_q = gt_measurement(4:7);
            q = obj.nominal_state(4:7);

            delta = zeros(6, 1);
            delta(1:3) = gt_p - obj.nominal_state(1:3);
            delta_q = quatmultiply(quatconj(q'), gt_q');
            if delta_q(1) < 0
                delta_q = -delta_q;
            end
            angle = asin(norm(delta_q(2:4)));
            if abs(angle) < 1e-6
                axis = zeros(3, 1);
            else
                axis = delta_q(2:4) / norm(delta_q(2:4));
            end
            delta(4:6) = angle * axis;

            % Compute state errors
            errors = K * delta;

            % Inject errors into nominal state
            obj.nominal_state(1:3) = obj.nominal_state(1:3) + errors(1:3); % Update position
            dq = axang2quat([errors(4:6)' norm(errors(4:6))]);
            obj.nominal_state(4:7) = quatmultiply(q', dq); % Update rotation
            obj.nominal_state(4:7) = obj.nominal_state(4:7) / norm(obj.nominal_state(4:7));
            if obj.nominal_state(4) < 0
                obj.nominal_state(4:7) = -obj.nominal_state(4:7);
            end
            obj.nominal_state(8:end) = obj.nominal_state(8:end) + errors(7:end); % Update the rest

            % Reset errors to zero and modify the error covariance matrix
            G = eye(18);
            G(4:6, 4:6) = eye(3) - skew_symmetric(0.5 * errors(4:6));
            obj.error_covar = G * obj.error_covar * G';
        end


        function obj = predictNominalState(obj, imu_measurement)
            p = obj.nominal_state(1:3);
            q = obj.nominal_state(4:7);
            v = obj.nominal_state(8:10);
            a_b = obj.nominal_state(11:13);
            w_b = obj.nominal_state(14:16);
            g = obj.nominal_state(17:19);

            w_m = imu_measurement(2:4)';
            a_m = imu_measurement(5:7)';
            dt = imu_measurement(1) - obj.last_predict_time;

            % Prediction equations
            w_m = w_m - w_b;
            a_m = a_m - a_b;

            % Integrate quaternion
            angle = norm(w_m);
            axis = w_m / angle;
            R_w = axang2rotm([ axis' 0.5 * dt * angle]);
            q_w = rotm2quat(R_w);
            q_half_next = quatmultiply(q', q_w);

            R_w = axang2rotm([axis' dt * angle]);
            q_w = rotm2quat(R_w);
            q_next = quatmultiply(q', q_w);
            if q_next(1) < 0
                q_next = -q_next; % Force positive real part
            end

            % RK4 integration for velocity and position
            R = quat2rotm(q');
            R_half_next = quat2rotm(q_half_next);
            R_next = quat2rotm(q_next);
            v_k1 = R(1:3, 1:3) * a_m + g;
            v_k2 = R_half_next(1:3, 1:3) * a_m + g;
            v_k3 = v_k2;
            v_k4 = R_next(1:3, 1:3) * a_m + g;
            v_next = v + dt * (v_k1 + 2 * v_k2 + 2 * v_k3 + v_k4) / 6;

            % Integrate position
            p_k1 = v;
            p_k2 = v + 0.5 * dt * v_k1;
            p_k3 = v + 0.5 * dt * v_k2;
            p_k4 = v + dt * v_k3;
            p_next = p + dt * (p_k1 + 2 * p_k2 + 2 * p_k3 + p_k4) / 6;

            obj.nominal_state(1:3) = p_next;
            obj.nominal_state(4:7) = q_next;
            obj.nominal_state(8:10) = v_next;
        end
        function obj = predictErrorCovar(obj, imu_measurement)
            w_m = imu_measurement(2:4)';
            a_m = imu_measurement(5:7)';
            a_b = obj.nominal_state(11:13);
            w_b = obj.nominal_state(14:16);
            q = obj.nominal_state(4:7);

            % Compute rotation matrix from quaternion
            R = quat2rotm(q');

            % Compute F matrix
            F = zeros(18,18);
            F(1:3, 7:9) = eye(3);
            F(4:6, 4:6) = -skew_symmetric(w_m - w_b);
            F(4:6, 13:15) = -eye(3);
            F(7:9, 4:6) = -R * skew_symmetric(a_m - a_b);
            F(7:9, 10:12) = -R;

            % Time step
            dt = imu_measurement(1) - obj.last_predict_time;

            % Compute transition matrix Phi using Taylor expansion
            Fdt = F * dt;
            Fdt2 = Fdt * Fdt;
            Fdt3 = Fdt2 * Fdt;
            Phi = eye(18) + Fdt + 0.5 * Fdt2 + (1/6) * Fdt3;

            % Compute process noise
            Qc_dt = 0.5 * dt * obj.noise_covar;
            obj.error_covar = Phi * (obj.error_covar + Qc_dt) * Phi' + Qc_dt;
        end
    end
end