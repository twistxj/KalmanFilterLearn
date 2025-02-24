classdef ImuParameters
    properties
        frequency = 200;
        sigma_a_n = 0.0;  % Accelerometer noise (m/s/sqrt(s))
        sigma_w_n = 0.0;  % Gyroscope noise (rad/s/sqrt(s))
        sigma_a_b = 0.0;  % Accelerometer bias noise (m/s/sqrt(s^5))
        sigma_w_b = 0.0;  % Gyroscope bias noise (rad/s/sqrt(s^3))
    end
end

