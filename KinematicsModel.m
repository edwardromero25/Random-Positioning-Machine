classdef KinematicsModel
    properties
        f_hz
        omega_rpm_min
        omega_rpm_max
        t_transition_s
        t_hold_s
        alpha_0_deg
        beta_0_deg
        x_cm
        y_cm
        z_cm
        duration_hours
    end

    methods
        function obj = KinematicsModel(f_hz, omega_rpm_min, omega_rpm_max, t_transition_s, t_hold_s, alpha_0_deg, beta_0_deg, x_cm, y_cm, z_cm, duration_hours)
            obj.f_hz = f_hz;
            obj.omega_rpm_min = omega_rpm_min;
            obj.omega_rpm_max = omega_rpm_max;
            obj.t_transition_s = t_transition_s;
            obj.t_hold_s = t_hold_s;
            obj.alpha_0_deg = alpha_0_deg;
            obj.beta_0_deg = beta_0_deg;
            obj.x_cm = x_cm;
            obj.y_cm = y_cm;
            obj.z_cm = z_cm;
            obj.duration_hours = duration_hours;
        end
        
        function random_positioning_machine_profile = generate_random_positioning_machine_profile(obj)
            dt = 1 / obj.f_hz;  
            time_array = 0:dt:(obj.duration_hours * 3600);  
            random_positioning_machine_profile = zeros(size(time_array));  
            
            current_rpm = 0; 
            random_positioning_machine_profile(1) = current_rpm;  
            
            i = 2;  
            while i <= length(time_array) 
                target_rpm = (obj.omega_rpm_max - obj.omega_rpm_min) * rand() + obj.omega_rpm_min;  
                target_rpm = target_rpm * (2 * randi([0 1]) - 1); 
                 
                steps_ramp = floor(obj.t_transition_s / dt);
                
                steps_hold = floor(obj.t_hold_s / dt);
        
                idx_end_ramp = min(i + steps_ramp - 1, length(random_positioning_machine_profile));
                ramp = linspace(current_rpm, target_rpm, idx_end_ramp - i + 1);
                random_positioning_machine_profile(i:idx_end_ramp) = ramp;
       
                idx_end_hold = min(idx_end_ramp + steps_hold, length(random_positioning_machine_profile));
                random_positioning_machine_profile((idx_end_ramp+1):idx_end_hold) = target_rpm;
                
                current_rpm = target_rpm;
                i = idx_end_hold + 1;
            end
        end

        function rad = rpm_to_rad_sec(~, rpm)
            rad = rpm * pi / 30;
        end

        function rad = deg_to_rad(~, degrees)
            rad = deg2rad(degrees);
        end

        function [time_array, omega_alpha, omega_beta, g_local_2, a_local_2, a_tot_local_2] = calculate_acceleration(obj)
            dt = 1 / obj.f_hz;
            time_array = 0:dt:(obj.duration_hours * 3600);

            omega_alpha = obj.rpm_to_rad_sec(obj.generate_random_positioning_machine_profile());
            omega_beta = obj.rpm_to_rad_sec(obj.generate_random_positioning_machine_profile());

            alpha_0 = obj.deg_to_rad(obj.alpha_0_deg);
            beta_0 = obj.deg_to_rad(obj.beta_0_deg);

            alpha_t = omega_alpha .* time_array + alpha_0;
            beta_t = omega_beta .* time_array + beta_0;

            x = obj.x_cm / 100;
            y = obj.y_cm / 100;
            z = obj.z_cm / 100;

            omega_tot = [
                omega_alpha .* ones(size(time_array));
                omega_beta .* cos(alpha_t);
                omega_beta .* sin(alpha_t)
            ];

            omega_tot_dot = [
                zeros(size(time_array));
                -omega_alpha .* omega_beta .* sin(alpha_t);
                omega_alpha .* omega_beta .* cos(alpha_t)
            ];

            r = [
                x * cos(beta_t) + z * sin(beta_t);
                y * cos(alpha_t) + x * sin(alpha_t).*sin(beta_t) - z * sin(alpha_t).*cos(beta_t);
                y * sin(alpha_t) - x * cos(alpha_t).*sin(beta_t) + z * cos(alpha_t).*cos(beta_t)
            ];

            omega_cross_r = cross(omega_tot', r')';
            omega_cross_omega_cross_r = cross(omega_tot', omega_cross_r')';
            omega_dot_cross_r = cross(omega_tot_dot', r')';
            a = -(omega_dot_cross_r + omega_cross_omega_cross_r);

            g = [0; 0; 1];

            R_y_T = zeros(3, 3, length(beta_t));
            R_x_T = zeros(3, 3, length(alpha_t));
            for i = 1:length(beta_t)
                R_y_T(:,:,i) = [
                    cos(beta_t(i)), 0, -sin(beta_t(i));
                    0, 1, 0;
                    sin(beta_t(i)), 0, cos(beta_t(i))
                ];
                R_x_T(:,:,i) = [
                    1, 0, 0;
                    0, cos(alpha_t(i)), sin(alpha_t(i));
                    0, -sin(alpha_t(i)), cos(alpha_t(i))
                ];
            end

            a_local_2 = zeros(3, length(time_array));
            g_local_2 = zeros(3, length(time_array));
            for i = 1:length(time_array)
                a_local_2(:,i) = R_y_T(:,:,i) * (R_x_T(:,:,i) * a(:,i)) / 9.8;
                g_local_2(:,i) = R_y_T(:,:,i) * (R_x_T(:,:,i) * g);
            end
            a_tot_local_2 = g_local_2 + a_local_2;
        end
    end
end
