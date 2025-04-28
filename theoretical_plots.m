function theoretical_plots()
    set(0, 'DefaultAxesFontName', 'Raleway');
    set(0, 'DefaultTextFontName', 'Raleway');
    set(0, 'DefaultAxesFontWeight', 'normal');
    set(0, 'DefaultTextFontWeight', 'normal');

    sampling_rate = 750; % hz
    angular_velocity_min = 10; % rpm
    angular_velocity_max = 20; % rpm
    transition_time = 1; % s
    hold_time = 4; % s
    outer_initial_angular_position = 0; % deg
    inner_initial_angular_position = 0; % deg
    distance_from_center = 0; % cm
    simulation_duration = 0.1; % h

    model = KinematicsModel(sampling_rate, angular_velocity_min, angular_velocity_max, transition_time, hold_time, ...
        outer_initial_angular_position, inner_initial_angular_position, ...
        distance_from_center, distance_from_center, distance_from_center, simulation_duration);
    [time_array, omega_alpha, omega_beta, g_local_2, a_local_2, a_tot_local_2] = model.calculate_acceleration();
    time_in_hours = time_array / 3600;
    
    omega_alpha_rpm = omega_alpha * 30 / pi;
    omega_beta_rpm = omega_beta * 30 / pi;

    figure;
    plot(time_in_hours, omega_alpha_rpm, ...
         time_in_hours, omega_beta_rpm);
    title("Angular Velocity vs. Time", 'FontWeight', 'normal');
    xlabel("Time (h)");
    ylabel("Angular Velocity (rpm)");
    legend("Outer", "Inner");

    g_x_avg = cumsum(g_local_2(1,:)) ./ (1:length(g_local_2(1,:)));
    g_y_avg = cumsum(g_local_2(2,:)) ./ (1:length(g_local_2(2,:)));
    g_z_avg = cumsum(g_local_2(3,:)) ./ (1:length(g_local_2(3,:)));
    g_magnitude_avg = sqrt(g_x_avg.^2 + g_y_avg.^2 + g_z_avg.^2);
    
    magnitude = mean(g_magnitude_avg);
    magnitude_label = sprintf('Magnitude: %.3g', magnitude);

    figure;
    hold on;
    plot(time_in_hours, g_magnitude_avg, 'b', 'DisplayName', magnitude_label);
    title("Time-Averaged Gravitational Acceleration", 'FontWeight', 'normal');
    xlabel("Time (h)");
    ylabel("Acceleration (g)");
    legend('show');

    figure;
    plot(time_in_hours, g_x_avg, 'm', ...
         time_in_hours, g_y_avg, 'g', ...
         time_in_hours, g_z_avg, 'k');
    title("Time-Averaged Gravitational Acceleration", 'FontWeight', 'normal');
    xlabel("Time (h)");
    ylabel("Acceleration (g)");
    legend("X", "Y", "Z");

    a_x_avg = cumsum(a_local_2(1,:)) ./ (1:length(a_local_2(1,:)));
    a_y_avg = cumsum(a_local_2(2,:)) ./ (1:length(a_local_2(2,:)));
    a_z_avg = cumsum(a_local_2(3,:)) ./ (1:length(a_local_2(3,:)));
    a_magnitude_avg = sqrt(a_x_avg.^2 + a_y_avg.^2 + a_z_avg.^2);

    magnitude = mean(a_magnitude_avg);
    magnitude_label = sprintf('Magnitude: %.3g', magnitude);

    figure;
    hold on;
    plot(time_in_hours, a_magnitude_avg, 'b', 'DisplayName', magnitude_label);
    title("Time-Averaged Non-Gravitational Acceleration", 'FontWeight', 'normal');
    xlabel("Time (h)");
    ylabel("Acceleration (g)");
    legend('show');

    figure;
    plot(time_in_hours, a_x_avg, 'm', ...
         time_in_hours, a_y_avg, 'g', ...
         time_in_hours, a_z_avg, 'k');
    title("Time-Averaged Non-Gravitational Acceleration", 'FontWeight', 'normal');
    xlabel("Time (h)");
    ylabel("Acceleration (g)");
    legend("X", "Y", "Z");

    distribution = FibonacciLattice("theoretical", a_tot_local_2(1,:), a_tot_local_2(2,:), a_tot_local_2(3,:)).getDistribution();
    distribution_label = sprintf('Distribution: %d', distribution);

    figure;
    hold on;
    [u, v] = meshgrid(linspace(0, 2*pi, 25), linspace(0, pi, 25));
    x = cos(u) .* sin(v);
    y = sin(u) .* sin(v);
    z = cos(v);
    mesh(x, y, z, 'EdgeColor', [0.5 0.5 0.5], 'FaceAlpha', 0, ...
         'EdgeAlpha', 0.5, 'HandleVisibility', 'off');
    plot3(a_tot_local_2(1,:), a_tot_local_2(2,:), a_tot_local_2(3,:), 'b', 'DisplayName', distribution_label);
    title("Orientation Distribution", 'FontWeight', 'normal');
    xlabel('X (g)');
    ylabel('Y (g)');
    zlabel('Z (g)');
    legend('show');
    xticks([-1 -0.5 0 0.5 1]);
    yticks([-1 -0.5 0 0.5 1]);
    zticks([-1 -0.5 0 0.5 1]);
    grid off;
    view(3);
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    axis equal;
    hold off;

    figure;
    hold on;
    title("Orientation Distribution", 'FontWeight', 'normal');
    xlabel('X (g)');
    ylabel('Y (g)');
    zlabel('Z (g)');
    legend('show');
    xticks([-1 -0.5 0 0.5 1]);
    yticks([-1 -0.5 0 0.5 1]);
    zticks([-1 -0.5 0 0.5 1]);
    grid off;
    view(3);
    xlim([-1 1]);
    ylim([-1 1]);
    zlim([-1 1]);
    axis equal;
    [u, v] = meshgrid(linspace(0, 2*pi, 25), linspace(0, pi, 25));
    x = cos(u) .* sin(v);
    y = sin(u) .* sin(v);
    z = cos(v);
    mesh(x, y, z, 'EdgeColor', [0.5 0.5 0.5], ...
         'FaceAlpha', 0, 'EdgeAlpha', 0.5, 'HandleVisibility', 'off');
    animated_line = plot3(NaN, NaN, NaN, 'r', 'LineWidth', 1, ...
                          'DisplayName', distribution_label);
    for k = 1:50:length(a_tot_local_2)
        set(animated_line, 'XData', a_tot_local_2(1,1:k), ...
                           'YData', a_tot_local_2(2,1:k), ...
                           'ZData', a_tot_local_2(3,1:k));
        drawnow;
    end
    hold off;
end
