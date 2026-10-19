%  HOHMANN TRANSFER SIMULATOR
%  Two-impulse transfer between coplanar circular orbits (LEO to GEO default)

clc; clear; close all;

%% 1. User Input
fprintf("\n======================================\n");
fprintf("  Hohmann Transfer Simulator\n");
fprintf("  Press Enter To Use Default Value.\n");
fprintf("======================================\n\n");

fprintf("-- Central Body ------------------------\n");
fprintf("   Presets: Earth mu=3.986004418e14  R=6378137\n");
fprintf("            Moon  mu=4.9048695e12    R=1737400\n");
fprintf("            Mars  mu=4.282837e13     R=3389500\n\n");
mu   = input_val("  Grav. Parameter  [m^3/s^2] (default 3.986004418e14) : ", 3.986004418e14);
Re   = input_val("  Body Radius      [m]       (default 6378137       ) : ", 6378137);

fprintf("\n-- Orbits ------------------------------\n");
fprintf("   Target Altitude Default 35785.863 km Is GEO For Earth\n\n");
alt1 = input_val("  Initial Altitude [km]      (default 400           ) : ", 400) * 1e3;
alt2 = input_val("  Target Altitude  [km]      (default 35785.863     ) : ", 35785.863) * 1e3;
fprintf("\n");

r1 = Re + alt1;           % Initial circular orbit radius [m]
r2 = Re + alt2;           % Target circular orbit radius [m]

if mu <= 0 || Re <= 0 || alt1 <= 0 || r2 <= r1
    fprintf("Invalid Input: Need Positive Mu, Radius, And Initial Altitude, With Target Above Initial Orbit.\n");
    return;
end

t_max = 2 * pi * sqrt((50 * r1)^3 / mu);   % Upper bound on propagation time [s]

%% 2. Analytic Hohmann Solution
a_t   = (r1 + r2) / 2;                    % Transfer orbit semi-major axis
v1c   = sqrt(mu / r1);                    % Circular speed at r1
v2c   = sqrt(mu / r2);                    % Circular speed at r2
vp_an = sqrt(mu * (2/r1 - 1/a_t));        % Transfer speed at periapsis
va_an = sqrt(mu * (2/r2 - 1/a_t));        % Transfer speed at apoapsis
dv1_an = vp_an - v1c;
dv2_an = v2c - va_an;
dv_an  = dv1_an + dv2_an;
tof_an = pi * sqrt(a_t^3 / mu);

%% 3. Numerical Solution (Shooting Method)
fprintf("\nSolving Departure Burn By Shooting...\n");
v_lo   = 1.001 * v1c;                     % Bracket: barely above circular
v_hi   = 0.99 * sqrt(2) * v1c;            % Bracket: just below escape
if apoapsis_radius(v_hi, r1, mu, t_max) < r2
    fprintf("Target Orbit Too Far For This Solver (Needs r2 Below About 49 r1).\n");
    return;
end
miss   = @(vp) apoapsis_radius(vp, r1, mu, t_max) - r2;
vp_num = fzero(miss, [v_lo, v_hi], optimset("TolX", 1e-12));

% Full trajectory for plotting (time grid only; the event ends the run)
t_grid = linspace(0, 1.1 * tof_an, 2000);
[t, Y] = propagate_to_apoapsis(vp_num, r1, mu, t_grid);

r_ap    = sqrt(Y(end,1)^2 + Y(end,2)^2);
v_ap    = sqrt(Y(end,3)^2 + Y(end,4)^2);
dv1_num = vp_num - v1c;
dv2_num = v2c - v_ap;
dv_num  = dv1_num + dv2_num;
tof_num = t(end);

%% 4. Conservation Checks
r   = sqrt(Y(:,1).^2 + Y(:,2).^2);
spd = sqrt(Y(:,3).^2 + Y(:,4).^2);
E   = 0.5 * spd.^2 - mu ./ r;                 % Specific energy
h   = Y(:,1) .* Y(:,4) - Y(:,2) .* Y(:,3);    % Specific angular momentum
E_err = abs((E - E(1)) / E(1));
h_err = abs((h - h(1)) / h(1));

%% 5. Results
pct = @(a, b) abs(a - b) / abs(b) * 100;

fprintf("\n======================================\n");
fprintf("  Hohmann Transfer Results\n");
fprintf("======================================\n");
fprintf("  Initial Orbit Radius   : %.1f km\n", r1 / 1e3);
fprintf("  Target Orbit Radius    : %.1f km\n", r2 / 1e3);
fprintf("  Apoapsis Reached       : %.3f km\n", r_ap / 1e3);
fprintf("--------------------------------------\n");
fprintf("  %-22s %12s %12s %10s\n", "Quantity", "Numerical", "Analytic", "Error [%]");
fprintf("  %-22s %12.4f %12.4f %10.2e\n", "Departure Dv [m/s]", dv1_num, dv1_an, pct(dv1_num, dv1_an));
fprintf("  %-22s %12.4f %12.4f %10.2e\n", "Arrival Dv [m/s]",   dv2_num, dv2_an, pct(dv2_num, dv2_an));
fprintf("  %-22s %12.4f %12.4f %10.2e\n", "Total Dv [m/s]",     dv_num,  dv_an,  pct(dv_num, dv_an));
fprintf("  %-22s %12.2f %12.2f %10.2e\n", "Transfer Time [s]",  tof_num, tof_an, pct(tof_num, tof_an));
fprintf("--------------------------------------\n");
fprintf("  Transfer Time          : %.3f hr\n", tof_num / 3600);
fprintf("  Max Energy Drift       : %.2e (Relative)\n", max(E_err));
fprintf("  Max Ang. Mom. Drift    : %.2e (Relative)\n", max(h_err));
fprintf("======================================\n\n");

%% 6. Plots
c = lines(4);
theta = linspace(0, 2*pi, 400);

figure("Name", "Hohmann Transfer Simulator", "Position", [100 80 1200 700]);
tl = tiledlayout(2, 2, "TileSpacing", "compact");

% Plot 1: Orbit Geometry
nexttile([2 1]);
fill(Re/1e3 * cos(theta), Re/1e3 * sin(theta), c(1,:), ...
     "FaceAlpha", 0.3, "EdgeColor", "none"); hold on;
plot(r1/1e3 * cos(theta), r1/1e3 * sin(theta), "--", "Color", c(2,:), "LineWidth", 1.2);
plot(r2/1e3 * cos(theta), r2/1e3 * sin(theta), "--", "Color", c(3,:), "LineWidth", 1.2);
plot(Y(:,1)/1e3, Y(:,2)/1e3, "-", "Color", c(4,:), "LineWidth", 2.2);
plot(Y(1,1)/1e3, Y(1,2)/1e3, "ko", "MarkerFaceColor", "k", "MarkerSize", 7);
plot(Y(end,1)/1e3, Y(end,2)/1e3, "ks", "MarkerFaceColor", "k", "MarkerSize", 7);
axis equal; grid on;
xlabel("X (km)"); ylabel("Y (km)");
title("Transfer Geometry");
legend(["Central Body", "Initial Orbit", "Target Orbit", "Transfer Orbit", ...
        "Departure Burn", "Arrival Burn"], "Location", "northeast");

% Plot 2: Speed vs Time
nexttile;
plot(t/3600, spd/1e3, "-", "Color", c(4,:), "LineWidth", 2); hold on;
yline(v1c/1e3, "--", "Initial Circular Speed");
yline(v2c/1e3, "--", "Target Circular Speed");
grid on;
xlabel("Time (hr)"); ylabel("Speed (km/s)");
title("Speed Along Transfer");

% Plot 3: Conservation Errors
nexttile;
semilogy(t/3600, max(E_err, eps), "-", "Color", c(2,:), "LineWidth", 1.8); hold on;
semilogy(t/3600, max(h_err, eps), "-", "Color", c(3,:), "LineWidth", 1.8);
grid on;
xlabel("Time (hr)"); ylabel("Relative Drift");
title("Conservation Check");
legend(["Specific Energy", "Angular Momentum"], "Location", "southeast");

title(tl, "Hohmann Transfer Simulator");

fprintf("Done.\n");

%  LOCAL FUNCTIONS
function dydt = two_body(y, mu)
    % Two-body equations of motion. State: [x; y; vx; vy]
    r    = sqrt(y(1)^2 + y(2)^2);
    dydt = [y(3); y(4); -mu * y(1) / r^3; -mu * y(2) / r^3];
end

function [value, isterminal, direction] = apoapsis_event(~, y)
    % Radial velocity (times r) crossing zero from positive to negative
    value      = y(1) * y(3) + y(2) * y(4);
    isterminal = 1;
    direction  = -1;
end

function [t, Y] = propagate_to_apoapsis(vp, r1, mu, tspan)
    % Start at periapsis (r1, 0) moving in +y, stop at apoapsis
    opts = odeset("RelTol", 1e-12, "AbsTol", 1e-9, "Events", @apoapsis_event);
    y0   = [r1; 0; 0; vp];
    [t, Y] = ode45(@(~, y) two_body(y, mu), tspan, y0, opts);
end

function val = input_val(prompt, default)
    % Prompt user for input; returns default if Enter is pressed
    raw = input(prompt, "s");
    if isempty(strtrim(raw))
        val = default;
    else
        val = str2double(raw);
        if isnan(val)
            fprintf("  Invalid Input, Using Default: %g\n", default);
            val = default;
        end
    end
end

function ra = apoapsis_radius(vp, r1, mu, t_max)
    % Apoapsis radius reached for a given periapsis speed
    [~, Y] = propagate_to_apoapsis(vp, r1, mu, [0, t_max]);
    ra = sqrt(Y(end,1)^2 + Y(end,2)^2);
end