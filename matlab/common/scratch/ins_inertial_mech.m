function nav_state = ins_inertial_mech(nav_state_prev, imu, sim_dt)


% inertial mechanisation of the ins equations (no bias model)

% previous time position, velocity and attitude
pos_prev = nav_state_prev(1:3);
vel_prev = nav_state_prev(4:6);
euler_prev = nav_state_prev(7:9);

% gravity vector
gn = [0; 0; 9.81];

% find euler rates
E = bodyrates_to_E(euler_prev);
euler_dot = E*imu(4:6);

% Integrate Euler Angles
est_euler = int_o1(euler_dot, euler_prev, sim_dt);

% heading modulation
est_euler(3) = head_mod(est_euler(3));

% compose body to LV-LH transformation matrix
cbn = euler_to_R(est_euler); % body to navigation

% Transform IMU acceleration
fn = cbn*imu(1:3) + gn;

% Integrate acceleration
est_velocity = int_o1(fn, vel_prev, sim_dt);

% Integrate velocity
est_position = int_o1(est_velocity, pos_prev, sim_dt);

nav_state = [est_position; est_velocity; est_euler];