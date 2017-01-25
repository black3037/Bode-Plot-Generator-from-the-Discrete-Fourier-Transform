omega = logspace(1,3,11);
%omega = [omega 2000 5000];   % add a couple of higher freq points
%omega = [1000 2000 5000];    % only do the high freq points

Tsim = 0.0001;         % Needed by simulink
Tdiscrete = 0.0001;    %        "
 
wn = 100; zeta = 0.1;
G = tf(wn^2,[1 2*zeta*wn wn^2]);        % Try both of these models
G = tf(wn^2,[1 2*zeta*wn wn^2 0]);
[m,p,w] = bode(G,{omega(1),omega(end)});
m = squeeze(m);
p = squeeze(p);
 
N = length(omega);
identified_mag = zeros(1,N);
identified_phase = zeros(1,N);
phase_shift = zeros(1,N);

for i = 1:N
    sim_termination = 0;
    sin_freq = omega(i);     % Needed by simulink
    wc = sin_freq/10;
    phase_shift(i) = 1.5*Tdiscrete*wc*360/2*pi;
    % run simulation for this frequency
    % sim_termination, sim_mag_ratio, and sim_phase are output by the sim
    sim('simulation');
    
    identified_mag(i) = sim_mag_ratio;
    identified_phase(i) = sim_phase;
    
    % print some things to command window
    msg = sprintf('termination = %f, freq = %f rad/s', sim_termination, sin_freq);
    disp(msg)
end
 
phase2 = unwrap(identified_phase);
figure(1);
subplot(2,1,1)
loglog(omega,identified_mag,'*',w,m);
ylabel('mag ratio')
subplot(2,1,2)
semilogx(omega,identified_phase*180/pi,'*',w,p,omega,phase2*180/pi + phase_shift,'O');
xlabel('freq (rad/s)'); ylabel('phase (deg)')

