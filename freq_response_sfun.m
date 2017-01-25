function [sys,x0,str,ts]=freq_response_sfun(t,x,u,flag,Tsim,Tdiscrete,sin_freq)
 
% Because we are using persistent variables, we have to be careful with
% multiple uses of this s-function. All instances must calculate these
% persistent variables same because they are stored with this function. 
persistent vars
if (flag == 0)
    
	%{
	*--*--**--*--**--*--**--*--**--*--**--*--*
	Set up delay variables
	*--*--**--*--**--*--**--*--**--*--**--*--*
	%}
		
    % calculate the number of sample delays and find actual delay
    vars.num_delays = round(Tdiscrete/Tsim);
    
    vars.Td=vars.num_delays*Tsim;
    
    vars.delay_count = 0; % counter for delay to create discrete period
    
    vars.w = sin_freq; % frequency of sin wave driving system
	
	%{
	*--*--**--*--**--*--**--*--**--*--**--*--*
	Setup wave state variables
	*--*--**--*--**--*--**--*--**--*--**--*--*
	%}
    
    vars.sin_out = 0;  % Sine wave output
    
    vars.wc = vars.w/10; % Cutoff frequency
    
    vars.sin_state = 0; % Low-pass filter cos correlation
    
    vars.cos_state = 0; % Low-pass filter sin correlation
    
	%{
	*--*--**--*--**--*--**--*--**--*--**--*--*
	Setup moving average (DFT) variables
	*--*--**--*--**--*--**--*--**--*--**--*--*
	%}
    
    num_of_periods = 3;
    
    T = 2*pi*(num_of_periods/vars.w);
    
    vars.chunk_size = round(T/vars.Td); % Size of chunk for running average
    
    vars.N = 0; % Number of points for running average
    
    vars.cos_sum = 0;
    
    vars.sin_sum = 0;
    
    vars.a = 0;
    
    vars.b = 0;
      
	%{
	*--*--**--*--**--*--**--*--**--*--**--*--*
	Setup state-space
	*--*--**--*--**--*--**--*--**--*--**--*--*
	%}
    
    vars.A = [0 1 ; -vars.wc^2 -1.414*vars.wc];  % A matrix
    
    vars.B = [0 ; vars.wc^2];% A matrix
    
    vars.xa = [0 ; 0]; % Euler integration
    
    vars.xb = [0 ; 0]; % Euler integration
    
    %{
	*--*--**--*--**--*--**--*--**--*--**--*--*
	Setup termination condition variables
	*--*--**--*--**--*--**--*--**--*--**--*--*
	%}
    
    vars.mag = [0;0;0;0];
    
    vars.num_DFTs = 0;

    vars.termination = 0;

end
 
switch flag
 
  %%%%%%%%%%%%%%%%%%%%%%%%
  % Initialization 
  %%%%%%%%%%%%%%%%%%%%%%%%
  case 0         
    [sys,x0,str,ts] = mdlInitializeSizes(Tsim);
 
  %%%%%%%%%%%%%%%%%%%%%%%%
  % Update
  %%%%%%%%%%%%%%%%%%%%%%%%
  case 2
    sys = [];                   % states are stored in vars
    vars = mdlUpdate(t,u,vars);
 
  %%%%%%%%%%%%%%%%%%%%%%%%
  % Output 
  %%%%%%%%%%%%%%%%%%%%%%%%
  case 3
    sys = mdlOutputs(t,u,vars); 
 
  %%%%%%%%%%%%%%%%%%%%%%%%
  % Terminate 
  %%%%%%%%%%%%%%%%%%%%%%%%
  case 9
    sys = []; % do nothing
 
  otherwise
    error(['unhandled flag = ',num2str(flag)]);
 
  end
 
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
function [sys,x0,str,ts] = mdlInitializeSizes(Tsim)
 
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;     % states are stored in persistent vars
    sizes.NumOutputs     = 6;
    sizes.NumInputs      = 1;
    sizes.DirFeedthrough = 1;     % input can directly affect output 
    sizes.NumSampleTimes = 1;
 
    sys = simsizes(sizes);
    str = [];
    x0  = [];
    ts  = [Tsim 0]; % sample time: [period, offset]
                    % run every continuous time step
                 
%=============================================================================
% mdlUpdate
% Compute discrete time update.
%=============================================================================
function vars = mdlUpdate(t,u,vars)
    
    vars.delay_count = vars.delay_count + 1;
    
    if ( vars.delay_count >= vars.num_delays )
 
        vars.delay_count = 0;           % reset the counter
        
        sin_val = sin(vars.w*t);
        cos_val = cos(vars.w*t);
        vars.sin_out = sin_val;
        
        % x(k+1) = x(k) + Td(Ax(k) + Bu(k))
        vars.xa = vars.xa + vars.Td*(vars.A*vars.xa + vars.B*sin_val*u);
        vars.xb = vars.xb + vars.Td*(vars.A*vars.xb + vars.B*cos_val*u);
        
        % Send xa and xb to sin/cos states for output
        vars.sin_state = vars.xa(1);
        vars.cos_state = vars.xb(1);
        
        %% Running average (DFT)
        %{
        
        Runnign average (DFT) works by summing together a correlation
        with the controller output and sin/cos values. After the algorithm
        has met a certain 'chunk' size, a conditional is met and the chunk
        size is averaged. Note: this works better with higher number of 
        periods.
        
        %}
        
        vars.N = vars.N + 1; % Keep track of # of samples
    
        % Running sum of controller input correlation
        vars.sin_sum = vars.sin_sum + u*sin_val;
        vars.cos_sum = vars.cos_sum + u*cos_val;
    
        if (vars.N >= vars.chunk_size)

            % Calculate sin and cos average values
            
            vars.a = (2*vars.cos_sum/vars.N);
            vars.b = (2*vars.sin_sum/vars.N);

            % Reset sum start
            vars.cos_sum = 0;
            vars.sin_sum = 0;

            % Reset parameters
            vars.N = 0;

            % Check for termination
            vars = check_for_termination(vars,t);

        end
  
    end 
    
%=============================================================================
% mdlOutputs
% Return the output vector for the S-function
%=============================================================================
function vars = check_for_termination(vars,t)
    % keep track of the last 4 calculations of the magnitude from DFT
    vars.mag(4) = vars.mag(3);
    vars.mag(3) = vars.mag(2);
    vars.mag(2) = vars.mag(1);
    vars.mag(1) = sqrt(vars.a^2 + vars.b^2);
    average_mag = (vars.mag(1) + vars.mag(2) + vars.mag(1) + vars.mag(2))/4;
    
    % find the fractional deviation of last 4 from the average
    delta(1) = abs((average_mag-vars.mag(1))/average_mag);
    delta(2) = abs((average_mag-vars.mag(2))/average_mag);
    delta(3) = abs((average_mag-vars.mag(3))/average_mag);
    delta(4) = abs((average_mag-vars.mag(4))/average_mag);
    
    % keep track of the number of times DFT has been calculated
    vars.num_DFTs = vars.num_DFTs + 1;
    
    vars.termination = 0;
    % if 4 last magnitudes are with 10% of the average then terminate normally
    if (delta(1)<0.01)&&(delta(2)<0.01)&&(delta(3)<0.01)&&(delta(4)<0.01)
        vars.termination = 1;
    end
    % if 1000 DFTs have been calculated then terminate with a warning
    if (vars.num_DFTs > 1000)
        vars.termination = 2;
    end
    % if time is longer than 90 seconds then terminate with a warning
    if (t>90)
        vars.termination = 3;
    end
        
    % if terminating then assign values in workspace
    if (vars.termination ~= 0)
        assignin('base','sim_termination',vars.termination);
        assignin('base','sim_mag_ratio',average_mag);
        assignin('base','sim_phase',atan2(vars.a,vars.b));
        assignin('base','a',vars.a);
        assignin('base','b',vars.b);
    end

function sys = mdlOutputs(t,u,vars)
 
    sys = [vars.sin_out;
           2*vars.sin_state(1);
           2*vars.cos_state(1);
           vars.a;
           vars.b;
           vars.termination]; 
