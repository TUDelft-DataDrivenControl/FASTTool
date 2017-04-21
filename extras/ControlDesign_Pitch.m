% Example: Controller design for NREL 5MW Model wind turbine
% S.P. Mulders (s.p.mulders@tudelft.nl), J.W. van Wingerden
% (j.w.vanwingerden@tudelft.nl)

%% Ask for file names
[FullFileName,FullPathName] = uigetfile('*.mat', 'Select full wind turbine model .mat-file');
[LinFileName,LinPathName] = uigetfile('*.mat', 'Select linearized wind turbine model .mat-file');
P_plotting = true;
FA_plotting = true;

% Load linearized model
load([FullPathName, FullFileName])
load([LinPathName, LinFileName])

% PITCH: design gain-scheduled PI pitch controller for wind speeds where
% pitch angle ~= fine pitch
P_FromTo = [(find((diff(Lin.Pitch) > 0) == 1, 1) + 1) length(Lin.V)];
P_Kp = -0.003*ones(1, diff(P_FromTo)+1);
P_Ki = -0.001*ones(1, diff(P_FromTo)+1);
FA_Gain = -0.01;

% Create continuous model from linearized state-space matrices
for i=1:length(sysm)
    SYSTURB.A(:,:,i) = sysm{i}.A;
    SYSTURB.B(:,:,i) = sysm{i}.B;
    SYSTURB.C(:,:,i) = sysm{i}.C;
    SYSTURB.D(:,:,i) = sysm{i}.D;
end
wtg = ss(...
    SYSTURB.A, ...
    SYSTURB.B, ...
    SYSTURB.C, ...
    SYSTURB.D, ...
    'InputName', sysm{1}.inputname, ...
    'OutputName', sysm{1}.outputname, ...
    'StateName', sysm{1}.statename);

% Convert generator speed to rad/s
C_InputNames = {'Horizontal wind speed [m/s]'; 'Generator torque [Nm]'; 'Collective blade pitch [rad]'};
C_OutputNames = {'Generator speed [rad/s]'; 'Nacelle acceleration x-axis [m/s^2]'};
C_InputSelect = [1 8 9];
C_OutputSelect = [33 12];

wtg = wtg(C_OutputSelect, C_InputSelect); % Select input/outputs
wtg.C(1,:) = wtg.C(1,:)*pi/30; % Convert speed from RPM to rad/s
wtg.B(:,3,:) = wtg.B(:,3,:); % Convert pitch angle from deg to rad

wtg.OutputName = C_OutputNames; % Set input names in state-space model
wtg.InputName = C_InputNames; % Set output names in state-space model

%% Filter design
% Design low-pass filter (to remove high-frequency dynamics)
LP = zpk([], -Control.LPFCutOff, Control.LPFCutOff); % Make transfer function of the low-pass filter

% Discretize the low-pass filter (for implementation reasons)
dLP = ss(c2d(LP, Control.DT)); % Discretization of the first order low-pass filter

%% PITCH: Plots
if P_plotting
    % Creating transfer functions of the pitch controller
    for i = 1:(diff(P_FromTo)+1)
        PIpitch(:,:,i) = tf([P_Kp(i) P_Ki(i)],[1 0]);
    end

    % Select pitch and wind speed (> rated) as inputs and rotor speed as output
    P_WTG = minreal(wtg(1, 3, P_FromTo(1):P_FromTo(2))); % Only pitch to generator speed
    P_WTG_wind = minreal(wtg(1, [1 3], P_FromTo(1):P_FromTo(2))); % Wind speed and pitch to generator speeed

    % Create open-loop systems
    P_OL = LP*PIpitch*P_WTG; % Open-loop with controller, only pitch to generator speed
    P_OL_wind = series(append(1,LP*PIpitch),[P_WTG_wind; [0 tf([1 0],1)]]); % Open-loop with controller, wind speed and pitch to generator speeed

    % Create closed-loop systems
    P_CL = feedback(P_OL,1);
    P_CL_wind = feedback(P_OL_wind,1,2,1,-1);

    % Set bode plotting options
    opts = bodeoptions('cstprefs');    
    opts.Grid = 'on';
    opts.Title.String = 'Open-loop WT, from pitch [rad] to generator speed [rad/s]';
    opts.PhaseWrapping = 'on';
    opts.XLimMode = 'manual';
    opts.XLim = [1e-2 1e2];
    
    % Check stability using the gain (> 2) and phase (> 45 deg) margins
    figure
    subplot(1,3,1)
    bode(P_WTG, opts)
    
    opts.Title.String = 'Open-loop WT*C, from pitch [rad] to generator speed [rad/s]';
    subplot(1,3,2)
    bode(P_OL, opts)
    [Gm,Pm] = margin(P_OL);
    fprintf('Stability check for PI pitch controller: \n')
    P_WindSpeeds = Lin.V(P_FromTo(1)):Lin.V(P_FromTo(2));
    for i = 1:(diff(P_FromTo)+1)
        fprintf('[%i m/s] Gain margin: %2.2f, Phase margin: %2.2f deg; \n', P_WindSpeeds(i), 20*log10(Gm(i)), Pm(i))
    end

    % Check the performance (overshoot < 2 rad/s) with wind gust from IEC 61400-1
    t = 0:0.005:10.5;
    vgust = 1; %0.33*(0.11*Vnom/(1+0.1*(0.7*H/D)));
    v = -0.37 * vgust .* sin(3*pi.*t./10.5) .* (1-cos(2*pi.*t./10.5)) + 0.4;
    t = 0:0.005:50;
    u = zeros(1,length(t)) + 0.4;
    u(1:length(v)) = v;
    subplot(1,3,3)
    lsim(P_CL_wind(1,1,:), u, t)
    title('Closed-loop system responses to gust in wind')
    hold on
    plot(t, [2*ones(1,length(t)); -2*ones(1,length(t))], 'r')
    legend('Closed-loop system responses', '2 RPM overshoot boundary')
end

%% FORE-AFT: Plots
if FA_plotting
    % Select correct input/outputs from the linear models
    FA_WTG = minreal(wtg(2,3,:));
    FA_WTG_wind = minreal(wtg(2,[1 3],:));
    
    % Fore-aft damping controller
    FA_C = tf(1,[1 0])*FA_Gain;
    
    % Open-loop plant*controller
    FA_OL = FA_WTG*FA_C;

    % Set bode plotting options
    opts = bodeoptions('cstprefs');    
    opts.Grid = 'on';
    opts.Title.String = 'Open-loop WT, from pitch [rad] to nacelle accelleration [m/s^2]';
    opts.PhaseWrapping = 'on';
    opts.XLimMode = 'manual';
    opts.XLim = [1e-2 1e2];
    
    % Create closed-loop and sensitivity function
    FA_CL = feedback(FA_WTG,FA_C);
    FA_S = 1/(1+FA_WTG*FA_C);

    % Bode plot
    figure
    subplot(1,2,1)
    bode(FA_WTG, opts)
    
    subplot(1,2,2)
    opts.Title.String = 'Sensitivity, from pitch [rad] to nacelle accelleration [m/s^2]';
    bode(FA_S, opts)
end