function ControllerDesign(Control)

% Filter design
% Design low-pass filter (to remove high-frequency dynamics)
LP = zpk([], -Control.LPFCutOff, Control.LPFCutOff); % Make transfer function of the low-pass filter

% Discretize the low-pass filter (for implementation reasons)
dLP = ss(c2d(LP, Control.DT)); % Discretization of the first order low-pass filter

assignin('base', 'LP', LP);
assignin('base', 'dLP', dLP);