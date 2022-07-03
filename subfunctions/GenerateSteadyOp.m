%% Initialization code
function varargout = GenerateSteadyOp(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GenerateSteadyOp_OpeningFcn, ...
                   'gui_OutputFcn',  @GenerateSteadyOp_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%% Opening function
function GenerateSteadyOp_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for GenerateSteadyOp
%handles.output = hObject;

% Get input
handles.Blade = varargin{1};
handles.Airfoil = varargin{2};
handles.Drivetrain = varargin{3};
handles.Control = varargin{4};
handles.AirDensity = varargin{5};

% Propose values for the pitch angle range
set(handles.PitchAngle_From, 'String', '0')
set(handles.PitchAngle_From, 'String', '25')
set(handles.PitchAngle_From, 'String', '5')

% Propose values for the wind speed range
set(handles.WindSpeed_From, 'String', '0')
set(handles.WindSpeed_To, 'String', int2str(5*ceil(handles.Control.WindSpeed.Cutout/5) + 5))
set(handles.WindSpeed_Step, 'String', '0.1')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GenerateSteadyOp wait for user response (see UIRESUME)
uiwait(handles.GenerateSteadyOp);

%% Output function
function varargout = GenerateSteadyOp_OutputFcn(hObject, eventdata, handles) 

%varargout{1} = handles.output;
delete(hObject)

%% Checkbox for calculation of rotor performance is toggled
function CalcPerformance_checkbox_Callback(hObject, eventdata, handles)

% Toggle status of editing pitch angle range
if get(hObject, 'Value') == 1
    set(handles.PitchAngle_From, 'Enable', 'on');
    set(handles.PitchAngle_To, 'Enable', 'on');
    set(handles.PitchAngle_Step, 'Enable', 'on');
else
    set(handles.PitchAngle_From, 'Enable', 'off');
    set(handles.PitchAngle_To, 'Enable', 'off');
    set(handles.PitchAngle_Step, 'Enable', 'off');
end

%% Checkbox for calculation of turbine operation is toggled
function CalcOperation_checkbox_Callback(hObject, eventdata, handles)

% Toggle status of editing wind speed range
if get(hObject, 'Value') == 1
    set(handles.WindSpeed_From, 'Enable', 'on');
    set(handles.WindSpeed_To, 'Enable', 'on');
    set(handles.WindSpeed_Step, 'Enable', 'on');
else
    set(handles.WindSpeed_From, 'Enable', 'off');
    set(handles.WindSpeed_To, 'Enable', 'off');
    set(handles.WindSpeed_Step, 'Enable', 'off');
end

%% Cancel without calculations
function Cancel_Callback(hObject, eventdata, handles)

uiresume(handles.GenerateSteadyOp);

%% Perform the calculations
function Calculate_Callback(hObject, eventdata, handles)
% hObject    handle to Calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disable window
buttons = findall(handles.GenerateSteadyOp, 'Type', 'UIControl');
for j = 1:length(buttons)
    set(buttons(j), 'Enable', 'off');
end
pause(0.1)

% Get geometry and air density from handles
Blade = handles.Blade;
Airfoil = handles.Airfoil;
Drivetrain = handles.Drivetrain;
Control = handles.Control;
AirDensity = handles.AirDensity;

Pitchi = [];
if get(handles.CalcPerformanceFine_checkbox, 'Value') == 1
    Pitchi = [Pitchi, Control.Pitch.Fine];
end

if get(handles.CalcPerformance_checkbox, 'Value') == 1
    % Pitch angle range
    Beta = str2double(get(handles.PitchAngle_From, 'String')):str2double(get(handles.PitchAngle_Step, 'String')):str2double(get(handles.PitchAngle_To, 'String'));
    if isnan(sum(Beta))
        errordlg('Invalid pitch angle range.', 'Error')
    else
        Pitchi = [Pitchi, Beta];
    end
end

% Find rotor performance coefficients, when requested
if ~isempty(Pitchi)

    disp('Calculating rotor performance...')

    TSRj = 0:0.1:20;
    CPij = zeros(length(Pitchi), length(TSRj));
    CTij = zeros(length(Pitchi), length(TSRj));
    CQij = zeros(length(Pitchi), length(TSRj));
    lgd_label = cell(1, length(Pitchi));
    
    for i = 1:length(Pitchi)
        disp(['Pitch angle = ', num2str(Pitchi(i))])
        lgd_label{i} = [num2str(Pitchi(i), '%.1f'), ' deg'];
        for j = 11:length(TSRj) % Tip speed ratios below 1 omitted to avoid problems in the performance calculations

            if (TSRj(j) - ceil(TSRj(j))) == 0
                disp(['TSR = ', num2str(TSRj(j))])
            end

            [CTij(i,j), CQij(i,j)] = PerformanceCoefficients(Blade, Airfoil, Pitchi(i), TSRj(j));
            CPij(i,j) = CQij(i,j)*TSRj(j);

            if CPij(i,j) < 0 && TSRj(j) > 1
                CPij(i,j) = 0;
                CQij(i,j) = 0;
                break
            end
        end        
    end
    
    % Plot
    disp('Drawing plots...')
    Plot = figure();
    set(Plot, 'Name', 'Power coefficient curve(s)')
    plot(TSRj,CPij(1,:))
    xlim([0 max(TSRj)])
    ylim([0 ceil(20*max(max(CPij)))/20+0.05])
    set(gca, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'Fontsize', 8);
    xlabel('Tip speed ratio [-]')
    ylabel('Power coefficient [-]')
    hold on
    for i = 2:length(Pitchi)
        plot(TSRj,CPij(i,:))
    end
    hold off
    grid on
    legend(lgd_label)
    pause(0.1)
    
    % Send data to Matlab workspace
    assignin('base', 'Rotor_Pitch', Pitchi(:));
    assignin('base', 'Rotor_Lamda', TSRj(:));
    assignin('base', 'Rotor_cP', CPij);
    assignin('base', 'Rotor_cT', CTij);
    assignin('base', 'Rotor_cQ', CQij);
end

% Find turbine operation conditions, when requested
if get(handles.CalcOperation_checkbox, 'Value') == 1
% Wind speed range
U = str2double(get(handles.WindSpeed_From, 'String')):str2double(get(handles.WindSpeed_Step, 'String')):str2double(get(handles.WindSpeed_To, 'String'));
if sum(U) == 0 || isnan(sum(U))
    errordlg('Invalid wind speed range.', 'Error')
else
    disp('Generating steady operating curves...')
    
    % Determine steady state curves (external function, which is also used in 'Linearization.m'
    [CT, CQ, OmegaU, PitchAngle] = SteadyState(Blade, Airfoil, Drivetrain, Control, U, AirDensity);

    TSR = OmegaU.*Blade.Radius(end)./U;
    CP = CQ.*TSR;
    TSR(1) = 0;
    CP(1) = 0;

    P = 0.5*AirDensity*pi*Blade.Radius(end)^2 * U.^3 .* CP .* Drivetrain.Gearbox.Efficiency .* Drivetrain.Generator.Efficiency;
    T = 0.5*AirDensity*pi*Blade.Radius(end)^2 * U.^2 .* CT;
    Q = 0.5*AirDensity*pi*Blade.Radius(end)^3 * U.^2 .* CQ;
    RPM = OmegaU * 60/(2*pi);

    Prated = max(P);

    % Plot
    disp('Drawing plots...')
    Plot = figure();
    set(Plot, 'Name', 'Steady power curve')
    plot(U,P/1e6)
    xlim([0 max(U)])
    ylim([0 ceil(Prated/1e6)+0.5])
    set(gca, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'Fontsize', 8);
    xlabel('Wind speed [m/s]')
    ylabel('Electrical power [MW]')

    % Send data to Matlab workspace
    assignin('base', 'WindSpeed', U(:));
    assignin('base', 'ElectricalPower', P(:));
    assignin('base', 'RotorThrust', T(:));
    assignin('base', 'RotorSpeed', RPM(:));
    assignin('base', 'RotorTorque', Q(:));
    assignin('base', 'PitchAngle', PitchAngle(:));
    assignin('base', 'TipSpeedRatio', TSR(:));
    assignin('base', 'PowerCoefficient', CP(:));
    assignin('base', 'ThrustCoefficient', CT(:));
    assignin('base', 'TorqueCoefficient', CQ(:));% Get geometry from handles
end
end

% Enable window
for j = 1:length(buttons)
    set(buttons(j), 'Enable', 'on');
end

uiresume(handles.GenerateSteadyOp);

%% Pitch angle steps
function PitchAngle_From_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function PitchAngle_To_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function PitchAngle_Step_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Wind speed steps
function WindSpeed_From_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function WindSpeed_To_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function WindSpeed_Step_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%{
% Unused functions
function CalcPerformanceFine_checkbox_Callback(hObject, eventdata, handles)
function PitchAngle_From_Callback(hObject, eventdata, handles)
function PitchAngle_To_Callback(hObject, eventdata, handles)
function PitchAngle_Step_Callback(hObject, eventdata, handles)
function WindSpeed_From_Callback(hObject, eventdata, handles)
function WindSpeed_To_Callback(hObject, eventdata, handles)
function WindSpeed_Step_Callback(hObject, eventdata, handles)
%}
