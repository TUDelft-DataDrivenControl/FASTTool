%% Initialization code
function varargout = Linearization(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Linearization_OpeningFcn, ...
                   'gui_OutputFcn',  @Linearization_OutputFcn, ...
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

%% Opening function
function Linearization_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Blade = varargin{1};
handles.Airfoil = varargin{2};
handles.Tower = varargin{3};
handles.Nacelle = varargin{4};
handles.Control = varargin{5};
handles.Drivetrain = varargin{6};
handles.AirDensity = varargin{7};

% Propose values for the wind speed range
set(handles.WindSpeed_From, 'String', int2str(ceil(handles.Control.WindSpeed.Cutin)))
set(handles.WindSpeed_To, 'String', int2str(floor(handles.Control.WindSpeed.Cutout)))
set(handles.WindSpeed_Step, 'String', '1')
set(handles.LinAmount, 'String', '1')
set(handles.LinRotations, 'String', '1')

% Check initial status of LinearizeAboveRated only checkbox
if handles.LinearizeAboveRatedOnly_checkbox.Value
    handles.LinearizationMode = 'LinearizeAboveRated';
else
    handles.LinearizationMode = 'Linearize';
end
    
% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.Linearization);

%% Closing function
function Linearization_CloseRequestFcn(hObject, eventdata, handles)

uiresume(hObject);

%% Close button
function Close_Callback(hObject, eventdata, handles)

uiresume(handles.Linearization);

%% Output function
function varargout = Linearization_OutputFcn(hObject, eventdata, handles) 

% Close figure
delete(hObject)

%% Set full load wind speed range
function SetFullLoad_Callback(hObject, eventdata, handles)

% Set input
Blade = handles.Blade;
Airfoil = handles.Airfoil;
Drivetrain = handles.Drivetrain;
Control = handles.Control;
AirDensity = handles.AirDensity;

% Find rated wind speed by comparing demanded torque at rotational speed C
disp('Searching for rated wind speed...')
OmegaC = Control.Torque.SpeedC*2*pi/60;
U1 = Control.WindSpeed.Cutin;
U2 = Control.WindSpeed.Cutout;
[~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, (OmegaC/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/U1);
Qr1 = CalculateAerodynamicTorque(handles, CQr, U1);
[~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, (OmegaC/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/U2);
Qr2 = CalculateAerodynamicTorque(handles, CQr, U2);

if Qr1 > Control.Torque.Demanded
    Urated = Control.WindSpeed.Cutin;
elseif Qr2 < Control.Torque.Demanded
    Urated = Control.WindSpeed.Cutout;
else
    success = false;
    for iter = 1:100
        U_new = U1 + abs((Control.Torque.Demanded-Qr1)/(Qr1-Qr2))*(U2-U1);
        [~, CQr] = PerformanceCoefficients(Blade, Airfoil, Control.Pitch.Fine, (OmegaC/Drivetrain.Gearbox.Ratio)*Blade.Radius(end)/U_new);
        Qr_new = 0.5*CQr*AirDensity*U_new^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;
        if abs((Qr_new - Control.Torque.Demanded)/Control.Torque.Demanded) < 0.005
            success = true;
            break
        end
        if Qr_new < Control.Torque.Demanded
           U1 = U_new;
           Qr1 = Qr_new;
        else
           U2 = U_new;
           Qr2 = Qr_new;
        end
    end
    Urated = U_new;
    if ~success
        warning('Tolerance not met during iteration of axial induction factor')
    else
        % Propose values for the wind speed range
        set(handles.WindSpeed_From, 'String', int2str(ceil(Urated)));
        set(handles.WindSpeed_To, 'String', int2str(floor(handles.Control.WindSpeed.Cutout)))
        disp(' ')
        disp('... full load wind speed range set')
    end
end

%% Set file names
function SetLinModel_Callback(hObject, eventdata, handles)

% Set file name
[FileName,PathName] = uiputfile('*.mat', 'Set file name for linearized model');
handles.LinModel = [PathName, FileName];

if FileName
    % Enable buttons
    set(handles.Linearize, 'Enable', 'on');

    % Display file name
    set(handles.FileName_text, 'String', FileName)

    % Store in handles
    guidata(hObject, handles);
end

%% Linearization
function Linearize_Callback(hObject, eventdata, handles)

% Wind speed range
WindSpeeds = str2double(get(handles.WindSpeed_From, 'String')):str2double(get(handles.WindSpeed_Step, 'String')):str2double(get(handles.WindSpeed_To, 'String'));
TSim = str2double(get(handles.LinTime, 'String'));
if sum(WindSpeeds) == 0 || isnan(sum(WindSpeeds))
    errordlg('Invalid wind speed range', 'Error')
elseif min(WindSpeeds) < handles.Control.WindSpeed.Cutin || max(WindSpeeds) > handles.Control.WindSpeed.Cutout
    errordlg('Wind speeds must be inside operational range', 'Error')
else
    
% Disable window
buttons = findall(handles.Linearization, 'Type', 'UIControl');
for j = 1:length(buttons)
    set(buttons(j), 'Enable', 'off');
end
pause(0.1)

LinAmount = str2double(get(handles.LinAmount, 'String'));
LinRotations = str2double(get(handles.LinRotations, 'String'));
LinMode = handles.LinearizationMode;

disp('Starting linearization...')
disp(' ')

% Get geometry from handles
Blade = handles.Blade;
Airfoil = handles.Airfoil;
Tower = handles.Tower;
Nacelle = handles.Nacelle;
Control = handles.Control;
Drivetrain = handles.Drivetrain;
AirDensity = handles.AirDensity;

% Run modal analysis for 0 rpm
disp('Finding non-rotating mode shapes...')

% Run BModes for tower
evalc([...
'[y11_shape, y11_coeff, y11_freq,', ...
' y12_shape, y12_coeff, y12_freq,', ...
' y21_shape, y21_coeff, y21_freq,', ...
' y22_shape, y22_coeff, y22_freq] = BModes(Blade,Tower,Nacelle,Control,2,0);']);

% Store in handles
Tower.ForeAft1_coeff = y21_coeff;
Tower.ForeAft2_coeff = y22_coeff;
Tower.SideSide1_coeff = y11_coeff;
Tower.SideSide2_coeff = y12_coeff;

% Run BModes for blade
evalc([...
'[y11_shape, y11_coeff, y11_freq,', ...
' y12_shape, y12_coeff, y12_freq,', ...
' y21_shape, y21_coeff, y21_freq,', ...
' y22_shape, y22_coeff, y22_freq] = BModes(Blade,Tower,Nacelle,Control,1,0);']);

% Store in handles
Blade.Flap1_coeff = y11_coeff;
Blade.Flap2_coeff = y12_coeff;
Blade.Edge1_coeff = y21_coeff;
Blade.Edge2_coeff = y22_coeff;

% Steady state curves
disp('Determining steady state rotational speeds and pitch angles...')
[~, ~, OmegaU, PitchAngle] = SteadyState(Blade, Airfoil, Drivetrain, Control, WindSpeeds, AirDensity);
RPM = OmegaU * 60/(2*pi);

% Avoid some common errors with linearization
% (1) Set the gearbox efficiency (to avoid error that ADAMS cannot handle
% nonideal gearboxes) and remember the true efficiency to reset later
TrueGearboxEfficiency = Drivetrain.Gearbox.Efficiency;
Drivetrain.Gearbox.Efficiency = 1;

%{
% (2) Scale the HSS inertia from the reference machine (to avoid issues
% with reaching very high rotor speeds during the linearization of some
% direct-drive machines). Assuming that P ~ Mass, Radius? ~ P^(2/3), we get
% a relation in the shape of I/Iref = (P/Pref)^(5/3) * (RPM_ref/RPM)^2.
Prated = Control.Torque.SpeedC*(2*pi/60) *  Control.Torque.Demanded * Drivetrain.Generator.Efficiency;
Drivetrain.Generator.HSSInertia = 534.116 * (Prated/(5e6))^(5/3) * (12.1*97/Control.Torque.SpeedC)^2;
%}

% Turbine input file
disp('Writing input files...')

% Simulink settings
FAST_InputFileName = [pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'FAST.fst'];

% Turbine input files
AeroDyn(Blade,Airfoil,Tower,LinMode,AirDensity);
ServoDyn(Drivetrain,Control,LinMode);

% Preload the OutList
load([pwd, filesep 'subfunctions' filesep 'OutList.mat'])
assignin('base', 'OutList', OutList);

% Run linearization
sysm = cell(length(WindSpeeds),1);
for j = 1:length(WindSpeeds)
        
    % Status update
    disp(['Linearizing at U = ', num2str(WindSpeeds(j), '%5.2f'), ' m/s, ', num2str(RPM(j), '%5.2f'), ' rpm, ', num2str(PitchAngle(j), '%5.2f'), ' deg pitch'])
    
    % Set initial RPM and pitch angle in ElastoDyn input file
    ElastoDyn(Blade,Tower,Nacelle,Drivetrain,Control,LinMode,RPM(j),PitchAngle(j));
    
    % Set linearization times for 10 deg azimuth step (after 30 s)
    LinAziPositions = linspace(0,360*LinRotations,LinAmount+1);
    LinTimes = TSim + Control.DT * round(LinAziPositions(2:end)/(RPM(j)*6) / Control.DT);
    TMax = max(LinTimes)+1.0;
    
    FASTinput(Control.DT, TMax, LinMode, LinTimes);

    % Wind input file
    Wind.Type = 1;
    Wind.T = TMax;
    InflowWind(Wind,WindSpeeds(j),Tower.HubHeight,Blade.Radius(end))

    % Run FAST and prevent console output
    assignin('base', 'TMax', TMax);
    assignin('base', 'FAST_InputFileName', FAST_InputFileName);
    evalc('sim(''OpenLoop'',TMax);');

% Extract steady state solution after 60 seconds and average over 36
% steps
    A = 0;
    B = 0;
    C = 0;
    D = 0;
    for i = 1:LinAmount
        LinName = [pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'FAST.SFunc.', int2str(i), '.lin'];
        data = ReadFASTLinear(LinName);
        A = A + 1/LinAmount * data.A;
        B = B + 1/LinAmount * data.B;
        C = C + 1/LinAmount * data.C;
        D = D + 1/LinAmount * data.D;
    end
    sysm{j} = ss(A, B, C, D, 'InputName', data.u_desc,'Outputname', data.y_desc, 'StateName', data.x_desc);
    
    Lin.V(j) = data.y_op{1};
    Lin.Torque(j) = data.y_op{5};
    Lin.Pitch(j) =  data.y_op{34}*pi/180;
    Lin.GSpeed(j) = data.y_op{33}*pi/30;
    Lin.RSpeed(j) = data.y_op{38}*pi/30;
	Lin.x_op{j} = cell2mat(data.x_op);
    Lin.y_op{j} = cell2mat(data.y_op);
    Lin.u_op{j} = cell2mat(data.u_op);
end

% Reset the gearbox efficiency
Drivetrain.Gearbox.Efficiency = TrueGearboxEfficiency;

save(handles.LinModel, 'sysm', 'Lin');
disp(' ')
disp(' ')
disp('... Linearization complete!')

% Enable window
for j = 1:length(buttons)
    set(buttons(j), 'Enable', 'on');
end

end

%% Aerodynamic torque
function Q = CalculateAerodynamicTorque(handles, CQr, windSpeed)
    % Set input
    Blade = handles.Blade;
    Drivetrain = handles.Drivetrain;    
    AirDensity = handles.AirDensity;
    
    % Calculate aerodynamic torque
    Q = 0.5*CQr*AirDensity*windSpeed^2*pi*Blade.Radius(end)^3*Drivetrain.Gearbox.Efficiency/Drivetrain.Gearbox.Ratio;

%% Wind speed steps
function WindSpeed_From_Callback(hObject, eventdata, handles)
function WindSpeed_From_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function WindSpeed_To_Callback(hObject, eventdata, handles)
function WindSpeed_To_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function WindSpeed_Step_Callback(hObject, eventdata, handles)
function WindSpeed_Step_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LinAmount_Callback(hObject, eventdata, handles)
function LinAmount_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LinRotations_Callback(hObject, eventdata, handles)
function LinRotations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LinearizeAboveRatedOnly_checkbox.
function LinearizeAboveRatedOnly_checkbox_Callback(hObject, eventdata, handles)
if get(hObject,'Value') 
    % Set above-rated wind speed region
    SetFullLoad_Callback(hObject, eventdata, handles)
    
    handles.LinearizationMode = 'LinearizeAboveRated';
else
    % Propose values for the wind speed range
    set(handles.WindSpeed_From, 'String', int2str(ceil(handles.Control.WindSpeed.Cutin)))
    set(handles.WindSpeed_To, 'String', int2str(floor(handles.Control.WindSpeed.Cutout)))
    
    handles.LinearizationMode = 'Linearize';
end

% Update handles structure
guidata(hObject, handles);



function LinTime_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function LinTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
