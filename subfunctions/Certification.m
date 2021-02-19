%% Initialization code
function varargout = Certification(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Certification_OpeningFcn, ...
                   'gui_OutputFcn',  @Certification_OutputFcn, ...
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
function Certification_OpeningFcn(hObject, eventdata, handles, varargin)

% Set background image
h = axes('Units', 'Normalized', 'position', [0 0 1 1]);
uistack(h, 'bottom');
img = imread(['graphics' filesep 'certification.png']);
imagesc(img);
set(h, 'HandleVisibility', 'off', 'visible','off')

% Get input
handles.Blade = varargin{1};
handles.Airfoil = varargin{2};
handles.Tower = varargin{3};
handles.Nacelle = varargin{4};
handles.Drivetrain = varargin{5};
handles.Control = varargin{6};
handles.CertificationSettings = varargin{7};
handles.AirDensity = varargin{8};

% Disable simulation time and wind speed if a stepped wind is requested
if handles.CertificationSettings.Wind.Type == 2
    set(handles.Runtime_textbox, 'Enable', 'off');
    set(handles.WindSpeed_textbox, 'Enable', 'off');
else
    set(handles.Runtime_textbox, 'Enable', 'on');
    set(handles.WindSpeed_textbox, 'Enable', 'on');
end

% Update input fields
set(handles.Runtime_textbox, 'String', num2str(handles.CertificationSettings.Run.Time));
handles.CertificationSettings.Wind.T = handles.CertificationSettings.Run.Time;
U = num2str(handles.CertificationSettings.Run.WindSpeed(1));
if length(handles.CertificationSettings.Run.WindSpeed) > 1
    for i = 2:length(handles.CertificationSettings.Run.WindSpeed)
        U = [U, ', ', num2str(handles.CertificationSettings.Run.WindSpeed(i))];
    end
end
set(handles.WindSpeed_textbox, 'String', U);
set(handles.NSeeds_textbox, 'String', num2str(handles.CertificationSettings.Run.Seeds));
set(handles.Edit_SeedNumber, 'String', num2str(handles.CertificationSettings.Run.SeedNumber));
handles.OutputFile = {};
set(handles.OutputFiles_text, 'String', '<empty>')

% Update wind icons
if handles.CertificationSettings.Wind.Type == 1       % Steady wind
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_steady.png']))
    set(handles.WindSettings, 'TooltipString', 'Steady wind')
elseif handles.CertificationSettings.Wind.Type == 2   % Stepped wind
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_step.png']))
    set(handles.WindSettings, 'TooltipString', 'Stepped wind')
elseif handles.CertificationSettings.Wind.Type == 3   % Normal wind profile (NWP)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_nwp.png']))
    set(handles.WindSettings, 'TooltipString', 'Normal wind profile (NWP)')
elseif handles.CertificationSettings.Wind.Type == 4   % Normal turbulence model (NTM)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_ntm.png']))
    set(handles.WindSettings, 'TooltipString', 'Normal turbulence model (NTM)')
elseif handles.CertificationSettings.Wind.Type == 5   % Annual extreme wind speed (EWM)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_ewm1.png']))
    set(handles.WindSettings, 'TooltipString', 'Annual extreme wind speed (EWM)')
elseif handles.CertificationSettings.Wind.Type == 6   % 50-year extreme wind speed (EWM)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_ewm50.png']))
    set(handles.WindSettings, 'TooltipString', '50-year extreme wind speed (EWM)')
elseif handles.CertificationSettings.Wind.Type == 7   % Extreme wind shear (EWS)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_ews.png']))
    set(handles.WindSettings, 'TooltipString', 'Extreme wind shear (EWS)')
elseif handles.CertificationSettings.Wind.Type == 8   % Extreme turbulence model (ETM)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_etm.png']))
    set(handles.WindSettings, 'TooltipString', 'Extreme turbulence model (ETM)')
elseif handles.CertificationSettings.Wind.Type == 9   % Extreme operating gust (EOG)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_eog.png']))
    set(handles.WindSettings, 'TooltipString', 'Extreme operating gust (EOG)')
elseif handles.CertificationSettings.Wind.Type == 10   % Extreme direction change (EDC)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_edc.png']))
    set(handles.WindSettings, 'TooltipString', 'Extreme direction change (EDC)')
elseif handles.CertificationSettings.Wind.Type == 11   % Extreme coherent gust (ECD)
    set(handles.WindSettings, 'CData', imread(['icons' filesep 'wind_ecg.png']))
    set(handles.WindSettings, 'TooltipString', 'Extreme coherent gust (ECD)')
end

% Update operations icons
if handles.CertificationSettings.Mode.Type == 1       % Power production
    set(handles.OperationSettings, 'CData', imread(['icons' filesep 'operation_prod.png']))
    set(handles.OperationSettings, 'TooltipString', 'Power production')
elseif handles.CertificationSettings.Mode.Type == 2   % Power production with fault
    set(handles.OperationSettings, 'CData', imread(['icons' filesep 'operation_prodfault.png']))
    set(handles.OperationSettings, 'TooltipString', 'Power production with fault')
elseif handles.CertificationSettings.Mode.Type == 3   % Startup
    set(handles.OperationSettings, 'CData', imread(['icons' filesep 'operation_start.png']))
    set(handles.OperationSettings, 'TooltipString', 'Startup')
elseif handles.CertificationSettings.Mode.Type == 4   % Normal shutdown
    set(handles.OperationSettings, 'CData', imread(['icons' filesep 'operation_shutdown.png']))
    set(handles.OperationSettings, 'TooltipString', 'Normal shutdown')
elseif handles.CertificationSettings.Mode.Type == 5   % Emergency shutdown
    set(handles.OperationSettings, 'CData', imread(['icons' filesep 'operation_emergency.png']))
    set(handles.OperationSettings, 'TooltipString', 'Emergency shutdown')
elseif handles.CertificationSettings.Mode.Type == 6   % Idling
    set(handles.OperationSettings, 'CData', imread(['icons' filesep 'operation_idle.png']))
    set(handles.OperationSettings, 'TooltipString', 'Idling')
elseif handles.CertificationSettings.Mode.Type == 7   % Parked
    set(handles.OperationSettings, 'CData', imread(['icons' filesep 'operation_idle.png']))
    set(handles.OperationSettings, 'TooltipString', 'Parked')
end

% Update random seed fields
if handles.CertificationSettings.Run.RandomSeed
    set(handles.Edit_SeedNumber, 'Enable', 'off')
    set(handles.Radio_RandomSeedYes, 'Value', 1)
    set(handles.Radio_RandomSeedNo, 'Value', 0)
else
    set(handles.Edit_SeedNumber, 'Enable', 'on')
    set(handles.Radio_RandomSeedYes, 'Value', 0)
    set(handles.Radio_RandomSeedNo, 'Value', 1)
end

% Disable start button
set(handles.Start, 'Enable', 'off')

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.Certification);

%% Output function
function varargout = Certification_OutputFcn(hObject, eventdata, handles) 

% Output
varargout{1} = handles.CertificationSettings;

% Close figure
delete(hObject)

%% Closing function
function Certification_CloseRequestFcn(hObject, eventdata, handles)
uiresume(hObject);

%% OK button
function Close_Callback(hObject, eventdata, handles)
uiresume(handles.Certification);

%% WindSettings settings
function WindSettings_Callback(hObject, eventdata, handles)

% Disable window
buttons = findall(handles.Certification, 'Type', 'UIControl');
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'off');
end

% Operation settings
handles.CertificationSettings = Wind(...
    handles.CertificationSettings, ...
    handles.Blade, ...
    handles.Tower);

% Enable window
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'on');
end

% Disable simulation time and wind speed if a stepped wind is requested
if handles.CertificationSettings.Wind.Type == 2
    set(handles.Runtime_textbox, 'Enable', 'off');
    set(handles.WindSpeed_textbox, 'Enable', 'off');
    set(handles.Runtime_textbox, 'String', handles.CertificationSettings.Run.Time);
    set(handles.WindSpeed_textbox, 'String', handles.CertificationSettings.Run.WindSpeed(1));
else
    set(handles.Runtime_textbox, 'Enable', 'on');
    set(handles.WindSpeed_textbox, 'Enable', 'on');
end

% Disable start button
if isempty(handles.OutputFile)
    set(handles.Start, 'Enable', 'off')
end

% Update icon
if handles.CertificationSettings.Wind.Type == 1       % Steady wind
    set(hObject, 'CData', imread(['icons' filesep 'wind_steady.png']))
    set(hObject, 'TooltipString', 'Steady wind')
elseif handles.CertificationSettings.Wind.Type == 2   % Stepped wind
    set(hObject, 'CData', imread(['icons' filesep 'wind_step.png']))
    set(hObject, 'TooltipString', 'Stepped wind')
elseif handles.CertificationSettings.Wind.Type == 3   % Normal wind profile (NWP)
    set(hObject, 'CData', imread(['icons' filesep 'wind_nwp.png']))
    set(hObject, 'TooltipString', 'Normal wind profile (NWP)')
elseif handles.CertificationSettings.Wind.Type == 4   % Normal turbulence model (NTM)
    set(hObject, 'CData', imread(['icons' filesep 'wind_ntm.png']))
    set(hObject, 'TooltipString', 'Normal turbulence model (NTM)')
elseif handles.CertificationSettings.Wind.Type == 5   % Annual extreme wind speed (EWM)
    set(hObject, 'CData', imread(['icons' filesep 'wind_ewm1.png']))
    set(hObject, 'TooltipString', 'Annual extreme wind speed (EWM)')
elseif handles.CertificationSettings.Wind.Type == 6   % 50-year extreme wind speed (EWM)
    set(hObject, 'CData', imread(['icons' filesep 'wind_ewm50.png']))
    set(hObject, 'TooltipString', '50-year extreme wind speed (EWM)')
elseif handles.CertificationSettings.Wind.Type == 7   % Extreme wind shear (EWS)
    set(hObject, 'CData', imread(['icons' filesep 'wind_ews.png']))
    set(hObject, 'TooltipString', 'Extreme wind shear (EWS)')
elseif handles.CertificationSettings.Wind.Type == 8   % Extreme turbulence model (ETM)
    set(hObject, 'CData', imread(['icons' filesep 'wind_etm.png']))
    set(hObject, 'TooltipString', 'Extreme turbulence model (ETM)')
elseif handles.CertificationSettings.Wind.Type == 9   % Extreme operating gust (EOG)
    set(hObject, 'CData', imread(['icons' filesep 'wind_eog.png']))
    set(hObject, 'TooltipString', 'Extreme operating gust (EOG)')
elseif handles.CertificationSettings.Wind.Type == 10   % Extreme direction change (EDC)
    set(hObject, 'CData', imread(['icons' filesep 'wind_edc.png']))
    set(hObject, 'TooltipString', 'Extreme direction change (EDC)')
elseif handles.CertificationSettings.Wind.Type == 11   % Extreme coherent gust (ECD)
    set(hObject, 'CData', imread(['icons' filesep 'wind_ecg.png']))
    set(hObject, 'TooltipString', 'Extreme coherent gust (ECD)')
end

% Update handles structure
guidata(hObject, handles);

%% Operation settings
function OperationSettings_Callback(hObject, eventdata, handles)

% Disable window
buttons = findall(handles.Certification, 'Type', 'UIControl');
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'off');
end

% Operation settings
handles.CertificationSettings = Operation(...
    handles.CertificationSettings);

% Enable window
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'on');
end

% Disable simulation time and wind speed if a stepped wind is requested
if handles.CertificationSettings.Wind.Type == 2
    set(handles.Runtime_textbox, 'Enable', 'off');
    set(handles.WindSpeed_textbox, 'Enable', 'off');
else
    set(handles.Runtime_textbox, 'Enable', 'on');
    set(handles.WindSpeed_textbox, 'Enable', 'on');
end

% Disable start button
if isempty(handles.OutputFile)
    set(handles.Start, 'Enable', 'off')
end

% Update icon
if handles.CertificationSettings.Mode.Type == 1       % Power production
    set(hObject, 'CData', imread(['icons' filesep 'operation_prod.png']))
    set(hObject, 'TooltipString', 'Power production')
elseif handles.CertificationSettings.Mode.Type == 2   % Power production with fault
    set(hObject, 'CData', imread(['icons' filesep 'operation_prodfault.png']))
    set(hObject, 'TooltipString', 'Power production with fault')
elseif handles.CertificationSettings.Mode.Type == 3   % Startup
    set(hObject, 'CData', imread(['icons' filesep 'operation_start.png']))
    set(hObject, 'TooltipString', 'Startup')
elseif handles.CertificationSettings.Mode.Type == 4   % Normal shutdown
    set(hObject, 'CData', imread(['icons' filesep 'operation_shutdown.png']))
    set(hObject, 'TooltipString', 'Normal shutdown')
elseif handles.CertificationSettings.Mode.Type == 5   % Emergency shutdown
    set(hObject, 'CData', imread(['icons' filesep 'operation_emergency.png']))
    set(hObject, 'TooltipString', 'Emergency shutdown')
elseif handles.CertificationSettings.Mode.Type == 6   % Idling
    set(hObject, 'CData', imread(['icons' filesep 'operation_idle.png']))
    set(hObject, 'TooltipString', 'Idling')
elseif handles.CertificationSettings.Mode.Type == 7   % Parked
    set(hObject, 'CData', imread(['icons' filesep 'operation_idle.png']))
    set(hObject, 'TooltipString', 'Parked')
end

% Update handles structure
guidata(hObject, handles);

%% Open Simulink interface
function Simulink_Callback(hObject, eventdata, handles)
open_system('FAST')

%% View output signals
function Output_Callback(hObject, eventdata, handles)
load_system('FAST')
open_system('FAST/Scope')

%% Wind speed range - text box
function WindSpeed_textbox_Callback(hObject, eventdata, handles)
try
    eval(['U = [', get(hObject, 'String'), '];']);
    if sum(isnan(U)) > 0 || isempty(isnan(U))
        U = num2str(handles.CertificationSettings.Run.WindSpeed(1));
        if length(handles.CertificationSettings.Run.WindSpeed) > 1
            for i = 2:length(handles.CertificationSettings.Run.WindSpeed)
                U = [U, ', ', num2str(handles.CertificationSettings.Run.WindSpeed(i))];
            end
        end
        set(handles.WindSpeed_textbox, 'String', U);
    end
catch
    U = num2str(handles.CertificationSettings.Run.WindSpeed(1));
    if length(handles.CertificationSettings.Run.WindSpeed) > 1
        for i = 2:length(handles.CertificationSettings.Run.WindSpeed)
            U = [U, ', ', num2str(handles.CertificationSettings.Run.WindSpeed(i))];
        end
    end
    set(handles.WindSpeed_textbox, 'String', U);
end

eval(['U = [', get(hObject, 'String'), '];']);
handles.CertificationSettings.Run.WindSpeed = sort(U);
U = num2str(handles.CertificationSettings.Run.WindSpeed(1));
if length(handles.CertificationSettings.Run.WindSpeed) > 1
    for i = 2:length(handles.CertificationSettings.Run.WindSpeed)
        U = [U, ', ', num2str(handles.CertificationSettings.Run.WindSpeed(i))];
    end
end
set(handles.WindSpeed_textbox, 'String', U);
guidata(hObject, handles);
UpdateOutputName(handles);

function WindSpeed_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Seed number - text box
function NSeeds_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Run.Seeds))
end
handles.CertificationSettings.Run.Seeds = ceil(str2double(get(hObject,'String')));
guidata(hObject, handles);
UpdateOutputName(handles);

function NSeeds_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Simulation runtime - text box
function Runtime_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Run.Time))
end
handles.CertificationSettings.Run.Time = str2double(get(hObject,'String'));
handles.CertificationSettings.Wind.T = handles.CertificationSettings.Run.Time;
guidata(hObject, handles);

function Runtime_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Set output file name
function SetOutputName_Callback(hObject, eventdata, handles)
[FileName,PathName] = uiputfile('*.mat', 'Set output file name(s)');
if FileName
    
    % Store in handles
    handles.OutputFile = {PathName, FileName(1:end-4)};
    guidata(hObject, handles);
    
    % Update output file names
    UpdateOutputName(handles);
    
    % Enable start button
    set(handles.Start, 'Enable', 'on')

end

%% Update output file names
function UpdateOutputName(handles)

if ~isempty(handles.OutputFile)

i = 1;
if length(handles.CertificationSettings.Run.WindSpeed) * handles.CertificationSettings.Run.Seeds <= 4
    for U = handles.CertificationSettings.Run.WindSpeed
        for seed  = 1:handles.CertificationSettings.Run.Seeds

            OutputFile = ['../', handles.OutputFile{2}];
            if length(handles.CertificationSettings.Run.WindSpeed) > 1
                OutputFile = [OutputFile, '_U=', num2str(U,'%2.2f')];
            end
            if handles.CertificationSettings.Run.Seeds > 1
                OutputFile = [OutputFile, '_seed=', int2str(seed)];
            end
            OutputFiles{i} = [OutputFile, '.mat'];
            i = i + 1;

        end
    end
else
    for U = handles.CertificationSettings.Run.WindSpeed
        for seed  = 1:handles.CertificationSettings.Run.Seeds

            OutputFile = ['../', handles.OutputFile{2}];
            if length(handles.CertificationSettings.Run.WindSpeed) > 1
                OutputFile = [OutputFile, '_U=', num2str(U,'%2.2f')];
            end
            if handles.CertificationSettings.Run.Seeds > 1
                OutputFile = [OutputFile, '_seed=', int2str(seed)];
            end
            OutputFiles{i} = [OutputFile, '.mat'];
            i = i + 1;

            if i == 3
                break
            end

        end

        if i == 3
            break
        end
    end

    OutputFiles{3} = '  ...';
    OutputFile = ['../', handles.OutputFile{2}];
    if length(handles.CertificationSettings.Run.WindSpeed) > 1
        OutputFile = [OutputFile, '_U=', num2str(handles.CertificationSettings.Run.WindSpeed(end),'%2.2f')];
    end
    if handles.CertificationSettings.Run.Seeds > 1
        OutputFile = [OutputFile, '_seed=', int2str(handles.CertificationSettings.Run.Seeds)];
    end
    OutputFiles{4} = [OutputFile, '.mat'];

end

set(handles.OutputFiles_text, 'String', OutputFiles)

end

%% Start simulation
function Start_Callback(hObject, eventdata, handles)

% Disable window - except for 'scope'
buttons = findall(handles.Certification, 'Type', 'UIControl');
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'off');
end
set(handles.Output, 'Enable', 'on');

% Set output filename
if isempty(handles.OutputFile)
    handles.OutputFile = [pwd, filesep 'output'];
end

% Temporarily turn off warnings
warning('off','all')

% Load parameters from handles
Blade = handles.Blade;
Airfoil = handles.Airfoil;
Tower = handles.Tower;
Nacelle = handles.Nacelle;
Drivetrain = handles.Drivetrain;
Control = handles.Control;
CertificationSettings = handles.CertificationSettings;
AirDensity = handles.AirDensity;
TMax = CertificationSettings.Run.Time;
T = CertificationSettings.Wind.T;
Ly = CertificationSettings.Wind.Ly;
Lz = CertificationSettings.Wind.Lz;
dt = CertificationSettings.Wind.dt;
Ny = CertificationSettings.Wind.Ny;
Nz = CertificationSettings.Wind.Nz;
windSeed = handles.CertificationSettings.Run.SeedNumber;

% Run modal analysis for 0 rpm
disp('Preparing input files...')

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
if CertificationSettings.Wind.Type == 2
    [~, ~, OmegaU, PitchAngle] = SteadyState(Blade, Airfoil, Drivetrain, Control, CertificationSettings.Wind.Step, AirDensity);
    OmegaU = OmegaU*ones(length(CertificationSettings.Run.WindSpeed));
    PitchAngle = PitchAngle*ones(length(CertificationSettings.Run.WindSpeed));
else
    [~, ~, OmegaU, PitchAngle] = SteadyState(Blade, Airfoil, Drivetrain, Control, CertificationSettings.Run.WindSpeed, AirDensity);
end
RPM = OmegaU * 60/(2*pi);

% Initialize controller
disp('Setting controller parameters...')
assignin('base', 'Drivetrain', Drivetrain)
assignin('base', 'Control', Control)

% Turbine input files
TMax = CertificationSettings.Run.Time;
FASTinput(Control.DT, TMax);
AeroDyn(Blade,Airfoil,Tower,string(CertificationSettings.Mode.Type), AirDensity);

% Send to base workspace and make structures available for Simulink and run the simulation
assignin('base', 'FAST_InputFileName', [pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'FAST.fst']);
assignin('base', 'TMax', TMax);
disp('')

% Loop over wind speeds and seeds
for j = 1:length(CertificationSettings.Run.WindSpeed)
    U = CertificationSettings.Run.WindSpeed(j);
    for seed = 1:CertificationSettings.Run.Seeds

        % Output file name
        OutputFile = [handles.OutputFile{1}, handles.OutputFile{2}];
        if length(handles.CertificationSettings.Run.WindSpeed) > 1
            OutputFile = [OutputFile, '_U=', num2str(U,'%2.2f')];
        end
        if handles.CertificationSettings.Run.Seeds > 1
            OutputFile = [OutputFile, '_seed=', int2str(seed)];
        end
        OutputFile = [OutputFile, '.mat'];

        % Find initial RPM and pitch angle
        if CertificationSettings.Wind.Type == 2
            Ui = CertificationSettings.Wind.Step;
        else
            Ui = U;
        end
        
        if Ui < Control.WindSpeed.Cutin || Ui > Control.WindSpeed.Cutout
            RPM_Init = 0;
            P_InitAngle = Control.Pitch.Max;
        else
            RPM_Init = RPM(j);
            P_InitAngle = PitchAngle(j);
        end
        
        if CertificationSettings.Mode.Type == 3     % Startup
            RPM_Init = 0;
            P_InitAngle = Control.Pitch.Max;
        elseif CertificationSettings.Mode.Type == 6 % Idling
            RPM_Init = 0;
            P_InitAngle = Control.Pitch.Max;
        elseif CertificationSettings.Mode.Type == 7	% Parked
            RPM_Init = 0;
            P_InitAngle = Control.Pitch.Max;
        end
        
        assignin('base', 'RPM_Init', RPM_Init);
        assignin('base', 'T_GenSpeedInit', RPM_Init*Drivetrain.Gearbox.Ratio);
        assignin('base', 'P_InitAngle', P_InitAngle);

        % Set operation mode in ElastoDyn file
        ElastoDyn(Blade,Tower,Nacelle,Drivetrain,Control,string(CertificationSettings.Mode.Type),RPM_Init,P_InitAngle);

        % Set operation mode in ServoDyn file
        ServoDyn(Drivetrain,Control,string(CertificationSettings.Mode.Type),CertificationSettings.Mode.Actiontime);

        % Wind input file
        disp('Generating wind file...')
        InflowWind(CertificationSettings.Wind,U,Tower.HubHeight,Blade.Radius(end),windSeed)

        % Preload the OutList
        load([pwd, filesep 'subfunctions' filesep 'OutList.mat'])
        assignin('base', 'OutList', OutList);
        assignin('base', 'CertificationSettings', CertificationSettings);

        % Call Simulink
        disp(['Running FAST (U = ', num2str(U,'%5.2f'), ' m/s, seed ', int2str(seed), '/', int2str(CertificationSettings.Run.Seeds), ')'])
        evalc('sim(''FAST'',TMax);');

        % Extract output
        filename = [pwd, filesep 'subfunctions' filesep 'inputfiles' filesep 'FAST.SFunc.out'];
        delimiter = '\t';
        startRow = 6;
        formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
        Output = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
        fclose(fileID);

        % Rename and save vectors
        save(OutputFile, 'Legend')
        for i = 1:length(OutList)
            eval([OutList{i}, ' = Output{i};']);
            eval(['save(OutputFile, ''', OutList{i}, ''', ''-append'');']);
        end

    end
end

% Enable window
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'on');
end

% Disable simulation time and wind speed if a stepped wind is requested
if handles.CertificationSettings.Wind.Type == 2
    set(handles.Runtime_textbox, 'Enable', 'off');
    set(handles.WindSpeed_textbox, 'Enable', 'off');
else
    set(handles.Runtime_textbox, 'Enable', 'on');
    set(handles.WindSpeed_textbox, 'Enable', 'on');
end

% Turn on warnings
warning('on','all')
disp('')
disp('Completed!')



function Edit_SeedNumber_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 1
    set(hObject, 'String', '1')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', '1')
end
handles.CertificationSettings.Run.SeedNumber = ceil(str2double(get(hObject,'String')));
guidata(hObject, handles);
UpdateOutputName(handles);


% --- Executes during object creation, after setting all properties.
function Edit_SeedNumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Radio_RandomSeedYes.
function Radio_RandomSeedYes_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.Edit_SeedNumber, 'Enable', 'off')
    handles.CertificationSettings.Run.RandomSeed = true;
end
guidata(hObject, handles);

% --- Executes on button press in Radio_RandomSeedNo.
function Radio_RandomSeedNo_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    set(handles.Edit_SeedNumber, 'Enable', 'on')
    handles.CertificationSettings.Run.RandomSeed = false;
end
guidata(hObject, handles);
