%% Initialization code
function varargout = ControlDesign(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ControlDesign_OpeningFcn, ...
                   'gui_OutputFcn',  @ControlDesign_OutputFcn, ...
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
function ControlDesign_OpeningFcn(hObject, eventdata, handles, varargin)

% Set background image
h = axes('Units', 'Normalized', 'position', [0 0 1 1]);
uistack(h, 'bottom');
img = imread('graphics\control.png');
imagesc(img);
set(h, 'HandleVisibility', 'off', 'visible','off')

% Get input
handles.Control = varargin{1};
handles.LinModel = varargin{2};
handles.Input = varargin(1);

% Update input fields
set(handles.WindSpeed_Cutin_textbox, 'String', num2str(handles.Control.WindSpeed.Cutin));
set(handles.WindSpeed_Cutout_textbox, 'String', num2str(handles.Control.WindSpeed.Cutout));
set(handles.Control_LPFCutOff_textbox, 'String', num2str(handles.Control.LPFCutOff));
set(handles.Control_DT_textbox, 'String', num2str(handles.Control.DT));
set(handles.Brake_Torque_textbox, 'String', num2str(handles.Control.Brake.Torque));
set(handles.Brake_Deploytime_textbox, 'String', num2str(handles.Control.Brake.Deploytime));
set(handles.Brake_Delay_textbox, 'String', num2str(handles.Control.Brake.Delay));
set(handles.Torque_Demanded_textbox, 'String', num2str(handles.Control.Torque.Demanded));
set(handles.Torque_Min_textbox, 'String', num2str(handles.Control.Torque.Min));
set(handles.Torque_Limit_textbox, 'String', num2str(handles.Control.Torque.Limit));
set(handles.Torque_Slewrate_textbox, 'String', num2str(handles.Control.Torque.Slewrate));
set(handles.Torque_SpeedA_textbox, 'String', num2str(handles.Control.Torque.SpeedA));
set(handles.Torque_SpeedB_textbox, 'String', num2str(handles.Control.Torque.SpeedB));
set(handles.Torque_SpeedB2_textbox, 'String', num2str(handles.Control.Torque.SpeedB2));
set(handles.Torque_SpeedC_textbox, 'String', num2str(handles.Control.Torque.SpeedC));
set(handles.Torque_OptGain_textbox, 'String', num2str(handles.Control.Torque.OptGain));
set(handles.Pitch_Fine_textbox, 'String', num2str(handles.Control.Pitch.Fine));
set(handles.Pitch_Max_textbox, 'String', num2str(handles.Control.Pitch.Max));
set(handles.Pitch_Min_textbox, 'String', num2str(handles.Control.Pitch.Min));
set(handles.Pitch_Maxrate_textbox, 'String', num2str(handles.Control.Pitch.Maxrate));
set(handles.Pitch_Minrate_textbox, 'String', num2str(handles.Control.Pitch.Minrate));
set(handles.Pitch_StartupPitch_textbox, 'String', num2str(handles.Control.Pitch.StartupPitch));
set(handles.Pitch_StartupPitchRate_textbox, 'String', num2str(handles.Control.Pitch.StartupPitchRate));
set(handles.Pitch_StartupSpeed_textbox, 'String', num2str(handles.Control.Pitch.StartupSpeed));
set(handles.ForeAft_MaxPitchAmplitude_textbox, 'String', num2str(handles.Control.ForeAft.MaxPitchAmplitude));
set(handles.ForeAft_Gain_textbox, 'String', num2str(handles.Control.ForeAft.Gain));

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.ControlDesign);

%% Closing function
function ControlDesign_CloseRequestFcn(hObject, eventdata, handles)
button = questdlg('Save changes?');
if strcmp(button, 'Yes')
    handles.Save = true;
    guidata(hObject, handles);
    uiresume(hObject);
elseif strcmp(button, 'No')
    handles.Save = false;
    guidata(hObject, handles);
    uiresume(hObject);
end

%% Save button
function Apply_Callback(hObject, eventdata, handles)
handles.Save = true;
guidata(hObject, handles);
uiresume(handles.ControlDesign);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.ControlDesign);

%% Output function
function varargout = ControlDesign_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save
    varargout{1} = handles.Control;
else
    varargout = handles.Input;
end

% Close figure
delete(hObject)

%% Cut-in wind speed - text box
function WindSpeed_Cutin_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.WindSpeed.Cutin))
end
handles.Control.WindSpeed.Cutin = str2double(get(hObject,'String'));
guidata(hObject, handles);
function WindSpeed_Cutin_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Cut-out wind speed - text box
function WindSpeed_Cutout_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.WindSpeed.Cutout))
end
handles.Control.WindSpeed.Cutout = str2double(get(hObject,'String'));
guidata(hObject, handles);
function WindSpeed_Cutout_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Low-pass filter cut-off frequency - text box
function Control_LPFCutOff_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.LPFCutOff))
end
handles.Control.LPFCutOff = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Control_LPFCutOff_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Sampling time of the controller and simulation - text box
function Control_DT_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.DT))
end
handles.Control.DT = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Control_DT_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Brake torque - text box
function Brake_Torque_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Brake.Torque))
end
handles.Control.Brake.Torque = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Brake_Torque_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Brake deploy time - text box
function Brake_Deploytime_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Brake.Deploytime))
end
handles.Control.Brake.Deploytime = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Brake_Deploytime_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Brake delay - text box
function Brake_Delay_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Brake.Delay))
end
handles.Control.Brake.Delay = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Brake_Delay_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% Demanded torque - text box
function Torque_Demanded_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.Demanded))
end
handles.Control.Torque.Demanded = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_Demanded_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Torque minimum - text box
function Torque_Min_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.Min))
end
handles.Control.Torque.Min = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_Min_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Torque limit - text box
function Torque_Limit_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.Limit))
end
handles.Control.Torque.Limit = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_Limit_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Torque slew rate - text box
function Torque_Slewrate_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.Slewrate))
end
handles.Control.Torque.Slewrate = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_Slewrate_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Omega A - text box
function Torque_SpeedA_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.SpeedA))
end
handles.Control.Torque.SpeedA = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_SpeedA_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Omega B - text box
function Torque_SpeedB_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.SpeedB))
end
handles.Control.Torque.SpeedB = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_SpeedB_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Omega B2 - text box
function Torque_SpeedB2_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.SpeedB2))
end
handles.Control.Torque.SpeedB2 = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_SpeedB2_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Omega C - text box
function Torque_SpeedC_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.SpeedC))
end
handles.Control.Torque.SpeedC = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_SpeedC_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Optimal mode gain - text box
function Torque_OptGain_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Torque.OptGain))
end
handles.Control.Torque.OptGain = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Torque_OptGain_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Fine pitch - text box
function Pitch_Fine_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Fine))
end
handles.Control.Pitch.Fine = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Pitch_Fine_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Maximum pitch angle - text box
function Pitch_Max_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Max))
end
handles.Control.Pitch.Max = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Pitch_Max_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Minimum pitch angle - text box
function Pitch_Min_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Min))
end
handles.Control.Pitch.Min = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Pitch_Min_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Maximum pitch rate - text box
function Pitch_Maxrate_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Maxrate))
end
handles.Control.Pitch.Maxrate = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Pitch_Maxrate_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Minimum pitch rate - text box
function Pitch_Minrate_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Minrate))
end
handles.Control.Pitch.Minrate = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Pitch_Minrate_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Startup pitch - text box
function Pitch_StartupPitch_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.StartupPitch))
end
handles.Control.Pitch.StartupPitch = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Pitch_StartupPitch_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Startup pitch rate - text box
function Pitch_StartupPitchRate_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.StartupPitchRate))
end
handles.Control.Pitch.StartupPitchRate = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Pitch_StartupPitchRate_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Startup speed for pitch to fine pitch - text box
function Pitch_StartupSpeed_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.StartupSpeed))
end
handles.Control.Pitch.StartupSpeed = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Pitch_StartupSpeed_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Set pitch gain
function Set_Pitch_Gain_Callback(hObject, eventdata, handles)
% Disable window
buttons = findall(handles.ControlDesign, 'Type', 'UIControl');
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'off');
end

% Operation settings
handles.Control = PitchGain(...
    handles.Control);

% Enable window
for i = 1:length(buttons)
    set(buttons(i), 'Enable', 'on');
end

% Update handles structure
guidata(hObject, handles);

%% Open Simulink interface
function Simulink_ForeAft_Callback(hObject, eventdata, handles)
load_system('FAST')
open_system('FAST/Controller/Fore-Aft Tower Control')

%% Open Simulink interface
function Simulink_Torque_Callback(hObject, eventdata, handles)
load_system('FAST')
open_system('FAST/Controller/Torque Control')

%% Open Simulink interface
function Simulink_Pitch_Callback(hObject, eventdata, handles)
load_system('FAST')
open_system('FAST/Controller/Pitch Control')

%% Fore-Aft controller gain - text box
function ForeAft_Gain_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.ForeAft.Gain))
end
handles.Control.ForeAft.Gain = str2double(get(hObject,'String'));
guidata(hObject, handles);
function ForeAft_Gain_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Fore-Aft maximum pitch amplitude contribution - text box
function ForeAft_MaxPitchAmplitude_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.ForeAft.MaxPitchAmplitude))
end
handles.Control.ForeAft.MaxPitchAmplitude = str2double(get(hObject,'String'));
guidata(hObject, handles);
function ForeAft_MaxPitchAmplitude_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Plot the (partial load) torque control curve
function plotbutton_Callback(hObject, eventdata, handles)
    % Prepare vectors for plotting
    OmegaA = handles.Control.Torque.SpeedA;
    OmegaB = handles.Control.Torque.SpeedB;
    OmegaB2 = handles.Control.Torque.SpeedB2;
    OmegaC = handles.Control.Torque.SpeedC;
    OmegaRegion2 = [OmegaB + 0:((OmegaB2-OmegaB)/50):OmegaB2, OmegaB2];
    TorqueRegion2 = handles.Control.Torque.OptGain * (2*pi*OmegaRegion2/60).^2;
    TorqueA = 0;
    TorqueC = handles.Control.Torque.Demanded;
    Omega = [OmegaA, OmegaRegion2, OmegaC];
    Torque = [TorqueA, TorqueRegion2, TorqueC];
    OmegaIdeal = [0:(max(Omega)/50):1.1*max(Omega), 1.1*max(Omega)];
    TorqueIdeal = handles.Control.Torque.OptGain * (2*pi*OmegaIdeal/60).^2;

    % Plot
    Plot = figure();
    set(Plot, 'Name', 'Control curve')
    plot(OmegaIdeal,TorqueIdeal,'--')
    xlim([0 max(OmegaIdeal)])
    ylim('auto')
    set(gca, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'Fontsize', 8);
    xlabel('HSS rotational speed [rad/s]')
    ylabel('HSS torque [Nm]')
    hold on
    plot(Omega,Torque)
    hold off
    pause(0.1)
