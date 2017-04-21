%% Initialization code
function varargout = DrivetrainDesign(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DrivetrainDesign_OpeningFcn, ...
                   'gui_OutputFcn',  @DrivetrainDesign_OutputFcn, ...
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
function DrivetrainDesign_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Drivetrain = varargin{1};
handles.Input = varargin;

% Set background imagfe
handles.Background = axes('Units', 'Normalized', 'position', [0 0 1 1]);
uistack(handles.Background, 'bottom');
if handles.Drivetrain.Gearbox.Ratio == 1
    img = imread('graphics\drivetrain_directdrive.png');
else
    img = imread('graphics\drivetrain_geared.png');
end
imagesc(img);
set(handles.Background, 'Visible', 'off')

% Update input fields
set(handles.Generator_Efficiency_textbox, 'String', num2str(100*handles.Drivetrain.Generator.Efficiency));
set(handles.HSS_Inertia_textbox, 'String', num2str(handles.Drivetrain.Generator.HSSInertia));
set(handles.Gearbox_Efficiency_textbox, 'String', num2str(100*handles.Drivetrain.Gearbox.Efficiency));
set(handles.Gearbox_Ratio_textbox, 'String', num2str(handles.Drivetrain.Gearbox.Ratio));

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.DrivetrainDesign);

%% Closing function
function DrivetrainDesign_CloseRequestFcn(hObject, eventdata, handles)
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
uiresume(handles.DrivetrainDesign);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.DrivetrainDesign);

%% Output function
function varargout = DrivetrainDesign_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save
    varargout{1} = handles.Drivetrain;
else
    varargout = handles.Input;
end

% Close figure
delete(hObject)

%% Generator efficiency - text box
function Generator_Efficiency_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Drivetrain.Generator.Efficiency))
end
handles.Drivetrain.Generator.Efficiency = 0.01*str2double(get(hObject,'String'));
guidata(hObject, handles);
function Generator_Efficiency_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Inertia around high-speed shaft - text box
function HSS_Inertia_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Drivetrain.Generator.HSSInertia))
end
handles.Drivetrain.Generator.HSSInertia = str2double(get(hObject,'String'));
guidata(hObject, handles);
function HSS_Inertia_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Gearbox efficiency - text box
function Gearbox_Efficiency_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Drivetrain.Gearbox.Efficiency))
end
handles.Drivetrain.Gearbox.Efficiency = 0.01*str2double(get(hObject,'String'));
guidata(hObject, handles);
function Gearbox_Efficiency_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Gearbox ratio - text box
function Gearbox_Ratio_textbox_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Drivetrain.Gearbox.Ratio))
end
handles.Drivetrain.Gearbox.Ratio = str2double(get(hObject,'String'));

% Update background image
axes(handles.Background);
if handles.Drivetrain.Gearbox.Ratio == 1
    img = imread('graphics\drivetrain_directdrive.png');
else
    img = imread('graphics\drivetrain_geared.png');
end
imagesc(img);
set(handles.Background, 'Visible', 'off')
guidata(hObject, handles);
function Gearbox_Ratio_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
