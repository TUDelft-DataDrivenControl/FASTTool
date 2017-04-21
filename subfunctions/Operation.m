%% Initialization code
function varargout = Operation(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Operation_OpeningFcn, ...
                   'gui_OutputFcn',  @Operation_OutputFcn, ...
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
function Operation_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.CertificationSettings = varargin{1};
handles.Input = varargin;

% Update input fields
set(handles.OperationMode, 'Value', handles.CertificationSettings.Mode.Type);
set(handles.OperationEvent, 'String', num2str(handles.CertificationSettings.Mode.Actiontime));

if get(handles.OperationMode, 'Value') == 1       % Power production
    set(handles.OperationEvent, 'Visible', 'off')
    set(handles.Stop_label, 'Visible', 'off')
    set(handles.Stop_unit, 'Visible', 'off')
elseif get(handles.OperationMode, 'Value') == 2   % Power production with fault
    set(handles.OperationEvent, 'Visible', 'on')
    set(handles.Stop_label, 'Visible', 'on')
    set(handles.Stop_unit, 'Visible', 'on')
    set(handles.Stop_label, 'String', 'Fault at:')
elseif get(handles.OperationMode, 'Value') == 3   % Startup
    set(handles.OperationEvent, 'Visible', 'on')
    set(handles.Stop_label, 'Visible', 'on')
    set(handles.Stop_unit, 'Visible', 'on')
    set(handles.Stop_label, 'String', 'Start at:')
elseif get(handles.OperationMode, 'Value') == 4   % Normal shutdown
    set(handles.OperationEvent, 'Visible', 'on')
    set(handles.Stop_label, 'Visible', 'on')
    set(handles.Stop_unit, 'Visible', 'on')
    set(handles.Stop_label, 'String', 'Shutdown at:')
elseif get(handles.OperationMode, 'Value') == 5   % Emergency shutdown
    set(handles.OperationEvent, 'Visible', 'on')
    set(handles.Stop_label, 'Visible', 'on')
    set(handles.Stop_unit, 'Visible', 'on')
    set(handles.Stop_label, 'String', 'Shutdown at:')
elseif get(handles.OperationMode, 'Value') == 6   % Idling
    set(handles.OperationEvent, 'Visible', 'off')
    set(handles.Stop_label, 'Visible', 'off')
    set(handles.Stop_unit, 'Visible', 'off')
elseif get(handles.OperationMode, 'Value') == 7   % Parked
    set(handles.OperationEvent, 'Visible', 'off')
    set(handles.Stop_label, 'Visible', 'off')
    set(handles.Stop_unit, 'Visible', 'off')
end

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.Operation);

%% Closing function
function Operation_CloseRequestFcn(hObject, eventdata, handles)
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
uiresume(handles.Operation);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.Operation);

%% Output function
function varargout = Operation_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save

    handles.CertificationSettings.Mode.Type = get(handles.OperationMode,'Value');
    varargout{1} = handles.CertificationSettings;
    
else
    
    varargout = handles.Input;
    
end

% Close figure
delete(hObject)

%% Operation mode
function OperationMode_Callback(hObject, eventdata, handles)
if get(hObject, 'Value') == 1       % Power production
    set(handles.OperationEvent, 'Visible', 'off')
    set(handles.Stop_label, 'Visible', 'off')
    set(handles.Stop_unit, 'Visible', 'off')
elseif get(hObject, 'Value') == 2   % Power production with fault
    set(handles.OperationEvent, 'Visible', 'on')
    set(handles.Stop_label, 'Visible', 'on')
    set(handles.Stop_unit, 'Visible', 'on')
    set(handles.Stop_label, 'String', 'Fault at:')
elseif get(hObject, 'Value') == 3   % Startup
    set(handles.OperationEvent, 'Visible', 'on')
    set(handles.Stop_label, 'Visible', 'on')
    set(handles.Stop_unit, 'Visible', 'on')
    set(handles.Stop_label, 'String', 'Start at:')
elseif get(hObject, 'Value') == 4   % Normal shutdown
    set(handles.OperationEvent, 'Visible', 'on')
    set(handles.Stop_label, 'Visible', 'on')
    set(handles.Stop_unit, 'Visible', 'on')
    set(handles.Stop_label, 'String', 'Shutdown at:')
elseif get(hObject, 'Value') == 5   % Emergency shutdown
    set(handles.OperationEvent, 'Visible', 'on')
    set(handles.Stop_label, 'Visible', 'on')
    set(handles.Stop_unit, 'Visible', 'on')
    set(handles.Stop_label, 'String', 'Shutdown at:')
elseif get(hObject, 'Value') == 6   % Idling
    set(handles.OperationEvent, 'Visible', 'off')
    set(handles.Stop_label, 'Visible', 'off')
    set(handles.Stop_unit, 'Visible', 'off')
elseif get(hObject, 'Value') == 7   % Parked
    set(handles.OperationEvent, 'Visible', 'off')
    set(handles.Stop_label, 'Visible', 'off')
    set(handles.Stop_unit, 'Visible', 'off')
end
function OperationMode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Time stamp for stop - text box
function OperationEvent_Callback(hObject, eventdata, handles)
if str2double(get(hObject,'String')) < 0
    set(hObject, 'String', '0')
elseif isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.CertificationSettings.Mode.Actiontime))
end
handles.CertificationSettings.Mode.Actiontime = str2double(get(hObject,'String'));
guidata(hObject, handles);
function OperationEvent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
