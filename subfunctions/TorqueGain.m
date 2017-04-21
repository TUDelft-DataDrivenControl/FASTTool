%% Initialization code
function varargout = TorqueGain(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TorqueGain_OpeningFcn, ...
                   'gui_OutputFcn',  @TorqueGain_OutputFcn, ...
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
function TorqueGain_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Kp = varargin{1};
handles.Ti = varargin{2};
handles.Scheduled = varargin{3};
handles.Input = varargin(1:3);

% Update input fields
handles.TableSize = length(handles.Kp) - 1;
for i = 1:handles.TableSize
    if isnan(handles.Ti(i))
        TableData(i,1) = num2cell([]);
    else
        TableData(i,1) = num2cell(handles.Ti(i));
    end
    if isnan(handles.Kp(i))
        TableData(i,2) = num2cell([]);
    else
        TableData(i,2) = num2cell(handles.Kp(i));
    end
end
set(handles.Constant_Ti_textbox, 'String', num2str(handles.Ti(end)));
set(handles.Constant_Kp_textbox, 'String', num2str(handles.Kp(end)));
set(handles.Table, 'Data', TableData);

% Check for gain scheduling
if handles.Scheduled
    set(handles.ConstantGain, 'Value', 0);
    set(handles.GainScheduled, 'Value', 1);
    set(handles.Constant_Ti_textbox, 'Enable', 'off')
    set(handles.Constant_Kp_textbox, 'Enable', 'off')
    set(handles.Table, 'Enable', 'on')
    set(handles.EditCells_checkbox, 'Enable', 'on')
else
    set(handles.ConstantGain, 'Value', 1);
    set(handles.GainScheduled, 'Value', 0);
    set(handles.Constant_Ti_textbox, 'Enable', 'on')
    set(handles.Constant_Kp_textbox, 'Enable', 'on')
    set(handles.Table, 'Enable', 'off')
    set(handles.EditCells_checkbox, 'Enable', 'off')
end

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.SetPIGain);

%% Closing function
function SetPIGain_CloseRequestFcn(hObject, eventdata, handles)
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
uiresume(handles.SetPIGain);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.SetPIGain);

%% Output function
function varargout = TorqueGain_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save

    % Get geometry from table
    Table = get(handles.Table, 'Data');
    
    % Empty gain vector
    Ti = nan(handles.TableSize + 1,1);
    Kp = nan(handles.TableSize + 1,1);

    % Find invalid cells
    for i = 1:size(Table,1)
        for j = 1:size(Table,2)
            invalid(i,j) = ...
                isempty(cell2mat(Table(i,j))) + ...
                sum(isnan(cell2mat(Table(i,j))));
        end
    end

    % Extract geometry from table
    for i = 1:size(Table,1)
        if invalid(i,1)
            Ti(i) = 0;
        else
            Ti(i) = Table{i,1};
        end
        if invalid(i,2)
            Kp(i) = 0;
        else
            Kp(i) = Table{i,2};
        end
    end
    Ti(end) = handles.Ti(end);
    Kp(end) = handles.Kp(end);

    varargout{1} = Kp;
    varargout{2} = Ti;
    varargout{3} = get(handles.GainScheduled, 'Value');
    
else
    
    varargout = handles.Input;
    
end

% Close figure
delete(hObject)

%% Edit cells - checkbox
function EditCells_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value') == 1
    set(handles.Table, 'ColumnEditable', [true true]);
else
    set(handles.Table, 'ColumnEditable', [false false]);
end

%% Table cell selection
function Table_CellSelectionCallback(hObject, eventdata, handles)
handles.Selection = eventdata.Indices;
guidata(hObject, handles); 

%% Table keyboard functions
function Table_KeyPressFcn(hObject, eventdata, handles)

% Paste functionality
if strcmpi(char(eventdata.Modifier),'control') && strcmp(eventdata.Key, 'v')
    
    % Get and reshape clipboard data
    Paste = clipboard('paste');
    rows = numel(strsplit(strtrim(Paste), '\n'));
    Paste = strsplit(strtrim(Paste));
    columns = numel(Paste)/rows;
    Paste = reshape(Paste,columns,rows)';
    
    % Convert numeric cells
    for i = 1:rows
        for j = 1:columns
        	Paste{i,j} = str2double(Paste{i,j});
        end
    end

    % Target cells
    Table = get(hObject, 'Data');
    Selection = handles.Selection;
    i1 = Selection(1,1);
    i2 = Selection(1,1) + (rows-1);
    j1 = Selection(1,2);
    j2 = Selection(1,2) + (columns-1);
    if i2 > size(Table,1)
        if i2 > handles.TableSize
            i2 = handles.TableSize;
        end
    end
    if j2 > size(Table,2)
        j2 = size(Table,2);
    end
    
    % Constrain data types
    Table(i1:i2,j1:j2) = Paste(1:(1+i2-i1),1:(1+j2-j1));
    for i = i1:i2
        for j = 1:size(Table,2)
            if ~isnumeric(Table{i,j})
                Table{i,j} = NaN;
            end
        end
    end

    % Update table
    set(hObject, 'Data', Table);

end

%% Constant gain - radio button
function ConstantGain_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.Constant_Ti_textbox, 'Enable', 'on')
    set(handles.Constant_Kp_textbox, 'Enable', 'on')
    set(handles.Table, 'Enable', 'off')
    set(handles.EditCells_checkbox, 'Enable', 'off')
end
    
%% Gain scheduling - radio button
function GainScheduled_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.Constant_Ti_textbox, 'Enable', 'off')
    set(handles.Constant_Kp_textbox, 'Enable', 'off')
    set(handles.Table, 'Enable', 'on')
    set(handles.EditCells_checkbox, 'Enable', 'on')
end

%% Constant integration constant - text box
function Constant_Ti_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Ti(end)))
end
handles.Ti(end) = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Ti_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Constant gain - text box
function Constant_Kp_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Kp(end)))
end
handles.Kp(end) = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Kp_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
