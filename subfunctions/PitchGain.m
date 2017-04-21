%% Initialization code
function varargout = PitchGain(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PitchGain_OpeningFcn, ...
                   'gui_OutputFcn',  @PitchGain_OutputFcn, ...
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
function PitchGain_OpeningFcn(hObject, eventdata, handles, varargin)

% Get input
handles.Control = varargin{1};
handles.Input = varargin;

% Update input fields
handles.TableSize = length(handles.Control.Pitch.KpGS);
for i = 1:handles.TableSize
    if isnan(handles.Control.Pitch.ScheduledPitchAngles(i))
        TableData(i,1) = num2cell([]);
    else
        TableData(i,1) = num2cell(handles.Control.Pitch.ScheduledPitchAngles(i));
    end
    if isnan(handles.Control.Pitch.KpGS(i))
        TableData(i,2) = num2cell([]);
    else
        TableData(i,2) = num2cell(handles.Control.Pitch.KpGS(i));
    end
    if isnan(handles.Control.Pitch.KiGS(i))
        TableData(i,3) = num2cell([]);
    else
        TableData(i,3) = num2cell(handles.Control.Pitch.KiGS(i));
    end
end
set(handles.Constant_Ki_textbox, 'String', num2str(handles.Control.Pitch.Ki));
set(handles.Constant_Kp_textbox, 'String', num2str(handles.Control.Pitch.Kp));
set(handles.TableSize_textbox, 'String', length(handles.Control.Pitch.KpGS));
set(handles.TableSize_slider, 'Value', length(handles.Control.Pitch.KpGS));
set(handles.Table, 'Data', TableData);

% Check for gain scheduling
if handles.Control.Pitch.Scheduled
    set(handles.ConstantGain, 'Value', 0);
    set(handles.GainScheduled, 'Value', 1);
    set(handles.Constant_Ki_textbox, 'Enable', 'off')
    set(handles.Constant_Kp_textbox, 'Enable', 'off')
    set(handles.Table, 'Enable', 'on')
    set(handles.EditCells_checkbox, 'Enable', 'on')
    set(handles.EditCells_checkbox, 'Value', 0)
    set(handles.TableSize_textbox, 'Enable', 'off')
    set(handles.TableSize_slider, 'Enable', 'off')
else
    set(handles.ConstantGain, 'Value', 1);
    set(handles.GainScheduled, 'Value', 0);
    set(handles.Constant_Ki_textbox, 'Enable', 'on')
    set(handles.Constant_Kp_textbox, 'Enable', 'on')
    set(handles.Table, 'Enable', 'off')
    set(handles.EditCells_checkbox, 'Enable', 'off')
    set(handles.EditCells_checkbox, 'Value', 0)
    set(handles.TableSize_textbox, 'Enable', 'off')
    set(handles.TableSize_slider, 'Enable', 'off')
end

% Update handles structure
guidata(hObject, handles);

% Halt window
uiwait(handles.PitchGain);

%% Closing function
function PitchGain_CloseRequestFcn(hObject, eventdata, handles)
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
uiresume(handles.PitchGain);

%% Cancel button
function Cancel_Callback(hObject, eventdata, handles)
handles.Save = false;
guidata(hObject, handles);
uiresume(handles.PitchGain);

%% Output function
function varargout = PitchGain_OutputFcn(hObject, eventdata, handles) 

% Set output
if handles.Save

    % Get geometry from table
    Table = get(handles.Table, 'Data');
    
    % Empty gain vector
    ScheduledPitchAngles = nan(handles.TableSize,1);
    KiGS = nan(handles.TableSize,1);
    KpGS = nan(handles.TableSize,1);

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
            ScheduledPitchAngles(i) = 0;
        else
            ScheduledPitchAngles(i) = Table{i,1};
        end
        if invalid(i,2)
            KpGS(i) = 0;
        else
            KpGS(i) = Table{i,2};
        end
        if invalid(i,3)
            KiGS(i) = 0;
        else
            KiGS(i) = Table{i,3};
        end
    end

    
    handles.Control.Pitch.Ki = str2double(get(handles.Constant_Ki_textbox,'String'));
    handles.Control.Pitch.Kp = str2double(get(handles.Constant_Kp_textbox,'String'));
    handles.Control.Pitch.ScheduledPitchAngles = ScheduledPitchAngles(~isnan(ScheduledPitchAngles));
    handles.Control.Pitch.KpGS = KpGS(~isnan(KpGS));
    handles.Control.Pitch.KiGS = KiGS(~isnan(KiGS));
    handles.Control.Pitch.Scheduled = get(handles.GainScheduled,'Value');
    varargout{1} = handles.Control;
else
    varargout = handles.Input;
end

% Close figure
delete(hObject)

%% Edit cells - checkbox
function EditCells_checkbox_Callback(hObject, eventdata, handles)
if get(hObject, 'Value') == 1
    set(handles.Table, 'ColumnEditable', [true true true]);
    set(handles.TableSize_textbox, 'Enable', 'on')
    set(handles.TableSize_slider, 'Enable', 'on')
else
    set(handles.Table, 'ColumnEditable', [false false false]);
    set(handles.TableSize_textbox, 'Enable', 'off')
    set(handles.TableSize_slider, 'Enable', 'off')
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
    set(handles.Constant_Ki_textbox, 'Enable', 'on')
    set(handles.Constant_Kp_textbox, 'Enable', 'on')
    set(handles.Table, 'Enable', 'off')
    set(handles.EditCells_checkbox, 'Enable', 'off')
end
    
%% Gain scheduling - radio button
function GainScheduled_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.Constant_Ki_textbox, 'Enable', 'off')
    set(handles.Constant_Kp_textbox, 'Enable', 'off')
    set(handles.Table, 'Enable', 'on')
    set(handles.EditCells_checkbox, 'Enable', 'on')
end

%% Constant integration constant - text box
function Constant_Ki_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Ki(end)))
end
handles.Control.Pitch.Ki(end) = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Ki_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Constant gain - text box
function Constant_Kp_textbox_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(handles.Control.Pitch.Kp(end)))
end
handles.Control.Pitch.Kp(end) = str2double(get(hObject,'String'));
guidata(hObject, handles);
function Constant_Kp_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Table size - text box
function TableSize_textbox_Callback(hObject, eventdata, handles)

rows = ceil(str2double(get(hObject, 'String')));
Table = get(handles.Table, 'Data');
if rows < 2
    rows = 2;
elseif rows > get(handles.TableSize_slider, 'Max')
    rows = get(handles.TableSize_slider, 'Max');
elseif isnan(rows)
    rows = size(Table,1);
end
if rows < size(Table,1)
    Table = Table(1:rows,:);
    set(handles.Table, 'Data', Table);
elseif rows > size(Table,1)
    Table{rows,size(Table,2)} = [];
    set(handles.Table, 'Data', Table);
end
set(hObject, 'String', int2str(rows));
set(handles.TableSize_slider, 'Value', rows);
function TableSize_textbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% Table size - slider
function TableSize_slider_Callback(hObject, eventdata, handles)
rows = get(hObject, 'Value');
Table = get(handles.Table, 'Data');
if rows < 2
    set(hObject, 'Value', 2);
else
    if rows < size(Table,1)
        Table = Table(1:rows,:);
        set(handles.Table, 'Data', Table);
    elseif rows > size(Table,1)
        Table{rows,size(Table,2)} = [];
        set(handles.Table, 'Data', Table);
    end
    set(handles.TableSize_textbox, 'String', int2str(rows));
end
function TableSize_slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
