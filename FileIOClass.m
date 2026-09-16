classdef FileIOClass < handle
% The FileIOClass takes care of the interactions that FASTTool has with
% data files.
%
% The FileIOClass manages and is aware of the formats of the FASTTool data
% files. This is used to check that the appropriate content exists in files
% that are opened, to update project files to what is used in the current
% version of FASTTool, if needed, and to save the appropriate parameters to
% project and results files.
% 
% For other data-file interaction (for FAST, bModes, Aerodyn, etc.),
% the class operates as agnostic as possible. Where possible, it passes on
% information between the program and the files as plain text. For some
% data files knowledge about the file structure is used for easier
% implementation of reading and writing functions.
%
% The FileIOClass creates a temporary directory in which files are saved
% that are used by BModes and FAST for modal analysis, linearisation and
% simulation. This directory is deleted when the FileIOClass object is
% deleted.
%
% The functions in this class have no or limited input and output
% variables, because the relevant data are passed on via the data
% structures provided by the handles that are exchanged between this class,
% WindTurbineClass and the FASTTool class. Some input variables are used to
% create reusability of the functions, by specifying which elements in the
% data structure need to be applied. Some output variables are used for
% temporary data, or for status information.
%
% The FileIOClass is a handle class, because the FileIO object that is
% created from this class definition in FASTTool.StartupFcn is passed on
% to and used by many components of the app. As a handle object, this is
% done by reference instead of by making copies of the object.

    properties
        % Handles to main app and WindTurbine
        MainApp
        WindTurbine

        % File information relevant to the open project
        projectFileName
        projectPathName
        projectInfo % When saving: Version of FASTTool used and creation time

        % File information relevant to simulation results
        simulationFile
        simulationPath
        dataDescription % List of parameter names and explanations

        % Header and footer for table exports
        % If the number of lines in the table header is changed, older
        % exported tables will not be imported correctly.
        headerTableExport = ["> Table with blade data, exported from FASTTool"; ...
            "> Can be edited, to be imported back into FASTTool" ; ...
            "> Do not swap, remove or add intermediate columns"; ...
            "> Do not remove this text, the table header or the cell with 'end-of-table' below the table"; ...
            "> Rows may be added or removed and cells may contain calculations"; ...
            "> Text and calculations may be added outside the dedicated area"; ...
            " "];
        footerTableExport = "end-of-table";
        linesTableHeader = 8; % Number of lines in the table header

        % Information about the folder used for temporary files
        % Temporary files are files that are written by FASTTool for
        % internal use, without instructions from the user. Examples are
        % the input files for bModes and for FAST, and output files of
        % bModes.
        tempDirectory = 'FASTtmp';
        pathTempDirectory;

        % Parameters to save to and load from files
        % The first variable of each pair provides the name of the variable
        % stored in the related file. Except for the variable Info, the
        % second variable of each pair provides a description of it.
        % These lists are used to store the relevant parameters from the
        % data structures and to check the validity of files to be loaded.
        % For the variables with a description, the descriptions
        % are stored in the variable Info and then also saved in the file.
        ProjectParameters = {'Info' ''; 'Blade' ''; 'Airfoil' ''; 'Tower' ''; 'Nacelle' ''; 'Drivetrain' ''; 'Control' ''; 'Appearance' ''; 'CertificationSettings' ''};
        LinearModelParameters = {'Info' ''; 'Lin' ''; 'sysm' ''};
        RotorParameters = {'Info' ''; ...
            'tsr' 'Tip speed ratio corresponding to tsr_index [-]'; ...
            'pitch' 'Collective blade pitch angle corresponding to pitch_index [deg]'; ...
            'cP' 'Power coefficient f(pitch_index, tsr_index) [-]'; ...
            'cT' 'Thrust coefficient f(pitch_index, tsr_index) [-]'; ...
            'cQ' 'Torque coefficient f(pitch_index, tsr_index) [-]'};
        SteadyOpParameters = {'Info' ''; ...
            'ratedPower' 'Rated electrical output power [W]'; ...
            'windSpeed' 'Wind speed corresponding to windSpeed_index [m/s]'; ...
            'electricalPower' 'Electrical output power f(windSpeed_index) [W]'; ...
            'thrust' 'Thrust force on the rotor f(windSpeed_index) [N]'; ...
            'torque' 'Torque in main (low-speed) shaft f(windSpeed_index) [Nm]'; ...
            'rotorSpeed' 'Rotational speed of the rotor f(windSpeed_index) [rpm]'; ...
            'pitch' 'Collective blade pitch angle f(windSpeed_index) [deg]'; ...
            'tsr' 'Tip speed ratio f(windSpeed_index) [-]'; ...
            'cP' 'Power coefficient f(windSpeed_index) [-]'; ...
            'cT' 'Thrust coefficient f(windSpeed_index) [-]'; ...
            'cQ' 'Torque coefficient f(windSpeed_index) [-]'};
        ModalParameters = {'Info' ''; ...
            'cutInRotorSpeed' 'Lower bound of operational rotor speed range [rpm]'; ...
            'ratedRotorSpeed' 'Upper bound of operational rotor speed range [rpm]'; ...
            'rotorSpeed' 'Rotational speed of the rotor corresponding to rotorSpeed_index [rpm]'; ...
            'Tower_ForeAft1_freq' 'Natural frequency of 1st tower fore-aft mode f(rotorSpeed_index) [Hz]'; ...
            'Tower_ForeAft2_freq' 'Natural frequency of 2nd tower fore-aft mode f(rotorSpeed_index) [Hz]'; ...
            'Tower_SideSide1_freq' 'Natural frequency of 1st tower side-to-side mode f(rotorSpeed_index) [Hz]'; ...
            'Tower_SideSide2_freq' 'Natural frequency of 2nd tower side-to-side mode f(rotorSpeed_index) [Hz]'; ...
            'Blade_Flap1_freq' 'Natural frequency of 1st blade flapwise mode f(rotorSpeed_index) [Hz]'; ...
            'Blade_Flap2_freq' 'Natural frequency of 2nd blade flapwise mode f(rotorSpeed_index) [Hz]'; ...
            'Blade_Edge1_freq' 'Natural frequency of 1st blade edgewise mode f(rotorSpeed_index) [Hz]'; ...
            'Blade_Edge2_freq' 'Natural frequency of 2nd blade edgewise mode f(rotorSpeed_index) [Hz]'};
        ControlParameters = {'Info' ''; ...
            'responseType' 'Types of frequency response corresponding to responseType_index - LPF = low-pass filter | PI = PI controller | Notch = notch filters | Plant = nominal system | LoopGain = loop gain (open loop)'; ...
            'pitch' 'Collective blade pitch angle corresponding to pitch_index (full-load conditions only) [rad]'; ...
            'windSpeed' 'Wind speed corresponding to windSpeed_index [m/s]'; ...
            'omega' 'Frequency corresponding to omega_index [rad/s]'; ...
            'Magnitude' 'Magnitude of the frequency response f{responseType_index}(omega_index,pitch_index) [dB]'; ...
            'Phase' 'Phase of the frequency response f{responseType_index}(omega_index,pitch_index) [deg]'; ...
            'GM' 'Gain margin f(pitch_index) [dB]'; ...
            'GMFreq' 'Frequency where gain margin is established (= first occurrence of phase reaching -180 deg) f(pitch_index) [Hz]'; ...
            'PM' 'Phase margin f(pitch_index) [deg]'; ...
            'PMFreq' 'Frequency where phase margin is established (= gain cross-over frequency at 0 dB gain) f(pitch_index) [Hz]'; ...
            'IsGMFreqWithinLimit' 'Boolean indicating whether GMFreq falls in the range of omega f(pitch_index)'; ...
            'IsPMFreqWithinLimit' 'Boolean indicating whether PMFreq falls in the range of omega  f(pitch_index)'};
        SimulationParameters = {'Info' ''; 'Time' ''; 'Azimuth' ''}; % Incomplete list, only used to check validity of simulation results file
    end

    methods

        % setHandles stores references to the WindTurbine and MainApp
        % handle objects in local properties for later use.
        %
        %   Input arguments
        %     WindTurbine - Handle to object of WindTurbineClass
        %     MainApp - Handle to main app (FASTTool)
        %
        %   Output
        %     [-]
        %
        %   Called by
        %     FASTTool.StartupFcn
        function setHandles(obj, WindTurbine, MainApp)
            obj.WindTurbine = WindTurbine;
            obj.MainApp = MainApp;
        end

        % openProject opens the specified project file and stores its
        % content in a variable if it is an appropriate file.
        %
        %   Input arguments
        %     fileName - Name of the project file to open
        %     pathName - Path to the project file to open
        %     mode - Mode of use of the project in the user interface:
        %     project (for editing) or compare (to compare two project
        %     files)
        %
        %   Output
        %     Data - All parameters loaded from the specified file if
        %     succesfull. Otherwise empty.
        %
        %   Behaviour
        %     The content of the file is checked, to ensure that it is a
        %     file with project data. Data from project files created by
        %     earlier versions of FASTTool are updated to the needs for the
        %     current version. When the mode of use is to edit the project,
        %     the filename, path and version are stored for later saves of
        %     the project.
        %
        %   Called by
        %     FASTTool.OpenProject
        %     MergeDialogue.updateData
        function [Data] = openProject(obj, fileName, pathName, mode)
            if ~properMatFile(obj, fileName, pathName, 'project')
                Data = [];
                return;
            end
            
            Data = load([pathName,fileName]);
            Data = parseData(obj, Data);

            switch mode
                case 'compare'
                case 'project'
                    obj.projectFileName = fileName;
                    obj.projectPathName = pathName;
                    obj.projectInfo.version = obj.MainApp.version;
            end
        end

        % saveProject saves the relevant data from the open project to the
        % open project file.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The current time is determined and added to the project
        %     information. Properties from WindTurbine that specify the
        %     wind turbine are saved to the open project file, overwriting
        %     previous content if the file already exists.
        %
        %   Called by
        %     FASTTool.SaveProject
        %     FASTTool.sureToClose
        %     FileIOClass.saveProjectAs
        function saveProject(obj)
            % All relevant properties from WindTurbine are copied to local
            % parameters, such that these local parameters contain all and
            % only relevant information to be saved. The current time is
            % added to the information about the project file.
            Info = obj.projectInfo;
            Info.creationTime = sprintf('%s', datetime);
            Blade = obj.WindTurbine.Blade;
            Airfoil = obj.WindTurbine.Airfoil;
            Tower = obj.WindTurbine.Tower;
            Nacelle = obj.WindTurbine.Nacelle;
            Drivetrain = obj.WindTurbine.Drivetrain;
            Control = obj.WindTurbine.Control;
            Appearance = obj.WindTurbine.Appearance;
            CertificationSettings = obj.WindTurbine.CertificationSettings;

            % The local parameters are saved to the open project file,
            % ensuring that the file contains all relevant project
            % information and nothing else.
            save([obj.projectPathName,obj.projectFileName], ...
                'Info', ...
                'Blade', ...
                'Airfoil', ...
                'Tower', ...
                'Nacelle', ...
                'Drivetrain', ...
                'Control', ...
                'Appearance', ...
                'CertificationSettings');
        end

        % saveProjectAs saves the relevant data from the open project to
        % a specified new project file.
        %
        %   Input arguments
        %     fileName - Name of the project file to save to
        %     pathName - Path to the project file to save to
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The file information for the open project is updated with the
        %     new project file specifications. After this, the function to
        %     save the open project to file is called, which will save it
        %     to the newly set project file. No check is performed on the
        %     specified file, since that should already be tested before
        %     calling saveProjectAs.
        %
        %   Called by
        %     FASTTool.SaveProjectAs
        %     FASTTool.sureToClose
        function saveProjectAs(obj,  fileName, pathName)
            obj.projectFileName = fileName;
            obj.projectPathName = pathName;
            saveProject(obj);
        end

        % saveMergedProject saves the provided project data to a specified
        % project file.
        %
        %   Input arguments
        %     fileName - Name of the project file to save to
        %     pathName - Path to the project file to save to
        %     Data - Data structure with turbine specifications containing
        %     all and no more parameters as needed for a project file
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The current time and the specified version of FASTTool are
        %     copied to the file information parameter. After this, the
        %     project data is saved in the specified file, overwriting
        %     previous content if the file already exists. No check is
        %     performed on the specified file, since that should already be
        %     tested before calling saveMergedProject. No check is
        %     performed on the data, since it is created by merging two
        %     data sets that where opened with the current version of
        %     FASTTool and have therefore been checked and updated before.
        %
        %   Called by
        %     MergeDialogue.SaveasButtonPushed
        function saveMergedProject(obj,  fileName, pathName, Data)
            Data.Info.version = obj.MainApp.version;
            Data.Info.creationTime = sprintf('%s', datetime);

            save([pathName, fileName], "-struct", "Data", ...
                'Info', ...
                'Blade', ...
                'Airfoil', ...
                'Tower', ...
                'Nacelle', ...
                'Drivetrain', ...
                'Control', ...
                'Appearance', ...
                'CertificationSettings');
        end

        % getRestoreData stores the content of the open project in a
        % variable.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     Data - All parameters loaded from the open project file
        %
        %   Behaviour
        %     Data from the open project is loaded. Data from project files
        %     created by earlier versions of FASTTool are updated to the
        %     needs for the current version. The data is not checked, since
        %     it is from the currently opened project. For the same reason
        %     file information for the open project is also not updated.
        %     Parts of the data will be used to restore data in the
        %     interface that have been edited by the user interface. The
        %     app component that calls this function determines which data
        %     is restored. Data is not directly and entirely copied to
        %     WindTurbine in getRestoreData, because not all data needs
        %     to be restored. Edits done in other app components remain.
        %
        %   Called by
        %     FASTTool.getRestoredata
        function [Data] = getRestoreData(obj)
            Data = load([obj.projectPathName,obj.projectFileName]);
            Data = parseData(obj, Data);
        end

        % writeExportTable writes the header lines, table content and
        % footer line to a specified Excel file
        %
        %   Input arguments
        %     pathName - Path to the Excel file to save to
        %     fileName - Name of the Excel file to save to
        %     DataTable - Table with data from the user interface, to be
        %     editable in Excel
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The header lines specified in the class properties are
        %     written to the specified file, overwriting previous content
        %     if the file already exists. The provided table data and the
        %     footer line from the class properties are appended to the
        %     file.
        %     No check is performed on the specified file, since that
        %     should already be tested before calling writeExportTable. The
        %     file should be forced to have have the extention .xlsx, to
        %     ensure that writeExportTable exports the data to Excel.
        %
        %   Called by
        %     FASTTool.exportTable
        function writeExportTable(obj, pathName, fileName, DataTable)
            writetable(table(obj.headerTableExport), [pathName, fileName], WriteMode = "overwritesheet", AutoFitWidth=false, WriteVariableNames=false);
            writetable(DataTable, [pathName, fileName], WriteMode= "append", WriteVariableNames=true);
            writetable(table(obj.footerTableExport), [pathName, fileName], WriteMode = "append", AutoFitWidth=false, WriteVariableNames=false);
        end

        % readImportTable 
        %
        %   Input arguments
        %     pathName - Path to the Excel file to import from
        %     fileName - Name of the Excel file to import from
        %     type - Type of table for which data needs to be imported
        %
        %   Output
        %     DataTable - Table with data from the Excel file, to be used
        %     in the user interface
        %     message - Error message
        %
        %   Behaviour
        %     The first cell after the header lines is interpreted as the
        %     variable name for the first table column. If it doesn't
        %     correspond with the expected name for the specified type, the
        %     file is considered invalid.
        %     For a valid file, the remainder of the first column is
        %     searched for the footer. If it cannot be found, the file is
        %     flagged to have no end-of-table marker.
        %     For invalid files, Data is empty and an error message is
        %     returned. For valid files, the table is copied into Data and
        %     the message is empty.
        %
        %   Called by 
        %     FASTTool.importTable
        function [DataTable, message] = readImportTable(obj, pathName, fileName, type)
            range = sprintf("A%i:A%i", obj.linesTableHeader, obj.linesTableHeader);
            xVariableName = readcell([pathName, fileName], Range = range);
            xVariableName = xVariableName{1};
            if ~properTableFile(obj, xVariableName, type)
                DataTable = [];
                message = ['<' fileName '> is not a valid table file'];
                return
            end

            columnA = readcell([pathName, fileName], Range = "A:A");
            endTable = find(strcmp(columnA, obj.footerTableExport)) - obj.linesTableHeader - 1;
            if isempty(endTable)
                DataTable = [];
                message = ['<' fileName '> has no valid end-of-table marker'];
            else
                T = readtable([pathName, fileName], NumHeaderLines = obj.linesTableHeader, ReadVariableNames=false);
                DataTable = T(1:endTable,:);
                message = '';
            end
        end        
        
        %  prepareTempDirectory creates a folder for temporary files.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     status - Coded integer indication of the success of the
        %     preparation of the temporary folder.
        %     overFlowLength - Number of characters exceeding allowed path
        %     length
        %     tempDirName - Name of temporary directory
        %
        %   Behaviour
        %     This function ensures that a temporary folder is available
        %     for temporary files to be created during a modal analysis,
        %     linearisation or simulation, or returns a status that
        %     indicates why this is not possible. A created temporary
        %     folder is removed when closing the app by the function
        %     removeTempDirectory.
        %
        %   Called by
        %     FASTTool.prepareFAST
        function [status, overflowLength, tempDirName] = prepareTempDirectory(obj)
            % coding of status:
            %   0 = Directory could not be created
            %   1 = Directory exists or was successfully created
            %   2 = Path of current working directory is too long
            %   3 = Directory contains a file that cannot be deleted
            overflowLength = 0;
            tempDirName = obj.pathTempDirectory;

            if isempty(obj.pathTempDirectory)
                % For the modal analysis, linearisation and simulation, BModes
                % is run. BModes writes to the file BModes.out, in the
                % temporary directory. The path and filename that BModes uses
                % for this output file is truncated at 100 characters
                % (excluding the dot and extension). This means that the
                % pathname of the temporary directory plus a file-separation
                % characters plus 'BModes'may not contain more than 100
                % characters.
                lengthTempDir = length([pwd filesep obj.tempDirectory]);
                overflowLength = lengthTempDir - (100 - 7);
                if overflowLength > 0
                    status = 2;
                    return;
                else
                    obj.pathTempDirectory = [pwd filesep obj.tempDirectory];
                    tempDirName = obj.pathTempDirectory;
                end
            end
            if exist(obj.pathTempDirectory, 'dir')
                state = rmdir(obj.pathTempDirectory, 's');
                if state == 0
                    status = 3;
                else
                    status = mkdir(obj.pathTempDirectory);
                    if status == 0
                        obj.pathTempDirectory = '';
                    end
                end
            else
                status = mkdir(obj.pathTempDirectory);
                if status == 0
                    obj.pathTempDirectory = '';
                end
            end
        end

        % removeTempDirectory removes the folder for temporary files
        % created by prepareTempDirectory. It is called when closing the
        % app.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     status - Coded integer indication of the success of the
        %     preparation of the temporary folder.
        %     tempDirName - Name of temporary directory
        %
        %   Behaviour
        %     The function checks whether a temporary file has been created
        %     by prepareTempDirectory. Then, a check of its existence is
        %     performed, since the user could have removed it during the
        %     execution of FASTTool. If both checks are past, the function
        %     tries to remove the folder. The returned status is 1, unless
        %     the function failed to remove the folder. This can happen
        %     when a file inside the folder is still open. There have been
        %     occurrences when files where still open and not released by
        %     Matlab after a Simulink run.
        %
        %   Called by
        %     FASTTool.FASTToolCloseRequest
        function [success, tempDirName] = removeTempDirectory(obj)
            success = 1;
            tempDirName = obj.pathTempDirectory;

            if isempty(obj.pathTempDirectory)
                return;
            end
            if exist(obj.pathTempDirectory, 'dir')
                % It seems that files inside the temporary directory are
                % sometimes not yet released for deletion after Simulink
                % has finished and then the file and directory cannot be
                % removed. Therefore, this is captured and returned with
                % the status parameter.
                success = rmdir(obj.pathTempDirectory, 's');
            end
        end

        % getWriteLocation returns a string with the full path to the
        % folder for temporary files added to the provided filename.
        %
        %   Input arguments
        %     filename - Filename to add to the full path string
        %
        %   Output
        %     fileLocation - String with full path and filename
        %
        %   Behaviour
        %     The strings for the path and the filename are concatenated,
        %     using the filesep keyword to ensure the appropriate file
        %     separation character for the used platform.
        %
        %   Called by
        %     WindTurbineClass.simulate
        function [fileLocation] = getWriteLocation(obj, filename)
            fileLocation = [obj.pathTempDirectory filesep filename];
        end

        % writeDataFile writes the provided text to the specified file in
        % the folder for temporary files. 
        %
        %   Input arguments
        %     filename - Name of the file to add to the folder for
        %     temporary files and to which to write
        %     text - String to write to the temporary file
        %
        %   Output
        %     fileLocation - String with full path and filename
        %
        %   Behaviour
        %     The strings for the path to the folder with temporary files
        %     and the specified filename are concatenated, using the
        %     filesep keyword to ensure the appropriate file separation
        %     character for the used platform.
        %     The resulting file is opened and the text string is written
        %     to it. If the file already exists, its content is
        %     overwritten.
        %     The filename is not checked, since it is hard coded in the
        %     function that calls writeDataFile.
        %
        %   Called by
        %     WindTurbineClass.performModalAnalysis
        %     WindTurbineClass.linearise
        %     WindTurbineClass.simulate
        function [fileLocation] = writeDataFile(obj, filename, text)
            fileID = fopen([obj.pathTempDirectory filesep filename], 'wt');
            fprintf(fileID, '%s', text);
            fclose(fileID);
            fileLocation = [obj.pathTempDirectory filesep filename];
        end

        % readBModesOutput reads the content of the bModes output file,
        % without interpreting it.
        %
        %   Input arguments
        %     filename - Name of the bModes output file from which to read
        %     the data
        %
        %   Output
        %     data - String array with content from the bModes output file
        %
        %   Behaviour
        %     The strings for the path to the folder with temporary files
        %     and the specified filename are concatenated, using the
        %     filesep keyword to ensure the appropriate file separation
        %     character for the used platform.
        %     The resulting file is opened and read, using a format
        %     specifier to store the tab-delimited data in a string array.
        %     The header information is skipped, when reading the file.
        %     The filename is not checked, since it is hard coded in the
        %     function that calls writeDataFile.
        %
        %   Called by
        %     WindTurbineClass.performModalAnalysis
        %     WindTurbineClass.linearise
        %     WindTurbineClass.simulate
        function data = readBModesOutput(obj, filename)
            fileID = fopen([obj.pathTempDirectory filesep filename], 'r');
            textscan(fileID, '%[^\n\r]', 5, 'ReturnOnError', false);
            data = textscan(fileID, '%s%s%s%s%s%s%[^\n\r]', 'Delimiter', '\t', 'ReturnOnError', false);
            fclose(fileID);
        end

        % readFASTSimOutput reads the content of the FAST output file after
        % simulation, without interpreting it.
        %
        %   Input arguments
        %     filename - Name of the FAST output file from which to read
        %     the data
        %
        %   Output
        %     data - String array with content from the FAST output file
        %
        %   Behaviour
        %     The strings for the path to the folder with temporary files
        %     and the specified filename are concatenated, using the
        %     filesep keyword to ensure the appropriate file separation
        %     character for the used platform.
        %     The resulting file is opened and read, using a format
        %     specifier to store the tab-delimited data in a string array.
        %     The header information is skipped, when reading the file.
        %     The filename is not checked, since it is hard coded in the
        %     function that calls writeDataFile.
        %
        %   Called by
        %     WindTurbineClass.simulate
        function data = readFASTSimOutput(obj, filename)
            fileID = fopen([obj.pathTempDirectory filesep filename], 'r');
            textscan(fileID, '%[^\n\r]', 5, 'ReturnOnError', false);
            formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
            data = textscan(fileID, formatSpec, 'Delimiter', '\t', 'EmptyValue' ,NaN,'ReturnOnError', false);
            fclose(fileID);
        end

        % readFASTLinOutput opens the FAST output file after linearisation,
        % and calls readFASTLinear to extract data from it.
        %
        %   Input arguments
        %     filename - Name of the FAST output file from which to read
        %     the data
        %
        %   Output
        %     data - String array with content from the FAST output file
        %
        %   Behaviour
        %     The strings for the path to the folder with temporary files
        %     and the specified filename are concatenated, using the
        %     filesep keyword to ensure the appropriate file separation
        %     character for the used platform.
        %     The resulting file is opened and readFASTLinear is called to
        %     extract the data, to keep the FileIOClass agnostic about the 
        %     more complicated format used for this output file.
        %     The filename is checked, since a linearisation can be
        %     completed successfully, while not storing results for a few
        %     failed wind speeds. readFASTLinOutput can be called for these
        %     wind speeds, for which no output file exists.
        %
        %   Called by
        %     WindTurbineClass.linearise
        function data = readFASTLinOutput(obj, filename)
            fileID = fopen([obj.pathTempDirectory filesep filename], 'r');
            if (fileID == -1)
                error(['Linearization file "',fileName,'" could not be opened.']);
            end
            data = readFASTLinear(fileID);
            fclose(fileID);
        end
        
        % writeBladedFile opens a file with the name wind.wnd in the folder
        % for temporary files and calls writeBladed to write the data to it
        % in the appropriate format.
        %
        %   Input arguments
        %     [u, v, w, x, y, z] - Wind data generated by bladedStyle.m, needed to
        %     create Bladed style wind file
        %     U - 10-Minute average wind speed
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The temporary file wind.inp is opened and writeBladed is
        %     called to write the data to it in the appropriate format, to
        %     keep the FileIOClass agnostic about the more complicated
        %     format used for this wind file.
        %
        %   Called by
        %     WindTurbineClass.simulate
        function writeBladedFile(obj, u, v, w, x, y, z, U)
            fileID = fopen([obj.pathTempDirectory filesep 'wind.wnd'], 'w');
            writeBladed(fileID, u, v, w, x, y, z, U);
            fclose(fileID);
        end

        % setSimulationFile stores a path and filename in the class
        % properties for later reference of where simulation results need
        % to be saved.
        %
        %   Input arguments
        %     fileName - Name of the file to save simulation results to
        %     pathName - Path to the file to save simulation results to
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The filename and pathname are stored in the related
        %     properties of the class. The filename will be modified in
        %     saveSimulationResults when it is used for a batch of
        %     simulations. See saveSimulationResults > Behaviour for more
        %     information about this.
        %
        %   Called by
        %     FASTTool.setSimulationfile
        function setSimulationFile(obj, fileName, pathName)
            obj.simulationFile = fileName;
            obj.simulationPath = pathName;
        end

        % saveSimulationResults
        %
        %   Input arguments
        %     Output - Output array from the simulation to be saved
        %     isMultipleWind - Boolean indicating batch runs for multiple
        %     wind speeds
        %     isMultipleSeed - Boolean indicating batch runs for multiple
        %     seeds
        %     windSpeed - 10-Minute average wind speed of the simulation
        %     (only used for filenames of batch runs for multiple wind
        %     speeds)
        %     seed - Seed number of the simulation (only used for filenames
        %     of batch runs for multiple seeds)
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     For batch runs the filename is modified to include the wind
        %     speed and/or the seed number in the filename, when relevant.
        %     A filename is requested from the user only once, at the start
        %     of the batch run, to avoid multiple queries that interrupt
        %     the batch. Simulation results are saved after each simulation
        %     of the batch, and therefore the single filename for the batch
        %     needs to be extended to identify the individual simulations
        %     in the batch.
        %     The filename of the open project, the current time and the
        %     data description of simulation output are stored in the file
        %     information, that will saved along with the simulation
        %     output.
        %     The output from the array is saved after assigning the value
        %     of each element to a variable with the corresponding name,
        %     which is obtained from the data description.
        %
        %   Called by
        %     WindTurbineClass.simulate
        %
        %   Note to developers
        %     Output results in the warning that the argument might not be
        %     used. This is caused by the use of the parameter in a
        %     character string in the command eval. Therefore, its use is
        %     not recognised.
        function saveSimulationResults(obj, Output, isMultipleWind, isMultipleSeed, windSpeed, seed)

            % Output file name
            OutputFile = [obj.simulationPath, obj.simulationFile(1:end-4)];
            if isMultipleWind
                OutputFile = [OutputFile, '_U=', num2str(windSpeed,'%2.2f')];
            end
            if isMultipleSeed
                OutputFile = [OutputFile, '_seed=', int2str(seed)];
            end
            OutputFile = [OutputFile, '.mat'];

            % Store simulation info
            Info.projectFile = obj.projectFileName;
            Info.creationTime = sprintf('%s', datetime);
            Info.dataDescription = obj.dataDescription;
            save(OutputFile, 'Info');

            % Name and save vectors
            OutList = obj.dataDescription(:,1);
            for i = 1:length(OutList)
                eval([OutList{i}, ' = Output{i};']);
                eval(['save(OutputFile, ''', OutList{i}, ''', ''-append'');']);
            end
        end

        % saveMatFile saves a selection of data of the specified type
        % from WindTurbine to a .mat file.
        %
        %   Input arguments
        %     fileName - Name of the file to save to
        %     pathName - Path to the file to save to
        %     type - Type of data to save
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The data specified by type is extracted from WindTurbine.
        %     The filename of the open project and the data description are
        %     stored in the file information, that will be saved along with
        %     the data. The time at which the results are generated is
        %     added in WindTurbineClass after the analysis. For a linear
        %     control analysis, the filename of the used linear model is
        %     added in WindTurbineClass.analyseControl.
        %     From the extracted data from WindTurbine, the parameters that
        %     are specified by the data specification for the given type
        %     are saved to the .mat file.
        %     No check is performed on the specified file, since that
        %     should already be tested before calling saveMatFile. The
        %     file should be forced to have have the extention .mat.
        %
        %   Called by
        %     FASTTool.saveMatfile
        function saveMatFile(obj, fileName, pathName, type)
            switch type
                case 'rotoranalysis'
                    Data = obj.WindTurbine.RotorResults;
                    Data.Info.dataDescription = obj.RotorParameters(2:end,:);
                    Parameters = obj.RotorParameters;
                case 'steadyopanalysis'
                    Data = obj.WindTurbine.SteadyOpResults;
                    Data.Info.dataDescription = obj.SteadyOpParameters(2:end,:);
                    Parameters = obj.SteadyOpParameters;
                case 'modalanalysis'
                    Data = obj.WindTurbine.ModalResults;
                    Data.Info.dataDescription = obj.ModalParameters(2:end,:);
                    Parameters = obj.ModalParameters;
                case 'linearmodel'
                    Data = obj.WindTurbine.LinearModel;
                    Data.Info.dataDescription = obj.LinearModelParameters(2:end,:);
                    Parameters = obj.LinearModelParameters;
                    obj.WindTurbine.LinearModel.fileName = fileName;
                case 'controlanalysis'
                    Data = obj.WindTurbine.ControlResults;
                    Data.Info.dataDescription = obj.ControlParameters(2:end,:);
                    Parameters = obj.ControlParameters;
            end

            Data.Info.projectFile = obj.projectFileName;
            Data.Info.FASTToolVersion = obj.MainApp.version;
            save([pathName, fileName], "-struct", "Data", Parameters{1,1});
            for i = 2 : length(Parameters)
                save([pathName, fileName], "-struct", "Data", Parameters{i,1}, "-append");
            end
        end

        % loadMatFile loads data from the specified file with results or a
        % linear model and stores it in WindTurbine in the parameter
        % specified by type.
        %
        %     fileName - Name of the file to load from
        %     pathName - Path to the file to load from
        %     type - Type of data to load
        %
        %   Output
        %     success - Boolean to indicate whether the specified file
        %     contained data of the specified type
        %
        %   Behaviour
        %     The content of the file is checked, to ensure that it is a
        %     file with data of the specified type. If that is the case,
        %     the data is loaded and stored in the appropriate parameter in
        %     WindTurbine. The filename is also stored, to be displayed in
        %     the results viewer.
        %
        %   Called by
        %     FASTTool.loadMatfile
        function success = loadMatFile(obj, fileName, pathName, type)
            success = true;

            if ~properMatFile(obj, fileName, pathName, type)
                success = false;
                return;
            end

            switch type
                case 'rotoranalysis'
                    obj.WindTurbine.RotorResults = load([pathName,fileName]);
                    obj.WindTurbine.RotorResults.fileName = fileName;
                case 'steadyopanalysis'
                    obj.WindTurbine.SteadyOpResults = load([pathName,fileName]);
                    obj.WindTurbine.SteadyOpResults.fileName = fileName;
                case 'modalanalysis'
                    obj.WindTurbine.ModalResults = load([pathName,fileName]);
                    obj.WindTurbine.ModalResults.fileName = fileName;
                case 'linearmodel'
                    obj.WindTurbine.LinearModel = load([pathName,fileName]);
                    obj.WindTurbine.LinearModel.fileName = fileName;
                case 'controlanalysis'
                    obj.WindTurbine.ControlResults = load([pathName,fileName]);
                    obj.WindTurbine.ControlResults.fileName = fileName;
                case 'simulationresults'
                    obj.WindTurbine.SimulationResults = load([pathName,fileName]);
                    obj.WindTurbine.SimulationResults.fileName = fileName;
            end
        end

        % setDataDescription stores the provided list of parameter names
        % and the provided list of data descriptions in a local class
        % property.
        %
        %   Input arguments
        %     OutList - List of parameter names
        %     Legend - List with descriptions of the parameters
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The provided data is copied to the local class property, such
        %     that it is available to add the data descriptions to the file
        %     information when saving the results of a linearisation or
        %     simulation.
        %
        %   Called by
        %     WindTurbineClass.simulate
        function setDataDescription(obj, OutList, Legend)
            obj.dataDescription = [OutList Legend];
        end

        % parseData modifies or adds to data that was loaded from a project
        % file created in an earlier version of FASTTool, to update it to
        % the needs of the current version.
        %
        %   Input arguments
        %     Data - Data loaded from the project file
        %
        %   Output
        %     Data - Data that is updated to the needs of the current
        %     version of FASTTool
        %
        %   Behaviour
        %     ToDo: Describe this behaviour after cleaning up this function
        %     to make it more versatile for future versions and version
        %     differences.
        %
        %   Called by
        %     FileIOClass.openProject
        %     FileIOClass.getRestoreData
        function Data = parseData(obj, Data)

            % Variable 'Info' was introduced in version 2.0, so if it
            % doesn't exist, the project file is from before that
            % Currently, only parsing from latest version before 2.0 to
            % version 2.0 or 2.1 is implemented, so no further check is
            % done on the version of the project file.
            if ~isfield(Data, "Info")
                    Data.CertificationSettings.Run.StartupTime = 0;
                    Data.CertificationSettings.Wind.AirDensity = Data.AirDensity;
                    Data.CertificationSettings.Run.DT = Data.Control.DT;
                    if Data.CertificationSettings.Wind.Type == 2 % Type 2: Stepped wind
                        Data.CertificationSettings.Wind.Speed = Data.CertificationSettings.Wind.Step:1:Data.CertificationSettings.Run.WindSpeed;
                    else
                        if isscalar(Data.CertificationSettings.Run.WindSpeed)
                            Data.CertificationSettings.Wind.Speed = Data.CertificationSettings.Run.WindSpeed;
                        else
                            % As of version 2.0 only arrays of wind speeds
                            % with equidistant steps are allowed.
                            % Therefore, the new array is filled
                            % automatically from the first to the last
                            % value of the wind-speed array.
                            Data.CertificationSettings.Wind.Speed = ...
                                Data.CertificationSettings.Run.WindSpeed(1): ...
                                (Data.CertificationSettings.Run.WindSpeed(end)-Data.CertificationSettings.Run.WindSpeed(1))/(length(Data.CertificationSettings.Run.WindSpeed)-1): ...
                                Data.CertificationSettings.Run.WindSpeed(end);
                        end
                    end

                    Data.CertificationSettings.Mode.EventTime = Data.CertificationSettings.Mode.Actiontime;
                    
                    switch Data.CertificationSettings.Wind.Type
                        case 7 % Type 7: Extreme wind shear (EWS)
                            Data.CertificationSettings.Wind.EventTime = Data.CertificationSettings.Wind.EWS;

                        case 9 % Type 9: Extreme operating gust (EOG)
                            Data.CertificationSettings.Wind.EventTime = Data.CertificationSettings.Wind.EOG;

                        case 10 % Type 10: Extreme direction change (EDC)
                            Data.CertificationSettings.Wind.EventTime = Data.CertificationSettings.Wind.EDC;

                        case 11 % Type 11: Extreme coherent gust (ECG)
                            Data.CertificationSettings.Wind.EventTime = Data.CertificationSettings.Wind.ECG;
                        otherwise
                            Data.CertificationSettings.Wind.EventTime = 30;
                    end
            %    end
            end
        end

        % properMatFile checks whether the specified file contains data of
        % the required type
        %
        %   Input arguments
        %     fileName - Name of the file to check
        %     pathName - Path to the file to check
        %     type - Type of data that should be in the file
        %
        %   Output
        %     success - Boolean indicating whether the file contains data
        %     of the required type
        %
        %   Behaviour
        %     The function checks whether all the names of parameters
        %     declared for the specified type in the class property that
        %     lists parameter names and specifications for that type of
        %     data appears in the names of the parameters stored in the
        %     file.
        %
        %   Called by
        %     FileIOClass.openProject
        %     FileIOClass.loadMatFile
        function success = properMatFile(obj, fileName, pathName, type)
            success = true;
            contents = whos('-file', [pathName, fileName]);
            switch type
                case 'project'
                    Parameters = obj.ProjectParameters;
                case 'rotoranalysis'
                    Parameters = obj.RotorParameters;
                    if length(contents) ~= size(Parameters,1) % All parameters in obj.RotorParameters also appear in obj.SteadyOpParameters
                        success = false;
                    end
                case 'steadyopanalysis'
                    Parameters = obj.SteadyOpParameters;
                case 'modalanalysis'
                    Parameters = obj.ModalParameters;
                case 'linearmodel'
                    Parameters = obj.LinearModelParameters;
                case 'controlanalysis'
                    Parameters = obj.ControlParameters;
                case 'simulationresults'
                    Parameters = obj.SimulationParameters;
            end
            for i = 2 : length(Parameters) % The first parameter, 'Info', is not checked, since it was not part of project and linear-model files before version 2.0
                if ~ismember(Parameters{i,1}, {contents.name})
                    success = false;
                    return;
                end
            end
        end

        % properTableFile checks whether the specified file contains data
        % of the required type
        %
        %   Input arguments
        %     xVariableName - Name of the variable to check
        %     type - Type of data that should be in the file
        %
        %   Output
        %     success - Boolean indicating whether the file contains data
        %     of the required type
        %
        %   Behaviour
        %     The function compares xVariableName with the name of the
        %     variable that should be in the first column of the file.
        %     Variable names of the other columns are not checked, because
        %     the naming of the first variable is already unique for the
        %     possible file types. It is assumed that the integrity of the
        %     file doesn't need to be validated.
        %
        %   Called by
        %     FileIOClass.readImportTable
        function success = properTableFile(obj, xVariableName, type)
            success = true;
            switch type
                case 'bladetable'
                    if ~strcmp(xVariableName, 'Radius [m]')
                        success = false;
                    end
                case 'polartable'
                    if ~strcmp(xVariableName, 'alpha [deg]')
                        success = false;
                    end
                case 'geometrytable'
                    if ~strcmp(xVariableName, 'x [-]')
                        success = false;
                    end
                case 'towertable'
                    if ~strcmp(xVariableName, 'Height [m]')
                        success = false;
                    end
                case 'controltable'
                    if ~strcmp(xVariableName, 'Pitch [deg]')
                        success = false;
                    end
            end
        end

    end
end