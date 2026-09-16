classdef WindTurbineClass < handle
% The WindTurbineClass manages the domain-related content and
% procedures.
%
% The properties of the WindTurbineClass provide the overarching data
% structure for the specifications of the wind turbine. They also contain
% the results from computations, for visualisation and storage purposes.
%
% Details of this data structure are not specified explicitly in the
% property declaration of the class. Turbine specifications (Project data)
% are further detailed by loading the project (mat) file. The FileIOClass
% updates the data structure to what is used in the current version of
% FASTTool, if needed. Details of the data structure for the results and
% linear model are specified by additions of variables to the main objects
% in the methods of the WindTurbineClass.
%
% Many of the functions in this class analyse certain aspects of the wind
% turbine, based on the provided specifications. This is partly done by
% internal functions (setRatedPower, performRotorAnalysis,
% performSteadyOpAnalysis, and analyseControl) and partly by using external
% tools, directly or in Simulink (performModalAnalysis, linearise, and
% simulate). Other functions extent and support the input of wind-turbine
% specifications, using knowledge about wind turbines to do some (simple)
% calculations for this. The function viewModel opens and shows the model
% in Simulink. It is included in this class, rather than in the FASTTool
% class because the Simulink model is considered to be a domain-specific
% part of the tool.
%
% The functions in this class have no or limited input and output
% variables, because the relevant data are passed on via the data
% structures provided by the handles that are exchanged between this class,
% FileIOClass and the FASTTool class.
%
% The WindTurbineClass is a handle class, because the WindTurbine object
% that is created from this class definition in FASTTool.StartupFcn is
% passed on to and used by many components of the app. As a handle object,
% this is done by reference instead of by making copies of the object.

    properties
        % Handles to main app and FileIO
        MainApp
        FileIO

        % Location of app and other files (executables, Simulink model,
        % etc.) after installation. This is needed for calls to run or use
        % these files, since that doesn't go via the Matlab search path.
        appLocation

        % Project data
        % The data to be stored under these names will come from the
        % project file, which sets the data structure. If project files are
        % opened that are created by an older version of FASTTool,
        % FileIOClass.parseData should correct and complete the data
        % structure where needed, to make it compatible with the current
        % needs. When saving a project, all data in these structures will
        % be saved to file.
        Airfoil
        Appearance
        Blade
        CertificationSettings
        Control
        Drivetrain
        Nacelle
        Tower

        % Results
        % The results from the different types of analyses that can be
        % performed will be stored under these names. The data structures
        % will be created inside the functions that perform the analysis,
        % which are given in the comments behind the name. The data
        % specification (info) is added to it in FileIOClass.saveMatFile.
        % If previously saved results are loaded from file, the data
        % structure and content are set by loading the file in
        % FileIOClass.loadMatFile.
        % The related properties in FileIOClass (called *Parameters,
        % instead of *Results, so RotorParameters relates to RotorResults)
        % determine which parameters stored in these data structures will
        % be saved to file when analysis results are saved.
        RotorResults        % performRotorAnalysis
        SteadyOpResults     % setRatedPower, performSteadyOpAnalysis
        ModalResults        % performModalAnalysis
        ControlResults      % analyseControl
        SimulationResults   % simulate (indirectly, by saving the data structure 'Output' to file via the call to FileIOClass.saveSimulationResults)

        % Linear model
        % The data structure for the linear model is determine by
        % assignments in the function Linearise in this class. The filename
        % to which it is saved and the data specification (info) are added
        % to it in FileIOClass.saveMatFile. If a linear model are loaded
        % from file, the data structure and content is set by loading the
        % file in FileIOClass.loadMatFile.
        LinearModel
    end

    methods

        % setHandles stores references to the FileIO and MainApp handle
        % objects in local properties for later use.
        %
        %   Input arguments
        %     FileIO - Handle to object of FileIOClass
        %     MainApp - Handle to main app (FASTTool)
        %
        %   Output
        %     [-]
        %
        %   Called by
        %     FASTTool.StartupFcn
        function setHandles(obj, FileIO, MainApp)
            obj.FileIO = FileIO;
            obj.MainApp = MainApp;

            [pathName, ~, ~] = fileparts(mfilename("fullpath"));
            obj.appLocation = [pathName filesep];
        end

        % setRatedPower computes the value of the rated power from the
        % torque control settings and assigns it to the SteadyOpResults.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The rated power is determined from the torque control
        %     settings. Therefore, this function is called each time these
        %     settings are changed (including when a new project is
        %     opened).
        %     This function doesn't check whether the rated power is
        %     achieved inside the range of cut-in and cut-out wind speeds,
        %     but this check will be performed when the power curve is
        %     shown to the user in SteadyOpAnalysis.plotSteadyOpAnalysis.
        %     The function was introduced to be able to display the rated
        %     power in SteadyOpAnalysis.plotSteadyOpAnalysis and to use it
        %     as a reference value there to determine the rated wind speed.
        %     The rated power cannot be reliably obtained from the array of
        %     powers determined in
        %     WindTurbineClass.performSteadyOpAnalysis. If the range of
        %     analysed wind speeds is too small, it may not contain the
        %     rated power. Furthermore, the power in full-load wind speeds
        %     is not entirely constant in this array, leading to small
        %     deviations from the theoretical rated power and potentially a
        %     misidentification of rated wind speed if the power at very
        %     high wind speeds is fractionally higher than elsewhere in
        %     full load.
        %
        %   Called by
        %     Control1Specs.applyEdits
        function setRatedPower(obj)
            obj.SteadyOpResults.ratedPower = obj.Control.Torque.SpeedC*(2*pi/60) *  obj.Control.Torque.Demanded * obj.Drivetrain.Generator.Efficiency;
        end

        % performRotorAnalysis determines the power, thrust, and torque
        % coefficients of the rotor as a function of tip speed ratio.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The performance coefficients are determined for each pitch
        %     angle as requested by the user, for tip speed ratios between
        %     1 and 20. It uses the BEM calculations in
        %     performanceCoefficients.m for this.
        %
        %   Called by
        %     RotorAnalysis.RunButtonPushed
        function performRotorAnalysis(obj)

            % Store analysis info
            obj.RotorResults.Info.creationTime = sprintf('%s', datetime);

            pitch = obj.RotorResults.pitch;

            % Find rotor performance coefficients
            tsr = 0:0.1:20;
            startTSR = 11; % Tip speed ratios below 1 omitted to avoid problems in the performance calculations

            cP = zeros(length(pitch), length(tsr));
            cT = zeros(length(pitch), length(tsr));
            cQ = zeros(length(pitch), length(tsr));

            for i = 1:length(pitch)
                for j = startTSR:length(tsr)
                    reportProgress(obj.MainApp, ['Pitch angle: ', num2str(pitch(i)), ' deg']);
                    reportProgress(obj.MainApp, (j - startTSR) / length(tsr));

                    [cT(i,j), cQ(i,j)] = performanceCoefficients(obj.Blade, obj.Airfoil, pitch(i), tsr(j));
                    cP(i,j) = cQ(i,j)*tsr(j);

                    if cP(i,j) < 0 && tsr(j) >= 1
                        cP(i,j) = 0;
                        cQ(i,j) = 0;
                        break;
                    end
                end
            end

            % Store results
            obj.RotorResults.tsr = tsr;
            obj.RotorResults.cP = cP;
            obj.RotorResults.cT = cT;
            obj.RotorResults.cQ = cQ;
        end

        % performSteadyOpAnalysis determines various performance parameters
        % as a function of wind speed.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The performance parameters are determined for each wind speed
        %     as requested by the user. It uses steadyState.m to determine
        %     values for the minimum set of parameters that specify the
        %     performance of the turbine. Additional performance parameters
        %     are determined by simple calculations.
        %
        %   Called by
        %     SteadyOpAnalysis.RunButtonPushed
        function performSteadyOpAnalysis(obj)

            % Store analysis info
            obj.SteadyOpResults.Info.creationTime = sprintf('%s', datetime);

            wind = obj.SteadyOpResults.windSpeed;

            % Determine steady state curves
            [cT, cQ, omega, pitchAngle] = steadyState(obj.Blade, obj.Airfoil, obj.Drivetrain, obj.Control, wind, obj.CertificationSettings.Wind.AirDensity);

            tsr = omega.*obj.Blade.Radius(end)./wind;
            cP = cQ .* tsr;

            P = 0.5*obj.CertificationSettings.Wind.AirDensity*pi*obj.Blade.Radius(end)^2 * wind.^3 .* cP .* obj.Drivetrain.Gearbox.Efficiency .* obj.Drivetrain.Generator.Efficiency;
            T = 0.5*obj.CertificationSettings.Wind.AirDensity*pi*obj.Blade.Radius(end)^2 * wind.^2 .* cT;
            Q = 0.5*obj.CertificationSettings.Wind.AirDensity*pi*obj.Blade.Radius(end)^3 * wind.^2 .* cQ;
            rpm = omega * 60/(2*pi);

            % Store results
            obj.SteadyOpResults.electricalPower = P;
            obj.SteadyOpResults.thrust = T;
            obj.SteadyOpResults.torque = Q;
            obj.SteadyOpResults.rotorSpeed = rpm;
            obj.SteadyOpResults.pitch = pitchAngle;
            obj.SteadyOpResults.tsr = tsr;
            obj.SteadyOpResults.cP = cP;
            obj.SteadyOpResults.cT = cT;
            obj.SteadyOpResults.cQ = cQ;
        end

        % performModalAnalysis determines the natural frequencies and mode
        % shapes for tower modes and blade modes.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The modal analysis is prepared and run for the tower once,
        %     since tower modes are unique. The same is done for the blade
        %     modes, for 6 rotational speeds that envelop the operational
        %     speed range and stand-still. The blade natural frequencies
        %     are a function of rotational speed due to centrifugal
        %     stiffening. The tower natural frequencies are copied in an
        %     array of the same length, to facilitate plotting these for
        %     the same rotational speeds as the blade natural frequencies.
        %     The mode shapes are expressed in both the coordinates of the
        %     deformed tower or blade and in coefficients of a polynomial
        %     representation of these coordinates.
        %
        %   Called by
        %     ModalAnalysis.RunButtonPushed
        function performModalAnalysis(obj)
            
            % Store analysis info
            obj.ModalResults.Info.creationTime = sprintf('%s', datetime);

            % Run BModes for tower
            reportProgress(obj.MainApp, 'Running BModes for tower');
            fileLocation = writeDataFile(obj.FileIO, 'BModes_tower.dat', bModesTower(obj.Tower));
            fileLocation = writeDataFile(obj.FileIO, 'BModes.bmi', bModesInput(obj.Blade,obj.Tower,obj.Nacelle,obj.Control,2,0,fileLocation));
            [~, ~] = system(['"' obj.appLocation 'BModes" "' fileLocation '"']);
            data = readBModesOutput(obj.FileIO, 'BModes.out');
            [y11_shape, y11_coeff, y11_freq, ...
                y12_shape, y12_coeff, y12_freq, ...
                y21_shape, y21_coeff, y21_freq, ...
                y22_shape, y22_coeff, y22_freq] = bModesOutput(data);

            % Set rotor speeds
            rangeEnd = ceil((obj.Control.Torque.SpeedC/obj.Drivetrain.Gearbox.Ratio)/5)*5;
            rangeStep = rangeEnd/5;
            rotorSpeed = 0:rangeStep:rangeEnd;
            samples = ones(size(rotorSpeed));

            % Store results for tower
            obj.ModalResults.Tower_ForeAft1_shape = y21_shape;
            obj.ModalResults.Tower_ForeAft2_shape = y22_shape;
            obj.ModalResults.Tower_ForeAft1_freq = y21_freq * samples;
            obj.ModalResults.Tower_ForeAft2_freq = y22_freq * samples;
            obj.ModalResults.Tower_ForeAft1_coeff = y21_coeff;
            obj.ModalResults.Tower_ForeAft2_coeff = y22_coeff;
            obj.ModalResults.Tower_SideSide1_shape = y11_shape;
            obj.ModalResults.Tower_SideSide2_shape = y12_shape;
            obj.ModalResults.Tower_SideSide1_freq = y11_freq * samples;
            obj.ModalResults.Tower_SideSide2_freq = y12_freq * samples;
            obj.ModalResults.Tower_SideSide1_coeff = y11_coeff;
            obj.ModalResults.Tower_SideSide2_coeff = y12_coeff;

            % Store results for rotor speeds
            obj.ModalResults.rotorSpeed = rotorSpeed;
            obj.ModalResults.cutInRotorSpeed = obj.Control.Torque.SpeedA/obj.Drivetrain.Gearbox.Ratio;
            obj.ModalResults.ratedRotorSpeed = obj.Control.Torque.SpeedC/obj.Drivetrain.Gearbox.Ratio;

            % Empty blade mode frequency vectors
            obj.ModalResults.Blade_Flap1_freq = nan(size(rotorSpeed));
            obj.ModalResults.Blade_Flap2_freq = nan(size(rotorSpeed));
            obj.ModalResults.Blade_Edge1_freq = nan(size(rotorSpeed));
            obj.ModalResults.Blade_Edge2_freq = nan(size(rotorSpeed));

            % Cycle through rotor speeds
            for i = 1:length(rotorSpeed)

                % Run BModes for blade
                reportProgress(obj.MainApp, ['Running BModes for blade at: ', num2str(rotorSpeed(i)), ' of ', num2str(rangeEnd),' PRM']);
                fileLocation = writeDataFile(obj.FileIO, 'BModes_blade.dat', bModesBlade(obj.Blade));
                fileLocation = writeDataFile(obj.FileIO, 'BModes.bmi', bModesInput(obj.Blade,obj.Tower,obj.Nacelle,obj.Control,1,rotorSpeed(i),fileLocation));
                [~, ~] = system(['"' obj.appLocation 'BModes" "' fileLocation '"']);
                data = readBModesOutput(obj.FileIO, 'BModes.out');
                [y11_shape, y11_coeff, y11_freq, ...
                    y12_shape, y12_coeff, y12_freq, ...
                    y21_shape, y21_coeff, y21_freq, ...
                    y22_shape, y22_coeff, y22_freq] = bModesOutput(data);

                % Store results for mode frequencies
                obj.ModalResults.Blade_Flap1_freq(i) = y11_freq;
                obj.ModalResults.Blade_Flap2_freq(i) = y12_freq;
                obj.ModalResults.Blade_Edge1_freq(i) = y21_freq;
                obj.ModalResults.Blade_Edge2_freq(i) = y22_freq;

                % Store results for stand-still modes
                if rotorSpeed(i) == 0

                    obj.ModalResults.Blade_Flap1_shape = y11_shape;
                    obj.ModalResults.Blade_Flap2_shape = y12_shape;
                    obj.ModalResults.Blade_Edge1_shape = y21_shape;
                    obj.ModalResults.Blade_Edge2_shape = y22_shape;
                    obj.ModalResults.Blade_Flap1_coeff = y11_coeff;
                    obj.ModalResults.Blade_Flap2_coeff = y12_coeff;
                    obj.ModalResults.Blade_Edge1_coeff = y21_coeff;
                    obj.ModalResults.Blade_Edge2_coeff = y22_coeff;

                end
            end
        end

        % RescheduleControl maps the gain scheduling table to the pitch
        % angles in the loaded linear model by interpolation of the
        % existing settings.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The controller parameters in the existing gain scheduling
        %     table are interpolated to the pitch angles in the open linear
        %     model using the same interpolation method that is used in
        %     simulink for the control behaviour. Therefore, the behaviour
        %     of the controller after this rescheduling should remain the
        %     same.
        %
        %   Called by
        %     ControlAnalysis.RescheduleButtonPushed
        function RescheduleControl(obj)
            x = obj.Control.Pitch.ScheduledPitchAngles;
            xq = obj.ControlResults.pitch;
            obj.Control.Pitch.LowPassCutOffFreqGS = interp1(x, obj.Control.Pitch.LowPassCutOffFreqGS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.KpGS = interp1(x, obj.Control.Pitch.KpGS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.KiGS = interp1(x, obj.Control.Pitch.KiGS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.Notch_beta1GS = interp1(x, obj.Control.Pitch.Notch_beta1GS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.Notch_beta2GS = interp1(x, obj.Control.Pitch.Notch_beta2GS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.Notch_wnGS = interp1(x, obj.Control.Pitch.Notch_wnGS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.Notch2_beta1GS = interp1(x, obj.Control.Pitch.Notch2_beta1GS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.Notch2_beta2GS = interp1(x, obj.Control.Pitch.Notch2_beta2GS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.Notch2_wnGS = interp1(x, obj.Control.Pitch.Notch2_wnGS, clip(xq, min(x), max(x)),'linear','extrap');
            obj.Control.Pitch.ScheduledPitchAngles = obj.ControlResults.pitch;
        end

        % analyseControl creates transfer functions for the
        % pitch-controller elements individually and in combinations, for
        % the plant, and for the loop gain and it assesses the gain and
        % phase margins. 
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The transfer functions for the plant and the loop gain are
        %     determined for the loaded linear model. First, the pitch
        %     angles in full load are determined since the pitch controller
        %     is only active in full load. Then the magnitudes and phases
        %     are determined for the transfer functions of the
        %     pitch-controller elements individually and in combination,
        %     for the plant, and for the loop gain. Finally, the gain and
        %     phase margins are determined, and the frequencies at which
        %     these occur. The function also stores the information whether
        %     the gain and phase margins occur inside the range of
        %     frequencies used for the magnitude and phase data.
        %
        %   Called by
        %     ControlAnalysis.LineariseButtonPushed
        %     ControlAnalysis.OpenButtonPushed
        function analyseControl(obj)

            % Store analysis info
            obj.ControlResults.Info.linearModelFile = obj.LinearModel.fileName;
            obj.ControlResults.Info.creationTime = sprintf('%s', datetime);

            % Identify pitch angles for the full load region, which must be
            % larger than the fine-pitch angle
            obj.ControlResults.fullLoadIndex = find((obj.LinearModel.Lin.Pitch > obj.Control.Pitch.Fine*pi/180));
            pitch = obj.LinearModel.Lin.Pitch(obj.ControlResults.fullLoadIndex);
            obj.ControlResults.pitch = pitch;
            obj.ControlResults.windSpeed = obj.LinearModel.Lin.V(obj.ControlResults.fullLoadIndex);

            obj.ControlResults.omega_llimit = -2;
            obj.ControlResults.omega_ulimit = 2;
            obj.ControlResults.omega = logspace(obj.ControlResults.omega_llimit, obj.ControlResults.omega_ulimit, 1000);

            % Create transfer function of filters, controller and plant in
            % series for all pitch angles
            obj.ControlResults.Magnitude = cell(1, 9);
            obj.ControlResults.Phase = cell(1,9);
            obj.ControlResults.responseType = ["LPF" "PI" "Notch" "PI_Notch" "LPF_Notch" "LPF_PI" "LPF_PI_Notch" "Plant" "LoopGain"];
            Controller = computeController(pitch, obj.Control, "LPF");
            [obj.ControlResults.Magnitude{1}, obj.ControlResults.Phase{1}] = ...
                computeFrequencyResponse(Controller, obj.ControlResults.omega);

            Controller = computeController(pitch, obj.Control, "PI");
            [obj.ControlResults.Magnitude{2}, obj.ControlResults.Phase{2}] = ...
                computeFrequencyResponse(Controller, obj.ControlResults.omega);

            Controller = computeController(pitch, obj.Control, "Notch");
            [obj.ControlResults.Magnitude{3}, obj.ControlResults.Phase{3}] = ...
                computeFrequencyResponse(Controller, obj.ControlResults.omega);

            Controller = computeController(pitch, obj.Control, ["PI", "Notch"]);
            [obj.ControlResults.Magnitude{4}, obj.ControlResults.Phase{4}] = ...
                computeFrequencyResponse(Controller, obj.ControlResults.omega);

            Controller = computeController(pitch, obj.Control, ["LPF", "Notch"]);
            [obj.ControlResults.Magnitude{5}, obj.ControlResults.Phase{5}] = ...
                computeFrequencyResponse(Controller, obj.ControlResults.omega);

            Controller = computeController(pitch, obj.Control, ["LPF", "PI"]);
            [obj.ControlResults.Magnitude{6}, obj.ControlResults.Phase{6}] = ...
                computeFrequencyResponse(Controller, obj.ControlResults.omega);

            %Controller with all elements computed last, to be used later
            %for loop gain later
            Controller = computeController(pitch, obj.Control, ["LPF", "PI", "Notch"]);
            [obj.ControlResults.Magnitude{7}, obj.ControlResults.Phase{7}] = ...
                computeFrequencyResponse(Controller, obj.ControlResults.omega);

            %Transfer function plant
            Plant = tf(1,1)*ones(1,length(pitch));
            for i = 1:length(pitch)
                iGenSpeed = contains(obj.LinearModel.sysm{obj.ControlResults.fullLoadIndex(i),1}.OutputName, 'ED GenSpeed');
                iBlPitchCPC = contains(obj.LinearModel.sysm{obj.ControlResults.fullLoadIndex(i),1}.InputName, 'collective blade-pitch');
                Plant(1,i) = pi/30*obj.LinearModel.sysm{obj.ControlResults.fullLoadIndex(i),1}(iGenSpeed,iBlPitchCPC);
            end
            [obj.ControlResults.Magnitude{8}, obj.ControlResults.Phase{8}] = ...
                computeFrequencyResponse(Plant, obj.ControlResults.omega);

            %Transfer function loop gain
            LoopGain = tf(1,1)*ones(1,length(pitch));
            for i = 1:length(pitch)
                LoopGain(1,i) = series(Controller(:,i), Plant(:,i));
            end
            [obj.ControlResults.Magnitude{9}, obj.ControlResults.Phase{9}] = ...
                computeFrequencyResponse(LoopGain, obj.ControlResults.omega);

            % Ensure the data is always stored column-wise for single line plots
            if size(obj.ControlResults.Magnitude{1}, 1) == 1
                for i = 1 : 9
                    obj.ControlResults.Magnitude{i} = obj.ControlResults.Magnitude{i}(:);
                    obj.ControlResults.Phase{i} = obj.ControlResults.Phase{i}(:);
                end
            end

            % Stability margins calculation
            cell_prealloc = cell(1,length(pitch));
            S = struct('GainMargin', cell_prealloc, ...
                'GMFrequency', cell_prealloc, ...
                'PhaseMargin', cell_prealloc, ...
                'PMFrequency', cell_prealloc, ...
                'DelayMargin', cell_prealloc, ...
                'DMFrequency', cell_prealloc, ...
                'Stable', cell_prealloc);
            obj.ControlResults.GM = nan(1, length(pitch));
            obj.ControlResults.PM = nan(1, length(pitch));
            obj.ControlResults.GMFreq = nan(1, length(pitch));
            obj.ControlResults.PMFreq = nan(1, length(pitch));
            obj.ControlResults.IsGMFreqWithinLimit = true(1, length(pitch));
            obj.ControlResults.IsPMFreqWithinLimit = true(1, length(pitch));

            LoopGainMagResponseAbs = db2mag(obj.ControlResults.Magnitude{9}); % for margins calculation
            for i = 1:length(pitch)
                S(i) = allmargin(LoopGainMagResponseAbs(:,i), obj.ControlResults.Phase{9}(:,i), obj.ControlResults.omega');
                try
                    obj.ControlResults.GM(i) = mag2db(S(i).GainMargin(1));
                    obj.ControlResults.PM(i) = S(i).PhaseMargin(1);
                catch
                    S(i) = allmargin(LoopGain(1,i));
                    try
                        obj.ControlResults.GM(i) = mag2db(S(i).GainMargin(1));
                        obj.ControlResults.PM(i) = S(i).PhaseMargin(1);
                    catch
                        continue;
                    end
                end

                obj.ControlResults.GMFreq(i) = S(i).GMFrequency(1);
                obj.ControlResults.IsGMFreqWithinLimit(i) = (log10(obj.ControlResults.GMFreq(i)) >= obj.ControlResults.omega_llimit) && (log10(obj.ControlResults.GMFreq(i)) <= obj.ControlResults.omega_ulimit);

                obj.ControlResults.PMFreq(i) = S(i).PMFrequency(1);
                obj.ControlResults.IsPMFreqWithinLimit(i) = (log10(obj.ControlResults.PMFreq(i)) >= obj.ControlResults.omega_llimit) && (log10(obj.ControlResults.PMFreq(i)) <= obj.ControlResults.omega_ulimit);
            end
        end

        % linearise linearises the behaviour of the wind turbine for wind
        % speeds in the full-load range.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     success - Coded integer indicating success of the linearisation or what went wrong
        %     exception - Exception object created when catching an error
        %     during linearisation
        %
        %   Behaviour
        %     First, the range of wind speeds in full load for which the
        %     linearisation needs to be performed is determined. If this
        %     range is not empty, preparations for the linearisation that
        %     are independent of wind speed are done. After this, for each
        %     wind speed the wind-speed dependent preparations are done, a
        %     simulation is performed that ends with a linearisation and
        %     the resulting linear model is stored in LinearModel. When a
        %     linearisation fails, no linearisation is done for the
        %     remaining wind speeds in the range.
        %     The indicator of success is coded as
        %      -1 - Linearisation completed succesfully for some, but not
        %           all wind speeds
        %      -2 - Linearisation failed
        %       0 - No wind speeds above rated were found
        %       1 - Linearisation completed successfully for all wind
        %           speeds
        %     When linearisation failed, the exception object that is
        %     created is passed on, so that it can be rethrown in
        %     ControlAnalysis.linearise after informing the user about the
        %     failure. Otherwise, exception is empty.
        %     Catching errors and passing the exception object is done
        %     inside this function, instead of wrapping it around the call
        %     to this function, to perform the 'try' per wind speed and to
        %     enable saving results if linearisation of only part of the
        %     wind speed range was successful.
        %
        %   Called by
        %     ControlAnalysis.LineariseButtonPushed
        function [success, exception] = linearise(obj)

            success = 1;
            exception = [];
            warnStruct = warning('off', 'backtrace'); % suppress tracing information of warnings

            reportProgress(obj.MainApp, 'Determining full load range');
            reportProgress(obj.MainApp, 0.05);

            % Determine full-load wind-speed range

            % Find rated wind speed by determining the wind speed at which the torque at rotational speed C equals the demanded torque
            OmegaC = obj.Control.Torque.SpeedC*2*pi/60;
            U1 = obj.Control.WindSpeed.Cutin;
            U2 = obj.Control.WindSpeed.Cutout;
            [~, CQr] = performanceCoefficients(obj.Blade, obj.Airfoil, obj.Control.Pitch.Fine, (OmegaC/obj.Drivetrain.Gearbox.Ratio)*obj.Blade.Radius(end)/U1);
            Qr1 = 0.5*CQr*obj.CertificationSettings.Wind.AirDensity*U1^2*pi*obj.Blade.Radius(end)^3*obj.Drivetrain.Gearbox.Efficiency/obj.Drivetrain.Gearbox.Ratio;
            [~, CQr] = performanceCoefficients(obj.Blade, obj.Airfoil, obj.Control.Pitch.Fine, (OmegaC/obj.Drivetrain.Gearbox.Ratio)*obj.Blade.Radius(end)/U2);
            Qr2 = 0.5*CQr*obj.CertificationSettings.Wind.AirDensity*U2^2*pi*obj.Blade.Radius(end)^3*obj.Drivetrain.Gearbox.Efficiency/obj.Drivetrain.Gearbox.Ratio;

            if Qr1 > obj.Control.Torque.Demanded
                Urated = U1;
            elseif Qr2 < obj.Control.Torque.Demanded
                Urated = U2+1; % Adding '1' is arbitrary, to indicate in the later check that the cut-out wind speed is below the rated wind speed.
            else
                success2 = false;
                for iter = 1:100
                    U_new = U1 + abs((obj.Control.Torque.Demanded-Qr1)/(Qr1-Qr2))*(U2-U1);
                    [~, CQr] = performanceCoefficients(obj.Blade, obj.Airfoil, obj.Control.Pitch.Fine, (OmegaC/obj.Drivetrain.Gearbox.Ratio)*obj.Blade.Radius(end)/U_new);
                    Qr_new = 0.5*CQr*obj.CertificationSettings.Wind.AirDensity*U_new^2*pi*obj.Blade.Radius(end)^3*obj.Drivetrain.Gearbox.Efficiency/obj.Drivetrain.Gearbox.Ratio;
                    if abs((Qr_new - obj.Control.Torque.Demanded)/obj.Control.Torque.Demanded) < 0.005
                        success2 = true;
                        break;
                    end
                    if Qr_new < obj.Control.Torque.Demanded
                        U1 = U_new;
                        Qr1 = Qr_new;
                    else
                        U2 = U_new;
                        Qr2 = Qr_new;
                    end
                end
                Urated = U_new;
                if ~success2
                    warning('Tolerance not met during iteration of rated wind speed');
                end
            end

            % Set wind speed range
            windSpeeds = (ceil(Urated/obj.LinearModel.stepSize)*obj.LinearModel.stepSize) : obj.LinearModel.stepSize : obj.Control.WindSpeed.Cutout;
            if sum(windSpeeds) == 0 || isnan(sum(windSpeeds))
                success = 0;
            else
                % The next three settings for the simulation are hard
                % coded, but can be changed if linearisation is not successful
                TSim = 150; % Simulation time after which linearisation is performed
                LinAmount = 1; % Number of azimuths for which linearisation is averaged
                LinRotations = 1; % Number of rotations for which linearisation is averaged
                LinMode = 'Linearize'; % Alternative is 'LinearizeWithNRELSettings'. See servoDyn.m for explanation why this is not (to be) used.

                reportProgress(obj.MainApp, 'Finding mode shapes (Preparation for FAST)');

                % Run BModes for tower
                fileLocation = writeDataFile(obj.FileIO, 'BModes_tower.dat', bModesTower(obj.Tower));
                fileLocation = writeDataFile(obj.FileIO, 'BModes.bmi', bModesInput(obj.Blade,obj.Tower,obj.Nacelle,obj.Control,2,0,fileLocation));
                [~, ~] = system(['"' obj.appLocation 'BModes" "' fileLocation '"']);
                data = readBModesOutput(obj.FileIO, 'BModes.out');
                [~, y11_coeff, ~, ...
                    ~, y12_coeff, ~, ...
                    ~, y21_coeff, ~, ...
                    ~, y22_coeff, ~] = bModesOutput(data);

                % Store
                obj.Tower.ForeAft1_coeff = y21_coeff;
                obj.Tower.ForeAft2_coeff = y22_coeff;
                obj.Tower.SideSide1_coeff = y11_coeff;
                obj.Tower.SideSide2_coeff = y12_coeff;

                % Run BModes for blade
                fileLocation = writeDataFile(obj.FileIO, 'BModes_blade.dat', bModesBlade(obj.Blade));
                fileLocation = writeDataFile(obj.FileIO, 'BModes.bmi', bModesInput(obj.Blade,obj.Tower,obj.Nacelle,obj.Control,1,0,fileLocation));
                [~, ~] = system(['"' obj.appLocation 'BModes" "' fileLocation '"']);
                data = readBModesOutput(obj.FileIO, 'BModes.out');
                [~, y11_coeff, ~, ...
                    ~, y12_coeff, ~, ...
                    ~, y21_coeff, ~, ...
                    ~, y22_coeff, ~] = bModesOutput(data);

                % Store
                obj.Blade.Flap1_coeff = y11_coeff;
                obj.Blade.Flap2_coeff = y12_coeff;
                obj.Blade.Edge1_coeff = y21_coeff;
                obj.Blade.Edge2_coeff = y22_coeff;

                % Steady state curves
                reportProgress(obj.MainApp, 'Determining steady state rotational speeds and pitch angles');
                [~, ~, OmegaU, PitchAngle] = steadyState(obj.Blade, obj.Airfoil, obj.Drivetrain, obj.Control, windSpeeds, obj.CertificationSettings.Wind.AirDensity);
                RPM = OmegaU * 60/(2*pi);

                % Set the gearbox efficiency (to avoid error during linearisation that ADAMS cannot handle
                % nonideal gearboxes) and remember the true efficiency to reset latersetDataDescription
                TrueGearboxEfficiency = obj.Drivetrain.Gearbox.Efficiency;
                obj.Drivetrain.Gearbox.Efficiency = 1;

                % Turbine input file
                reportProgress(obj.MainApp, 'Preparing FAST input files');

                % Turbine input files
                writeDataFile(obj.FileIO, 'AeroDyn.dat', aeroDyn(obj.Blade,obj.Airfoil,obj.Tower,LinMode,obj.CertificationSettings.Wind.AirDensity));
                writeDataFile(obj.FileIO, 'AeroDyn_blade.dat', aeroDynBlade(obj.Blade));
                for i = 1:length(obj.Blade.IFoil)
                    writeDataFile(obj.FileIO, ['AeroDyn_' obj.Airfoil.Name{obj.Blade.IFoil(i)} '.dat'], aeroDynAirfoil(obj.Blade.IFoil(i),obj.Airfoil));
                end
                writeDataFile(obj.FileIO, 'ServoDyn.dat', servoDyn(obj.Drivetrain,obj.Control,LinMode));

                load([obj.appLocation 'OutList.mat'], 'OutList', 'Legend');
                assignin('base', 'OutList', OutList);
                assignin('base', 'Legend', Legend);

                % Run linearization
                sysm = cell(length(windSpeeds),1);
                for j = 1:length(windSpeeds)

                    % Status update
                    reportProgress(obj.MainApp, ['Linearising at U = ', num2str(windSpeeds(j), '%5.2f'), ' m/s, ', num2str(RPM(j), '%5.2f'), ' rpm, ', num2str(PitchAngle(j), '%5.2f'), ' deg pitch']);
                    reportProgress(obj.MainApp, (j / length(windSpeeds)) - 0.05);

                    % Set initial RPM and pitch angle in ElastoDyn input file
                    writeDataFile(obj.FileIO, 'ElastoDyn.dat', elastoDyn(obj.Blade,obj.Tower,obj.Nacelle,obj.Drivetrain,obj.Control,LinMode,RPM(j),PitchAngle(j)));
                    writeDataFile(obj.FileIO, 'ElastoDyn_blade.dat', elastoDynBlade(obj.Blade));
                    writeDataFile(obj.FileIO, 'ElastoDyn_tower.dat', elastoDynTower(obj.Tower));

                    % Set linearization times for 10 deg azimuth step (after 30 s)
                    LinAziPositions = linspace(0,360*LinRotations,LinAmount+1);
                    LinTimes = TSim + obj.CertificationSettings.Run.DT * round(LinAziPositions(2:end)/(RPM(j)*6) / obj.CertificationSettings.Run.DT);
                    TMax = max(LinTimes)+1.0;

                    FAST_InputFileName = writeDataFile(obj.FileIO, 'FAST.fst', FASTinput(obj.CertificationSettings.Run.DT, TMax, LinMode, LinTimes));


                    % Wind input file
                    type = 1; % Steady wind
                    writeDataFile(obj.FileIO, 'InflowWind.dat', inflowWind(type,windSpeeds(j),obj.Tower.HubHeight));

                    % Run FAST and prevent console output
                    assignin('base', 'TMax', TMax);
                    assignin('base', 'FAST_InputFileName', FAST_InputFileName);

                    % Capture errors during linearisation
                    try
                        evalc('sim(''OpenLoop'',TMax);');
                    catch exception1
                        exception = exception1;
                        if j == 1
                            success = -2;
                            break;
                        else
                            success = -1;
                            break;
                        end
                    end

                    % Extract steady state solution
                    A = 0;
                    B = 0;
                    C = 0;
                    D = 0;
                    for i = 1:LinAmount
                        data = readFASTLinOutput(obj.FileIO, ['FAST.SFunc.', int2str(i), '.lin']);
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
                obj.Drivetrain.Gearbox.Efficiency = TrueGearboxEfficiency;

                % Store analysis info
                if success == -1 || success == 1
                    obj.LinearModel.Info.creationTime = sprintf('%s', datetime);
                    obj.LinearModel.Lin = Lin;
                    obj.LinearModel.sysm = sysm;
                end

                evalin("base", "clear OutList Legend TMax FAST_InputFileName DT");
            end
            warning(warnStruct); % Restore backtracing of warnings to original state
        end

        % simulate simulates the wind turbine at the wind speeds provided
        % by the user, repeating the simulation for the requested number of
        % random realisations of the wind field.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     First, the simulations are prepared by doing modal analysis
        %     and writing the data files for FAST that apply to all wind
        %     speeds and seeds. Then, loops over the required wind speeds
        %     and number of seeds are run. In these loops the initial
        %     rotational speed and pitch angle are determined, the
        %     wind-speed and seed specific input files for FAST are
        %     written, the input files for the wind field are generated,
        %     the simulation is done and the results are extracted and
        %     saved.
        %
        %   Called by
        %     Simulation.RunButtonPushed
        function simulate(obj)

            % Set operational event time in simulation time (including
            % startup time that will be discarded)
            obj.CertificationSettings.Mode.Actiontime = obj.CertificationSettings.Run.StartupTime + obj.CertificationSettings.Mode.EventTime;

            % Set spatial and temporal resolution for wind files used by
            % bladedStyle and turbSim
            obj.CertificationSettings.Wind.dt = 0.1; % (s)
            obj.CertificationSettings.Wind.Ly = 2*obj.Blade.Radius(end); % (m)
            obj.CertificationSettings.Wind.Lz = 2*obj.Blade.Radius(end); % (m)
            obj.CertificationSettings.Wind.Ny = 40; % (-)
            obj.CertificationSettings.Wind.Nz = 40; % (-)

            reportProgress(obj.MainApp, 'Finding mode shapes and preparing FAST input files');
            reportProgress(obj.MainApp, 0.05);

            % Run BModes for tower
            fileLocation = writeDataFile(obj.FileIO, 'BModes_tower.dat', bModesTower(obj.Tower));
            fileLocation = writeDataFile(obj.FileIO, 'BModes.bmi', bModesInput(obj.Blade,obj.Tower,obj.Nacelle,obj.Control,2,0,fileLocation));
            [~, ~] = system(['"' obj.appLocation 'BModes" "' fileLocation '"']);
            data = readBModesOutput(obj.FileIO, 'BModes.out');
            [~, y11_coeff, ~, ...
                ~, y12_coeff, ~, ...
                ~, y21_coeff, ~, ...
                ~, y22_coeff, ~] = bModesOutput(data);

            % Store
            obj.Tower.ForeAft1_coeff = y21_coeff;
            obj.Tower.ForeAft2_coeff = y22_coeff;
            obj.Tower.SideSide1_coeff = y11_coeff;
            obj.Tower.SideSide2_coeff = y12_coeff;

            % Run BModes for blade
            fileLocation = writeDataFile(obj.FileIO, 'BModes_blade.dat', bModesBlade(obj.Blade));
            fileLocation = writeDataFile(obj.FileIO, 'BModes.bmi', bModesInput(obj.Blade,obj.Tower,obj.Nacelle,obj.Control,1,0,fileLocation));
            [~, ~] = system(['"' obj.appLocation 'BModes" "' fileLocation '"']);
            data = readBModesOutput(obj.FileIO, 'BModes.out');
            [~, y11_coeff, ~, ...
                ~, y12_coeff, ~, ...
                ~, y21_coeff, ~, ...
                ~, y22_coeff, ~] = bModesOutput(data);

            % Store
            obj.Blade.Flap1_coeff = y11_coeff;
            obj.Blade.Flap2_coeff = y12_coeff;
            obj.Blade.Edge1_coeff = y21_coeff;
            obj.Blade.Edge2_coeff = y22_coeff;

            % Turbine input files
            TMax = obj.CertificationSettings.Run.StartupTime + obj.CertificationSettings.Run.Time;
            FAST_InputFileName = writeDataFile(obj.FileIO, 'FAST.fst', FASTinput(obj.CertificationSettings.Run.DT, TMax));
            writeDataFile(obj.FileIO, 'AeroDyn.dat', aeroDyn(obj.Blade,obj.Airfoil,obj.Tower,string(obj.CertificationSettings.Mode.Type),obj.CertificationSettings.Wind.AirDensity));
            writeDataFile(obj.FileIO, 'AeroDyn_blade.dat', aeroDynBlade(obj.Blade));
            for i = 1:length(obj.Blade.IFoil)
                writeDataFile(obj.FileIO, ['AeroDyn_' obj.Airfoil.Name{obj.Blade.IFoil(i)} '.dat'], aeroDynAirfoil(obj.Blade.IFoil(i),obj.Airfoil));
            end

            % Send the required parameters to the base workspace and load
            % the Simulink model
            reportProgress(obj.MainApp, 'Preparing Simulink');
            reportProgress(obj.MainApp, 0.1);
            load([obj.appLocation 'OutList.mat'], 'OutList', 'Legend');
            setDataDescription(obj.FileIO, OutList, Legend)
            assignin('base', 'Drivetrain', obj.Drivetrain);
            assignin('base', 'Control', obj.Control);
            assignin('base', 'FAST_InputFileName', FAST_InputFileName);
            assignin('base', 'TMax', TMax);
            assignin('base', 'OutList', OutList);
            assignin('base', 'Legend', Legend);
            assignin('base', 'CertificationSettings', obj.CertificationSettings);
            load_system('FAST');
            open_system('FAST/Scope');

            % Determine whether batches of wind speeds or seeds are set
            isMultipleWind = false;
            isMultipleSeed = false;
            seeds = 1;
            if length(obj.CertificationSettings.Wind.Speed) > 1 && ...
                    obj.CertificationSettings.Wind.Type ~= 2 % Type 2 is stepped wind, which uses multiple wind speeds for one simulation
                isMultipleWind = true;
            end
            if obj.CertificationSettings.Run.Seeds > 1 && ismember(obj.CertificationSettings.Wind.Type, [4 5 6 8]) % Wind types 4, 5, 6 and 8 use stochastic wind fields. For other wind types, the setting for the seeds is ignored
                isMultipleSeed = true;
                seeds = obj.CertificationSettings.Run.Seeds;
            end

            % Loop over wind speeds and seeds
            for i = 1:length(obj.CertificationSettings.Wind.Speed)
                U = obj.CertificationSettings.Wind.Speed(i);
                for seed = 1:seeds


                    if obj.CertificationSettings.Wind.Type == 2 % Stepped wind uses multiple wind speeds for one simulation and is deterministic (so 'seed' is not used)
                        textRunID = sprintf('U = %2.2f - %2.2f m/s - ', obj.CertificationSettings.Wind.Speed(1), obj.CertificationSettings.Wind.Speed(end));
                    elseif obj.CertificationSettings.Wind.Type == 5 % For the annual extreme the wind speed is fixed (so 'U' is not used)
                        textRunID = sprintf('EWM1 | seed %i/%i - ', seed, seeds);
                    elseif obj.CertificationSettings.Wind.Type == 6 % For the 50-year extreme the wind speed is fixed (so 'U' is not used)
                        textRunID = sprintf('EWM50 | seed %i/%i - ', seed, seeds);
                    elseif ismember(obj.CertificationSettings.Wind.Type, [4 8]) % NTM and ETM use stochastic wind fields
                        textRunID = sprintf('U = %2.2f m/s | seed %i/%i - ', U, seed, seeds);
                    else
                        textRunID = sprintf('U = %2.2f m/s - ', U);
                    end
                    progress = 0.1 + 0.9 * ((i - 1) * seeds + seed) / (length(obj.CertificationSettings.Wind.Speed) * seeds + 1);
                    reportProgress(obj.MainApp, progress);
                    textProgress = [textRunID, 'Generating wind file'];
                    reportProgress(obj.MainApp, textProgress);

                    % Find initial RPM and pitch angle
                    [~, ~, OmegaU, P_InitAngle] = steadyState(obj.Blade, obj.Airfoil, obj.Drivetrain, obj.Control, U, obj.CertificationSettings.Wind.AirDensity);
                    RPM_Init = OmegaU * 60/(2*pi);
                    if obj.CertificationSettings.Mode.Type == 3     % Startup
                        RPM_Init = 0;
                        P_InitAngle = obj.Control.Pitch.Max;
                    elseif obj.CertificationSettings.Mode.Type == 6 % Idling
                        RPM_Init = 0;
                        P_InitAngle = obj.Control.Pitch.Max;
                    elseif obj.CertificationSettings.Mode.Type == 7	% Parked
                        RPM_Init = 0;
                        P_InitAngle = obj.Control.Pitch.Max;
                    end

                    % Set operation mode in ElastoDyn file
                    writeDataFile(obj.FileIO, 'ElastoDyn.dat', elastoDyn(obj.Blade,obj.Tower,obj.Nacelle,obj.Drivetrain,obj.Control,string(obj.CertificationSettings.Mode.Type),RPM_Init,P_InitAngle));
                    writeDataFile(obj.FileIO, 'ElastoDyn_blade.dat', elastoDynBlade(obj.Blade));
                    writeDataFile(obj.FileIO, 'ElastoDyn_tower.dat', elastoDynTower(obj.Tower));

                    % Set operation mode in ServoDyn file
                    writeDataFile(obj.FileIO, 'ServoDyn.dat', servoDyn(obj.Drivetrain,obj.Control,string(obj.CertificationSettings.Mode.Type),obj.CertificationSettings.Mode.Actiontime));

                    % Wind input file
                    writeDataFile(obj.FileIO, 'InflowWind.dat', inflowWind(obj.CertificationSettings.Wind.Type,U,obj.Tower.HubHeight));

                    if obj.CertificationSettings.Run.RandomSeed
                        rng('shuffle');
                        seedNumber = randi([-2147483648, 2147483647]);
                    else
                        seedNumber = obj.CertificationSettings.Run.SeedNumber + seed - 1;
                    end

                    switch obj.CertificationSettings.Wind.Type
                        case 1 % Type 1: Steady wind
                            % Nothing else needed

                        case 2 % Type 2: Stepped wind
                            writeDataFile(obj.FileIO, 'wind.sum', windSum(U,obj.Tower.HubHeight));
                            [x, y, z, u, v, w] = bladedStyle(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight,obj.Blade.Radius(end), obj.CertificationSettings.Run.StartupTime);
                            writeBladedFile(obj.FileIO,(u-U)/U,v/U,w/U,x,y,z,U);

                        case 3% Type 3: Normal wind profile (NWP)
                            % Nothing else needed

                        case 4 % Type 4: Normal turbulence model (NTM)
                            writeDataFile(obj.FileIO, 'wind.inp', turbSim(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight, seedNumber));
                            [fileLocation] = getWriteLocation(obj.FileIO, 'wind.inp');
                            [~, ~] = system(['"' obj.appLocation 'TurbSim" [/h] "' fileLocation '"']); % No console output
                              
                        case 5 % Type 5: Annual extreme wind speed (EWM1)
                            writeDataFile(obj.FileIO, 'wind.inp', turbSim(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight, seedNumber)); % 'U' is passed on and written in the turbSim input file, but will be ignored when running TurbSim
                            [fileLocation] = getWriteLocation(obj.FileIO, 'wind.inp');
                            [~, ~] = system(['"' obj.appLocation 'TurbSim" [/h] "' fileLocation '"']); % No console output

                        case 6 % Type 6: 50-year extreme wind speed (EWM50)
                            writeDataFile(obj.FileIO, 'wind.inp', turbSim(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight, seedNumber)); % 'U' is passed on and written in the turbSim input file, but will be ignored when running TurbSim
                            [fileLocation] = getWriteLocation(obj.FileIO, 'wind.inp');
                            [~, ~] = system(['"' obj.appLocation 'TurbSim" [/h] "' fileLocation '"']); % No console output

                        case 7 % Type 7: Extreme wind shear (EWS)
                            writeDataFile(obj.FileIO, 'wind.sum', windSum(U,obj.Tower.HubHeight));
                            [x, y, z, u, v, w] = bladedStyle(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight,obj.Blade.Radius(end), obj.CertificationSettings.Run.StartupTime);
                            writeBladedFile(obj.FileIO,(u-U)/U,v/U,w/U,x,y,z,U);

                        case 8 % Type 8: Extreme turbulence model (ETM)
                            writeDataFile(obj.FileIO, 'wind.inp', turbSim(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight, seedNumber));
                            [fileLocation] = getWriteLocation(obj.FileIO, 'wind.inp');
                            [~, ~] = system(['"' obj.appLocation 'TurbSim" [/h] "' fileLocation '"']); % No console output

                        case 9 % Type 9: Extreme operating gust (EOG)
                            writeDataFile(obj.FileIO, 'wind.sum', windSum(U,obj.Tower.HubHeight));
                            [x, y, z, u, v, w] = bladedStyle(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight,obj.Blade.Radius(end), obj.CertificationSettings.Run.StartupTime);
                            writeBladedFile(obj.FileIO,(u-U)/U,v/U,w/U,x,y,z,U);

                        case 10 % Type 10: Extreme direction change (EDC)
                            writeDataFile(obj.FileIO, 'wind.sum', windSum(U,obj.Tower.HubHeight));
                            [x, y, z, u, v, w] = bladedStyle(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight,obj.Blade.Radius(end), obj.CertificationSettings.Run.StartupTime);
                            writeBladedFile(obj.FileIO,(u-U)/U,v/U,w/U,x,y,z,U);

                        case 11 % Type 11: Extreme coherent gust (ECG)
                            writeDataFile(obj.FileIO, 'wind.sum', windSum(U,obj.Tower.HubHeight));
                            [x, y, z, u, v, w] = bladedStyle(obj.CertificationSettings.Wind,U,obj.Tower.HubHeight,obj.Blade.Radius(end), obj.CertificationSettings.Run.StartupTime);
                            writeBladedFile(obj.FileIO,(u-U)/U,v/U,w/U,x,y,z,U);
                    end

                    textProgress = [textRunID, 'Running Simulink/FAST'];
                    reportProgress(obj.MainApp, textProgress);

                    % Send the required parameters to the base workspace
                    % and run Simulink
                    assignin('base', 'RPM_Init', RPM_Init);
                    assignin('base', 'P_InitAngle', P_InitAngle);
                    assignin('base', 'T_GenSpeedInit', RPM_Init*obj.Drivetrain.Gearbox.Ratio);
                    evalc('sim(''FAST'',TMax);');

                    % Extract output
                    Output = readFASTSimOutput(obj.FileIO, 'FAST.SFunc.out');

                    % Discard results for startup period and let time
                    % array start at zero
                    startIndex = ceil(obj.CertificationSettings.Run.StartupTime / obj.CertificationSettings.Run.DT) + 1;
                    for j = 1:length(Output)
                        Output{j} = Output{j}(startIndex:end);
                    end
                    Output{1} = Output{1} - Output{1}(1);

                    saveSimulationResults(obj.FileIO, Output, isMultipleWind, isMultipleSeed, U, seed);
                end

                if obj.CertificationSettings.Wind.Type == 2 % Stepped wind, which uses multiple wind speeds for one simulation, or the annual or 50-year extreme for which the wind speed is fixed
                    evalin("base", "clear Drivetrain Control FAST_InputFileName TMax RPM_Init T_GenSpeedInit P_InitAngle OutList Legend CertificationSettings DT");
                    return;
                end
            end

            evalin("base", "clear Drivetrain Control FAST_InputFileName TMax RPM_Init T_GenSpeedInit P_InitAngle OutList Legend CertificationSettings DT");
        end

        % setRemainingBladeProperties supplements the blade data provided
        % by the user, to complete the required data for FAST.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     Some details of the blade data that are required for FAST are
        %     outside the scope of the intentional use of FASTTool. These
        %     are therefore set using simple calculations and knowledge
        %     based rules for blade designs.
        %
        %   Called by
        %     RotorSpecs.applyBladeTable
        function setRemainingBladeProperties(obj)
            % These are either directly and uniquely computable from user
            % inputs, or more challenging to determine but with limited
            % impact

            % Smooth blade Thickness
            t = zeros(size(obj.Blade.Radius));
            for i = 1:length(obj.Blade.NFoil)
                t_u = max(obj.Airfoil.Geometry{obj.Blade.IFoil(obj.Blade.NFoil(i))}(2,1:200) * obj.Blade.Chord(i));
                t_l = min(obj.Airfoil.Geometry{obj.Blade.IFoil(obj.Blade.NFoil(i))}(2,200:end) * obj.Blade.Chord(i));
                t(i) = t_u - t_l;
            end
            obj.Blade.Thickness = t;
            s = 2.*round((length(obj.Blade.NFoil)/5+1)/2)-1;
            t = conv(t(:),ones(1,s)/sum(s),'same')./t(:);
            t(1:(s-1)/2) = 1;
            t(end-(s-1)/2:end) = 1;
            obj.Blade.Thickness = obj.Blade.Thickness .* t;

            % Offsets
            i = find(obj.Blade.NFoil > 2);
            x = [obj.Blade.Radius(1), obj.Blade.Radius(i(1)), obj.Blade.Radius(end)];
            cg = [0, 0.2, 0.2];
            sc = [0, -0.03, 0.1];
            obj.Blade.cg = interp1(x,cg,obj.Blade.Radius) .* obj.Blade.Chord;
            obj.Blade.sc = interp1(x,sc,obj.Blade.Radius) .* obj.Blade.Chord;

            % Factors for inertia and torsional stiffness per airfoil
            tref =  [1.000; 0.700; 0.405; 0.350; 0.300; 0.250; 0.210; 0.180];
            Cflap = [0.446; 0.260; 0.147; 0.035; 0.027; 0.014; 0.010; 0.004];
            Cedge = [0.034; 0.028; 0.022; 0.019; 0.017; 0.014; 0.015; 0.018];
            Ctor =  [0.170; 0.136; 0.094; 0.020; 0.021; 0.025; 0.025; 0.022];

            % Historic: parameters were used in earlier data files for FAST
            % Ccgo =  [0.010; 0.018; 0.030; 0.060; 0.047; 0.037; 0.045; 0.060];
            % obj.Blade.ac = nan(size(obj.Blade.Radius));
            % obj.Blade.ac(obj.Blade.NFoil == 1) = 0.25;
            % obj.Blade.ac(obj.Blade.Thickness./obj.Blade.Chord < 0.350) = 0.125;
            % obj.Blade.ac(isnan(obj.Blade.ac)) = interp1(...
            %     obj.Blade.Radius(~isnan(obj.Blade.ac)), ...
            %     obj.Blade.ac(~isnan(obj.Blade.ac)), ...
            %     obj.Blade.Radius(isnan(obj.Blade.ac)));
            % obj.Blade.eo = interp1(tref,Ccgo,obj.Blade.Thickness./obj.Blade.Chord, 'pchip') .* obj.Blade.Chord;

            % obj.Blade structural properties estimated from thickness
            obj.Blade.FlapIner = interp1(tref,Cflap,obj.Blade.Thickness./obj.Blade.Chord, 'pchip').*(obj.Blade.Mass.*obj.Blade.Thickness.^2);
            obj.Blade.EdgeIner = interp1(tref,Cedge,obj.Blade.Thickness./obj.Blade.Chord, 'pchip').*(obj.Blade.Mass.*obj.Blade.Chord.^2);
            obj.Blade.GJ = interp1(tref,Ctor,obj.Blade.Thickness./obj.Blade.Chord, 'pchip').*(obj.Blade.EIflap+obj.Blade.EIedge);
            obj.Blade.EA = 1.3e7 * obj.Blade.Mass;
            obj.Blade.PitchAxis = interp1(tref,[0.5; 0.45; 0.40; 0.375; 0.375; 0.375; 0.375; 0.375], obj.Blade.Thickness./obj.Blade.Chord, 'pchip');
        end

        % setRemainingTowerProperties supplements the towr data provided
        % by the user, to complete the required data for FAST.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     Some details of the tower data that are required for FAST are
        %     outside the scope of the intentional use of FASTTool. Most
        %     are calculated from the provided geometrical inputs.
        %
        %   Called by
        %     TowerSpecs.setTowerTable
        function setRemainingTowerProperties(obj)
            % These are directly and uniquely computable from user inputs
            obj.Tower.ShearModulus = 80.8e9;
            obj.Tower.YoungsModulus = 210e9;
            obj.Tower.Mass = obj.Tower.EffectiveDensity * pi* (obj.Tower.Diameter.^2 - (obj.Tower.Diameter-2*obj.Tower.WallThickness).^2) ./ 4;
            obj.Tower.EI = obj.Tower.YoungsModulus * pi/64*(obj.Tower.Diameter.^4 - (obj.Tower.Diameter-2*obj.Tower.WallThickness).^4);
            obj.Tower.GJ = obj.Tower.ShearModulus * pi/32*(obj.Tower.Diameter.^4 - (obj.Tower.Diameter-2*obj.Tower.WallThickness).^4);
            obj.Tower.EA = obj.Tower.YoungsModulus * pi*obj.Tower.Diameter.*obj.Tower.WallThickness;
            obj.Tower.Iner = 0.5 * obj.Tower.Mass .* (obj.Tower.Diameter/2 - obj.Tower.WallThickness).^2;
        end

        % setNewAerofoilData sets data for a dummy aerofoil.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     The aerofoil data for a newly created aerofoil are filled
        %     with preliminary values, which can be replaced by the user.
        %     The data do not represent a true aerofoil, but they fulfil
        %     the constraints for usage in the tool, such as providing a
        %     geometry that can be used in the 3D visualisation.
        %
        %   Called by
        %     AerofoilSpecs.NewButtonPushed
        %     RotorSpecs.applyBladeTable
        function setNewAerofoilData(obj)
            i = length(obj.Airfoil.Name);
            obj.Airfoil.Geometry{i} = [1,0,1;0,0,0];
            obj.Airfoil.Alpha{i} = [-180;0;180];
            obj.Airfoil.Cl{i} = [0;0;0];
            obj.Airfoil.Cd{i} = [0;0;0];
            obj.Airfoil.Cm{i} = [0;0;0];
            obj.Airfoil.CnSlope(i) = 0;
            obj.Airfoil.StallAngle1(i) = 0;
            obj.Airfoil.StallAngle2(i) = 0;
            obj.Airfoil.CritCn1(i) = 0;
            obj.Airfoil.CritCn2(i) = 0;
            parseAerofoil(obj, obj.Airfoil.Name{i});
        end

        % parseAerofoil maps and/or extends basic aerofoil data to meet the
        % format requirements for this data type.
        %
        %   Input arguments
        %     [-]
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     Various aspects of the dataset are checked and repaired, it
        %     they don't meet the requirements. Duplication of coordinates
        %     of the aerofoil shape is not allowed, and these are removed.
        %     Aerofoil geometry should be stored in 399 coordinates,
        %     starting at the trailing edge, moving over the upper surface
        %     to the leading edge and back to the trailing edge via the
        %     lower surface. Existing data is reordered and resampled to
        %     achieve this. Lift and drag coefficients must be present for
        %     angles of attack from -180 to 180 degrees. Missing data is
        %     supplemented from the data of the NACA 64-618 aerofoil.
        %
        %   Called by
        %     WindTurbineClass.setNewAerofoilData
        %     AerofoilSpecs.processImport
        function parseAerofoil(obj, name)
            i = find(strcmp(name, obj.Airfoil.Name));

            % Aerofoil coordinates
            x = obj.Airfoil.Geometry{i}(1,:);
            y = obj.Airfoil.Geometry{i}(2,:);

            % Find upper and lower surface
            j = find(x == min(x));
            j = j(round(length(j)/2));
            if mean(y(1:j)) > 0
                x_u = x(1:j);
                x_l = x(j:end);
                y_u = y(1:j);
                y_l = y(j:end);
            else
                x_u = x(j:end);
                x_l = x(1:j);
                y_u = y(j:end);
                y_l = y(1:j);
            end

            % Remove duplicates
            [x_u,j] = unique(x_u);
            y_u = y_u(j);
            [x_l,j] = unique(x_l);
            y_l = y_l(j);

            % Normalize chord length
            x_u = x_u - min(x_u);
            x_u = x_u/max(x_u);
            x_l = x_l - min(x_l);
            x_l = x_l/max(x_l);

            % Interpolate
            x = 0.5 + 0.5*cos(linspace(0,2*pi,399));
            y_u = interp1(x_u, y_u, x(1:200), 'pchip');
            y_l = interp1(x_l, y_l, x(200:end), 'pchip');
            y = [y_u(1:200), y_l(2:end)];

            % Update geometry
            obj.Airfoil.Geometry{i} = zeros(2,399);
            obj.Airfoil.Geometry{i}(1,:) = x;
            obj.Airfoil.Geometry{i}(2,:) = y;

            % Complete aerofoil data with behavior of the NACA 64-618
            if i > 2
                if min(obj.Airfoil.Alpha{i}) > -180
                    obj.Airfoil.Cl{i} = [obj.Airfoil.Cl{8}(obj.Airfoil.Alpha{8} < min(obj.Airfoil.Alpha{i})); obj.Airfoil.Cl{i}];
                    obj.Airfoil.Cd{i} = [obj.Airfoil.Cd{8}(obj.Airfoil.Alpha{8} < min(obj.Airfoil.Alpha{i})); obj.Airfoil.Cd{i}];
                    obj.Airfoil.Cm{i} = [obj.Airfoil.Cm{8}(obj.Airfoil.Alpha{8} < min(obj.Airfoil.Alpha{i})); obj.Airfoil.Cm{i}];
                    obj.Airfoil.Alpha{i} = [obj.Airfoil.Alpha{8}(obj.Airfoil.Alpha{8} < min(obj.Airfoil.Alpha{i})); obj.Airfoil.Alpha{i}];
                end
                if max(obj.Airfoil.Alpha{i}) < 180
                    obj.Airfoil.Cl{i} = [obj.Airfoil.Cl{i}; obj.Airfoil.Cl{8}(obj.Airfoil.Alpha{8} > max(obj.Airfoil.Alpha{i}))];
                    obj.Airfoil.Cd{i} = [obj.Airfoil.Cd{i}; obj.Airfoil.Cd{8}(obj.Airfoil.Alpha{8} > max(obj.Airfoil.Alpha{i}))];
                    obj.Airfoil.Cm{i} = [obj.Airfoil.Cm{i}; obj.Airfoil.Cm{8}(obj.Airfoil.Alpha{8} > max(obj.Airfoil.Alpha{i}))];
                    obj.Airfoil.Alpha{i} = [obj.Airfoil.Alpha{i}; obj.Airfoil.Alpha{8}(obj.Airfoil.Alpha{8} > max(obj.Airfoil.Alpha{i}))];
                end
            end
        end

        % viewModel opens the Simulink model and shows it at the requested
        % level or component.
        %
        %   Input arguments
        %     part - The part of the simulink model that the user wants to
        %     view
        %
        %   Output
        %     [-]
        %
        %   Behaviour
        %     This function opens the model in Simulink, showing the
        %     desired part:
        %       full   - The entire Simulink model
        %       torque - The model of the torque controller
        %       pitch  - The model of the torque controller
        %       damper - The model of the fore-aft damping controller
        %
        %   Called by
        %     FASTTool.ViewFullModel
        %     FASTTool.ViewTorqueController
        %     FASTTool.ViewPitchController
        %     FASTTool.ViewForeAftDamper
        function viewModel(~, part)
            if ~bdIsLoaded('FAST')
                load_system('FAST');
            end
            switch part
                case 'full'
                    open_system('FAST');
                    % Call open_system a second time, because the first
                    % time a simulink editor window is opened it shows the
                    % pitch controller
                    open_system('FAST');
                case 'torque'
                    open_system('FAST/Controller/Torque Control');
                case 'pitch'
                    open_system('FAST/Controller/Pitch Control');
                case 'damper'
                    open_system('FAST/Controller/Fore-Aft Tower Control');
            end
        end

    end
end
