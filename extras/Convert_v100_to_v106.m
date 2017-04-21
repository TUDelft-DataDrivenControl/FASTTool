
%% Ask for file names
[OldFileName,OldPathName] = uigetfile('*.mat', 'Select old project file');
[NewFileName,NewPathName] = uiputfile('*.mat', 'Name for new project file');

%% Open old project file
load([OldPathName, OldFileName])

%% Add baseline values for new fields
try
    Control.DT = 0.008;
    Control.LPFCutOff = 3;
    Control.ForeAft.MaxPitchAmplitude = 5;
    Control.Pitch.Scheduled = 0;
    Control.Pitch.ScheduledPitchAngles = transpose(linspace(1,25,14)*pi/180);
    Control.Pitch.KpGS = zeros(14,1);
    Control.Pitch.KiGS = zeros(14,1);
    Control.Pitch.Kp = -0.003;
    Control.Pitch.Ki = -0.001;
    Control.Pitch.StartupPitch = 45;
    Control.Pitch.StartupPitchRate = 1;
    Control.Pitch.StartupSpeed = 0.25;
    Control.Torque.SpeedB2 = 1073;
    Control.Torque.Min = 200;
    Control.Brake.Delay = 1.5;
    CertificationSettings.Mode.Actiontime = 30;
catch
    disp('An error occured in the "Add baseline values for new fields" section')
end

%% Rename fields
try
    Airfoil.StallAngle1 = Airfoil.StallAngle;
    Airfoil.StallAngle2 = -Airfoil.StallAngle;
    Airfoil.CritCn1 = Airfoil.StallPosCn;
    Airfoil.CritCn2 = Airfoil.StallNegCn;
    CertificationSettings.Wind.ECG = CertificationSettings.Wind.ECD;
catch
    disp('An error occured in the "Rename fields" section')
end

%% Remove old fields
try
    Control.Pitch = rmfield(Control.Pitch,'Ti');
    Drivetrain.Generator = rmfield(Drivetrain.Generator,'Speed');
    Drivetrain.Generator = rmfield(Drivetrain.Generator,'Torque');
    Drivetrain.Generator = rmfield(Drivetrain.Generator,'Power');
    Airfoil = rmfield(Airfoil,'StallAngle');
    Airfoil = rmfield(Airfoil,'StallPosCn');
    Airfoil = rmfield(Airfoil,'StallNegCn');
    Airfoil = rmfield(Airfoil,'ZeroCnAngle');
    CertificationSettings.Mode = rmfield(CertificationSettings.Mode,'Stop');
    CertificationSettings.Wind = rmfield(CertificationSettings.Wind,'ECD');
catch
    disp('An error occured in the "Remove old fields" section')
end

%% Correct for possible duplicate angles of attack in the original file
for i = 1:length(Airfoil.Alpha)
    if length(Airfoil.Alpha{i}) - length(unique(Airfoil.Alpha{i})) ~= 0
        [Airfoil.Alpha{i},j] = unique(Airfoil.Alpha{i});
        Airfoil.Cl{i} = Airfoil.Cl{i}(j);
        Airfoil.Cd{i} = Airfoil.Cd{i}(j);
        Airfoil.Cm{i} = Airfoil.Cm{i}(j);
    end
end

%% Save new project file
save([NewPathName, NewFileName], ...
    'Blade', ...
    'Airfoil', ...
    'Tower', ...
    'Nacelle', ...
    'Drivetrain', ...
    'Control', ...
    'CertificationSettings', ...
    'Appearance')