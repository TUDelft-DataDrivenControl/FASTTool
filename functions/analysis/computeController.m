function Controller = computeController(pitch, Control, controllerElements)
% computeController determines transfer functions of the pitch controller.
%
%   Syntax
%     Controller = computeController(pitch, Control, controllerElements)
%
%   Input arguments
%     pitch - Array with pitch angles for which the transfer function needs
%     to be determined and for which control settings are available
%     Control - Structure with control settings
%     controllerElements - String array, indicating which elements of the
%     controller need to be included in the transfer function:
%       "LPF" - Low pass filter
%       "PI" - Gains of the PI controller
%       "Notch" - The notch filters
%
%   Output
%     Controller - Transfer functions for all pitch angles, obtained by
%     multiplication of the transfer functions of the included elements
%
%   Behaviour
%     When the controller is set to be gain scheduled, the control
%     settings for each of the input pitch angles will be obtained by
%     interpolation of the gain-scheduled control settings.
%     When the controller is set to be constant, the constant control
%     settings will be used for all of the input pitch angles.
%
%   Called by
%     WindTurbineClass.analyseControl

    % This conditional code sets the controller settings. If the condition
    % is met, the gain-scheduled control settings will be interpolated to
    % the pitch angles at which the transfer functions are required.
    % Otherwise, the constant controller settings are used for all pitch
    % angles.
    %
    % The pitch angles of the gain scheduling are chosen by the user, and
    % the pitch angles at which the transfer functions are needed are
    % determined by the provided linear model. These two sets of pitch
    % angles don't have to be the same, and therefore the interpolation is
    % used.
    %
    % The interpolation is the same as the one used by the actual
    % controller in Simulink.
    if Control.Pitch.Scheduled
        x = Control.Pitch.ScheduledPitchAngles;
        xq = pitch;
        LowPassCutOffFreq = interp1(x, Control.Pitch.LowPassCutOffFreqGS, clip(xq, min(x), max(x)),'linear','extrap');
        Kp = interp1(x, Control.Pitch.KpGS, clip(xq, min(x), max(x)),'linear','extrap');
        Ki = interp1(x, Control.Pitch.KiGS, clip(xq, min(x), max(x)),'linear','extrap');
        Notch_beta1 = interp1(x, Control.Pitch.Notch_beta1GS, clip(xq, min(x), max(x)),'linear','extrap');
        Notch_beta2 = interp1(x, Control.Pitch.Notch_beta2GS, clip(xq, min(x), max(x)),'linear','extrap');
        Notch_wn = interp1(x, Control.Pitch.Notch_wnGS, clip(xq, min(x), max(x)),'linear','extrap');
        Notch2_beta1 = interp1(x, Control.Pitch.Notch2_beta1GS, clip(xq, min(x), max(x)),'linear','extrap');
        Notch2_beta2 = interp1(x, Control.Pitch.Notch2_beta2GS, clip(xq, min(x), max(x)),'linear','extrap');
        Notch2_wn = interp1(x, Control.Pitch.Notch2_wnGS, clip(xq, min(x), max(x)),'linear','extrap');
    else
        onesPitch = ones(length(pitch));
        LowPassCutOffFreq = Control.Pitch.LowPassCutOffFreq * onesPitch;
        Kp = Control.Pitch.Kp * onesPitch;
        Ki = Control.Pitch.Ki * onesPitch;
        Notch_beta1 = Control.Pitch.Notch_beta1 * onesPitch;
        Notch_beta2 = Control.Pitch.Notch_beta2 * onesPitch;
        Notch_wn = Control.Pitch.Notch_wn * onesPitch;
        Notch2_beta1 = Control.Pitch.Notch2_beta1 * onesPitch;
        Notch2_beta2 = Control.Pitch.Notch2_beta2 * onesPitch;
        Notch2_wn = Control.Pitch.Notch2_wn * onesPitch;
    end

    % The next code first creates an identity transfer function and then
    % multiplies it with the transfer functions of each of the controller
    % elements that need to be included. This multiplication of transfer
    % functions represents the total transfer function, because it is used
    % for a linear-system analysis.
    Controller = tf(1,1)*ones(1,length(pitch));
    for i = 1:length(pitch)
        if any(contains(controllerElements, "LPF"))
            if Control.Pitch.LowPassOrder == 1
                Controller(1,i) = Controller(1,i)...
                    *tf(LowPassCutOffFreq(i),...
                    [1 LowPassCutOffFreq(i)]);
            else
                Controller(1,i) = Controller(1,i)...
                    *tf(LowPassCutOffFreq(i)^2,...
                    [1 2/sqrt(2)*LowPassCutOffFreq(i) LowPassCutOffFreq(i)^2]);
            end
        end    
        if any(contains(controllerElements, "PI"))
            if all([Kp(i) Ki(i)] == 0)
                Controller(1,i) = Controller(1,i);
            else
                Controller(1,i) = Controller(1,i)...
                    *tf([Kp(i) Ki(i)], [1 0]);
            end
        end
        if any(contains(controllerElements, "Notch"))
            if any([Notch_beta1(i) Notch_beta2(i) Notch_wn(i)] == 0) ...
                    && any([Notch2_beta1(i) Notch2_beta2(i) Notch2_wn(i)] == 0)
                Controller(1,i) = Controller(1,i);
            else
                Controller(1,i) = Controller(1,i)...
                    *tf([1 2*Notch_beta1(i)*Notch_wn(i) Notch_wn(i)^2],...
                    [1 2*Notch_beta2(i)*Notch_wn(i) Notch_wn(i)^2])...
                    *tf([1 2*Notch2_beta1(i)*Notch2_wn(i) Notch2_wn(i)^2],...
                    [1 2*Notch2_beta2(i)*Notch2_wn(i) Notch2_wn(i)^2]);
            end
        end
    end
