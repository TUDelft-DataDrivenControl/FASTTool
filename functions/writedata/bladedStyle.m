function [x, y, z, u, v, w] = bladedStyle(Wind, U, H, R, startupTime)
% bladedStyle generates wind data needed to create a wind file in Bladed
% style. Bladed style wind files are used for deterministic wind fields.
% Wind data can be generated for different types of wind conditions.
% 
%
%   Syntax
%     [x, y, z, u, v, w] = bladedStyle(Wind, U, H, R, startupTime)
%
%   Input arguments
%     Wind - Data structure with wind data. Wind.Type determines for which
%     wind conditions the wind data needs to be generated.
%     U - 10-Minute average wind speed. For stepped wind U is the first and
%     lowest wind speed in the series of steps
%     H - Hub height
%     R - Rotor radius
%     startupTime - Time period for the transient of the startup of the
%     simulation for which the results will be discarded.
%
%   Output
%     [x, y, z, u, v, w] - Data needed to create a Bladed style wind file
%
%   Behaviour
%     For most wind types, the wind data will be generated for the average
%     wind speed U. However, for stepped wind the steps in the wind speed
%     are taken from Wind.Speed.
%     The wind conditions in the startup period will be the same as those
%     that are generated for the start of the simulation for which the
%     results will be kept.
%
%   Called by
%     WindTurbineClass.simulate

% ToDo: This function needs further commenting. The person who wrote it is
% no longer available for the development of FASTTool.

switch Wind.Type
    case 2 % Type 2: Stepped wind

        % Wind profile
        t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
        if rem(length(t),2) ~= 0
            t = [t, max(t)+Wind.dt];
        end
        x = U*t;
        y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
        z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
        windIndex = 1 + max(0,floor(length(Wind.Speed)*(t-startupTime)/(Wind.T-startupTime)));
        windIndex = min(windIndex, length(Wind.Speed));
        u = Wind.Speed(windIndex);

        u = repmat(u(:),[1,Wind.Ny,Wind.Nz]);
        v = zeros(size(u));
        w = zeros(size(u));

    case 7 % Type 7: Extreme wind shear (EWS)

        % IEC class
        if Wind.Class(2) == 1
            Iref = 0.16;
        elseif Wind.Class(2) == 2
            Iref = 0.14;
        elseif Wind.Class(2) == 3
            Iref = 0.12;
        end
        if H >= 60
            Lambda1 = 0.7*H;
        else
            Lambda1 = 42;
        end

        % NTM values
        sigma1 = Iref*(0.75*U+5.6);
        alpha = 0.2;

        % EWS parameters
        beta = 6.4;
        Tg = 12;
        T0 = startupTime + Wind.EventTime;

        % Wind profile
        t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
        if rem(length(t),2) ~= 0
            t = [t, t+Wind.dt];
        end
        x = U*t;
        y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
        z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
        [T,~,Z] = ndgrid(t,y,z);
        u = U*(Z/H).^alpha + (Z-H)/(2*R) * (2.5 + 0.2*beta*sigma1*(2*R/Lambda1).^0.25) .* (1 -cos(2*pi*(T-T0)/Tg));
        u(T<T0) = U*(Z(T<T0)/H).^alpha;
        u(T>T0+Tg) = U*(Z(T>T0+Tg)/H).^alpha;
        v = zeros(size(u));
        w = zeros(size(u));

    case 9 % Type 9: Extreme operating gust (EOG)

        % IEC class
        if Wind.Class(1) == 1
            Vref = 50;
        elseif Wind.Class(1) == 2
            Vref = 42.5;
        elseif Wind.Class(1) == 3
            Vref = 37.5;
        end
        if Wind.Class(2) == 1
            Iref = 0.16;
        elseif Wind.Class(2) == 2
            Iref = 0.14;
        elseif Wind.Class(2) == 3
            Iref = 0.12;
        end
        if H >= 60
            Lambda1 = 0.7*H;
        else
            Lambda1 = 42;
        end

        % EWM values
        Ve50 = 1.4*Vref;
        Ve1 = 0.8*Ve50;

        % NTM values
        alpha = 0.2;
        sigma1 = Iref*(0.75*U+5.6);

        % EOG parameters
        Vgust = min([1.35*(Ve1-U), 3.3*sigma1/(1+0.1*2*R/Lambda1)]);
        Tg = 10.5;
        T0 =  startupTime + Wind.EventTime;

        % Wind profile
        t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
        if rem(length(t),2) ~= 0
            t = [t, t+Wind.dt];
        end
        x = U*t;
        y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
        z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
        [T,~,Z] = ndgrid(t,y,z);
        u = U*(Z/H).^alpha - 0.37*Vgust*sin(3*pi*(T-T0)/Tg) .* (1-cos(2*pi*(T-T0)/Tg));
        u(T<T0) = U*(Z(T<T0)/H).^alpha;
        u(T>T0+Tg) = U*(Z(T>T0+Tg)/H).^alpha;
        v = zeros(size(u));
        w = zeros(size(u));

    case 10 % Type 10: Extreme direction change (EDC)

        % IEC class
        if Wind.Class(2) == 1
            Iref = 0.16;
        elseif Wind.Class(2) == 2
            Iref = 0.14;
        elseif Wind.Class(2) == 3
            Iref = 0.12;
        end
        if H >= 60
            Lambda1 = 0.7*H;
        else
            Lambda1 = 42;
        end

        % NTM values
        alpha = 0.2;
        sigma1 = Iref*(0.75*U+5.6);

        % EDC parameters
        T0 =  startupTime + Wind.EventTime;
        Tg = 6;

        % Wind profile
        t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
        if rem(length(t),2) ~= 0
            t = [t, t+Wind.dt];
        end
        x = U*t;
        y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
        z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
        [T,~,Z] = ndgrid(t,y,z);
        thetae = 4*atan(sigma1/(U*(1+0.1*2*R/Lambda1)));
        theta = 0.5*thetae * (1-cos(pi*(T-T0)/Tg));
        theta(T < T0) = 0;
        theta(T > T0+Tg) = thetae;
        u = U*(Z/H).^alpha .* cos(theta);
        v = U*(Z/H).^alpha .* sin(theta);
        w = zeros(size(u));

    case 11 % Type 11: Extreme coherent gust (ECG)

        % NTM values
        alpha = 0.2;

        % ECG parameters
        T0 =  startupTime + Wind.EventTime;
        Tg = 10;
        Vcg = 15;

        % Wind profile
        t = -Wind.Ly/U/2:Wind.dt:(Wind.T + Wind.Ly/U/2);
        if rem(length(t),2) ~= 0
            t = [t, t+Wind.dt];
        end
        x = U*t;
        y = linspace(Wind.Ly/2,-Wind.Ly/2,Wind.Ny);
        z = H + linspace(-Wind.Lz/2,Wind.Lz/2,Wind.Nz);
        [T,~,Z] = ndgrid(t,y,z);
        V = U*(Z/H).^alpha + 0.5*Vcg*(1-cos(pi*(T-T0)/Tg));
        V(T < T0) = U*(Z(T<T0)/H).^alpha;
        V(T > T0+Tg) = U*(Z(T>T0+Tg)/H).^alpha + Vcg;
        if U < 4
            thetacg = pi;
        else
            thetacg = 4*pi/U;
        end
        theta = 0.5*thetacg*(1-cos(pi*(T-T0)/Tg));
        theta(T < T0) = 0;
        theta(T > T0+Tg) = thetacg;
        u = V .* cos(theta);
        v = V .* sin(theta);
        w = zeros(size(u));
end