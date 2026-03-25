%% analyze_calibration_workspace.m
% Uses workspace variable "calib" (recommended) or "params" (fallback).
% Produces:
%  1) Reprojection error plots (built-in + per-image mean)
%  2) Radial distortion model plots
%  3) Tangential distortion model vector field + magnitude map
%  4) Combined distortion field (optional)
%
% Outputs figures and optionally saves PNGs.


%% ====== INPUTS ======
% If you saved MAT from previous scr11ipt:
load("C:\Users\PARSYSTEM\Desktop\EO photo\Results\M4_6th_Tan_Skew_calibration.mat");  % gives calib

% If calib exists, use it. Otherwise try params.
if exist("calib","var")
    params = calib.params;
    cfgName = calib.cfg.Name;
    if resultsDir == ""
        resultsDir = fullfile(pwd, "results_analysis");
    end
else
    if ~exist("params","var")
        error("Workspace needs either 'calib' or 'params'. Load your *_calibration.mat or run calibration first.");
    end
    cfgName = "CalibAnalysis";
    if resultsDir == ""
        resultsDir = fullfile(pwd, "results_analysis");
    end
end

%% ====== 1) REPROJECTION ERRORS ======
fig1 = figure("Name","Reprojection Errors (Built-in)","NumberTitle","off");
try
    % MATLAB's standard plot (if available)
    showReprojectionErrors(params);
catch
    % Fallback if showReprojectionErrors not available
    warning("showReprojectionErrors not available. Using custom scatter fallback.");
    customReprojScatter(params);
end

fig2 = figure("Name","Mean Reprojection Error per Image","NumberTitle","off");
plotPerImageMeanReproj(params);

%% ====== 2) DISTORTION MODELS (RADIAL + TANGENTIAL) ======
% We will visualize distortion in NORMALIZED image coordinates (x,y),
% where r = sqrt(x^2 + y^2).
%
% Radial model:  x' = x * L(r), y' = y * L(r)
% Tangential:    dx_t = 2*p1*x*y + p2*(r^2 + 2*x^2)
%               dy_t = p1*(r^2 + 2*y^2) + 2*p2*x*y

fig3 = figure("Name","Radial Distortion Model","NumberTitle","off");
plotRadialModel(params);

fig4 = figure("Name","Tangential Distortion Model","NumberTitle","off");
plotTangentialModel(params);

fig5 = figure("Name","Combined Distortion Field","NumberTitle","off");
plotCombinedDistortionField(params);
disp("Analysis done.");


%% ===================== Helper Functions =====================

function customReprojScatter(params)
    re = params.ReprojectionErrors; % M-by-2-by-N
    x = re(:,1,:); y = re(:,2,:);
    x = x(:); y = y(:);
    scatter(x, y, 6, "filled");
    grid on;
    axis equal;
    xlabel("Reprojection error X (px)");
    ylabel("Reprojection error Y (px)");
    title("Reprojection error scatter (fallback)");
end

function plotPerImageMeanReproj(params)
    re = params.ReprojectionErrors; % M-by-2-by-N
    mag = sqrt(re(:,1,:).^2 + re(:,2,:).^2);
    perImageMean = squeeze(mean(mag,1,"omitnan"));
    perImageMean = perImageMean(:);
    plot(perImageMean, "-o");
    grid on;
    xlabel("Image index");
    ylabel("Mean reprojection error magnitude (px)");
    title("Mean reprojection error per image");
end

function plotRadialModel(params)
    [fxfy, imageSize, kRad] = getIntrinsicsAndImage(params);

    % Build k1,k2,k3 (k3=0 if not estimated)
    k1 = getIdxOrZero(kRad,1);
    k2 = getIdxOrZero(kRad,2);
    k3 = getIdxOrZero(kRad,3);

    % Estimate r_max in normalized coords roughly at image corners
    fx = fxfy(1); fy = fxfy(2);
    w = imageSize(2); h = imageSize(1);
    % normalized corner ~ ((w/2)/fx, (h/2)/fy)
    rMax = sqrt(((w/2)/fx)^2 + ((h/2)/fy)^2);

    r = linspace(0, rMax, 400);
    L = 1 + k1*r.^2 + k2*r.^4 + k3*r.^6;

    % Radial displacement: r_d - r  where r_d = r*L
    rd = r .* L;
    dr = rd - r;

    % Plot L(r) and dr
    tiledlayout(2,1,"Padding","compact","TileSpacing","compact");

    nexttile;
    plot(r, L, "LineWidth", 1.2);
    grid on;
    xlabel("r (normalized)");
    ylabel("L(r)");
    title("Radial scale factor L(r) = 1 + k1 r^2 + k2 r^4 + k3 r^6");

    nexttile;
    plot(r, dr, "LineWidth", 1.2);
    grid on;
    xlabel("r (normalized)");
    ylabel("r_d - r (normalized)");
    title("Radial displacement (r_d - r)");
end

function plotTangentialModel(params)
    [fxfy, imageSize, ~, pTan] = getIntrinsicsAndImage(params);

    p1 = getIdxOrZero(pTan,1);
    p2 = getIdxOrZero(pTan,2);

    % === Grid sınırı (normalized coords) ===
    fx = fxfy(1); fy = fxfy(2);
    w = imageSize(2); h = imageSize(1);
    rMax = sqrt(((w/2)/fx)^2 + ((h/2)/fy)^2);

    n = 25;
    lim = 0.9*rMax;
    xv = linspace(-lim, lim, n);
    yv = linspace(-lim, lim, n);
    [X,Y] = meshgrid(xv,yv);
    R2 = X.^2 + Y.^2;

    % === Tangential distortion displacement ===
    dxt = 2*p1.*X.*Y + p2.*(R2 + 2*X.^2);
    dyt = p1.*(R2 + 2*Y.^2) + 2*p2.*X.*Y;

    mag = sqrt(dxt.^2 + dyt.^2);
    maxMag = max(mag(:));

    % === Görsel ölçekleme (SADECE görselleştirme) ===
    % hedef: okların max uzunluğu eksen genişliğinin ~%8'i olsun
    targetArrow = 0.08 * (2*lim);
    if maxMag < 1e-18
        scale = 1;  % zaten 0
    else
        scale = targetArrow / maxMag;
    end

    tiledlayout(1,2,"Padding","compact","TileSpacing","compact");

    nexttile;
    quiver(X, Y, dxt*scale, dyt*scale, 0); % 0 => autoscale kapalı, biz scale ediyoruz
    axis equal; grid on;
    xlabel("x (normalized)"); ylabel("y (normalized)");
    title(sprintf("Tangential vector field (scaled)  p1=%.3g, p2=%.3g", p1, p2));
    subtitle(sprintf("visual scale = %.3g  |d|_max=%.3g (normalized)", scale, maxMag));

    nexttile;
    imagesc(xv, yv, mag);
    axis image; set(gca,"YDir","normal");
    colorbar;
    xlabel("x (normalized)"); ylabel("y (normalized)");
    title("|Tangential distortion| magnitude");

    % CLim'i düzgün ayarla (mag >= 0 olmalı)
    if maxMag < 1e-18
        clim([0 1]);   % tamamen sıfırsa da boş görünmesin
    else
        clim([0 maxMag]);
    end
end


function plotCombinedDistortionField(params)
    [fxfy, imageSize, kRad, pTan] = getIntrinsicsAndImage(params);

    k1 = getIdxOrZero(kRad,1);
    k2 = getIdxOrZero(kRad,2);
    k3 = getIdxOrZero(kRad,3);
    p1 = getIdxOrZero(pTan,1);
    p2 = getIdxOrZero(pTan,2);

    fx = fxfy(1); fy = fxfy(2);
    w = imageSize(2); h = imageSize(1);
    rMax = sqrt(((w/2)/fx)^2 + ((h/2)/fy)^2);

    n = 25;
    lim = 0.9*rMax;
    xv = linspace(-lim, lim, n);
    yv = linspace(-lim, lim, n);
    [X,Y] = meshgrid(xv,yv);
    R2 = X.^2 + Y.^2;

    % Radial
    L = 1 + k1*R2 + k2*(R2.^2) + k3*(R2.^3);
    Xr = X.*L; Yr = Y.*L;
    dxr = Xr - X; dyr = Yr - Y;

    % Tangential
    dxt = 2*p1.*X.*Y + p2.*(R2 + 2*X.^2);
    dyt = p1.*(R2 + 2*Y.^2) + 2*p2.*X.*Y;

    % Combined displacement
    dx = dxr + dxt;
    dy = dyr + dyt;

    mag = sqrt(dx.^2 + dy.^2);

    tiledlayout(1,2,"Padding","compact","TileSpacing","compact");

    nexttile;
    quiver(X, Y, dx, dy, "AutoScale","on");
    axis equal; grid on;
    xlabel("x (normalized)"); ylabel("y (normalized)");
    title("Combined distortion field (radial + tangential)");

    nexttile;
    imagesc(xv, yv, mag);
    axis image; set(gca,"YDir","normal");
    colorbar;
    xlabel("x (normalized)"); ylabel("y (normalized)");
    title("|Combined distortion| magnitude");
end

function [fxfy, imageSize, kRad, pTan] = getIntrinsicsAndImage(params)
    if isprop(params,"FocalLength")
        fxfy = params.FocalLength(:).';
    else
        error("params.FocalLength not found");
    end
    if isprop(params,"ImageSize")
        imageSize = params.ImageSize;
    else
        % If ImageSize isn't stored, infer from principal point (fallback)
        if isprop(params,"PrincipalPoint")
            cc = params.PrincipalPoint;
            imageSize = [round(2*cc(2)), round(2*cc(1))];
        else
            error("params.ImageSize not found and cannot infer.");
        end
    end

    if isprop(params,"RadialDistortion")
        kRad = params.RadialDistortion(:).';
    else
        kRad = [];
    end
    if isprop(params,"TangentialDistortion")
        pTan = params.TangentialDistortion(:).';
    else
        pTan = [];
    end
end

function val = getIdxOrZero(vec, idx)
    if isempty(vec) || numel(vec) < idx
        val = 0;
    else
        val = vec(idx);
    end
end
