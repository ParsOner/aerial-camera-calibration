%% run_calibration_model.m
% Script-based camera calibration + reporting/analysis in workspace.
% Compatible with your 4-model workflow:
%   - Radial 2 coeffs (4th order) or 3 coeffs (6th order)
%   - Tangential on/off
%   - Skew on/off
%
% Outputs:
%   calib  : struct with params/errors/analysis/score + helper tables
%   saves  : MAT file into resultsDir

clear; clc;

%% ========= USER INPUTS =========
imagesDir   = "C:\Users\PARSYSTEM\Desktop\EO photo";   % <-- CHANGE THIS
resultsDir  = "C:\Users\PARSYSTEM\Desktop\EO photo\Results";                  % <-- CHANGE THIS
squareSize  = 25;                                     % mm (or your unit)
worldUnits  = "mm";                                   % "mm" default per docs

% Choose model config here (change cfg for other 3 models)
% Model 1: 4th-order, NO tangential, NO skew
%cfg = makeModelCfg("M1_4th_noTan_noSkew", 2, false, false);

% Model 2: 4th-order, YES tangential, NO skew
%cfg = makeModelCfg("M2_4th_Tan_noSkew", 2, true,  false);

% Model 3: 6th-order, YES tangential, NO skew
%cfg = makeModelCfg("M3_6th_Tan_noSkew", 3, true,  false);

% Model 4: 6th-order, YES tangential, YES skew
 cfg = makeModelCfg("M4_6th_Tan_Skew",   3, true,  true);

%% ========= RUN =========
calib = calibrateAndAnalyze(imagesDir, squareSize, resultsDir, cfg, worldUnits);

% Convenience: push key outputs to workspace with readable names
params  = calib.params;
errors  = calib.errors;     % cameraCalibrationErrors object
report  = calib.report;     % struct of tables + summary
analysis = calib.analysis;  % reprojection statistics etc.
score   = calib.score;      % scalar

disp("Done. Workspace variables: calib, params, errors, report, analysis, score");


%% ===================== FUNCTIONS =====================

function cfg = makeModelCfg(name, numRadialCoeffs, estimateTangential, estimateSkew)
    cfg = struct();
    cfg.Name = string(name);
    cfg.NumRadialDistortionCoefficients = numRadialCoeffs;  % 2 or 3
    cfg.EstimateTangentialDistortion = logical(estimateTangential);
    cfg.EstimateSkew = logical(estimateSkew);
end

function calib = calibrateAndAnalyze(imagesDir, squareSize, resultsDir, cfg, worldUnits)

    arguments
        imagesDir (1,1) string
        squareSize (1,1) double {mustBePositive}
        resultsDir (1,1) string
        cfg (1,1) struct
        worldUnits (1,1) string = "mm"
    end

    if ~isfolder(imagesDir)
        error("imagesDir not found: %s", imagesDir);
    end
    if ~isfolder(resultsDir)
        mkdir(resultsDir);
    end

    % Load images
    imds = imageDatastore(imagesDir, ...
        "IncludeSubfolders", false, ...
        "FileExtensions", [".jpg",".jpeg",".png",".bmp",".tif",".tiff"]);

    if numel(imds.Files) < 2
        error("Need >=2 calibration images. Found: %d", numel(imds.Files));
    end

    % Detect checkerboard corners
    [imagePoints, boardSize, imagesUsed] = detectCheckerboardPoints(imds.Files);

    usedIdx = find(imagesUsed);
    if numel(usedIdx) < 2
        error("Checkerboard detected in <2 images. Improve images or pattern visibility.");
    end

    % Read one used image to get ImageSize
    I = imread(imds.Files{usedIdx(1)});
    imageSize = [size(I,1), size(I,2)];

    % World points of checkerboard corners
    % In R2025b docs, patternWorldPoints is used for checkerboard/circle grids.
    worldPoints = patternWorldPoints("checkerboard", boardSize, squareSize);

    % Estimate camera parameters with model options
    % Name-Value args are documented: EstimateSkew, NumRadialDistortionCoefficients,
    % EstimateTangentialDistortion, ImageSize, WorldUnits.  :contentReference[oaicite:2]{index=2}
    [params, imagesUsed2, errors] = estimateCameraParameters( ...
        imagePoints, worldPoints, ...
        ImageSize=imageSize, ...
        WorldUnits=worldUnits, ...
        EstimateSkew=cfg.EstimateSkew, ...
        NumRadialDistortionCoefficients=cfg.NumRadialDistortionCoefficients, ...
        EstimateTangentialDistortion=cfg.EstimateTangentialDistortion);

    % cameraCalibrationErrors contains IntrinsicsErrors and ExtrinsicsErrors. :contentReference[oaicite:3]{index=3}
    % IntrinsicsErrors properties include FocalLengthError, PrincipalPointError, RadialDistortionError,
    % TangentialDistortionError, SkewError. :contentReference[oaicite:4]{index=4}
    % ExtrinsicsErrors include RotationVectorsError, TranslationVectorsError. :contentReference[oaicite:5]{index=5}

    % Build report tables + analysis
    report  = buildReportTables(params, errors, cfg);
    analysis = analyzeReprojection(params);
    score   = computeScore(params, errors, analysis);

    % Pack everything
    calib = struct();
    calib.cfg = cfg;
    calib.imagesDir = imagesDir;
    calib.squareSize = squareSize;
    calib.worldUnits = worldUnits;

    calib.imageFiles = imds.Files;
    calib.imagesUsed_detect = imagesUsed;
    calib.imagesUsed_estimation = imagesUsed2;

    calib.params  = params;
    calib.errors  = errors;
    calib.report  = report;
    calib.analysis = analysis;
    calib.score   = score;

    % Print a compact console summary (workspace still has full structs)
    fprintf("\n=== %s ===\n", cfg.Name);
    fprintf("Used images (detected): %d / %d\n", nnz(imagesUsed), numel(imagesUsed));
    fprintf("Used images (estimation): %d\n", nnz(imagesUsed2));
    fprintf("Overall mean reprojection error (px): %.4f\n", analysis.overallMeanErrorPx);
    fprintf("Score (lower is better): %.6f\n\n", score);

    % Save results
    outMat = fullfile(resultsDir, cfg.Name + "_calibration.mat");
    save(outMat, "calib");
    fprintf("Saved: %s\n", outMat);
end

function report = buildReportTables(params, errors, cfg)
    % Extract intrinsics values from cameraParameters object safely
    [fxfy, cxcy, skewVal, kRad, pTan] = getIntrinsics(params);

    % Extract intrinsics uncertainties from errors.IntrinsicsErrors
    ie = errors.IntrinsicsErrors;
    fxFyErr = toRow(ie.FocalLengthError);
    cxyErr  = toRow(ie.PrincipalPointError);

    % These exist even if not estimated; if not estimated they may be 0 or absent.
    skewErr = safeProp(ie, "SkewError", 0);
    kRadErr = toRow(safeProp(ie, "RadialDistortionError", zeros(1,numel(kRad))));
    pTanErr = toRow(safeProp(ie, "TangentialDistortionError", zeros(1,numel(pTan))));

    % Intrinsics table
    intrTbl = table();
    intrTbl.Parameter = ["fx";"fy";"cx";"cy";"skew"];
    intrTbl.Value     = [fxfy(1); fxfy(2); cxcy(1); cxcy(2); skewVal];
    intrTbl.StdError  = [fxFyErr(1); fxFyErr(2); cxyErr(1); cxyErr(2); skewErr];

    % Distortion table (radial/tangential depending on model)
    % Radial length is 2 or 3 based on cfg
    distNames = strings(0,1);
    distVal   = [];
    distErr   = [];

    if ~isempty(kRad)
        for i=1:numel(kRad)
            distNames(end+1,1) = "k" + i;
            distVal(end+1,1)   = kRad(i);
            distErr(end+1,1)   = getIdxOrZero(kRadErr,i);
        end
    end

    if cfg.EstimateTangentialDistortion
        for i=1:numel(pTan)
            distNames(end+1,1) = "p" + i;
            distVal(end+1,1)   = pTan(i);
            distErr(end+1,1)   = getIdxOrZero(pTanErr,i);
        end
    else
        % still list p1 p2 as 0 (to keep your score formula consistent across models)
        distNames(end+1,1) = "p1"; distVal(end+1,1) = 0; distErr(end+1,1) = 0;
        distNames(end+1,1) = "p2"; distVal(end+1,1) = 0; distErr(end+1,1) = 0;
    end

    distTbl = table(distNames, distVal, distErr, ...
        'VariableNames', ["Parameter","Value","StdError"]);

    % Extrinsics tables (per-image)
    ex = errors.ExtrinsicsErrors;
    rvErr = safeProp(ex, "RotationVectorsError", []);
    tvErr = safeProp(ex, "TranslationVectorsError", []);

    % Params extrinsics:
    [rotVecs, transVecs] = getExtrinsics(params);

    extrTbl = table();
    extrTbl.ImageIndex = (1:size(rotVecs,1))';
    extrTbl.RotationVector = rotVecs;
    extrTbl.TranslationVector = transVecs;

    % Add uncertainties if sizes match
    if ~isempty(rvErr) && all(size(rvErr)==size(rotVecs))
        extrTbl.RotationVectorStdError = rvErr;
    else
        extrTbl.RotationVectorStdError = nan(size(rotVecs));
    end
    if ~isempty(tvErr) && all(size(tvErr)==size(transVecs))
        extrTbl.TranslationVectorStdError = tvErr;
    else
        extrTbl.TranslationVectorStdError = nan(size(transVecs));
    end

    report = struct();
    report.intrinsicsTable = intrTbl;
    report.distortionTable = distTbl;
    report.extrinsicsTable = extrTbl;

    % Bonus: print pretty standard error report (same style as docs)
    % displayErrors exists for cameraCalibrationErrors and its sub-objects. :contentReference[oaicite:6]{index=6}
    report.prettyTextIntrinsics = evalc("displayErrors(errors.IntrinsicsErrors, params);");
    report.prettyTextExtrinsics = evalc("displayErrors(errors.ExtrinsicsErrors, params);");
end

function analysis = analyzeReprojection(params)
    % ReprojectionErrors: M-by-2-by-N (points, xy, images)
    if isprop(params, "ReprojectionErrors")
        re = params.ReprojectionErrors;
    else
        error("params does not have ReprojectionErrors property in this MATLAB version.");
    end

    % Per-point magnitude
    mag = sqrt(re(:,1,:).^2 + re(:,2,:).^2);      % Mx1xN

    % Mean error per image (scalar per image)
    perImageMean = squeeze(mean(mag, 1, "omitnan"));   % 1xN -> Nx1 after squeeze
    perImageMean = perImageMean(:);

    % Overall mean reprojection error (scalar)
    overallMean = mean(perImageMean, "omitnan");

    % Overall mean component errors (x and y) (scalar each)
    xAbsMean = mean(abs(re(:,1,:)), "all", "omitnan");
    yAbsMean = mean(abs(re(:,2,:)), "all", "omitnan");

    analysis = struct();
    analysis.perImageMeanErrorPx = perImageMean;
    analysis.overallMeanErrorPx  = overallMean;
    analysis.meanAbsErrorX_Px    = xAbsMean;
    analysis.meanAbsErrorY_Px    = yAbsMean;

    % Optional: also compute RMS (more "statistical")
    analysis.rmsErrorX_Px = sqrt(mean(re(:,1,:).^2, "all", "omitnan"));
    analysis.rmsErrorY_Px = sqrt(mean(re(:,2,:).^2, "all", "omitnan"));
    analysis.rmsErrorMag_Px = sqrt(mean(mag.^2, "all", "omitnan"));
end

function score = computeScore(params, errors, analysis)
    % Implements your score idea in a reusable way.
    % We map MATLAB params -> fc, cc, kc(1..5), alpha_c + uncertainties.
    % Note: MATLAB uses RadialDistortion (k1..k3), TangentialDistortion (p1,p2), Skew.
    % We'll create a unified kc = [k1 k2 p1 p2 k3] style vector (Caltech-like ordering).

    epsDen = 1e-12;

    [fc, cc, alpha_c, kRad, pTan] = getIntrinsics(params);

    % build kc = [k1 k2 p1 p2 k3] with missing terms as 0
    k1 = getIdxOrZero(kRad,1);
    k2 = getIdxOrZero(kRad,2);
    k3 = getIdxOrZero(kRad,3);
    p1 = getIdxOrZero(pTan,1);
    p2 = getIdxOrZero(pTan,2);
    kc = [k1 k2 p1 p2 k3];

    ie = errors.IntrinsicsErrors;
    fc_err = toRow(ie.FocalLengthError);
    cc_err = toRow(ie.PrincipalPointError);
    alpha_err = safeProp(ie, "SkewError", 0);

    kRadErr = toRow(safeProp(ie, "RadialDistortionError", zeros(1,numel(kRad))));
    pTanErr = toRow(safeProp(ie, "TangentialDistortionError", zeros(1,numel(pTan))));

    % Align kc_err in same ordering
    k1e = getIdxOrZero(kRadErr,1);
    k2e = getIdxOrZero(kRadErr,2);
    k3e = getIdxOrZero(kRadErr,3);
    p1e = getIdxOrZero(pTanErr,1);
    p2e = getIdxOrZero(pTanErr,2);
    kc_err = [k1e k2e p1e p2e k3e];

    % Reprojection component term in your formula: sqrt(err(1)^2 + err(2)^2)
    % We'll use mean absolute x/y errors from analysis.
    err_xy = sqrt(analysis.meanAbsErrorX_Px^2 + analysis.meanAbsErrorY_Px^2);

    % Score sum
    score = abs(fc_err(1))/max(abs(fc(1)),epsDen) + ...
            abs(fc_err(2))/max(abs(fc(2)),epsDen) + ...
            abs(cc_err(1))/max(abs(cc(1)),epsDen) + ...
            abs(cc_err(2))/max(abs(cc(2)),epsDen);

    % kc terms (avoid division by zero)
    for i = 1:5
        score = score + abs(kc_err(i))/max(abs(kc(i)),epsDen);
    end

    % alpha/skew term
    score = score + abs(alpha_err)/max(abs(alpha_c),epsDen);

    % reprojection term
    score = score + err_xy;
end

function [fc, cc, skewVal, kRad, pTan] = getIntrinsics(params)
    % Robust extraction from cameraParameters object
    fc = [nan nan];
    cc = [nan nan];
    skewVal = 0;
    kRad = [];
    pTan = [];

    if isprop(params, "FocalLength")
        fc = toRow(params.FocalLength);
    end
    if isprop(params, "PrincipalPoint")
        cc = toRow(params.PrincipalPoint);
    end
    if isprop(params, "Skew")
        skewVal = params.Skew;
    end
    if isprop(params, "RadialDistortion")
        kRad = toRow(params.RadialDistortion);
    end
    if isprop(params, "TangentialDistortion")
        pTan = toRow(params.TangentialDistortion);
    end
end

function [rotVecs, transVecs] = getExtrinsics(params)
    % cameraParameters commonly provides RotationVectors and TranslationVectors.
    rotVecs = [];
    transVecs = [];

    if isprop(params, "RotationVectors")
        rotVecs = params.RotationVectors;
    else
        % fallback if only matrices are present
        if isprop(params, "RotationMatrices")
            R = params.RotationMatrices; % 3x3xN
            n = size(R,3);
            rotVecs = zeros(n,3);
            for i=1:n
                rotVecs(i,:) = rotationMatrixToVector(R(:,:,i));
            end
        end
    end

    if isprop(params, "TranslationVectors")
        transVecs = params.TranslationVectors;
    end
end

function v = rotationMatrixToVector(R)
    % Minimal helper: convert rotation matrix to rotation vector (Rodrigues)
    axang = rotm2axang(R);
    v = axang(1:3) * axang(4);
end

function out = safeProp(obj, propName, defaultVal)
    if isprop(obj, propName)
        out = obj.(propName);
    else
        out = defaultVal;
    end
end

function r = toRow(x)
    r = x(:).';
end

function val = getIdxOrZero(vec, idx)
    if isempty(vec) || numel(vec) < idx
        val = 0;
    else
        val = vec(idx);
    end
end
