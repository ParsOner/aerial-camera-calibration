%% sweep_calibration_image_count_rev.m
% Sweep number of calibration images N for ONE selected camera model.
% Runs K random subsets for each N, estimates camera parameters, and computes:
%   - Score (assignment-style, but numerically robust near-zero parameters)
%   - Mean reprojection error magnitude (pixels)
%   - RMS reprojection components err1/err2 (pixels)
%   - Intrinsics & distortion parameters (optional)
%
% REQUIREMENTS:
%   - Computer Vision Toolbox (estimateCameraParameters, detectCheckerboardPoints)
%
% Written to be drop-in runnable: configure the USER SETTINGS section.

clear; clc;

%% ========= USER SETTINGS =========
imagesDir   = "C:\Users\PARSYSTEM\Desktop\EO photo";            % <-- change
resultsDir  = "C:\Users\PARSYSTEM\Desktop\EO photo\Results";    % <-- change

squareSize  = 25;     % checkerboard square size (mm or your unit)
worldUnits  = "mm";

% Choose N values to test (must be <= number of valid images)
Ns          = [5 8 10 15 20];

K           = 10;      % number of random subsets per N
rngSeed     = 42;      % fixed seed for reproducibility (change for different sweeps)

% ----- MODEL SWITCHES (your current case: sixth order + tan + skew) -----
radialOrder     = 3;       % 2 -> [k1 k2], 3 -> [k1 k2 k3] (k3 = r^6 term)
estimateTangential = true; % yes/no
estimateSkew       = true; % yes/no
numRadial = radialOrder;   % MATLAB uses 2 or 3

% ----- Score robustness settings (prevents blow-ups near zero) -----
scoreCfg = struct();
scoreCfg.alpha_floor = 1.0;                         % pixels
scoreCfg.kc_floor    = [1e-3 1e-3 5e-4 5e-4 1e-3];   % [k1 k2 p1 p2 k3]
scoreCfg.epsDen      = 1e-12;                        % tiny epsilon for intrinsics denominators

% Saving
saveFigures = true;

%% ========= INPUT VALIDATION / SETUP =========
if ~isfolder(imagesDir)
    error("imagesDir does not exist: %s", imagesDir);
end
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

% Discover images
exts = ["*.png","*.jpg","*.jpeg","*.tif","*.tiff","*.bmp"];
files = [];
for e = exts
    files = [files; dir(fullfile(imagesDir, e))]; %#ok<AGROW>
end
if isempty(files)
    error("No images found in %s", imagesDir);
end

filePaths = string(fullfile({files.folder}, {files.name}))';
numFiles  = numel(filePaths);

fprintf("Found %d image files.\n", numFiles);

% Detect checkerboard points in all images
imagePointsAll = {};
imagesUsed = false(numFiles,1);
boardSize = [];
for i = 1:numFiles
    I = imread(filePaths(i));
    [imagePoints, bsz] = detectCheckerboardPoints(I);
    if ~isempty(imagePoints)
        imagePointsAll{i} = imagePoints; %#ok<SAGROW>
        imagesUsed(i) = true;
        if isempty(boardSize)
            boardSize = bsz;
        else
            % Keep only images with consistent board size
            if any(bsz ~= boardSize)
                imagesUsed(i) = false;
            end
        end
    end
end

usedIdx = find(imagesUsed);
if isempty(usedIdx)
    error("No valid checkerboard detections found.");
end

fprintf("Valid detections in %d/%d images.\n", numel(usedIdx), numFiles);

% World points for detected checkerboard
worldPoints = generateCheckerboardPoints(boardSize, squareSize);

% Clamp Ns to valid range
maxN = numel(usedIdx);
Ns = Ns(Ns <= maxN);
if isempty(Ns)
    error("All Ns exceed number of usable images (%d).", maxN);
end

%% ========= SWEEP =========
rng(rngSeed);

runRows = [];  % table rows (we'll build incrementally)
rowCount = 0;

for n = Ns
    fprintf("\n=== N = %d (K=%d) ===\n", n, K);

    for k = 1:K
        rowCount = rowCount + 1;

        % Random subset of usable images
        if n == numel(usedIdx)
            subset = usedIdx;  % deterministic (all images)
        else
            subset = usedIdx(randperm(numel(usedIdx), n));
        end

        % Stack image points into Mx2xN format for estimateCameraParameters
        imagePoints = zeros(size(imagePointsAll{subset(1)},1), 2, n);
        for j = 1:n
            imagePoints(:,:,j) = imagePointsAll{subset(j)};
        end

        % Estimate camera parameters
        opts = {'EstimateSkew', logical(estimateSkew), ...
                'EstimateTangentialDistortion', logical(estimateTangential), ...
                'NumRadialDistortionCoefficients', numRadial, ...
                'WorldUnits', worldUnits};

        try
            [params, imagesUsedByEstimator, estimationErrors] = estimateCameraParameters( ...
                imagePoints, worldPoints, opts{:});
        catch ME
            warning("estimateCameraParameters failed at N=%d, k=%d: %s", n, k, ME.message);

            runRows = [runRows; makeFailRow(rowCount, n, k, subset, numRadial, estimateTangential, estimateSkew)]; %#ok<AGROW>
            continue;
        end

        % Compute reprojection errors
        % reprojectionErrors: Mx2xN
        reprojectionErrors = params.ReprojectionErrors;
        [meanReprojMag, stdReprojMag, err1_rms, err2_rms] = reprojMetrics(reprojectionErrors);

        % Compute Score (robust assignment-style)
        [score, parts, alpha, kc, kc_err, alpha_err] = computeScoreLikeAssignmentRobust( ...
            params, estimationErrors, err1_rms, err2_rms, scoreCfg);

        % Save row
        runRows = [runRows; makeRunRow(rowCount, n, k, subset, ...
            score, meanReprojMag, stdReprojMag, err1_rms, err2_rms, ...
            parts, alpha, alpha_err, kc, kc_err, ...
            numRadial, estimateTangential, estimateSkew, imagesUsedByEstimator)]; %#ok<AGROW>
    end
end

runsTbl = struct2table(runRows);

% Summary by N
summaryTbl = groupsummary(runsTbl, "N", ["mean","std"], ...
    ["Score","reprojMeanMag","err1","err2"]);

% Save
stamp = datestr(now, "yyyymmdd_HHMMSS");
modelTag = sprintf("M4_%s_%s_%s", ...
    ternary(radialOrder==3,"6th","4th"), ...
    ternary(estimateTangential,"Tan","NoTan"), ...
    ternary(estimateSkew,"Skew","NoSkew"));

runsCSV    = fullfile(resultsDir, modelTag + "_sweepN_runs_" + stamp + ".csv");
summaryCSV = fullfile(resultsDir, modelTag + "_sweepN_summary_" + stamp + ".csv");
writetable(runsTbl, runsCSV);
writetable(summaryTbl, summaryCSV);

fprintf("\nSaved:\n  %s\n  %s\n", runsCSV, summaryCSV);

%% ========= PLOTS =========
if saveFigures
    % Score vs N (mean±std)
    fig1 = figure('Name','Score vs N');
    hold on;
    errorbar(summaryTbl.N, summaryTbl.mean_Score, summaryTbl.std_Score, 'o-');
    xlabel('N (number of images)');
    ylabel('Score (robust assignment metric)');
    title("Score vs N  (" + modelTag + ")");
    grid on;
    saveas(fig1, fullfile(resultsDir, modelTag + "_Score_vs_N_" + stamp + ".png"));

    % Reprojection magnitude vs N
    fig2 = figure('Name','Reprojection mean magnitude vs N');
    hold on;
    errorbar(summaryTbl.N, summaryTbl.mean_reprojMeanMag, summaryTbl.std_reprojMeanMag, 'o-');
    xlabel('N (number of images)');
    ylabel('Mean reprojection error magnitude (px)');
    title("Mean Reprojection Error vs N  (" + modelTag + ")");
    grid on;
    saveas(fig2, fullfile(resultsDir, modelTag + "_Reproj_vs_N_" + stamp + ".png"));

    % err1/err2 vs N
    fig3 = figure('Name','err1 err2 vs N');
    hold on;
    errorbar(summaryTbl.N, summaryTbl.mean_err1, summaryTbl.std_err1, 'o-');
    errorbar(summaryTbl.N, summaryTbl.mean_err2, summaryTbl.std_err2, 'o-');
    xlabel('N (number of images)');
    ylabel('RMS reprojection component (px)');
    legend('err1 (RMS x)','err2 (RMS y)','Location','best');
    title("RMS Reprojection Components vs N  (" + modelTag + ")");
    grid on;
    saveas(fig3, fullfile(resultsDir, modelTag + "_ErrComponents_vs_N_" + stamp + ".png"));
end

%% ========= DONE =========
disp("Done.");

%% ===================== LOCAL FUNCTIONS =====================

function row = makeFailRow(rowCount, n, k, subset, numRadial, estTan, estSkew)
    row = struct();
    row.rowId = rowCount;
    row.N = n;
    row.k = k;
    row.subsetIdx = join(string(subset), ",");
    row.NumRadial = numRadial;
    row.EstTan = logical(estTan);
    row.EstSkew = logical(estSkew);

    row.Score = NaN;
    row.reprojMeanMag = NaN;
    row.reprojStdMag  = NaN;
    row.err1 = NaN;
    row.err2 = NaN;

    row.score_fx = NaN; row.score_fy = NaN; row.score_cx = NaN; row.score_cy = NaN;
    row.score_k_sum = NaN;
    row.score_alpha = NaN;
    row.score_err   = NaN;

    row.alpha = NaN; row.alpha_err = NaN;
    row.k1 = NaN; row.k2 = NaN; row.p1 = NaN; row.p2 = NaN; row.k3 = NaN;
    row.k1_err = NaN; row.k2_err = NaN; row.p1_err = NaN; row.p2_err = NaN; row.k3_err = NaN;

    row.alpha_near_zero = false;
    row.p1_near_zero = false;
    row.p2_near_zero = false;

    row.numImagesUsedByEstimator = NaN;
end

function row = makeRunRow(rowCount, n, k, subset, score, meanReprojMag, stdReprojMag, err1, err2, ...
    parts, alpha, alpha_err, kc, kc_err, numRadial, estTan, estSkew, imagesUsedByEstimator)

    row = struct();
    row.rowId = rowCount;
    row.N = n;
    row.k = k;
    row.subsetIdx = join(string(subset), ",");
    row.NumRadial = numRadial;
    row.EstTan = logical(estTan);
    row.EstSkew = logical(estSkew);

    row.Score = score;
    row.reprojMeanMag = meanReprojMag;
    row.reprojStdMag  = stdReprojMag;
    row.err1 = err1;
    row.err2 = err2;

    row.score_fx = parts.fx;
    row.score_fy = parts.fy;
    row.score_cx = parts.cx;
    row.score_cy = parts.cy;
    row.score_k_sum = parts.k_sum;
    row.score_alpha = parts.alpha;
    row.score_err   = parts.err;

    row.alpha = alpha;
    row.alpha_err = alpha_err;

    row.k1 = kc(1); row.k2 = kc(2); row.p1 = kc(3); row.p2 = kc(4); row.k3 = kc(5);
    row.k1_err = kc_err(1); row.k2_err = kc_err(2); row.p1_err = kc_err(3); row.p2_err = kc_err(4); row.k3_err = kc_err(5);

    row.alpha_near_zero = parts.flags.alpha_near_zero;
    row.p1_near_zero = parts.flags.p1_near_zero;
    row.p2_near_zero = parts.flags.p2_near_zero;

    row.numImagesUsedByEstimator = sum(imagesUsedByEstimator);
end

function [meanMag, stdMag, err1_rms, err2_rms] = reprojMetrics(reprojErr)
    % reprojErr: Mx2xN (pixels)
    e = reprojErr;
    mag = sqrt(e(:,1,:).^2 + e(:,2,:).^2);  % Mx1xN
    mag = mag(:);
    meanMag = mean(mag, 'omitnan');
    stdMag  = std(mag,  'omitnan');

    ex = e(:,1,:); ex = ex(:);
    ey = e(:,2,:); ey = ey(:);

    err1_rms = sqrt(mean(ex.^2, 'omitnan'));
    err2_rms = sqrt(mean(ey.^2, 'omitnan'));
end

function [score, parts, alpha, kc, kc_err, alpha_err] = computeScoreLikeAssignmentRobust(params, errors, err1, err2, cfg)
    % Robust assignment-style score:
    %   sum(|param_err|/max(|param|, floor)) + sqrt(err1^2+err2^2)
    % Floors prevent explosions when param ~ 0 (typical for skew or tangential distortion).

    if nargin < 5 || isempty(cfg)
        cfg = struct();
    end
    if ~isfield(cfg,'alpha_floor'), cfg.alpha_floor = 1.0; end
    if ~isfield(cfg,'kc_floor'),    cfg.kc_floor    = [1e-3 1e-3 5e-4 5e-4 1e-3]; end
    if ~isfield(cfg,'epsDen'),      cfg.epsDen      = 1e-12; end

    epsDen = cfg.epsDen;

    % Intrinsics
    fc = params.FocalLength(:).';          % [fx fy]
    cc = params.PrincipalPoint(:).';       % [cx cy]
    alpha = safeProp(params, "Skew", 0);

    % Distortion (MATLAB: RadialDistortion [k1 k2 k3?], Tangential [p1 p2])
    kRad = safeProp(params, "RadialDistortion", []);
    pTan = safeProp(params, "TangentialDistortion", []);

    k1 = getIdxOrZero(kRad,1);
    k2 = getIdxOrZero(kRad,2);
    k3 = getIdxOrZero(kRad,3);
    p1 = getIdxOrZero(pTan,1);
    p2 = getIdxOrZero(pTan,2);

    kc = [k1 k2 p1 p2 k3];

    % Errors from estimationErrors
    ie = errors.IntrinsicsErrors;

    fc_err = ie.FocalLengthError(:).';
    cc_err = ie.PrincipalPointError(:).';
    alpha_err = safeProp(ie, "SkewError", 0);

    kRadErr = safeProp(ie, "RadialDistortionError", zeros(size(kRad)));
    pTanErr = safeProp(ie, "TangentialDistortionError", zeros(size(pTan)));

    k1_err = getIdxOrZero(kRadErr,1);
    k2_err = getIdxOrZero(kRadErr,2);
    k3_err = getIdxOrZero(kRadErr,3);
    p1_err = getIdxOrZero(pTanErr,1);
    p2_err = getIdxOrZero(pTanErr,2);

    kc_err = [k1_err k2_err p1_err p2_err k3_err];

    % Components
    parts = struct();
    parts.fx = abs(fc_err(1)) / max(abs(fc(1)), epsDen);
    parts.fy = abs(fc_err(2)) / max(abs(fc(2)), epsDen);
    parts.cx = abs(cc_err(1)) / max(abs(cc(1)), epsDen);
    parts.cy = abs(cc_err(2)) / max(abs(cc(2)), epsDen);

    denom_kc = max(abs(kc), cfg.kc_floor);
    parts.k_terms = abs(kc_err) ./ denom_kc;
    parts.k_sum   = sum(parts.k_terms);

    denom_alpha = max(abs(alpha), cfg.alpha_floor);
    parts.alpha = abs(alpha_err) / denom_alpha;

    parts.err = sqrt(err1^2 + err2^2);

    score = parts.fx + parts.fy + parts.cx + parts.cy + parts.k_sum + parts.alpha + parts.err;

    parts.flags.alpha_near_zero = abs(alpha) < cfg.alpha_floor;
    parts.flags.p1_near_zero    = abs(p1) < cfg.kc_floor(3);
    parts.flags.p2_near_zero    = abs(p2) < cfg.kc_floor(4);
end

function out = safeProp(obj, propName, defaultVal)
    if isprop(obj, propName)
        out = obj.(propName);
    else
        out = defaultVal;
    end
end

function val = getIdxOrZero(vec, idx)
    if isempty(vec) || numel(vec) < idx
        val = 0;
    else
        val = vec(idx);
    end
end

function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end