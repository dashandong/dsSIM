% Image Reconstruction for Structured Illumination Microscopy(SIM)
% Copyright (C) 2023 by Dr. Dashan Dong (dongdashan@icloud.com)
%
% This program is free software:you can redistribute it and / or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see < https: // www.gnu.org / licenses /> .

% last edited @ 2023/09/09 by Dashan Dong

function dsSIM
    clc;
    clear;

    fprintf('\tdsSIM  Copyright (C) 2023  by Dashan Dong\n');
    fprintf('\tThis program comes with ABSOLUTELY NO WARRANTY.\n');
    fprintf('\tThis is free software, and you are welcome to redistribute it\n');
    fprintf('\tunder certain conditions.\n\n');

    %% 基本参数
    % 每个方向条纹移相次数
    numPhase = 3;
    % 照明条纹方向变化次数
    numDirection = 3;
    % Wiener滤波参数基数
    paraWiener = 0.1;
    % 成像系统放大率
    paraMag = 100.0;
    % 成像物镜数值孔径
    paraNA = 1.49;
    % 图像实际像素尺寸
    paraPixelSz = 6.5 / paraMag;
    % 高频增益幅度
    paraAttAmp = 0;
    % 高频增益范围
    paraAttSigma = 2 * pi * 1.5;

    %% 选择SIM原始数据TIFF文件(标准Tiff格式<4G)
    [filename_raw, pathname_raw] = uigetfile( ...
        '*.tif', 'Select SIM Raw Tiff Stack');

    if isequal(filename_raw, 0)
        error('SIM raw tiff selection cancelled');
    else
        disp('Open SIM raw tiff file:');
        disp(fullfile(pathname_raw, filename_raw));
    end

    %% 开始读取文件获取图像尺寸信息
    warning off;
    tiff_raw = Tiff([pathname_raw, filename_raw], 'r');
    szWidth = tiff_raw.getTag('ImageWidth');
    szHeight = tiff_raw.getTag('ImageLength');
    szStack = 1;

    while (~tiff_raw.lastDirectory())
        tiff_raw.nextDirectory();
        szStack = szStack + 1;
    end

    warning on;
    disp(['Tiff Size: ', num2str(szWidth), ' x ', num2str(szHeight), ' x ', num2str(szStack)]);

    if (szStack < numPhase * numDirection)
        error('Raw Tiff stack is too short to be SIM rawdata!');
    end

    %% 用户输入求参时图像堆栈范围
    disp('Start to solve the parameters...');
    numStackGroup = input([ ...
                               'How many time points to stack together?\n', ...
                               '(Range: 1-', ...
                               num2str(floor(szStack / (numPhase * numDirection))), ...
                           ', default 10):']);

    if isempty(numStackGroup)
        numStackGroup = 10;
    end

    numStackGroup = round(numStackGroup);

    if (numStackGroup < 1) || (numStackGroup > floor(szStack / (numPhase * numDirection)))
        error('Wrong number input!');
    end

    disp(['Stack ', num2str(numStackGroup), ' time points to stack.']);

    numStackStart = 1;

    if (szStack > numStackGroup * numPhase * numDirection)
        numStackStart = input([ ...
                                   'Start stacking from frame #?\n', ...
                                   '(Range: 1-', ...
                                   num2str(szStack - numStackGroup * numPhase * numDirection + 1), ...
                               ', default 1):']);

        if isempty(numStackStart)
            numStackStart = 1;
        end

        numStackStart = round(numStackStart);

        if (numStackStart < 1) || (numStackStart > szStack - numStackGroup * numPhase * numDirection + 1)
            error('Wrong number input!');
        end

    end

    imgRaw = zeros(szWidth, szHeight, numPhase * numDirection, 'single');
    warning off;

    for i = 1:numStackGroup

        for k = 1:(numPhase * numDirection)
            tiff_raw.setDirectory(numStackStart + (i - 1) * (numPhase * numDirection) + (k - 1));
            imgRaw(:, :, k) = imgRaw(:, :, k) + single(tiff_raw.read());
        end

    end

    imgRaw = imgRaw ./ numStackGroup;
    warning on;
    tiff_raw.close();

    %% 读取OTF
    [filename_otf, pathname_otf] = uigetfile( ...
        '*.tif', 'Select OTF');

    if isequal(filename_otf, 0)
        error('OTF Selection Cancelled!');
    end

    OTF = im2double(imread([pathname_otf, filename_otf]));
    OTF = imresize(OTF, [szWidth, szHeight], 'bilinear');
    OTF = OTF / max(OTF(:));
    OTFx2 = zeros(size(OTF) * 2);
    OTFx2((floor((szWidth - 1) / 2) + 2):(floor((szWidth - 1) / 2) + 1 + szWidth), ...
        (floor((szHeight - 1) / 2) + 2):(floor((szHeight - 1) / 2) + 1 + szHeight)) = OTF;
    OTF = fftshift(OTF);
    OTFx2 = fftshift(OTFx2);

    %% 读取Background
    [filename_bg, pathname_bg] = uigetfile( ...
        '*.tif', 'Select Background');

    if isequal(filename_bg, 0)
        warning('Background selection cancelled, NO background correction.');
        BG = zeros(szWidth, szHeight);
    end

    BG = im2double(imread([pathname_bg, filename_bg]));

    if ~isequal(size(BG), [szWidth, szHeight])
        warning('Worng background size, NO background correction.');
        BG = zeros(szWidth, szHeight);
    end

    % 扣除背景
    for i = 1:(numPhase * numDirection)
        imgRaw(:, :, i) = imgRaw(:, :, i) - BG;
    end

    %% 中间变量计算
    paraKm = 2 * pi * paraNA / 0.52;
    paraKWSz = 2 * pi / paraPixelSz / szWidth;
    paraKHSz = 2 * pi / paraPixelSz / szHeight;
    paraPhaseVector = ((1:numPhase) - 1) * 2 * pi / numPhase;
    paraPhaseVectorFull = ones(numPhase, numDirection);
    paraPhaseMatrix = ones(numPhase, 3);
    paraPhaseMatrixFull = ones(numPhase, 3, numDirection);
    paraIMatrix = ones(3, numDirection);

    % 构建相位矩阵
    for p = 1:numPhase
        paraPhaseVectorFull(p, :) = paraPhaseVector(p);
        paraPhaseMatrix(p, :) = paraPhaseMatrix(p, :) .* ...
            exp(1i * [-paraPhaseVector(p), 0, paraPhaseVector(p)]);
    end

    % 坐标架生成
    axisW = ((1:szWidth) - 1 - floor(szWidth / 2)) / szWidth;
    axisH = ((1:szHeight) - 1 - floor(szHeight / 2)) / szHeight;
    [meshW, meshH] = meshgrid(axisW, axisH);

    axisKW = ((1:szWidth) - 1 - floor(szWidth / 2)) * paraKWSz;
    axisKH = ((1:szHeight) - 1 - floor(szHeight / 2)) * paraKHSz;

    axisWx2 = ((1:(szWidth * 2)) - 1 - szWidth) / (szWidth * 2);
    axisHx2 = ((1:(szHeight * 2)) - 1 - szHeight) / (szHeight * 2);
    [meshWx2, meshHx2] = meshgrid(axisWx2, axisHx2);

    axisKWx2 = ((1:(szWidth * 2)) - 1 - szWidth) * paraKWSz;
    axisKHx2 = ((1:(szHeight * 2)) - 1 - szHeight) * paraKHSz;

    [meshKW, meshKH] = meshgrid (axisKW, axisKH);
    [~, meshKR] = cart2pol(meshKW, meshKH);
    meshKR_shift = ifftshift(meshKR);

    [meshKWx2, meshKHx2] = meshgrid(axisKWx2, axisKHx2);
    [~, meshKRx2] = cart2pol(meshKWx2, meshKHx2);
    meshKRx2_shift = ifftshift(meshKRx2);

    % 高频增益函数
    attFun = 1 - paraAttAmp * exp(-meshKRx2_shift .^ 2 / (2 * paraAttSigma ^ 2));

    % 波矢求解滤波器
    OTFMask = (meshKR_shift < 2 * paraKm);
    OTFRingMask = (meshKR_shift < 2 * paraKm) & (meshKR_shift > 1.3 * paraKm);
    OTFMaskx2 = (meshKRx2_shift < 2 * paraKm);

    paraIllVector = zeros(numDirection, 2);

    %% 傅里叶变换原始图像
    spectrumRaw = complex(zeros(szHeight, szWidth, ...
        numPhase * numDirection));

    for i = 1:(numPhase * numDirection)
        spectrumRaw(:, :, i) = fft2(ifftshift(imgRaw(:, :, i)));
    end

    %% 粗分离各方向三个频谱成分
    spectrumSep = complex(zeros(3, szHeight, szWidth, numDirection));

    for d = 1:numDirection
        spectrumSepEachP = complex(zeros(numPhase, szHeight * szWidth));

        for p = 1:numPhase
            spectrumTemp = spectrumRaw(:, :, p + (d - 1) * numPhase);
            spectrumSepEachP(p, :) = spectrumTemp(:);
        end

        spectrumSepEach = paraPhaseMatrix \ spectrumSepEachP;
        spectrumSep((((d - 1) * 3 * szHeight * szWidth) + 1) ...
            :(d * 3 * szHeight * szWidth)) = spectrumSepEach(:);
    end

    % 绘制粗分离频谱，波矢粗定位
    f = figure( ...
        'MenuBar', 'none', ...
        'ToolBar', 'none', ...
        'Name', ['SIM parameters for ', filename_raw], ...
        'Units', 'normalized', ...
        'Position', [0.1 0.1 0.8 0.8]);
    tiledlayout(f, 3, 4);
    spectrumCrossCorr = zeros(szHeight, szWidth, numDirection);
    disp('Locating vectors...');

    for d = 1:numDirection
        spectrum1 = squeeze(spectrumSep(1, :, :, d));
        spectrum0 = squeeze(spectrumSep(2, :, :, d));
        spectrum2 = squeeze(spectrumSep(3, :, :, d));
        spectrumCrossCorr1 = abs(fft2(conj(ifft2(spectrum1)) .* ifft2(spectrum0)));
        spectrumCrossCorr2 = abs(fft2(conj(ifft2(spectrum2)) .* ifft2(spectrum0)));
        spectrumCrossCorr(:, :, d) = spectrumCrossCorr1 + spectrumCrossCorr2;

        nexttile((d - 1) * 4 + 1);
        pbg = imagesc( ...
            [-floor(szWidth / 2), szWidth - 1 - floor(szWidth / 2)], ...
            [-floor(szHeight / 2), szHeight - 1 - floor(szHeight / 2)], ...
            ifftshift(spectrumCrossCorr(:, :, d) .* OTFRingMask));
        colormap(contrast(pbg.CData));
        title(['Direction #', num2str(d)]);
        axis image off;
        hold on;

        spectrumLoc = ifftshift(spectrumCrossCorr1 .* OTFRingMask);
        [~, indMax] = max(spectrumLoc(:));
        [indH, indW] = ind2sub(size(spectrumLoc), indMax);
        indW = indW - floor(szWidth / 2) - 1;
        indH = indH - floor(szHeight / 2) - 1;
        plot(indW, indH, 'go');
        paraIllVector(d, :) = [indW, indH];

        spectrumLoc = ifftshift(spectrumCrossCorr2 .* OTFRingMask);
        [~, indMax] = max(spectrumLoc(:));
        [indH, indW] = ind2sub(size(spectrumLoc), indMax);
        indW = indW - floor(szWidth / 2) - 1;
        indH = indH - floor(szHeight / 2) - 1;
        plot(indW, indH, 'ro');

        paraIllVector(d, :) = (paraIllVector(d, :) - [indW, indH]) / 2;
    end

    % 绘制精细波矢位置图
    for d = 1:numDirection
        nexttile((d - 1) * 4 + 2);
        pbg = imagesc( ...
            [-floor(szWidth / 2), szWidth - 1 - floor(szWidth / 2)], ...
            [-floor(szHeight / 2), szHeight - 1 - floor(szHeight / 2)], ...
            ifftshift(spectrumCrossCorr(:, :, d) .* OTFRingMask));
        colormap(gray);
        axis image;
        xlim([-2 2] + paraIllVector(d, 1))
        ylim([-2 2] + paraIllVector(d, 2))
        title(['Direction #', num2str(d), ' Zoom In']);
        hold on;
        plot(paraIllVector(d, 1), paraIllVector(d, 2), 'go');
    end

    drawnow;

    %% 迭代求解波矢
    for d = 1:numDirection
        disp(['Vector iteration for direction #', num2str(d), ' ...']);
        iterGradStep = 1e-7;
        iterBeta = 0.8;
        vector = paraIllVector(d, :);
        corValw1 = getCorrVal(spectrumSep(:, :, :, d), ...
            vector - [0.5 * iterGradStep, 0]);
        corValh1 = getCorrVal(spectrumSep(:, :, :, d), ...
            vector - [0, 0.5 * iterGradStep]);
        corValw2 = getCorrVal(spectrumSep(:, :, :, d), ...
            vector + [0.5 * iterGradStep, 0]);
        corValh2 = getCorrVal(spectrumSep(:, :, :, d), ...
            vector + [0, 0.5 * iterGradStep]);
        corGrad = [corValw2 - corValw1, corValh2 - corValh1] ...
            / iterGradStep;
        iterStep = 0.1 / max(corGrad);
        corStep = iterStep * corGrad;

        nexttile((d - 1) * 4 + 2);
        p = plot(vector(1), vector(2), 'rx');
        drawnow;

        while norm(corStep) > iterGradStep
            vector = vector + corStep;
            corValw1 = getCorrVal(spectrumSep(:, :, :, d), ...
                vector - [0.5 * iterGradStep, 0]);
            corValh1 = getCorrVal(spectrumSep(:, :, :, d), ...
                vector - [0, 0.5 * iterGradStep]);
            corValw2 = getCorrVal(spectrumSep(:, :, :, d), ...
                vector + [0.5 * iterGradStep, 0]);
            corValh2 = getCorrVal(spectrumSep(:, :, :, d), ...
                vector + [0, 0.5 * iterGradStep]);
            corGrad = [corValw2 - corValw1, corValh2 - corValh1] ...
                / iterGradStep;
            corStep = iterBeta * iterStep * corGrad + ...
                (1 - iterBeta) * corStep;

            p.XData = [p.XData, vector(1)];
            p.YData = [p.YData, vector(2)];
            xlabel(['k_x=', num2str(vector(1))]);
            ylabel(['k_y=', num2str(vector(2))]);
            drawnow;
        end

        vector = vector + corStep;
        paraIllVector(d, :) = vector;
    end

    %% 求解初始相位和调制深度
    for d = 1:numDirection
        disp(['Solving initial phase and modulation depth for direction #', num2str(d), ' ...']);
        [gamma, phi] = getParameters(spectrumSep(:, :, :, d), ...
            paraIllVector(d, :), d);
        paraPhaseVectorFull(:, d) = wrapToPi( ...
            paraPhaseVectorFull(:, d) - phi);
        paraIMatrix(:, d) = [2 / gamma, 1, 2 / gamma];

        for p = 1:numPhase
            paraPhaseMatrixFull(p, :, d) = ...
                paraPhaseMatrixFull(p, :, d) .* exp(1i * ...
                [-paraPhaseVectorFull(p, d), ...
                 0, ...
                 paraPhaseVectorFull(p, d)]);
        end

    end

    %% 存储求解报告
    exportgraphics(f, [pathname_raw, filesep, filename_raw, '.pdf'], 'ContentType', 'vector');

    %% 按照精细参数重新分离频谱
    spectrumSep2x = complex( ...
        zeros(3, szHeight * 2, szWidth * 2, numDirection));

    for d = 1:numDirection
        spectrumSepEachP = complex(zeros(numPhase, szHeight * szWidth));

        for p = 1:numPhase
            spectrumTemp = spectrumRaw(:, :, p + (d - 1) * numPhase);
            spectrumSepEachP(p, :) = spectrumTemp(:);
        end

        spectrumSepEach = squeeze(paraPhaseMatrixFull(:, :, d)) ...
            \ spectrumSepEachP;
        spectrumSep((((d - 1) * 3 * szHeight * szWidth) + 1) ...
            :(d * 3 * szHeight * szWidth)) = spectrumSepEach(:);

        for i = 1:3
            spectrumSep(i, :, :, d) = paraIMatrix(i, d) .* ...
                spectrumSep(i, :, :, d);
            spectrumTemp = fftshift(squeeze(spectrumSep(i, :, :, d)));
            spectrumSep2x(i, ...
                (floor((szHeight - 1) / 2) + 2):(floor((szHeight - 1) / 2) + 1 + szHeight), ...
                (floor((szHeight - 1) / 2) + 2):(floor((szHeight - 1) / 2) + 1 + szHeight), ...
                d) = spectrumTemp * 4;
        end

    end

    %% Wiener逆滤波
    spectrumSep2x = ifftshift(spectrumSep2x, 2);
    spectrumSep2x = ifftshift(spectrumSep2x, 3);

    spectrumSum = complex(zeros(szHeight * 2, szWidth * 2));
    spectrumWiener = complex(zeros(szHeight * 2, szWidth * 2));

    for d = 1:numDirection
        spectrumSum = spectrumSum + ...
            specCombine(spectrumSep2x(:, :, :, d), paraIllVector(d, :));
        spectrumWiener = spectrumWiener + ...
            wienerCombine(paraIllVector(d, :));
    end

    %% 预览结果，调整Wiener参数
    close(f);
    uiResult = uifigure('Name', 'SIM Result', ...
        'Units', 'normalized', ...
        'Position', [0.1, 0.1, 0.8, 0.8]);
    uiLayout = uigridlayout(uiResult, [1, 3]);
    uiLayout.RowHeight = {'1x', 50, 20};
    uiLayout.ColumnWidth = {'1x'};
    uiAxes = uiaxes(uiLayout);
    uiAxes.Layout.Row = 1;
    uiAxes.Layout.Column = 1;
    spectrumResult = spectrumSum ./ (spectrumWiener + paraWiener ^ 2);
    simResult = abs(real(fftshift(ifft2(spectrumResult))));
    rimg = imagesc(uiAxes, simResult);
    axis(uiAxes, 'image');
    axis(uiAxes, 'off');
    colormap(uiAxes, gray);
    clim(uiAxes, [0, 0.5 * max(simResult(:))]);
    uiWienerSlider = uislider (uiLayout, 'Limits', [-5, 5], 'ValueChangedFcn', @drawResult);
    uiWienerSlider.Layout.Row = 2;
    uiWienerSlider.Layout.Column = 1;
    uiProcessButton = uibutton(uiLayout, 'push', ...
        'Text', 'Start Batch Process', ...
        'ButtonPushedFcn', @batchProcess);
    uiProcessButton.Layout.Row = 3;
    uiProcessButton.Layout.Column = 1;

    %% 批量处理TIff
    function batchProcess(src, event)
        disp('Start batch processing ...');
        paraWiener = paraWiener * 10 ^ (uiWienerSlider.Value);
        uiProcessButton.Enable = 'off';
        uiWienerSlider.Enable = 'off';
        numStackStart = mod(numStackStart - 1, 9) + 1;
        numStackGroup = (szStack - numStackStart + 1) / (numPhase * numDirection);
        uiProgress = uiprogressdlg(uiResult, 'Title', 'Batch processing...', ...
            'Message', ['0 / ', num2str(numStackGroup)]);

        % 打开文件
        tiff_outDescription = {'ImageJ', ...
                                   'unit=\u00B5m', ...
                                   'loop=false', ...
                                   'min=0', ...
                                   ['max=', num2str(round(0.5 * max(simResult(:))))]};
        warning off;
        tiff_raw = Tiff([pathname_raw, filename_raw], 'r');
        tiff_out = Tiff([pathname_raw, 'SIM_', filename_raw], 'w');
        tiff_outTag = struct('ImageWidth', szWidth * 2, ...
            'ImageLength', szHeight * 2, ...
            'Photometric', Tiff.Photometric.MinIsBlack, ...
            'BitsPerSample', 32, ...
            'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
            'SamplesPerPixel', 1, ...
            'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        tiff_out.setTag('ImageDescription', sprintf('%s\n', tiff_outDescription{:}));
        tiff_out.setTag('XResolution', 2 / paraPixelSz);
        tiff_out.setTag('YResolution', 2 / paraPixelSz);
        tiff_out.setTag('ResolutionUnit', 1);

        tiff_WF = Tiff([pathname_raw, 'WF_', filename_raw], 'w');
        tiff_WFTag = struct('ImageWidth', szWidth, ...
            'ImageLength', szHeight, ...
            'Photometric', Tiff.Photometric.MinIsBlack, ...
            'BitsPerSample', 32, ...
            'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
            'SamplesPerPixel', 1, ...
            'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        tiff_WF.setTag('ImageDescription', sprintf('%s\n', tiff_outDescription{:}));
        tiff_WF.setTag('XResolution', 1 / paraPixelSz);
        tiff_WF.setTag('YResolution', 1 / paraPixelSz);
        tiff_WF.setTag('ResolutionUnit', 1);
        warning on;

        for g = 1:numStackGroup

            for f = 1:(numPhase * numDirection)
                warning off;
                tiff_raw.setDirectory(numStackStart + (g - 1) * (numPhase * numDirection) + (f - 1));
                imgRaw(:, :, f) = single(tiff_raw.read());
                warning on;
                spectrumRaw(:, :, f) = fft2(ifftshift(imgRaw(:, :, f)));
            end

            for d = 1:numDirection
                spectrumSepEachP = complex(zeros(numPhase, szHeight * szWidth));

                for p = 1:numPhase
                    spectrumTemp = spectrumRaw(:, :, p + (d - 1) * numPhase);
                    spectrumSepEachP(p, :) = spectrumTemp(:);
                end

                spectrumSepEach = squeeze(paraPhaseMatrixFull(:, :, d)) ...
                    \ spectrumSepEachP;
                spectrumSep((((d - 1) * 3 * szHeight * szWidth) + 1) ...
                    :(d * 3 * szHeight * szWidth)) = spectrumSepEach(:);

                for r = 1:3
                    spectrumSep(r, :, :, d) = paraIMatrix(r, d) .* ...
                        spectrumSep(r, :, :, d);
                    spectrumTemp = fftshift(squeeze(spectrumSep(r, :, :, d)));
                    spectrumSep2x(r, ...
                        (floor((szHeight - 1) / 2) + 2):(floor((szHeight - 1) / 2) + 1 + szHeight), ...
                        (floor((szHeight - 1) / 2) + 2):(floor((szHeight - 1) / 2) + 1 + szHeight), ...
                        d) = spectrumTemp * 4;
                end

            end

            spectrumSep2x = ifftshift(spectrumSep2x, 2);
            spectrumSep2x = ifftshift(spectrumSep2x, 3);

            spectrumSum = complex(zeros(szHeight * 2, szWidth * 2));
            spectrumWiener = complex(zeros(szHeight * 2, szWidth * 2));

            for d = 1:numDirection
                spectrumSum = spectrumSum + ...
                    specCombine(spectrumSep2x(:, :, :, d), paraIllVector(d, :));
                spectrumWiener = spectrumWiener + ...
                    wienerCombine(paraIllVector(d, :));
            end

            spectrumResult = spectrumSum ./ (spectrumWiener + paraWiener ^ 2);
            simResult = abs(real(fftshift(ifft2(spectrumResult))));
            tiff_out.setTag(tiff_outTag);
            tiff_out.write(single(simResult));
            tiff_out.writeDirectory();
            tiff_WF.setTag(tiff_WFTag);
            tiff_WF.write(sum(single(imgRaw), 3));
            tiff_WF.writeDirectory();
            uiProgress.Message = [num2str(g), ' / ', num2str(numStackGroup)];
        end

        tiff_raw.close();
        tiff_out.close();
        tiff_WF.close();
    end

    %% 重绘预览结果回调函数
    function drawResult(src, event)
        paraWienerTemp = paraWiener * 10 ^ (src.Value);
        spectrumResult = spectrumSum ./ (spectrumWiener + paraWienerTemp ^ 2);
        simResult = real(fftshift(ifft2(spectrumResult)));
        rimg.CData = simResult;
        clim(uiAxes, [0, 0.5 * max(simResult(:))]);
    end

    %% 获取特定频移波矢下的相关值
    function corrValve = getCorrVal(spectrumSep3, shiftVector)
        shiftMask = exp(1i * 2 * pi * ...
            (shiftVector(1) * meshW + shiftVector(2) * meshH));
        mask0 = OTFMask;
        mask1 = real(fft2(ifftshift(fftshift(ifft2(OTFMask)) ...
            .* shiftMask)));
        mask2 = real(fft2(ifftshift(fftshift(ifft2(OTFMask)) ...
            .* conj(shiftMask))));
        mask1(mask1 < 0.5) = 0;
        mask2(mask2 < 0.5) = 0;
        mask1(mask1 ~= 0) = 1;
        mask2(mask2 ~= 0) = 1;
        OTF1 = real(fft2(ifftshift(fftshift(ifft2(OTF)) ...
            .* shiftMask)));
        OTF2 = real(fft2(ifftshift(fftshift(ifft2(OTF)) ...
            .* conj(shiftMask))));
        spec0 = squeeze(spectrumSep3(2, :, :));
        spec1 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(1, :, :)))) .* shiftMask));
        spec2 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(3, :, :)))) .* conj(shiftMask)));
        maskCom1 = mask0 .* mask1;
        maskCom2 = mask0 .* mask2;
        specCom1 = maskCom1 .* OTF .* OTF1 .* ...
            spec0 .* conj(spec1);
        specCom2 = maskCom2 .* OTF .* OTF2 .* ...
            spec0 .* conj(spec2);
        corrValve = (abs(sum(specCom1(:))) / sum(maskCom1(:)) + ...
            abs(sum(specCom2(:))) / sum(maskCom2(:))) / 2;
    end

    %% 组合拼接特定波矢下的三个频谱
    function spectrumNew = specCombine(spectrumSep3, shiftVector)
        shiftMask = exp(1i * 2 * pi * ...
            (shiftVector(1) * meshWx2 + shiftVector(2) * meshHx2));
        mask0 = OTFMaskx2;
        mask1 = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* shiftMask)));
        mask2 = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* conj(shiftMask))));
        mask1(mask1 < 0.5) = 0;
        mask2(mask2 < 0.5) = 0;
        mask1(mask1 ~= 0) = 1;
        mask2(mask2 ~= 0) = 1;
        OTF1 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
            .* shiftMask)));
        OTF2 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
            .* conj(shiftMask))));
        spec0 = squeeze(spectrumSep3(2, :, :));
        spec1 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(1, :, :)))) .* shiftMask));
        spec2 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(3, :, :)))) .* conj(shiftMask)));
        attFun1 = fft2(ifftshift(fftshift(ifft2( ...
            attFun)) .* shiftMask));
        attFun2 = fft2(ifftshift(fftshift(ifft2( ...
            attFun)) .* conj(shiftMask)));
        sum0 = mask0 .* attFun .* conj(OTFx2) .* spec0;
        sum1 = mask1 .* attFun1 .* conj(OTF1) .* spec1;
        sum2 = mask2 .* attFun2 .* conj(OTF2) .* spec2;
        spectrumNew = sum0 + sum1 +sum2;
    end

    %% 求解组合拼接特定波矢下的OTF
    function wienerNew = wienerCombine(shiftVector)
        shiftMask = exp(1i * 2 * pi * ...
            (shiftVector(1) * meshWx2 + shiftVector(2) * meshHx2));
        mask0 = OTFMaskx2;
        mask1 = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* shiftMask)));
        mask2 = real(fft2(ifftshift(fftshift(ifft2(OTFMaskx2)) ...
            .* conj(shiftMask))));
        mask1(mask1 < 0.5) = 0;
        mask2(mask2 < 0.5) = 0;
        mask1(mask1 ~= 0) = 1;
        mask2(mask2 ~= 0) = 1;
        OTF1 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
            .* shiftMask)));
        OTF2 = real(fft2(ifftshift(fftshift(ifft2(OTFx2)) ...
            .* conj(shiftMask))));
        attFun1 = fft2(ifftshift(fftshift(ifft2( ...
            attFun)) .* shiftMask));
        attFun2 = fft2(ifftshift(fftshift(ifft2( ...
            attFun)) .* conj(shiftMask)));
        sum0 = mask0 .* attFun .* conj(OTFx2) .* OTFx2;
        sum1 = mask1 .* attFun1 .* conj(OTF1) .* OTF1;
        sum2 = mask2 .* attFun2 .* conj(OTF2) .* OTF2;
        wienerNew = sum0 + sum1 +sum2;
    end

    %% 求解给定波矢下的调制深度和初始相位
    function [gamma, phi] = getParameters(spectrumSep3, shiftVector, ind)
        shiftMask = exp(1i * 2 * pi * ...
            (shiftVector(1) * meshW + shiftVector(2) * meshH));
        mask0 = OTFMask;
        mask1 = real(fft2(ifftshift(fftshift(ifft2(OTFMask)) ...
            .* shiftMask)));
        mask2 = real(fft2(ifftshift(fftshift(ifft2(OTFMask)) ...
            .* conj(shiftMask))));
        mask1(mask1 < 0.5) = 0;
        mask2(mask2 < 0.5) = 0;
        mask1(mask1 ~= 0) = 1;
        mask2(mask2 ~= 0) = 1;
        OTF1 = real(fft2(ifftshift(fftshift(ifft2(OTF)) ...
            .* shiftMask)));
        OTF2 = real(fft2(ifftshift(fftshift(ifft2(OTF)) ...
            .* conj(shiftMask))));
        spec0 = squeeze(spectrumSep3(2, :, :));
        spec1 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(1, :, :)))) .* shiftMask));
        spec2 = fft2(ifftshift(fftshift(ifft2( ...
            squeeze(spectrumSep3(3, :, :)))) .* conj(shiftMask)));
        maskcom1 = mask0 .* mask1;
        maskcom2 = mask0 .* mask2;
        test1 = (maskcom1 .* OTF .* spec1) ./ ...
            (maskcom1 .* OTF1 .* spec0);
        test2 = (maskcom2 .* OTF .* spec2) ./ ...
            (maskcom2 .* OTF2 .* spec0);
        testData1 = test1(maskcom1 ~= 0);
        testData2 = test2(maskcom2 ~= 0);
        theta = ((-180:1:179) + 0.5) / 180 * pi;
        rho1 = histcounts(angle(testData1), linspace(-pi, pi, 361));
        rho2 = histcounts(angle(testData2), linspace(-pi, pi, 361));
        [alpha1, fitdata1] = getElipseAlpha(theta, rho1);
        [alpha2, fitdata2] = getElipseAlpha(theta, rho2);
        phi = sign(alpha1) * (abs(alpha1) + abs(alpha2)) / 2;

        nexttile((ind - 1) * 4 + 3);
        polarplot(theta, rho1, 'b.', theta, rho2, 'r.');
        hold on
        polarplot(fitdata1(1, :), fitdata1(2, :), 'b-');
        polarplot(fitdata2(1, :), fitdata2(2, :), 'r-');
        title(['\phi=', num2str(rad2deg(phi)), ' deg']);
        drawnow;

        contrast = 0:0.01:2;
        contrastX = 0.005:0.01:1.995;
        contrastData1 = abs(testData1);
        contrastData2 = abs(testData2);
        contrastData1 = sqrt(contrastData1 .* contrastData1(end:-1:1));
        contrastData2 = sqrt(contrastData2 .* contrastData2(end:-1:1));
        contrastCount1 = histcounts(contrastData1, contrast);
        contrastCount1 = contrastCount1 ./ numel(contrastData1);
        contrastCount1 = contrastCount1 / 0.01;
        contrastCount2 = histcounts(contrastData2, contrast);
        contrastCount2 = contrastCount2 ./ numel(contrastData2);
        contrastCount2 = contrastCount2 / 0.01;
        contrastDist1 = fitdist(contrastData1, 'gev');
        contrastDist2 = fitdist(contrastData2, 'gev');
        contrastFit1 = pdf(contrastDist1, contrastX);
        contrastFit2 = pdf(contrastDist2, contrastX);
        gamma = (contrastDist1.mu + contrastDist2.mu);

        nexttile((ind - 1) * 4 + 4);
        plot(contrastX, contrastCount1, 'b.');
        hold on;
        plot(contrastX, contrastCount2, 'r.');
        plot(contrastX, contrastFit1, 'b-');
        plot(contrastX, contrastFit2, 'r-');
        title(['\gamma =', num2str(gamma)]);
        drawnow;
    end

    %% 极坐标下最小二乘拟合椭圆
    function [alpha_rad, ellipse] = getElipseAlpha(Theta, Rho)
        [x, y] = pol2cart(Theta(:), Rho(:));
        mean_x = mean(x(:));
        mean_y = mean(y(:));
        x = x - mean_x;
        y = y - mean_y;
        X = [x .^ 2, x .* y, y .^ 2, x, y];
        cm = sum(X) / (X' * X);
        [A, B, C, D, E] = deal(cm(1), cm(2), cm(3), cm(4), cm(5));
        orientation_rad = 1/2 * atan(B / (C - A));
        cos_phi = cos(orientation_rad);
        sin_phi = sin(orientation_rad);
        [A, B, C, D, E] = deal( ...
            A * cos_phi ^ 2 - B * cos_phi * sin_phi + C * sin_phi ^ 2, ...
            0, ...
            A * sin_phi ^ 2 + B * cos_phi * sin_phi + C * cos_phi ^ 2, ...
            D * cos_phi - E * sin_phi, ...
            D * sin_phi + E * cos_phi);
        [mean_x, mean_y] = deal( ...
            cos_phi * mean_x - sin_phi * mean_y, ...
            sin_phi * mean_x + cos_phi * mean_y);

        X0 = mean_x - D / 2 / A;
        Y0 = mean_y - E / 2 / C;
        F = 1 + (D ^ 2) / (4 * A) + (E ^ 2) / (4 * C);
        [A, B] = deal(sqrt(F / A), sqrt(F / C));

        R = [cos_phi sin_phi; -sin_phi cos_phi];

        ellipse_x_r = X0 + A * cos(Theta);
        ellipse_y_r = Y0 + B * sin(Theta);
        rotated_ellipse = R * [ellipse_x_r; ellipse_y_r];
        [ellipse(1, :), ellipse(2, :)] = cart2pol( ...
            rotated_ellipse(1, :), rotated_ellipse(2, :));
        alpha_rad = atan2(-sin_phi * X0 + cos_phi * Y0, ...
            cos_phi * X0 + sin_phi * Y0);
    end

end
