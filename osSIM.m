% Optical Sectioned Structured Illumination Microscopy Image Processing

function osSIM()
    clc;
    %% Basic Parameters
    % Number of SIM phases
    numPhase = 3;
    % Number of SIM directions
    numDirection = 3;
    % Winener filter parameter
    paraWiener = 0.5;
    % Magnification
    paraMag = 50.0 * 250 / 200;
    % Numerical aperture
    paraNA = 0.8;
    % Pixel size
    paraPixelSz = 6.5 / paraMag;
    % Gain for high frequency
    paraAttAmp = 0;
    % High gain range for high frequency
    paraAttSigma = 2 * pi * 1.5;
    % Wavelength for emission
    paraWavelength = 0.52;

    %% Select SIM raw data (standard Tiff format <4G)
    [filename_raw, pathname_raw] = uigetfile( ...
        '*.tif', 'Select SIM Raw Tiff Stack');

    if isequal(filename_raw, 0)
        error('SIM raw tiff selection cancelled');
    else
        disp('Open SIM raw tiff file:');
        disp(fullfile(pathname_raw, filename_raw));
    end

    %% Reading size of raw data
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

    %% Stack and average raw data
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

    %% Average the stack
    imgRaw = imgRaw ./ numStackGroup;
    warning on;
    tiff_raw.close();

    %% Background substraction with optical sectioned SIM
    %  osImg = sqrt(2) * var(imgRaw);
    %  or    = sqrt(((A-B).^2+(A-C).^2+(B-C).^2)*2) / 3;
    osImg = zeros(szWidth, szHeight, numDirection, 'single');
    osBG = zeros(szWidth, szHeight, numDirection, 'single');
    for d = 1:numDirection
        for p = 1:numPhase
            for p2 = 1:numPhase
                if p2 ~= p
                    osImg(:, :, d) = osImg(:, :, d) + (imgRaw(:, :, p + (d - 1) * numPhase) - imgRaw(:, :, p2 + (d - 1) * numPhase)).^2;
                end
            end
            osBG(:, :, d) = osBG(:, :, d) + imgRaw(:, :, p + (d - 1) * numPhase) / numPhase;
        end
        osImg(:, :, d) = sqrt(osImg(:, :, d)) / numPhase;
        osBG(:, :, d) = osBG(:, :, d) - osImg(:, :, d);
        for p = 1:numPhase
            imgRaw(:, :, p + (d - 1) * numPhase) = imgRaw(:, :, p + (d - 1) * numPhase) - osBG(:, :, d);
        end
    end

    tiff_os = Tiff([pathname_raw, 'os_', filename_raw], 'w');
    tiff_osTag = struct('ImageWidth', szWidth, ...
        'ImageLength', szHeight, ...
        'Photometric', Tiff.Photometric.MinIsBlack, ...
        'BitsPerSample', 32, ...
        'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
        'SamplesPerPixel', 1, ...
        'PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    tiff_os.setTag(tiff_osTag);
    tiff_os.write(mean(osImg, 3));
    tiff_os.writeDirectory();
    tiff_os.close();
end
