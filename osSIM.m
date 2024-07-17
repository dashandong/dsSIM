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
    paraMag = 100.0;
    % Numerical aperture
    paraNA = 1.49;
    % Pixel size
    paraPixelSz = 6.5 / paraMag;
    % Gain for high frequency
    paraAttAmp = 1;
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

    osImg1 = sqrt((imgRaw(:,:,1)-imgRaw(:,:,2)).^2+(imgRaw(:,:,1)-imgRaw(:,:,3)).^2+(imgRaw(:,:,2)-imgRaw(:,:,3)).^2);
    osImg2 = sqrt((imgRaw(:,:,4)-imgRaw(:,:,5)).^2+(imgRaw(:,:,4)-imgRaw(:,:,6)).^2+(imgRaw(:,:,5)-imgRaw(:,:,6)).^2);
    osImg3 = sqrt((imgRaw(:,:,7)-imgRaw(:,:,8)).^2+(imgRaw(:,:,7)-imgRaw(:,:,9)).^2+(imgRaw(:,:,8)-imgRaw(:,:,9)).^2);

    BG1 = mean(imgRaw(:,:,1:3), 3) - osImg1;
    BG2 = mean(imgRaw(:,:,4:6), 3) - osImg2;
    BG3 = mean(imgRaw(:,:,7:9), 3) - osImg3;

    imagesc(imgRaw(:,:,1) - BG1);
end
