function [spikesParams, spikesData] = detect_spikes_multiunit (vVec, siMs, varargin)
%% Detects spikes from a multiunit recording
% Usage: [spikesParams, spikesData] = detect_spikes_multiunit (vVec, siMs, varargin)
% Explanation:
%       Detects spikes from a multiunit recording via the following algorithm:
%           1. (optional) bandpass filters signal at 'FiltFreq' Hz
%           2. computes a slope threshold:
%               i. a baseline slope noise is computed from 
%                   the slope vector within 'BaseWindow'
%               ii. a signal-to-noise ratio is either provided ('Signal2Noise')
%                   or computed relative to the ratio of maximum slope 
%                       over baseline slope
%               iii. the slope threshold is signal-to-noise * baseline slope
%           3. finds all the local maxima of the slope vector 
%               (inflection points of the original signal) and 
%               calls the index a spike if it matches 3 criteria:
%                   i. occurs at least minDelayMs after stimulation start
%                   ii. occurs no more than maxDelayMs after stimulation start
%                   iii. slope value is >= the slope threshold 
% Example(s):
%       [spikesParams, spikesData] = detect_spikes_multiunit([0, 0, 0, 1, 0], 1, 'StimVec', [0; 0; 1; 0; 0]);
%
% Outputs:
%       spikesParams- Used and detected parameters, with fields:
%                       siMs
%                       idxStimStart
%                       minDelayMs
%                       maxDelayMs
%                       signal2Noise
%                       baseWindow
%                       minDelaySamples
%                       idxDetectStart
%                       detectStartMs
%                       baseSlopeNoise
%                       slopeThreshold
%                       nSpikesTotal
%                       idxFirstSpike
%                       firstSpikeMs
%                       vMin
%                       vMax
%                       vRange
%                       slopeMin
%                       slopeMax
%                       slopeRange
%                   specified as a scalar structure
%       spikesData  - Detected spikes data, with fields:
%                       slopes
%                       isPeakSlope
%                       isSpike
%                       idxSpikes
%                       spikeTimesMs
%                   specified as a scalar structure
% Arguments:
%       vVec        - voltage vector
%                   must be a TODO
%       siMs        - sampling interval in ms
%                   must be a TODO
%       varargin    - 'FiltFreq': cutoff frequency(ies) (Hz or normalized) 
%                                   for a bandpass filter
%                   must be a numeric a two-element vector
%                   default == none
%                   - 'tVec': original time vector
%                   must be a numeric array or a cell array of numeric arrays
%                   default == created from siMs and number of samples
%                   - 'StimVec': stimulation vector
%                   must be a numeric array or a cell array of numeric arrays
%                   default == []
%                   - 'IdxStimStart': index of stimulation start
%                   must be a positive integer scalar
%                   default == detected from StimVec
%                   - 'StimStartMs': index of stimulation start
%                   must be a positive integer scalar
%                   default == tVec(idxStimStart)
%                   - 'BaseWindow': baseline window for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == [tVec(1), stimStartMs]
%                   - 'MinDelayMs': minimum delay after stim start (ms)
%                   must be a nonnegative scalar
%                   default == 25 ms
%                   - 'MaxDelayMs': maximum delay after stim start (ms)
%                   must be a nonnegative scalar
%                   default == none
%                   - 'IdxDetectStart': index of detection start
%                   must be a positive integer scalar
%                   default == idxStimStart + minDelayMs / siMs
%                   - 'IdxDetectEnd': index of detection start
%                   must be a positive integer scalar
%                   default == idxStimStart + maxDelayMs / siMs
%                   - 'Signal2Noise': signal-to-noise ratio
%                   must be a positive scalar
%                   default == 0.5 * (max(slopes(idxDetectStart:end)) 
%                                       / baseSlopeNoise)
%
% Requires:
%       cd/compute_baseline_noise.m
%       cd/compute_time_window.m
%       cd/create_error_for_nargin.m
%       cd/create_time_vectors.m
%       cd/create_logical_array.m
%       cd/find_closest.m
%       cd/ispositiveintegerscalar.m
%       cd/parse_stim.m
%
% Used by:
%       cd/compute_oscillation_duration.m
%       cd/parse_multiunit.m

% File History:
% 2019-05-03 Moved from parse_multiunit.m
% 2019-05-04 Added input parser
% 2019-05-14 Added 'FiltFreq' as an optional argument
% 2019-05-14 Added 'MaxDelayMs' as an optional argument
% 2019-05-30 Changed signal-to-noise default to be dependent on maximum slope
% 2019-07-22 Added maxRangeOfInterestMs and fixed MaxDelayMs
% 2020-01-06 MB - Now saves tVec in spikesData
% 2020-01-07 Fixed usage of create_time_vectors.m
% 2020-01-07 Now detects stimulation start if 'StimVec' is passed in
%               and idxStimStart not provided
% 2020-01-07 Changed default baseWindow to be [tVec(1), stimStartMs]
% 2020-01-07 Now returns error if idxStimStart is 1 and 
%               baseWindow not provided
% 2020-01-07 Added 'StimStartMs' as an optional argument
% TODO: Finish documentation
% 

%% Hard-coded parameters
MS_PER_S = 1000;
relSnrThres2Max = 0.1;
maxRangeOfInterestMs = 10000;   % 1000 ms or 10 seconds
idxEndOfInterest = [];          % set later

%% Default values for optional arguments
filtFreqDefault = NaN;          % set later
tVecDefault = [];               % set later
stimVecDefault = [];            % set later
baseWindowDefault = [];         % set later
idxStimStartDefault = [];       % set later
stimStartMsDefault = [];        % set later
minDelayMsDefault = 25;         % 25 ms
maxDelayMsDefault = [];         % set later
idxDetectStartDefault = [];     % set later
idxDetectEndDefault = [];       % set later
signal2NoiseDefault = [];       % set later


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vVec', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vVec must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'siMs', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FiltFreq', filtFreqDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'tVec', tVecDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['tVec must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'StimVec', stimVecDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['StimVec must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'BaseWindow', baseWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['BaseWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'IdxStimStart', idxStimStartDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['IdxStimStart must be either empty ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'StimStartMs', stimStartMsDefault, ...
    @(x) assert(isempty(x) || isscalar(x), ...
                ['StimStartMs must be either empty ', ...
                    'or a numeric scalar!']));
addParameter(iP, 'MinDelayMs', minDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'MaxDelayMs', maxDelayMsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'IdxDetectStart', idxDetectStartDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'IdxDetectEnd', idxDetectEndDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'Signal2Noise', signal2NoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, vVec, siMs, varargin{:});
filtFreq = iP.Results.FiltFreq;
tVec = iP.Results.tVec;
stimVec = iP.Results.StimVec;
baseWindow = iP.Results.BaseWindow;
idxStimStart = iP.Results.IdxStimStart;
stimStartMs = iP.Results.StimStartMs;
minDelayMs = iP.Results.MinDelayMs;
maxDelayMs = iP.Results.MaxDelayMs;
idxDetectStart = iP.Results.IdxDetectStart;
idxDetectEnd = iP.Results.IdxDetectEnd;
signal2Noise = iP.Results.Signal2Noise;

%% Preparation
% Create time vectors
if isempty(tVec)
    % Count the number of samples
    nSamples = count_samples(vVec);
       
    % Create time vector(s)
    tVec = create_time_vectors(nSamples, 'SamplingIntervalMs', siMs, ...
                                'TimeUnits', 'ms');
end

% Decide on the stimulation start index and time if not provided
if isempty(idxStimStart) && isempty(stimStartMs)
    if ~isempty(stimVec)
        % Detect the stimulation start index from the provided stim vector
        stimParams = parse_stim(stimVec, 'SamplingIntervalMs', siMs);
        idxStimStart = stimParams.idxStimStart;
    else
        % Assume stimulation started at index 1
        idxStimStart = 1;
    end
    stimStartMs = tVec(idxStimStart);
elseif ~isempty(idxStimStart) && isempty(stimStartMs)
    stimStartMs = tVec(idxStimStart);
elseif isempty(idxStimStart) && ~isempty(stimStartMs)
    idxStimStart = find_closest(tVec, stimStartMs);
else
    if ~isequal(stimStartMs, tVec(idxStimStart))
        error('stimStartMs and tVec(idxStimStart) does not match!');
    end
end

% Decide on the baseline window
if isempty(baseWindow)
    if idxStimStart == 1
        error(['Please pass in a baseline window if stimulation starts ', ...
                'at the beginning of the vector!']);
    else
        baseWindow = compute_time_window(tVec, 'TimeEnd', stimStartMs);
    end
end

% Find the starting index for detecting a spike
if isempty(idxDetectStart)
    % Compute the minimum delay in samples
    minDelaySamples = round(minDelayMs ./ siMs);

    % Find the starting index for detecting a spike
    idxDetectStart = max(1, idxStimStart + minDelaySamples);
end

% Compute the number of samples
nSamples = numel(tVec);

% Find the ending index for detecting a spike
if isempty(idxDetectEnd) && ~isempty(maxDelayMs)
    % Compute the maximum delay in samples
    maxDelaySamples = round(maxDelayMs ./ siMs);

    % Find the ending index for computing slope and value ranges
    idxDetectEnd = min(numel(tVec), idxStimStart + maxDelaySamples);
elseif ~isempty(idxDetectEnd) && isempty(maxDelayMs)
    % Compute the maximum delay in samples and in ms
    maxDelaySamples = idxDetectEnd - idxStimStart;
    maxDelayMs = maxDelaySamples .* siMs;
else
    % The last index is idxDetectEnd
    idxDetectEnd = nSamples;

    % Compute the maximum delay in samples and in ms
    maxDelaySamples = idxDetectEnd - idxStimStart;
    maxDelayMs = maxDelaySamples .* siMs;
end

% Find the ending index for computing slope and value ranges
if isempty(idxEndOfInterest)
    % Compute the maximum range of interest in samples
    maxRangeOfInterestSamples = round(maxRangeOfInterestMs ./ siMs);

    % Find the ending index for computing slope and value ranges
    idxEndOfInterest = min(max(idxStimStart + maxRangeOfInterestSamples, 1), nSamples);
end

% Find the corresponding times
detectStartMs = tVec(idxDetectStart);
detectEndMs = tVec(idxDetectEnd);
rangeOfInterestEndMs = tVec(idxEndOfInterest);

%% Do the job
% Bandpass filter if requested
if ~isnan(filtFreq)
    siSeconds = siMs / MS_PER_S;    
    vVecFilt = freqfilter(vVec, filtFreq, siSeconds, 'FilterType', 'band');
else
    vVecFilt = vVec;
end

% Compute all instantaneous slopes in uV/ms == mV/s
slopes = diff(vVecFilt) ./ siMs;

% Compute a baseline slope noise in mV/s
baseSlopeNoise = compute_baseline_noise(slopes, tVec(1:(end-1)), baseWindow);

% Compute a default signal-to-noise ratio
if isempty(signal2Noise)
    signal2Noise = 1 + relSnrThres2Max * ...
                        (max(slopes(idxDetectStart:end)) ./ baseSlopeNoise - 1);
end

% Compute a slope threshold in mV/s
slopeThreshold = baseSlopeNoise * signal2Noise;

% Determine whether each slope is a local maximum
[~, indPeakSlopes] = findpeaks(slopes);
isPeakSlope = create_logical_array(indPeakSlopes, [nSamples - 1, 1]);

% Create all indices minus 1
allIndices = transpose(1:nSamples);

% Detect spikes after idxDetectStart
isSpike = [false; slopes >= slopeThreshold] & [false; isPeakSlope] & ...
            allIndices >= idxDetectStart & allIndices <= idxDetectEnd;
idxSpikes = find(isSpike);

% Compute the overall spike count
nSpikesTotal = numel(idxSpikes);

% Index of first spike
if nSpikesTotal == 0
    idxFirstSpike = NaN;
else
    idxFirstSpike = idxSpikes(1);
end

% Store spike times
if nSpikesTotal == 0
    spikeTimesMs = [];
    firstSpikeMs = NaN;
else
    spikeTimesMs = tVec(idxSpikes);
    firstSpikeMs = spikeTimesMs(1);    
end

% Find the last index in the range of interest
idxEndOfInterest1 = min(numel(vVecFilt), idxEndOfInterest);

% Query the maximum and range of vVec after detectStartMs
vVecTrunc = vVecFilt(idxDetectStart:idxEndOfInterest1);
vMin = min(vVecTrunc);
vMax = max(vVecTrunc);
vRange = vMax - vMin;

% Find the last index in the range of interest
idxEndOfInterest2 = min(numel(slopes), idxEndOfInterest);

% Query the maximum and range of slope after detectStartMs
slopesTrunc = slopes(idxDetectStart:idxEndOfInterest2);
slopeMin = min(slopesTrunc);
slopeMax = max(slopesTrunc);
slopeRange = slopeMax - slopeMin;

%% Save in output
spikesParams.siMs = siMs;
spikesParams.idxStimStart = idxStimStart;
spikesParams.minDelayMs = minDelayMs;
spikesParams.maxDelayMs = maxDelayMs;
spikesParams.filtFreq = filtFreq;
spikesParams.signal2Noise = signal2Noise;
spikesParams.baseWindow = baseWindow;
spikesParams.idxDetectStart = idxDetectStart;
spikesParams.idxDetectEnd = idxDetectEnd;
spikesParams.idxEndOfInterest = idxEndOfInterest;
spikesParams.detectStartMs = detectStartMs;
spikesParams.detectEndMs = detectEndMs;
spikesParams.rangeOfInterestEndMs = rangeOfInterestEndMs;
spikesParams.baseSlopeNoise = baseSlopeNoise;
spikesParams.slopeThreshold = slopeThreshold;
spikesParams.nSpikesTotal = nSpikesTotal;
spikesParams.idxFirstSpike = idxFirstSpike;
spikesParams.firstSpikeMs = firstSpikeMs;
spikesParams.vMin = vMin;
spikesParams.vMax = vMax;
spikesParams.vRange = vRange;
spikesParams.slopeMin = slopeMin;
spikesParams.slopeMax = slopeMax;
spikesParams.slopeRange = slopeRange;

spikesData.vVecFilt = vVecFilt;
spikesData.tVec = tVec;
spikesData.slopes = slopes;
spikesData.isPeakSlope = isPeakSlope;
spikesData.isSpike = isSpike;
spikesData.idxSpikes = idxSpikes;
spikesData.spikeTimesMs = spikeTimesMs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

signal2NoiseDefault = 3;        

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
