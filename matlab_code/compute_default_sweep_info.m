function [baseWindow, fitWindow, baseNoise, sweepWeights] = ...
                compute_default_sweep_info (tVecs, data, varargin)
%% Computes default windows, noise, weights and errors
% Usage: [baseWindow, fitWindow, baseNoise, sweepWeights] = ...
%               compute_default_sweep_info (tVecs, data, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       baseWindow  - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       tVecs       - time vector(s) for plotting
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       data        - data vectors(s)
%                   Note: If a cell array, each element must be a vector
%                         If a non-vector array, each column is a vector
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'BaseWindow': baseline window for each trace
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == first half of the trace
%                   - 'FitWindow': time window to fit for each trace
%                   must be a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == second half of the trace
%                   - 'BaseNoise': baseline noise value(s)
%                   must be a numeric vector
%                   default == if data is not empty, 
%                               apply compute_baseline_noise.m 
%                               otherwise, all sweeps have noise == 1
%                   - 'SweepWeights': sweep weights for averaging
%                   must be empty or a numeric vector with length == nSweeps
%                   default == 1 ./ baseNoise
%
% Requires:
%       cd/argfun.m
%       cd/compute_baseline_noise.m
%       cd/force_column_vector.m
%       cd/iscellnumeric.m
%       cd/isnumericvector.m
%       cd/match_format_vector_sets.m
%
% Used by:
%       cd/compute_single_neuron_errors.m
%       cd/m3ha_import_raw_traces.m
%       cd/plot_fitted_traces.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-11-01 Moved from plot_fitted_traces.m
% 2018-11-28 Default window(s) for fitting must now be after baseline window(s)
% 

%% Default values for optional arguments
baseWindowDefault = [];         % set later
fitWindowDefault = [];          % set later
baseNoiseDefault = [];          % set later
sweepWeightsDefault = [];       % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'tVecs', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'data', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vec1s must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'BaseWindow', baseWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['BaseWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'FitWindow', fitWindowDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['FitWindow must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'BaseNoise', baseNoiseDefault, ...
    @(x) assert(isnumericvector(x), 'BaseNoise must be a numeric vector!'));
addParameter(iP, 'SweepWeights', sweepWeightsDefault, ...
    @(x) assert(isnumericvector(x), 'SweepWeights must be a numeric vector!'));

% Read from the Input Parser
parse(iP, tVecs, data, varargin{:});
baseWindow = iP.Results.BaseWindow;
fitWindow = iP.Results.FitWindow;
baseNoise = iP.Results.BaseNoise;
sweepWeights = iP.Results.SweepWeights;

%% Preparation
% If both tVecs and data are empty, return empty vectors
if isempty(tVecs) && isempty(data)
    return
end

% Makes sure any numeric vector is a column
[baseWindow, fitWindow, baseNoise, sweepWeights] = ...
    argfun(@(x) force_column_vector(x, 'IgnoreNonVectors', true), ...
            baseWindow, fitWindow, baseNoise, sweepWeights);

% Matches tVecs and data so that they are both cell arrays 
%   of the same number of column vectors
[tVecs, data] = match_format_vector_sets(tVecs, data, 'ForceCellOutputs', true);

%% Do the job
% Find the minimum and maximum times and center times
if isempty(baseWindow) || isempty(fitWindow)
    % Find the minimum times as a column vector
    minTimes = cellfun(@min, tVecs);

    % Find the maximum times as a column vector
    maxTimes = cellfun(@max, tVecs);

    % Find the center times as a column vector
    centerTimes = (minTimes + maxTimes) / 2;
end

% Set default baseline window(s)
if isempty(baseWindow)
    % See if fitWindow is provided
    if isempty(fitWindow)
        baseWindow = transpose([minTimes, centerTimes]);
    else
        % Get the minimum time for each fitting window
        if iscell(fitWindow)
            minFitTimes = cellfun(@min, fitWindow);
        else
            minFitTimes = transpose(min(fitWindow, [], 1));
        end

        % Match minFitTimes with minTimes
        minFitTimes = match_row_count(minFitTimes, length(minTimes));

        % Get the default baseline window(s)
        baseWindow = transpose([minTimes, minFitTimes]);
    end
end

% Set default window(s) for fitting 
%   Note: assumes baseline window already set before
if isempty(fitWindow)
    % Get the maximum time for each baseline window
    if iscell(baseWindow)
        maxBaseTimes = cellfun(@max, baseWindow);
    else
        maxBaseTimes = transpose(max(baseWindow, [], 1));
    end

    % Match maxBaseTimes with maxTimes
    maxBaseTimes = match_row_count(maxBaseTimes, length(maxTimes));

    % Get the default fit window(s)
    fitWindow = transpose([maxBaseTimes, maxTimes]);
end

% Re-compute baseline noise if not provided
if isempty(baseNoise)
    if ~isempty(data)
        baseNoise = compute_baseline_noise(data, tVecs, baseWindow);
    else
        % Compute the number of sweeps
        nSweeps = numel(tVecs);

        % Assume a noise of 1
        baseNoise = ones(nSweeps, 1);
    end
end

% Re-compute sweep weights if not provided
if isempty(sweepWeights)
    % Compute sweep weights
    sweepWeights = 1 ./ baseNoise;

    % Normalize so that it sums to one
    sweepWeights = sweepWeights / sum(sweepWeights);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fitWindow = transpose([centerTimes, maxTimes]);

% Set default window(s) for fitting (must be after baseline window)

% Set default baseline window(s)
if isempty(baseWindow)
    baseWindow = transpose([minTimes, centerTimes]);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
