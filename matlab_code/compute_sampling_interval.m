function samplingIntervals = compute_sampling_interval (timeVecs, varargin)
%% Computes sampling intervals from time vectors
% Usage: samplingIntervals = compute_sampling_interval (timeVecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compute_sampling_interval([1; 3; 4; 6])
%       compute_sampling_interval([1; 3; 4; 6], 'IsRegular', false)
%       compute_sampling_interval(repmat((1:5)', 1, 5))
%       compute_sampling_interval(repmat({(1:5)'}, 1, 5))
%
% Outputs:
%       samplingIntervals   - sampling intervals
%                           specified as a numeric vector
% Arguments:
%       timeVecs    - time vectors
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'IsRegular': whether the sampling intervals are regular
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/isnum.m
%       cd/vecfun.m
%
% Used by:
%       cd/compute_all_pulse_responses.m
%       cd/compute_average_pulse_response.m
%       cd/create_average_time_vector.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/match_time_info.m

% File History:
% 2018-11-28 Created by Adam Lu
% 2019-09-19 Added 'IsRegular' as an optional argument
% 

%% Hard-coded parameters

%% Default values for optional arguments
isRegularDefault = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'timeVecs', ...
    @(x) assert(isnum(x) || iscellnumericvector(x), ...
                ['vecs must be either a numeric vector ', ...
                    'or a cell array of numeric vectors!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'IsRegular', isRegularDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, timeVecs, varargin{:});
isRegular = iP.Results.IsRegular;

%% Do the job
samplingIntervals = ...
    vecfun(@(x) compute_sampling_interval_helper(x, isRegular), timeVecs, ...
            'UniformOutput', true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function samplingInterval = compute_sampling_interval_helper (timeVec, isRegular)
%% Computes the sampling interval from one vector

if iscell(timeVec)
    % Apply to each cell
    samplingIntervalAll = ...
        cellfun(@(x) compute_sampling_interval_helper (x, isRegular), timeVec);
    
    % Compute the man sampling interval
    samplingInterval = nanmean(samplingIntervalAll);
    return
end

if numel(timeVec) >= 2
    if isRegular
        samplingInterval = timeVec(2) - timeVec(1);
    else
        samplingInterval = nanmean(diff(timeVec));
    end
else
    samplingInterval = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if iscell(timeVecs) 
    samplingIntervals = ...
        cellfun(@(x) compute_sampling_interval_helper(x, isRegular), timeVecs);
elseif size(timeVecs, 2) > 1
    samplingIntervals = ...
        vecfun(@(x) compute_sampling_interval_helper(x, isRegular), timeVecs);
else
    samplingIntervals = compute_sampling_interval_helper(timeVecs, isRegular);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%