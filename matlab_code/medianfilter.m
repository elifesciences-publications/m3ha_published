function vecsFilt = medianfilter (vecs, varargin)
%% Applies a median filter to vectors
% Usage: vecsFilt = medianfilter (vecs, filtWidth, si, varargin)
% Explanation:
%       Same as medfilt1() but with option of using a window in time units
%
% Example(s):
%       medianfilter(magic(4))
%       medianfilter(magic(4), 10, 10)
%       medianfilter(magic(4), 10, 2)
%       medianfilter(magic(4), 10, 1)
%
% Outputs:
%       vecsFilt    - filtered vector(s)
%                   specified as a numeric array
% Arguments:
%       vecs        - vector(s) to median filter
%                   must be a numeric array
%       filtWidth   - (opt) filter window width in time units
%                           must be a numeric scalar
%       si          - (opt) sampling interval in time units
%                           must be a numeric scalar
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the medfilt1() function
%
% Requires:
%       cd/apply_to_nonnan_part.m
%       cd/create_error_for_nargin.m
%       cd/find_nearest_odd.m
%       cd/isnum.m
%       cd/iscellnumericvector.m
%       cd/struct2arglist.m
%       cd/vecfun.m
%
% Used by:
%       cd/find_passive_params.m
%       cd/identify_CI_protocol.m
%       cd/parse_lts.m

% File History:
% 2019-01-14 Created by Adam Lu
% 2019-11-18 Now sets nanflag to be 'omitnan'
% 2019-11-18 Now sets padding to be 'truncate'
% 

%% Hard-coded parameters

%% Default values for optional arguments
filtWidthDefault = 3;       % 3 sample points by default
siDefault = 1;              % treat filtWidth as samples by default
% param1Default = [];       % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vec', ...
    @(x) assert(isnum(x) || iscellnumericvector(x), ...
                ['Vectors must be a numeric array ', ...
                    'or a cell array of numeric vectors!']));

% Add optional inputs to the Input Parser
addOptional(iP, 'filtWidth', filtWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addOptional(iP, 'si', siDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, vecs, varargin{:});
filtWidth = iP.Results.filtWidth;
si = iP.Results.si;
% param1 = iP.Results.param1;

% Keep unmatched arguments for the medfilt1() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Calculate the median filter window width in samples
%   Note: Round down to the nearest odd integer to preserve values!!
%           However, must be >= 1
filtWidthSamples = find_nearest_odd(filtWidth / si, 'Direction', 'down');

% Median filter vectors
vecsFilt = vecfun(@(x) medianfilter_helper(x, filtWidthSamples, ...
                                                otherArguments), ...
                    vecs, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecFilt = medianfilter_helper(vec, filtWidthSamples, otherArguments)
%% Applies a median filter to one vector

% Create a custom medfilt1 function
myFun = @(x) medfilt1(x, filtWidthSamples, otherArguments{:}, ...
                        'omitnan', 'truncate');

% Apply the custom smooth function to just the non-NaN part
%   Note: If this filter is applied to leading and trailing NaNs,
%           values become extrapolated, which is most likely not desired
vecFilt = apply_to_nonnan_part(myFun, vec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%