function varargout = find_first_match (candidates, array, varargin)
%% Returns the first matching index and match in an array for candidate(s)
% Usage: [index, matched] = find_first_match (candidates, array, varargin)
% Explanation:
%       TODO
%   Use find_in_strings.m to treat each element of cand
%       as a part of a candidate
%
% Example(s):
%       [index, matched] = find_first_match([2, 5, 1], 5:-1:1)
%       [index, matched] = find_first_match({'dog'; 'cat'}, ["dog"; "fly"; "cat"])
%       [index, matched] = find_first_match(["dog", "fly", "cat"], {'dog'; 'cat'})
%       [index, matched] = find_first_match(2:5, {1:3, 2:5, 6:8})
%       [index, matched] = find_first_match({6:8, 1:3}, {1:3, 2:5, 6:8})
%
% Outputs:
%       index       - the index with matching element in the array
%                   specified as a positive integer array (may contain NaN)
%       matched     - matching element in the array
%                   specified as an array of the same type as elements of array
%
% Arguments:
%       candidates  - candidates to be matched
%       array       - an array
%       varargin    - 'MatchMode': the matching mode
%                   must be an unambiguous, case-insensitive match to one of:
%                       'exact'     - cand must be identical to the members
%                       'parts'     - cand can be parts of the members
%                       'prefix'    - cand is the prefix of the members
%                       'keyword'   - cand is a part of the members
%                       'suffix'    - cand is the suffix of the members
%                       'regexp'    - cand is a regular expression
%                   default == 'parts'
%                   - 'IgnoreCase': whether to ignore differences in letter case
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/ismember_custom.m
%
% Used by:
%       cd/compute_activation_profile.m
%       cd/create_plot_movie.m
%       cd/find_matching_files.m
%       cd/match_positions.m
%       cd/m3ha_decide_on_plot_vars.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/m3ha_xolotl_plot.m
%       cd/parse_spike2_mat.m
%       cd/plot_calcium_imaging_traces.m
%       cd/plot_table_parallel.m

% File History:
% 2019-01-09 Created by Adam Lu
% 2019-01-10 Now returns matched
% 2020-01-01 Added 'prefix', 'keyword', 'suffix' as match modes
% 2020-05-13 Now accepts array as any cell array

%% Hard-coded parameters
validMatchModes = {'exact', 'parts', 'prefix', 'keyword', 'suffix', 'regexp'};

%% Default values for optional arguments
matchModeDefault = 'parts';         % can be parts by default
ignoreCaseDefault = false;          % don't ignore case by default

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
addRequired(iP, 'candidates');
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MatchMode', matchModeDefault, ...   % the search mode
    @(x) any(validatestring(x, validMatchModes)));
addParameter(iP, 'IgnoreCase', ignoreCaseDefault, ...   % whether to ignore case
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, candidates, array, varargin{:});
matchMode = validatestring(iP.Results.MatchMode, validMatchModes);
ignoreCase = iP.Results.IgnoreCase;

%% Do the job
% Find the index in array for each element in candidates
%   Note: If not found, NaN will be returned by ismember_custom.m
[~, index] = ismember_custom(candidates, array, ...
                            'MatchMode', matchMode, 'IgnoreCase', ignoreCase);

% Get all matched elements if requested
if nargout >= 2
    if iscell(array) && numel(index) == 1
        matched = array{index};
    else
        matched = array(index);
    end
end

%% Deal with outputs
varargout{1} = index;
if nargout >= 2
    varargout{2} = matched;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
