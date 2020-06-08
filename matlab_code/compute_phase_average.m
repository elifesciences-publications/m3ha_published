function [phaseAverage, indSelected] = compute_phase_average (values, varargin)
%% Computes the average of values over the last of a phase
% Usage: [phaseAverage, indSelected] = compute_phase_average (values, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       compute_phase_average(randn(30, 1), 'PhaseBoundaries', [15.5, 25.5])
%
% Outputs:
%       phaseAverage    - average value over the last of the phase
%                       specified as a numeric scalar
%
% Arguments:
%       values      - values for a phase 
%                       or all values if phaseBoundaries provided
%                   must be a numeric vector
%       varargin    - 'ReturnLastTrial': whether to return last attempt
%                                           instead of NaNs if criteria not met
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Indices': indices for the subvectors to extract 
%                       Note: if provided, would override 'EndPoints'
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == set in select_similar_values.m
%                   - 'EndPoints': endpoints for the phase
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == set in select_similar_values.m
%                   - 'PhaseBoundaries': vector of phase boundaries
%                   must be a numeric vector
%                   default == none
%                   - 'PhaseNumber': phase number to average
%                   must be a positive integer scalar
%                   default == 1
%                   - 'NLastOfPhase': number of values at the last of a phase
%                   must be empty or a positive integer scalar
%                   default == 10
%                   - 'NToAverage': number of values to average
%                   must be empty or a positive integer scalar
%                   default == set in select_similar_values.m
%                   - 'SelectionMethod': the selection method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'auto'          - set in select_similar_values.m
%                       'notNaN'        - select any non-NaN value
%                       'maxRange2Mean' - select vales so that the maximum 
%                                           range is within a percentage 
%                                           of the mean
%                   default == 'auto'
%                   - 'MaxRange2Mean': maximum percentage of range versus mean
%                   must be empty or a nonnegative scalar
%                   default == set in select_similar_values.m
%                   - Any other parameter-value pair for 
%                       the select_similar_values() function
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/select_similar_values.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/parse_phase_info.m
%       cd/plot_measures.m

% File History:
% 2019-05-13 Created by Adam Lu
% 2019-05-15 Add 'SelectionMethod' as an optional argument
% 2019-05-16 Add 'ReturnLastTrial' as an optional argument
% 2019-05-16 Now uses nanmean() instead of mean()
% 2019-07-25 Added 'PhaseBoundaries' and 'PhaseNumber' as optional arguments
% TODO: Make 'AverageWindows' as an optional argument and use compute_average_windows.m
% TODO: Return value no matter what?
% 

%% Hard-coded parameters
validSelectionMethods = {'auto', 'notNaN', 'maxRange2Mean'};
defaultNLastOfPhase = 10;

%% Default values for optional arguments
returnLastTrialDefault = false; % return NaN if criteria not met by default
indicesDefault = [];            % set in extract_subvectors.m
phaseBoundariesDefault = [];    % no phase boundaries by default
phaseNumberDefault = 1;         % the first phase by default
endPointsDefault = [];          % set in select_similar_values.m
nLastOfPhaseDefault = [];       % set later
nToAverageDefault = [];         % set in select_similar_values.m
selectionMethodDefault = 'auto';% set in select_similar_values.m
maxRange2MeanDefault = [];      % set in select_similar_values.m

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
addRequired(iP, 'values', ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ReturnLastTrial', returnLastTrialDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Indices', indicesDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Indices must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'PhaseBoundaries', phaseBoundariesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PhaseNumber', phaseNumberDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['EndPoints must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'NLastOfPhase', nLastOfPhaseDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['NLastOfPhase must be either empty ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'NToAverage', nToAverageDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                ['NToAverage must be either empty ', ...
                    'or a positive integer scalar!']));
addParameter(iP, 'SelectionMethod', selectionMethodDefault, ...
    @(x) any(validatestring(x, validSelectionMethods)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) assert(isempty(x) || isscalar(x) && x >= 0, ...
                ['MaxRange2Mean must be either empty ', ...
                    'or a nonnegative scalar!']));

% Read from the Input Parser
parse(iP, values, varargin{:});
returnLastTrial = iP.Results.ReturnLastTrial;
indices = iP.Results.Indices;
phaseBoundaries = iP.Results.PhaseBoundaries;
phaseNumber = iP.Results.PhaseNumber;
endPoints = iP.Results.EndPoints;
nLastOfPhase = iP.Results.NLastOfPhase;
nToAverage = iP.Results.NToAverage;
selectionMethod = validatestring(iP.Results.SelectionMethod, ...
                                    validSelectionMethods);
maxRange2Mean = iP.Results.MaxRange2Mean;

% Keep unmatched arguments for the select_similar_values() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on nLastOfPhase
if isempty(nLastOfPhase)
    nLastOfPhase = defaultNLastOfPhase;
end

%% Compute end points for the last of phase
if isempty(endPoints) && ~isempty(phaseBoundaries)
    % Number of entries
    nEntries = numel(values);

    % Create indices
    pValues = 1:nEntries;

    % Find the left boundaries for each phase
    phaseLeftBounds = [1; phaseBoundaries(:)];

    % Find the right boundaries for each phase
    phaseRightBounds = [phaseBoundaries(:); nEntries];

    % Find the first index for the selected phase
    idxFirstOfPhase = find(pValues >= phaseLeftBounds(phaseNumber), 1, 'first');

    % Find the last index for the selected phase
    idxLastOfPhase = find(pValues <= phaseRightBounds(phaseNumber), 1, 'last');

    % Find the last index to average for the selected phase
    idxLast = idxLastOfPhase;

    % Find the first index to average for the selected phase
    idxFirst = max(idxFirstOfPhase, idxLast - nLastOfPhase + 1);

    % Find the end points for the last of phase
    endPoints = [idxFirst, idxLast];
end

%% Compute average of values over the last of phase
% Select values similar to the last phase value
[valSelected, indSelected] = ...
    select_similar_values(values, 'ReturnLastTrial', returnLastTrial, ...
                        'EndPoints', endPoints, 'Indices', indices, ...
                        'NToSelect', nToAverage, 'Direction', 'backward', ...
                        'SelectionMethod', selectionMethod, ...
                        'MaxRange2Mean', maxRange2Mean, ...
                        otherArguments{:});

% Compute the phase average
phaseAverage = nanmean(valSelected);
% phaseAverage = compute_stats(values, 'mean', 'Indices', indSelected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%