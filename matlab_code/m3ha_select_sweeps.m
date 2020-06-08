function [swpInfo, fileBasesToUse] = m3ha_select_sweeps (varargin)
%% Selects file bases and row indices in swpInfo that will be used
% Usage: [swpInfo, fileBasesToUse] = m3ha_select_sweeps (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       swpInfo         - same as input table but with updated variables:
%                           toUse   - whether the sweep is to be fitted
%       fileBasesToUse  - file bases of sweeps to use
%                       specified as a cell array
%
% Arguments:
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'SwpInfo': a table of sweep info, with each row named by 
%                               the matfile base containing the raw data
%                   must a 2D table with row names being file bases
%                       and with the fields:
%                       cellidrow   - cell ID
%                       prow        - pharmacological condition
%                       grow        - conductance amplitude scaling
%                       toUse       - whether the sweep is to be used (optional)
%                   default == m3ha_load_sweep_info
%                   - 'RowsToUse' - row indices or row names in swpInfo 
%                                       of sweeps to fit
%                   must be a positive integer vector, a string array 
%                       or a cell array of character vectors
%                   default == []
%                   - 'DataMode': data mode
%                   must be a one of:
%                       0 - all data
%                       1 - all of g incr = 100%, 200%, 400%
%                       2 - all of g incr = 100%, 200%, 400% 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                       3 - all data 
%                               but exclude cell-pharm-g_incr sets 
%                               containing problematic sweeps
%                   default == 2
%                   - 'CasesDir' - the directory that contains 
%                                   'TAKE_OUT_*' folders with special cases
%                   must be a directory
%                   default == ~/m3ha/data_dclamp/take4/special_cases
%                   - 'PharmConditions': pharmacological condition(s)
%                                           to restrict to
%                   must be empty or some of:
%                       1 - control
%                       2 - GAT1 blockade
%                       3 - GAT3 blockade
%                       4 - dual blockade
%                   default == no restrictions
%                   - 'GIncrConditions': conductance amplitude condition(s) (%)
%                                           to restrict to
%                   must be empty or some of: 25, 50, 100, 200, 400, 800
%                   default == no restrictions
%                   - 'VHoldConditions': holding potential condition(s) (mV)
%                                           to restrict to
%                   must be empty or some of: -60, -65, -70
%                   default == no restrictions
%                   - 'CellIds': original cell ID(s) to restrict to
%                   must be empty or integer(s) between 1 and 49
%                   default == no restrictions
%                   - 'CellNames': original cell name(s) to restrict to
%                   must be empty or a character vector
%                       or a cell array of character vectors
%                   default == no restrictions
%                   - 'RepNums': repetition number(s) to restrict to
%                   must be empty or integer(s) between 1 and 5
%                   default == no restrictions
%                   
% Requires:
%       cd/extract_fileparts.m
%       cd/has_same_attributes.m
%       cd/is_var_in_table.m
%       cd/m3ha_find_files_to_take_out.m
%       cd/m3ha_load_sweep_info.m
%
% Used by:
%       cd/m3ha_compute_and_plot_statistics.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure04.m
%       cd/m3ha_select_cells.m
%       cd/m3ha_select_raw_traces.m

% TODO:
%       /media/adamX/m3ha/data_dclamp/dclampPassiveFitter.m
%       /media/adamX/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%       /media/adamX/m3ha/data_dclamp/PlotCorrelations.m
%       /media/adamX/m3ha/data_dclamp/dclampdatalog_analyze.m

% File History:
% 2018-11-18 Adapted from m3ha_find_ind_to_fit()
% 2018-12-06 Now adds a toUse column to sweep info
% 2019-11-26 Added dataMode == 0
% 2019-11-27 Added 'PharmConditions', 'GIncrConditions', 'VHoldConditions', ...
%               & 'CellIds' & 'RepNums' as optional arguments
% 2019-12-18 Added dataMode == 3
% 2019-12-18 Added 'CellNames' as an optional argument
% 2019-12-21 Added 'RowsToUse' as an optional argument
% 2019-12-21 Changed default dataMode to 0

%% Hard-coded parameters
pharmStr = 'prow';
gIncrStr = 'grow';
vHoldStr = 'vrow';
cellIdStr = 'cellidrow';
repNumStr = 'swpnrow';
toUseStr = 'toUse';
attributesToMatch = {cellIdStr, pharmStr, gIncrStr};
defaultCasesDir = '~/m3ha/data_dclamp/take4/special_cases';

%% Default values for optional arguments
verboseDefault = true;             % print to standard output by default
swpInfoDefault = [];
rowsToUseDefault = [];
dataModeDefault = 0;
casesDirDefault = '';               % set later
pharmConditionsDefault = [];
gIncrConditionsDefault = [];
vHoldConditionsDefault = [];
cellIdsDefault = [];
cellNamesDefault = {};
repNumsDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SwpInfo', swpInfoDefault, ...
    @(x) validateattributes(x, {'table'}, {'2d'}));
addParameter(iP, 'RowsToUse', rowsToUseDefault, ...
    @(x) assert(isempty(x) || ispositiveintegervector(x) || ...
                    iscellstr(x) || isstring(x), ...
                ['RowsToUse must be either a positive integer vector, ', ...
                    'a string array or a cell array of character arrays!']));
addParameter(iP, 'DataMode', dataModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'}));
addParameter(iP, 'CasesDir', casesDirDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    
addParameter(iP, 'PharmConditions', pharmConditionsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GIncrConditions', gIncrConditionsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'VHoldConditions', vHoldConditionsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'CellIds', cellIdsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'CellNames', cellNamesDefault, ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
        ['CellNames must be a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'RepNums', repNumsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Read from the Input Parser
parse(iP, varargin{:});
verbose = iP.Results.Verbose;
swpInfo = iP.Results.SwpInfo;
rowsToUse = iP.Results.RowsToUse;
dataMode = iP.Results.DataMode;
casesDir = iP.Results.CasesDir;
pharmConditions = iP.Results.PharmConditions;
gIncrConditions = iP.Results.GIncrConditions;
vHoldConditions = iP.Results.VHoldConditions;
cellIds = iP.Results.CellIds;
cellNames = iP.Results.CellNames;
repNums = iP.Results.RepNums;

%% Preparation
% Read in swpInfo if not provided
if isempty(swpInfo)
    swpInfo = m3ha_load_sweep_info;
end

% Extract the conductance amplitude scaling percentages
gIncrRow = swpInfo.(gIncrStr);

% Extract all the file bases
fileBases = swpInfo.Properties.RowNames;

% Set default special cases directory
if isempty(casesDir)
    casesDir = defaultCasesDir;
end

%% Do the job
% Print message
if verbose
    fprintf('Selecting the sweeps to use ... \n');
end

% Make sure a toUse column exists in the table
if ~is_var_in_table(toUseStr, swpInfo)
    % Record whether toUse column exists before
    toUseExistsBefore = false;

    % Initialize all sweeps to be used
    toUse = true(height(swpInfo), 1);

    % Place in swpInfo
    swpInfo = addvars(swpInfo, toUse);
else
    toUseExistsBefore = true;
end

% Replace with user's choice of rows to use
if ~isempty(rowsToUse)
    if toUseExistsBefore
        % Print message
        fprintf('Warning: toUse column already exists in swpInfo.\n');
        fprintf('RowsToUse will replace it!\n');
    end
    swpInfo{rowsToUse, toUseStr} = true;
    swpInfo{~rowsToUse, toUseStr} = false;
end

% Read the initial state of the toUse column
toUse = swpInfo.(toUseStr);

% Determine whether each sweep has
%   conductance amplitudes with 100%, 200% or 400% scaling 
%   (these are present in all experiments)
if dataMode == 1 || dataMode == 2
    isGIncrToUse = gIncrRow == 100 | gIncrRow == 200 | gIncrRow == 400;
end

% Get the file names of files to take out from specialCasesDir
%   Note: these are labelled with 'TAKE_OUT_*' and were
%           the result of voting by blab
if dataMode == 2 || dataMode == 3
    % Find all files to take out
    fileNamesToTakeOut = m3ha_find_files_to_take_out('CasesDir', casesDir);

    % If there are no files, return with message
    if isempty(fileNamesToTakeOut)
        fprintf('There are no files to take out under %s!\n', casesDir);
        fileBasesToUse = {};
        return
    end

    % Extract file bases
    fileBasesToTakeOut = extract_fileparts(fileNamesToTakeOut, 'base');

    % Determine whether each sweep as the same cell ID, pharm condition
    %   and conductance amplitude scaling as any of the files to take out
    isNotToUse = has_same_attributes(swpInfo, fileBasesToTakeOut, ...
                                        attributesToMatch);
end

% Find the sweep indices to fit
if dataMode == 1
    toUse = toUse & isGIncrToUse;
elseif dataMode == 2
    toUse = toUse & isGIncrToUse & ~isNotToUse;
elseif dataMode == 3
    toUse = toUse & ~isNotToUse;
end

%% Update toUse according to the conditions requested
toUse = restrict_to_conditions(toUse, swpInfo, pharmStr, pharmConditions);
toUse = restrict_to_conditions(toUse, swpInfo, gIncrStr, gIncrConditions);
toUse = restrict_to_conditions(toUse, swpInfo, vHoldStr, vHoldConditions);
toUse = restrict_to_conditions(toUse, swpInfo, cellIdStr, cellIds);
toUse = restrict_to_conditions(toUse, swpInfo, repNumStr, repNums);

%% Update toUse if cellNames provided
if ~isempty(cellNames)
    toUse = toUse & contains(fileBases, cellNames);
end

%% Output results
% Update the toUse column
swpInfo.(toUseStr) = toUse;

% Find the file bases to fit
fileBasesToUse = fileBases(toUse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function toUse = restrict_to_conditions (toUse, swpInfo, varStr, values)
%% Updates toUse according to the provided values for a condition

% Update toUse according to the condition
if ~isempty(values)
    % Extract all values for the condition
    varAll = swpInfo.(varStr);

    % Restrict to the requested values
    toUse = toUse & ismember(varAll, values);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Print message
fprintf('Warning: toUse column already exists in swpInfo.\n');
fprintf('RowsToUse will be ignored!\n');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
