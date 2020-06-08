function paramsTable = update_param_values (paramsTable, varargin)
%% Updates a parameters table with new values
% Usage: paramsTable = update_param_values (paramsTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       paramsTable - updated parameters table
%                   specified as a table
%
% Arguments:
%       paramsTable - a parameters table with the row names being the 
%                       parameter names and a 'Value' column
%                   must be a table
%       varargin    - name-value pairs that matches row names in the table
%                   must be a list of string-numeric pairs or a structure
%                   - 'IgnoreRange': whether to ignore parameter range
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/is_contained_in.m
%       cd/rmfield_custom.m
%
% Used by:    
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2018-10-31 Created by Adam Lu
% 2018-11-14 Now checks if each parameter exists in rownames
% 2018-11-14 Now check bounds if 'UpperBound' and 'LowerBound' fields exist
% 2020-03-12 Added 'IgnoreRange' as an optional argument

%% Hard-coded parameters
upperBoundStr = 'UpperBound';
lowerBoundStr = 'LowerBound';

%% Default values for optional arguments
% param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Allow extraneous options: these will be parameters to be updated
iP.KeepUnmatched = true;

% Add required inputs to the Input Parser
addRequired(iP, 'paramsTable', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, paramsTable, varargin{:});
% param1 = iP.Results.param1;

% Read the name-value pairs unmatched by the input parser
%   These will be name-value pairs for the parameters
allParamValuePairs = iP.Unmatched;

% Deal with parameter-value pairs for this function
% TODO: Use parse_and_remove_from_struct.m
if isfield(allParamValuePairs, 'IgnoreRange')
    ignoreRange = allParamValuePairs.IgnoreRange;
    allParamValuePairs = rmfield_custom(allParamValuePairs, 'IgnoreRange');
else
    ignoreRange = false;
end

%% Preparation
% Get all parameter names from the table
paramsTableRowNames = paramsTable.Properties.RowNames;

% Get all variable names from the table
paramsTableVariableNames = paramsTable.Properties.VariableNames;

%% Do the job
% Get all parameter names in a cell array
paramNames = fieldnames(allParamValuePairs);

% Check if all parameter names exist in the table
if ~is_contained_in(paramNames, paramsTableRowNames)
    fprintf('Original parameters table returned!\n');
    return
end

% Get all parameter values in a cell array
paramValuesCell = struct2cell(allParamValuePairs);

% If upper bounds and lower bounds exist, 
%   check whether the values are within bounds
if ~ignoreRange && is_contained_in({upperBoundStr, lowerBoundStr}, ...
                    paramsTableVariableNames, 'SuppressOutput', true)
    % Convert to a numeric array if poss
    paramValues = cell2mat(paramValuesCell);

    % Extract the upper and lower bounds for each parameter
    upperBounds = paramsTable{paramNames, upperBoundStr};
    lowerBounds = paramsTable{paramNames, lowerBoundStr};

    % Check if each value is within bounds
    [isAllWithinBounds, isWithinBound] = ...
        check_within_bounds(paramValues, lowerBounds, upperBounds);

    % Only update the parameter values within bounds
    if ~isAllWithinBounds
        paramNames = paramNames(isWithinBound);
        paramValuesCell = paramValuesCell(isWithinBound);
    end
end

% Replace all parameters with new values
paramsTable(paramNames, {'Value'}) = paramValuesCell;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if ~check_within_bounds(paramValues, lowerBounds, upperBounds)
    fprintf('Original parameters table returned!\n');
    return
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
