function varargout = save_params (paramsTable, varargin)
%% Saves parameters to a file
% Usage: fullPath = save_params (paramsTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       fullPath    - full path to file
%                   must be a string scalar or a character vector
% Arguments:    
%       paramsTable - a table for all the parameters
%                   specified as a 2-d table
%       varargin    - 'FileName' - file name to use
%                   must be a string scalar or a character vector
%                   default == strcat(create_time_stamp, '_params.csv')
%                   - 'OutFolder': directory to place parameters file
%                   must be a string scalar or a character vector
%                   default == pwd
%
% Requires:
%       cd/create_time_stamp.m
%       cd/construct_fullpath.m
%       cd/issheettype.m
%
% Used by:
%       cd/m3ha_fminsearch3.m
%       cd/m3ha_log_errors_params.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_pfiles2csv.m

% File History:
% 2018-10-16 Created by Adam Lu
% 2018-10-21 The first argument is now a parameters table
% 2018-10-31 Added 'OutFolder' as an optional parameter
% TODO: Make 'Suffix' an optional parameter
% 

%% Hard-coded parameters
suffixDefault = '_params.csv';

% TODO: Make optional parameter
verbose = false;

%% Default values for optional arguments
fileNameDefault = strcat(create_time_stamp, suffixDefault);
outFolderDefault = '';      % set in construct_fullpath.m

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
addRequired(iP, 'paramsTable', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FileName', fileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    

% Read from the Input Parser
parse(iP, paramsTable, varargin{:});
fileName = iP.Results.FileName;
outFolder = iP.Results.OutFolder;

%% Preparation
% Construct the full path
%   TODO: Expand to accept optional Suffix', etc.
fullPath = construct_fullpath(fileName, 'Directory', outFolder);

% Get the file extension
[~, ~, fileExt] = fileparts(fullPath);

%% Do the job
if verbose
    fprintf('Saving file %s ... \n', fullPath);
end
if issheettype(fileExt)
    % Write the table to the spreadsheet file
    writetable(paramsTable, fullPath, ...
               'WriteVariableNames', true, ...
               'WriteRowNames', true);
else
    % Save as a matfile
    save(fullPath, 'paramsTable', '-v7.3');
end

%% Outputs
varargout{1} = fullPath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

function save_params (fileName, paramNames, paramValues, ...
                                paramLowerBounds, paramUpperBounds, varargin)
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       paramNames  - parameter names
%                   must be a cell vector of character vectors
%       paramValues - parameter values
%                   must be a numeric vector
%       paramLBs    - parameter lower bounds
%                   must be a numeric vector
%       paramUBs    - parameter upper bounds
%                   must be a numeric vector

addRequired(iP, 'paramNames', ...
    @(x) iscellstr(x));
addRequired(iP, 'paramValues', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'paramLowerBounds', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'paramUpperBounds', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

parse(iP, fileName, paramNames, paramValues, ...
        paramLowerBounds, paramUpperBounds, varargin{:});

% Check relationships between arguments
%       TODO cd/isequallength.m
% TODO: Make a function isequallength.m
if length(unique([numel(paramNames), length(paramValues), ...
                  length(paramLowerBounds), length(paramUpperBounds)])) ~= 1
    fprintf(['Cannot save parameters because names, values, ', ...
             'lower bounds and upper bounds are not all the same length!!\n']);
    paramsTable = [];
    return
end

% Force as a column and rename the variables for the header
Name = force_column_cell(paramNames);
Value = force_column_vector(paramValues);
LowerBound = force_column_vector(paramLowerBounds);
UpperBound = force_column_vector(paramUpperBounds);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
