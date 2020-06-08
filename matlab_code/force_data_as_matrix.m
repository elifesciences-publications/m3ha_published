function [dataValues, varargout] = force_data_as_matrix (data, varargin)
%% Forces data values as a numeric matrix where each group is a column and each row is a sample
% Usage: [dataValues, groupLabels, sampleLabels] = ...
%               force_data_as_matrix (data, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples;
%       [v, g, s] = force_data_as_matrix(myTableNumeric)
%       [v, g, s] = force_data_as_matrix(transpose_table(myTableNumeric))
%
% Outputs:
%       dataValues      - data matrix where each column is a group
%                           and each row is a sample
%                       specified as a numeric matrix
%       groupLabels     - label for each group
%                       specified as a cell array of character vectors
%       sampleLabels    - label for each sample
%                       specified as a cell array of character vectors
%
% Arguments:
%       data        - data vectors or data table or path to data table
%                   Note: The dimension with fewer elements is taken as 
%                           the parameter
%                   must be a table or a numeric array
%                       or a cell array of numeric vectors
%                       or a string scalar or a character vector
%       varargin    - 'ForceVectorAsColumns': whether to convert row vectors 
%                                               to column vectors
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/create_trace_plot_movie.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%
% Used by:
%       cd/create_trace_plot_movie.m
%       cd/plot_chevron.m
%       cd/plot_chevron_bar_inset.m
%       cd/plot_violin.m
%       cd/test_difference.m

% File History:
% 2019-12-30 Moved from plot_chevron.m
% 2020-05-04 Added 'ForceVectorAsColumns' as an optional argument

%% Hard-coded parameters

%% Default values for optional arguments
forceVectorAsColumnsDefault = false;

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
addRequired(iP, 'data', ...
    @(x) validateattributes(x, {'numeric', 'cell', 'table', ...
                            'char', 'string'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceVectorAsColumns', forceVectorAsColumnsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, data, varargin{:});
forceVectorAsColumns = iP.Results.ForceVectorAsColumns;

%% Do the job
% Read from file
if ischar(data) || isstring(data)
    data = readtable(data);
end

% Read from the table or matrix
if istable(data)
    % Extract values
    dataValues = table2array(data);

    % Make sure each column is a group and transpose if necessary
    %   Note: This assumes that the number of groups are not greater than
    %           the number of samples
    nRowsOrig = size(dataValues, 1);
    nColumnsOrig = size(dataValues, 2);
    if nRowsOrig < nColumnsOrig
        dataValues = transpose(dataValues);
        tableTransposed = true;
    else
        tableTransposed = false;
    end
else
    % Force as a matrix
    dataValues = force_matrix(data, 'AlignMethod', 'leftadjustpad');
end

% Force as a column vector if requested
if forceVectorAsColumns && isrow(dataValues)
    dataValues = force_column_vector(dataValues);
end

% Decide on group labels
if nargout >= 2
    if istable(data) && tableTransposed && ~isempty(data.Properties.RowNames)
        % Use the row names if transposed
        groupLabels = force_column_vector(data.Properties.RowNames);
    elseif istable(data)
        % Use the variable names
        groupLabels = force_column_vector(data.Properties.VariableNames);
    else
        % Count the number of groups
        nGroups = size(dataValues, 2);

        % Label as group1, group2, ...
        groupLabels = create_labels_from_numbers(1:nGroups, 'Prefix', 'group');
    end
end

% Decide on sample labels
if nargout >= 3
    % Decide on sample labels
    if istable(data) && tableTransposed
        % Use the variable names if transposed
        sampleLabels = force_column_vector(data.Properties.VariableNames);
    elseif istable(data) && ~isempty(data.Properties.RowNames)
        % Use the row names if exist
        sampleLabels = force_column_vector(data.Properties.RowNames);
    else
        % Count the number of samples
        nSamples = size(dataValues, 1);

        % Create labels
        sampleLabels = create_labels_from_numbers(1:nSamples, 'Prefix', 'sample');
    end
end

%% Outputs
if nargout >= 2
    varargout{1} = groupLabels;
end
if nargout >= 3
    varargout{2} = sampleLabels;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%