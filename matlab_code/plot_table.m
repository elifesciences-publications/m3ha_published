function handles = plot_table (myTable, varargin)
%% Plots variables (columns) in a table
% Usage: handles = plot_table (myTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       load_examples
%       plot_table(myTableNumeric)
%       plot_table(myTableNumeric, 'PlotMode', 'separate')
%       TODO: plot_table(myTableNumeric, 'PlotMode', 'parallel')
%
% Outputs:
%       handles     - structure with fields:
%                       fig - figure handle(s) for the created figure(s)
%                   specified as a scalar structure
%
% Arguments:
%       myTable     - a table with variables to plot
%                   must be a table
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'tuning'    - circles
%                       'bar'       - horizontal bars
%                   default == 'tuning'
%                   - 'PlotMode': plotting mode for multiple columns
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'overlapped'    - overlapped in a single plot
%                       'parallel'      - in parallel in subPlots
%                       'separate'      - in separate figures
%                   default == 'overlapped'
%                   - 'VarsToPlot': variable (column) names of the table
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == plot all variables
%                   - 'VarLabels': variable labels
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == use distinct parts of variable names
%                   - 'DistinctParts': whether to extract distinct parts
%                                       or variable labels
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'PhaseVariables': variable (column) names for phases
%                   must be empty or a character vector or a string vector
%                       or a cell array of character vectors
%                   default == none
%                   - 'Delimiter': delimiter used for extracting distinct parts
%                   must be a string scalar or a character vector
%                   default == '_'
%                   - 'ReadoutLabel': label for the readout
%                   must be a string scalar or a character vector
%                   default == set by plot_tuning_curve.m
%                   - 'TableLabel': label for the table
%                   must be a string scalar or a character vector
%                   default == either a common prefix from variable names
%                               or the input table variable name
%                   - 'PLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == none ('suppress')
%                   - 'PTicks': x tick values for the parameter values
%                   must be a numeric vector
%                   default == 1:numel(pTickLabels)
%                   - 'PTickLabels': x tick labels
%                   must be a cell array of character vectors/strings
%                   default == row names or times if provided
%                   - 'OutFolder': output folder if FigNames not set
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigName': figure name for saving
%                   must be a string scalar or a character vector
%                   default == ''
%                   - Any other parameter-value pair for the plot_struct() 
%                       or the plot_tuning_curve() function
%
% Requires:
%       cd/char2rgb.m
%       cd/create_error_for_nargin.m
%       cd/extract_common_directory.m
%       cd/extract_common_prefix.m
%       cd/extract_fileparts.m
%       cd/plot_struct.m
%       cd/plot_tuning_curve.m
%
% Used by:
%       cd/plot_measures.m
%       cd/parse_multiunit.m
%       cd/plot_repetitive_protocols.m

% File History:
% 2018-12-18 Moved from plot_repetitive_protocols.m
% 2018-12-18 Now uses iP.KeepUnmatched
% 2018-12-18 Now uses extract_common_directory.m
% 2018-12-18 Now uses row names without processing if not file names
% 2019-03-17 Deal with timetables differently for RowNames
% 2019-03-17 Implemented plotting together (use plot_tuning_curve directly)
% 2019-03-17 Added 'PlotSeparately' as an optional argument
% 2019-03-25 Added 'PhaseVariables' as an optional argument
% 2019-05-08 Added 'PlotType' as an optional argument
% 2019-08-07 Added 'PTickLabels' as an optional argument
% 2019-08-07 Added 'PTicks' as an optional argument
% 2019-12-30 Changed 'PlotSeparately' to 'PlotMode'
% TODO: Merge with plot_table_parallel.m
% TODO: Return handles to plots
% TODO: Pass in figNames or figNumbers when plotting separately
% 

%% Hard-coded parameters
validPlotTypes = {'tuning', 'bar'};
validPlotModes = {'overlapped', 'parallel', 'separate'};

lineSpecOverlapped = '-';
lineWidthOverlapped = 1;
markerEdgeColorOverlapped = char2rgb('DarkOrchid');
markerFaceColorOverlapped = char2rgb('LightSkyBlue');

lineSpecParallel = [];
lineWidthParallel = [];
markerEdgeColorParallel = [];
markerFaceColorParallel = [];

lineSpecSeparate = 'o';
lineWidthSeparate = 1;
markerEdgeColorSeparate = char2rgb('DarkOrchid');
markerFaceColorSeparate = char2rgb('LightSkyBlue');

%% Default values for optional arguments
plotTypeDefault = 'tuning';
plotModeDefault = 'overlapped'; % plot columns overlapped by default
lineSpecDefault = '';
lineWidthDefault = [];
markerEdgeColorDefault = [];
markerFaceColorDefault = [];
varsToPlotDefault = {};         % plot all variables by default
varLabelsDefault = {};          % set later
distinctPartsDefault = true;    % extract distinct parts of variable names
                                %   by default
phaseVariablesDefault = {};     % no phases by default
delimiterDefault = '_';         % use '_' as delimiter by default
readoutLabelDefault = '';       % set later
tableLabelDefault = '';         % set later
pLabelDefault = 'suppress';     % No x label by default
pTicksDefault = [];
pTickLabelsDefault = {};
outFolderDefault = pwd;
figNameDefault = '';

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
addRequired(iP, 'myTable', ...
    @(x) validateattributes(x, {'table', 'timetable'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));
addParameter(iP, 'PlotMode', plotModeDefault, ...
    @(x) any(validatestring(x, validPlotModes)));
addParameter(iP, 'LineSpec', lineSpecDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'LineWidth', lineWidthDefault);
addParameter(iP, 'MarkerEdgeColor', markerEdgeColorDefault);
addParameter(iP, 'MarkerFaceColor', markerFaceColorDefault);
addParameter(iP, 'VarsToPlot', varsToPlotDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VarsToPlot must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'VarLabels', varLabelsDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VarLabels must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'DistinctParts', distinctPartsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PhaseVariables', phaseVariablesDefault, ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['VarsToPlot must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ReadoutLabel', readoutLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'TableLabel', tableLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PLabel', pLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'PTicks', pTicksDefault, ...
    @(x) isempty(x) || isnumericvector(x));
addParameter(iP, 'PTickLabels', pTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigName', figNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, myTable, varargin{:});
plotType = validatestring(iP.Results.PlotType, validPlotTypes);
plotMode = validatestring(iP.Results.PlotMode, validPlotModes);
lineSpec = iP.Results.LineSpec;
lineWidth = iP.Results.LineWidth;
markerEdgeColor = iP.Results.MarkerEdgeColor;
markerFaceColor = iP.Results.MarkerFaceColor;
varsToPlot = iP.Results.VarsToPlot;
varLabels = iP.Results.VarLabels;
distinctParts = iP.Results.DistinctParts;
phaseVariables = iP.Results.PhaseVariables;
delimiter = iP.Results.Delimiter;
readoutLabel = iP.Results.ReadoutLabel;
tableLabel = iP.Results.TableLabel;
pLabel = iP.Results.PLabel;
pTicks = iP.Results.PTicks;
pTickLabels = iP.Results.PTickLabels;
outFolder = iP.Results.OutFolder;
figName = iP.Results.FigName;

% Keep unmatched arguments for the line() function
otherArguments = iP.Unmatched;

%% Preparation
% Check if output directory exists
check_dir(outFolder);

% Restrict to variables to plot or extract the variable names
if ~isempty(varsToPlot)
    tableToPlot = myTable(:, varsToPlot);
else
    tableToPlot = myTable;
    varsToPlot = myTable.Properties.VarsToPlot;
end

% If provided, make sure there are an equal number of phase variables
%   and extract the phase vectors for each variable
if ~isempty(phaseVariables)
    % Force as a column cell array
    phaseVariables = force_column_cell(phaseVariables);

    % Match the number of phase variables to the number of variables to plot
    phaseVariables = match_row_count(phaseVariables, numel(varsToPlot));

    % Extract the phase vectors
    phaseVectors = cellfun(@(x) myTable{:, x}, phaseVariables, ...
                            'UniformOutput', false);
else
    phaseVectors = {};
end

% Extract distinct parts if requested
if distinctParts
    varLabels = extract_fileparts(varsToPlot, 'distinct', 'Delimiter', delimiter);
else
    varLabels = varsToPlot;
end

% Decide on table label
if isempty(tableLabel)
    % First try to extract a common prefix from the variables to plot
    tableLabel = extract_common_prefix(varsToPlot, 'Delimiter', delimiter);

    % If no such prefix exists, use the table variable name
    if isempty(tableLabel)
        tableLabel = inputname(1);
    end
end

% Decide on pTickLabels
if isempty(pTickLabels)
    if isfield(myTable.Properties, 'RowNames') && ...
            iscell(myTable.Properties.RowNames)
        % Get the row names
        rowNames = myTable.Properties.RowNames;

        % If all row names are file names, process them
        %   Otherwise, just use the row names as the x tick labels
        if all(isfile(rowNames))
            % Extract the distinct file bases
            pTickLabels = extract_fileparts(rowNames, 'distinct');

            % Replace all instances of '_' with '\_'
            pTickLabels = replace(pTickLabels, '_', '\_');
        else
            % Just use the row names
            pTickLabels = rowNames;
        end
    else
        % Use default x tick labels
        pTickLabels = {};
    end
end

% Decide on lineSpec
if isempty(lineSpec)
    switch plotMode
    end
end

% Decide on lineWidth
if isempty(lineWidth)
    switch plotMode
    case 'overlapped'
        lineWidth = lineWidthOverlapped;
    case 'parallel'
        lineWidth = lineWidthParallel;
    case 'separate'
        lineWidth = lineWidthSeparate;
    end
end

% Decide on markerEdgeColor
if isempty(markerEdgeColor)
    switch plotMode
    case 'overlapped'
        markerEdgeColor = markerEdgeColorOverlapped;
    case 'parallel'
        markerEdgeColor = markerEdgeColorParallel;
    case 'separate'
        markerEdgeColor = markerEdgeColorSeparate;
    end
end

% Decide on markerFaceColor
if isempty(markerFaceColor)
    switch plotMode
    case 'overlapped'
        markerFaceColor = markerFaceColorOverlapped;
    case 'parallel'
        markerFaceColor = markerFaceColorParallel;
    case 'separate'
        markerFaceColor = markerFaceColorSeparate;
    end
end

%% Do the job
switch plotMode
case 'overlapped'
    % Create a figure name if empty
    if isempty(figName)
        figName = fullfile(outFolder, [tableLabel, '.png']);
    end

    % Convert to an array
    if istimetable(tableToPlot)
        % Extract variables
        myArray = tableToPlot.Variables;
    else
        % Use table2array
        myArray = table2array(tableToPlot);
    end

    % Decide on x values
    if istimetable(tableToPlot)
        % Extract time
        xValues = tableToPlot.Properties.RowTimes;
    else
        % Count rows
        nRows = height(tableToPlot);

        % Use row numbers
        xValues = transpose(1:nRows);
    end

    % Decide on readout label
    if isempty(readoutLabel)
        readoutLabel = replace(tableLabel, '_', ' ');
    end

    % Decide on figure title
    if istimetable(tableToPlot)
        figTitle = replace([tableLabel, ' over time'], '_', ' ');
    else
        figTitle = '';
    end

    % Clear the current figure
    clf;

    % Plot a tuning curve
    switch plotType
    case 'tuning'
        handles = plot_tuning_curve(xValues, myArray, 'FigName', figName, ...
                        'PhaseVectors', phaseVectors, ...
                        'PTicks', pTicks, 'PTickLabels', pTickLabels, ...
                        'PLabel', pLabel, ...
                        'ReadoutLabel', readoutLabel, ...
                        'ColumnLabels', varLabels, ...
                        'FigTitle', figTitle, ...
                        'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', markerFaceColor, ...
                        otherArguments);
    case 'bar'
        % TODO
    otherwise
        error('plotType unrecognized!')
    end
case 'parallel'
    % TODO: Use plot_table_parallel.m
case 'separate'
    % Convert to a structure array
    myStruct = table2struct(tableToPlot);

    % Plot fields
    fig = plot_struct(myStruct, 'OutFolder', outFolder, ...
                        'PhaseVectors', phaseVectors, ...
                        'PlotType', plotType, ...
                        'FieldLabels', varLabels, ...
                        'PTicks', pTicks, 'PTickLabels', pTickLabels, ...
                        'PLabel', pLabel, ...
                        'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', markerFaceColor, ...
                         otherArguments);
    handles.fig = fig;
otherwise
    error('plotMode unrecognized!');
end

%% Outputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get the file bases
[~, fileBases, ~] = ...
    cellfun(@(x) fileparts(x), newPaths, 'UniformOutput', false);

% Create x tick labels
pTickLabels = cellfun(@(x) strrep(x, '_', '\_'), fileBases, ...
                        'UniformOutput', false);

% Will not work for durations data
if isempty(pTicks) && ~isempty(pTickLabels)
    pTicks = (1:numel(pTickLabels))';
end

% Does not work if pTicks is not also set
elseif isfield(myTable.Properties, 'RowTimes')
    % Convert time to minutes
    timeVec = minutes(myTable.Properties.RowTimes);

    % Convert to a cell array of character vectors
    pTickLabels = convert_to_char(timeVec);

%                   - 'PlotSeparately': whether to plot each column separately
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
plotSeparatelyDefault = false;  % plot variables together by default
addParameter(iP, 'PlotSeparately', plotSeparatelyDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
plotSeparately = iP.Results.PlotSeparately;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
