function figs = plot_struct (structArray, varargin)
%% Plot all fields in a structure array as tuning curves
% Usage: figs = plot_struct (structArray, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       figs        - figure handle(s) for the created figure(s)
%                   specified as a figure object handle column vector
%
% Arguments:    
%       structArray - a structure array containing scalar fields
%                   must be a 2-D structure array
%       varargin    - 'PlotType': type of plot
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'tuning'    - circles
%                       'bar'       - horizontal bars
%                   default == 'tuning'
%                   - 'LineSpec': line specification
%                   must be a character array
%                   default == '-'
%                   - 'PIsLog': whether parameter values are to be plotted 
%                               log-scaled
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false;
%                   - 'PTicks': x tick values for the parameter values
%                   must be a numeric vector
%                   default == []
%                   - 'PTickLabels': x tick labels in place of parameter values
%                   must be a cell array of character vectors/strings
%                   default == {}
%                   - 'PTickAngle': angle for parameter tick labels
%                   must be a numeric scalar
%                   default == 0
%                   - 'PLabel': label for the parameter
%                   must be a string scalar or a character vector
%                   default == 'Parameter'
%                   - 'FieldLabels': label for the field
%                   must be a cell array of character vectors/strings
%                   default == field name
%                   - 'ColorMap' - color map used when nColumnsToPlot > 1
%                   must be a 2-D numeric array with 3 columns
%                   default == set in plot_tuning_curve.m
%                   - 'FigTitles': titles for each figure
%                   must be a cell array of character vectors/strings
%                   default == [fieldLabel, ' vs. ', pLabel]
%                   - 'FigNumber': figure number for creating figure(s)
%                   must be a positive integer vector
%                   default == []
%                   - 'OutFolder': output folder if FigNames not set
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FigNames': figure names for saving
%                   must be a cell array of character vectors/strings
%                   default == {}
%                   - 'FigTypes': figure type(s) for saving; 
%                               e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by 
%                       the built-in saveas() function
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   - Any other parameter-value pair for the 
%                       plot_tuning_curve() or plot_bar() function
%
% Requires:
%       ~/Downloaded_Functions/rgb.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/force_column_cell.m
%       cd/isfigtype.m
%       cd/ispositiveintegervector.m
%       cd/match_row_count.m
%       cd/plot_bar.m
%       cd/plot_tuning_curve.m
%       cd/save_all_figtypes.m
%
% Used by:    
%       cd/plot_table.m

% File History:
% 2018-09-26 Created by Adam Lu
% 2018-12-15 Updated PTicks so that it is dependent on nEntries
% 2018-12-17 Now uses create_labels_from_numbers.m
% 2018-12-18 Now uses iP.KeepUnmatched
% 2018-12-18 Changed lineSpec default to o and singleColorDefault to SkyBlue
% 2019-03-14 Now saves the plots here
% 2019-05-08 Added 'PlotType' as an optional argument
% 2019-05-11 Added 'RBoundaries' as an optional argument
% 2019-06-11 Moved boundary plotting code to plot_bar.m and plot_tuning_curve.m
% 2019-11-24 Moved phase average computing code to plot_bar.m
% TODO: Return handles to plots
% TODO: Pass in figNames or figNumbers when plotting separately
% 

%% Hard-coded parameters
validPlotTypes = {'tuning', 'bar'};
maxNPTicks = 10;
barDirection = 'horizontal';
barReverseOrder = true;

%% Default values for optional arguments
plotTypeDefault = 'tuning';
lineSpecDefault = 'o';
lineWidthDefault = [];
markerEdgeColorDefault = [];
markerFaceColorDefault = [];
pIsLogDefault = false;
pTicksDefault = [];
pTickLabelsDefault = {};
pTickAngleDefault = 0;
pLabelDefault = 'Parameter';
fieldLabelsDefault = {};
colorMapDefault = [];           % set later
figTitlesDefault = {};          % set later
figNumberDefault = [];          % use current figure by default
outFolderDefault = pwd;
figNamesDefault = {};
figTypesDefault = 'png';

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
addRequired(iP, 'structArray', ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotType', plotTypeDefault, ...
    @(x) any(validatestring(x, validPlotTypes)));
addParameter(iP, 'LineSpec', lineSpecDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'LineWidth', lineWidthDefault);
addParameter(iP, 'MarkerEdgeColor', markerEdgeColorDefault);
addParameter(iP, 'MarkerFaceColor', markerFaceColorDefault);
addParameter(iP, 'PIsLog', pIsLogDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PTicks', pTicksDefault, ...
    @(x) isempty(x) || isnumericvector(x));
addParameter(iP, 'PTickLabels', pTickLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'PTickAngle', pTickAngleDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PLabel', pLabelDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FieldLabels', fieldLabelsDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'ColorMap', colorMapDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d', 'ncols', 3}));
addParameter(iP, 'FigTitles', figTitlesDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FigNumber', figNumberDefault, ...
    @(x) isempty(x) || ispositiveintegervector(x));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'FigNames', figNamesDefault, ...
    @(x) isempty(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'FigTypes', figTypesDefault, ...
    @(x) all(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, structArray, varargin{:});
plotType = validatestring(iP.Results.PlotType, validPlotTypes);
lineSpec = iP.Results.LineSpec;
lineWidth = iP.Results.LineWidth;
markerEdgeColor = iP.Results.MarkerEdgeColor;
markerFaceColor = iP.Results.MarkerFaceColor;
pIsLog = iP.Results.PIsLog;
pTicks = iP.Results.PTicks;
pTickLabels = iP.Results.PTickLabels;
pTickAngle = iP.Results.PTickAngle;
pLabel = iP.Results.PLabel;
fieldLabels = iP.Results.FieldLabels;
colorMap = iP.Results.ColorMap;
figTitles = iP.Results.FigTitles;
figNumber = iP.Results.FigNumber;
outFolder = iP.Results.OutFolder;
figNames = iP.Results.FigNames;
[~, figTypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);

% Keep unmatched arguments for the plot_tuning_curve() or plot_bar() function
otherArguments = iP.Unmatched;

% Check relationships between arguments
if ~isempty(pTicks) && ~isempty(pTickLabels) && ...
    numel(pTicks) ~= numel(pTickLabels)
    fprintf(['PTicks and PTickLabels must have ', ...
                'the same number of elements!\n']);
    figs = gobjects(0);
    return
end

%% Preparation
% Count the number of entries
nEntries = length(structArray);

% Return if there are no entries
if nEntries == 0
    figs = gobjects(0);
    return;
end

% Create a vector for the parameter values
pValues = transpose(1:nEntries);

% Decide on the number of parameter values to actually show
if isempty(pTicks)
    % Decide on the number of parameter values to show
    nPTicks = min(maxNPTicks, nEntries);

    % Evenly space them out starting with the first parameter
    pTicks = transpose(1:nPTicks) .* floor(nEntries/nPTicks);
else
    nPTicks = length(pTicks);
end

% Generate corresponding parameter value labels
if ~isempty(pTickLabels) 
    % Force as a column cell array
    pTickLabels = force_column_cell(pTickLabels);
    
    % Match the row counts
    pTickLabels = match_row_count(pTickLabels, nPTicks);
elseif isempty(pTickLabels)
    % Generate pTickLabels from pTicks
    pTickLabels = create_labels_from_numbers(pTicks);
end

% Get all the fields of the structArray as a cell array
allFields = fieldnames(structArray);

% Create figure names if not provided
if isempty(figNames)
    figNames = cellfun(@(x) fullfile(outFolder, [x, '_vs_', pLabel]), ...
                        allFields, 'UniformOutput', false);
end

% Count the number of fields
nFieldsOrig = numel(allFields);

% Take only the fields of the structArray that are numeric scalars
scalarStructArray = structArray;
for iField = 1:nFieldsOrig
    % Get the field name
    thisFieldName = allFields{iField};

    % Get the first instance of this field value
    thisFieldValue = structArray(1).(thisFieldName);

    % Remove this field if not a numeric scalar
    if ~isscalar(thisFieldValue) || ~isnumeric(thisFieldValue)
        fprintf(['Warning: the field %s is not a ', ...
                    'numeric scalar so will be removed!!\n'], ...
                    thisFieldName);
        scalarStructArray = rmfield(scalarStructArray, thisFieldName);
    end
end

% Get all the fields of the scalarStructArray as a cell array
allScalarFields = fieldnames(scalarStructArray);

% Count the number of fields
nFields = numel(allScalarFields);

% Return if there are no more fields
if nFields == 0
    figs = gobjects(0);
    return;
end

% Convert the data to a homogeneous array, with each column being a field
fieldData = table2array(struct2table(scalarStructArray));

%% Plot all fields
figs = gobjects(nFields, 1);
for iField = 1:nFields
    % Get the field value vector for this field
    fieldVals = fieldData(:, iField);

    % Set the field label for this field
    if ~isempty(fieldLabels)
        % Use the user-provided field label
        fieldLabel = fieldLabels{iField};
    else
        % Use the field name
        fieldLabel = allScalarFields{iField};
    end

    % Set the figure title
    if ~isempty(figTitles)
        % Use the user-provided figure title
        figTitle = figTitles{iField};
    else
        % Use the default
        if ~strcmpi(fieldLabel, 'suppress') && ~strcmpi(pLabel, 'suppress')
            figTitle = strrep([fieldLabel, ' vs. ', pLabel], '_', '\_');
        elseif ~strcmpi(fieldLabel, 'suppress')
            figTitle = strrep([fieldLabel, ' vs. parameter'], '_', '\_');
        else
            figTitle = 'Readout vs. parameter';
        end
    end
    
    % Set the figure name
    if ~isempty(figNames)
        figName = figNames{iField};
    else
        figName = '';
    end
    
    % Create a new figure
    figThis = set_figure_properties('FigNumber', figNumber);

    % Clear the figure
    clf(figThis);

    switch plotType
    case 'tuning'
        % Plot the tuning curve
        handles = ...
            plot_tuning_curve(pValues, fieldVals, 'PIsLog', pIsLog, ...
                        'PTicks', pTicks, 'PTickLabels', pTickLabels, ...
                        'PTickAngle', pTickAngle, ...
                        'PLabel', pLabel, 'ReadoutLabel', fieldLabel, ...
                        'ColorMap', colorMap, ...
                        'FigTitle', figTitle, 'FigHandle', figThis, ...
                        'LineSpec', lineSpec, 'LineWidth', lineWidth, ...
                        'MarkerEdgeColor', markerEdgeColor, ...
                        'MarkerFaceColor', markerFaceColor, ...
                        otherArguments);
        figThis = handles.fig;
    case 'bar'
        % Plot horizontal bars
        % TODO: Deal with pIsLog
        % TODO: Implement singlecolor
        handles = ...
            plot_bar(fieldVals, 'ForceVectorAsRow', false, ...
                        'ReverseOrder', barReverseOrder, ...
                        'BarDirection', barDirection, ...
                        'PValues', pValues, ...
                        'PTicks', pTicks, 'PTickLabels', pTickLabels, ...
                        'PTickAngle', pTickAngle, ...
                        'PLabel', pLabel, 'ReadoutLabel', fieldLabel, ...
                        'FigTitle', figTitle, 'FigHandle', figThis, ...
                        otherArguments);
        figThis = handles.fig;
    otherwise
        error('plotType unrecognized!')
    end

    if ~isempty(figName)
        save_all_figtypes(figThis, figName, figTypes);
    end

    figs(iField, 1) = figThis;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

