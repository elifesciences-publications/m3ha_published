function m3ha_compare_neuronparams2 (paramValues, paramNames, suffixes, varargin)
%% compare graphs across different sets of NEURON parameters
% Usage: m3ha_compare_neuronparams2 (paramValues, paramNames, suffixes, varargin)
% Arguments:
%       paramValues - a cell array of neuronparams to compare
%                   must be a cell array of numeric arrays
%       paramNames  - a cell array of corresponding neuronparamnames
%                   must be a cell array of cell arrays
%       suffixes    - a cell array of corresponding suffixes
%                   must be a cell array of character arrays/strings
%       varargin    - 'Celsius': temperature of simulation [degC]
%                   must be a numeric scalar
%                   default == 33 degC
%                   - 'Ek': K+ reversal potential [mV]
%                   must be a numeric scalar
%                   default == -100 mV
%                   - 'Ena': Na+ reversal potential [mV]
%                   must be a numeric scalar
%                   default == 88 mV
%                   - 'Eh': reversal potential of Ih [mV]
%                   must be a numeric scalar
%                   default == -43 mV
%                   - 'CaOut': extracellular [Ca++] [mM]
%                   must be a numeric scalar
%                   default == 2 mM
%                   - 'CaIn': intracellular [Ca++] [mM]
%                   must be a numeric scalar
%                   default == 2.4e-4 mM
%                   - 'PcabarIT': maximum Ca++ permeability [cm/s]
%                   must be a numeric scalar
%                   default == 0.2e-3 cm/s
%                   - 'ShiftmIT': depolarizing shift of activation curves [mV]
%                   must be a numeric scalar
%                   default == 1 mV
%                   - 'ShifthIT': depolarizing shift of inactivation curves [mV]
%                   must be a numeric scalar
%                   default == 1 mV
%                   - 'SlopemIT': scaling factor for slope of activation curves
%                   must be a numeric scalar
%                   default == 1
%                   - 'SlopehIT': scaling factor for slope of inactivation curves
%                   must be a numeric scalar
%                   default == 1
%                   - 'GhbarIh': maximum conductance of Ih [S/cm2]
%                   must be a numeric scalar
%                   default == 2.2e-5 S/cm2
%                   - 'ShiftmIh': depolarizing shift of activation curves [mV]
%                   must be a numeric scalar
%                   default == 0 mV
%                   - 'GkbarIKir': maximum conductance of IKir [S/cm2]
%                   must be a numeric scalar
%                   default == 2.0e-5 S/cm2
%                   - 'GkbarIA': maximum conductance of IA [S/cm2]
%                   must be a numeric scalar
%                   default == 5.5e-3 S/cm2
%                   - 'GnabarINaP': maximum conductance of INaP [S/cm2]
%                   must be a numeric scalar
%                   default == 5.5e-6 S/cm2
%                   - 'OutFolder': output directory
%                   must be a string scalar or a character vector
%                   default == pwd
%
% Requires:
%       cd/m3ha_compute_and_plot_geometry2.m
%       cd/m3ha_compute_and_plot_IT2.m
%       cd/m3ha_compute_and_plot_Ih2.m
%       cd/m3ha_compute_and_plot_IKir2.m
%       cd/m3ha_compute_and_plot_IA2.m
%       cd/m3ha_compute_and_plot_INaP2.m
%       cd/m3ha_compute_and_plot_all_IV2.m
%       /home/Matlab/Adams_Functions/find_in_strings.m
%       /home/Matlab/Adams_Functions/intersect_over_cells.m
%
% Used by:
%       cd/m3ha_optimizer_4compgabab.m

% File History:
% 2017-08-09 - Created
% 2017-08-10 - Finished
% 2017-08-10 - Now creates outFolder if doesn't exist
% 2017-08-26 - Now plots geometry if it changes
% 2017-08-29 - Removed 'LDend2', 'diamDend2To1'
% 2017-08-29 - Added plot+geometry
% 2017-09-22 - BT - Calls to compute_and_plot_x once
% 2018-01-24 - Added isdeployed
% 2018-08-15 Updated geometry parameters
% 2018-08-17 Now uses intersect_over_cells.m
% TODO accept just one suffix
% TODO: Plot the different conditions (suffixes) against each other in the same graphs
%       May need to make version 2s for all the compute_and_plot functions that will
%       do this when suffixes is passed as a cell array
%

% Parameters relevant for each set of plots
paramNamesGeo = {'diamSoma', 'LDend1', 'LDend2'};
paramPartsGeo = {'diamSoma', {'L', 'dend1'}, {'L', 'dend2'}};
paramNamesIT = {'caOut', 'caIn', ...
                'pcabarITSoma', 'pcabarITDend1', 'pcabarITDend2', ...
                'shiftmIT', 'shifthIT', ...
                'slopemIT', 'slopehIT'};
paramPartsIT = {{'c', 'out'}, {'c', 'in'}, ...
                {'cabar', 'soma'}, {'cabar', 'dend1'}, {'cabar', 'dend2'}, ...               
                {'shiftm', 'it'}, {'shifth', 'it'}, ...
                {'slopem', 'it'}, {'slopeh', 'it'}};
paramNamesIh = {'ehIh', 'ghbarIhSoma', 'ghbarIhDend1', ...
                'ghbarIhDend2', 'shiftmIh'};
paramPartsIh = {{'eh', 'ih'}, {'ghbar', 'soma'}, {'ghbar', 'dend1'}, ...
                {'ghbar', 'dend2'}, {'shift', 'ih'}};
paramNamesIKir = {'ek', 'gkbarIKirSoma', 'gkbarIKirDend1', 'gkbarIKirDend2'};
paramPartsIKir = {'ek', {'gkbar', 'ikir', 'soma'}, ...
                    {'gkbar', 'ikir', 'dend1'}, {'gkbar', 'ikir', 'dend2'}};
paramNamesIA = {'ek', 'gkbarIASoma', 'gkbarIADend1', 'gkbarIADend2'};
paramPartsIA = {'ek', {'gkbar', 'ia', 'soma'}, ...
                    {'gkbar', 'ia', 'dend1'}, {'gkbar', 'ia', 'dend2'}};
paramNamesINaP = {'ena', 'gnabarINaPSoma', 'gnabarINaPDend1', 'gnabarINaPDend2'};
paramPartsINaP = {'ena', {'gnabar', 'soma'}, ...
                    {'gnabar', 'dend1'}, {'gnabar', 'dend2'}};
paramNamesIV = {'ek', 'ena', 'ehIh', 'caOut', 'caIn', ...
                'pcabarITSoma', 'pcabarITDend1', 'pcabarITDend2', ...
                'shiftmIT', 'shifthIT', ...
                'slopemIT', 'slopehIT', ...
                'ghbarIhSoma', 'ghbarIhDend1', 'ghbarIhDend2', ...
                'shiftmIh', 'gkbarIKirSoma', ...
                'gkbarIKirDend1', 'gkbarIKirDend2', ...
                'gkbarIASoma', 'gkbarIADend1', 'gkbarIADend2', ...
                'gnabarINaPSoma', 'gnabarINaPDend1', 'gnabarINaPDend2'};
paramPartsIV = {'ek', 'ena', {'eh', 'ih'}, {'c', 'out'}, {'c', 'in'}, ...
                {'cabar', 'soma'}, {'cabar', 'dend1'}, {'cabar', 'dend2'}, ... 
                {'shiftm', 'it'}, {'shifth', 'it'}, ... 
                {'slopem', 'it'}, {'slopeh', 'it'}, ...
                {'ghbar', 'soma'}, {'ghbar', 'dend1'}, {'ghbar', 'dend2'}, ...
                {'shiftm', 'Ih'}, {'gkbar', 'ikir', 'soma'}, ...
                {'gkbar', 'ikir', 'dend1'}, {'gkbar', 'ikir', 'dend2'}, ...
                {'gkbar', 'ia', 'soma'},  {'gkbar', 'ia', 'dend1'}, ...
                {'gkbar', 'ia', 'dend2'}, {'gnabar', 'soma'}, ...
                {'gnabar', 'dend1'}, {'gnabar', 'dend2'}};

%% Default voltage vector limits and steps
vMinDefault = -120;         % lowest voltage to plot [mV]
vMaxDefault = 30;           % highest voltage to plot [mV]
vStepDefault = 0.05;        % voltage step for plotting [mV]

%% Default parameters
celsiusDefault = 33;        % default temperature of simulation [degC]
eKDefault = -100;           % default K+ reversal potential [mV]
eNaDefault = 88;            % default Na+ reversal potential [mV]
ehDefault = -43;            % default reversal potential of Ih [mV]
caOutDefault = 2;           % default extracellular [Ca++] [mM]
caInDefault = 2.4e-4;       % default intracellular [Ca++] [mM]
pcabarITDefault = 0.2e-3;   % default maximum Ca++ permeability [cm/s]
shiftmITDefault = 1;  % default depolarizing shift of activation curves [mV]
shifthITDefault = 1;  % default depolarizing shift of inactivation curves [mV]
slopemITDefault = 1;  % default scaling factor for slope of activation curves
slopehITDefault = 1;  % default scaling factor for slope of inactivation curves
ghbarIhDefault = 2.2e-5;    % default maximum conductance of Ih [S/cm2]
shiftmIhDefault = 0;    % default depolarizing shift of activation curves [mV]
gkbarIKirDefault = 2.0e-5;  % default maximum conductance of IKir [S/cm2]
gkbarIADefault = 5.5e-3;    % default maximum conductance of IA [S/cm2]
gnabarINaPDefault = 5.5e-6; % default maximum conductance of INaP [S/cm2]
%outFolderDefault = pwd;     % default output directory
%outFolderDefault = '/media/shareX/share/Brian/compare_neuronparams2/'; % XXX: Debugging
outFolderDefault = '/media/adamX/m3ha/optimizer4gabab/compare_neuronparams2DEBUG/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'm3ha_compare_neuronparams2';

% Add required inputs to the Input Parser
addRequired(iP, 'paramValues', ...      % neuronparams to compare
    @(x) assert(iscell(x) && min(cellfun(@isnumeric, x)), ...
        'paramValues must be a cell array of numeric arrays!'));
addRequired(iP, 'paramNames', ...       % corresponding neuronparamnames
    @(x) assert(iscell(x) && min(cellfun(@iscell, x)), ...
        'paramNames must be a cell array of cell arrays!'));
addRequired(iP, 'suffixes', ...         % corresponding suffixes
    @(x) assert(iscell(x) && (min(cellfun(@ischar, x)) ...
                                || min(cellfun(@isstring, x))), ...
        'suffixes must be a cell array of strings/character arrays!'));

% Add optional inputs to the Input Parser
addOptional(iP, 'v', [], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Celsius', celsiusDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Ek', eKDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Ena', eNaDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'Eh', ehDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'CaOut', caOutDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'CaIn', caInDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'PcabarIT', pcabarITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ShiftmIT', shiftmITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ShifthIT', shifthITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'SlopemIT', slopemITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'SlopehIT', slopehITDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'GhbarIh', ghbarIhDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'ShiftmIh', shiftmIhDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'GkbarIKir', gkbarIKirDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'GkbarIA', gkbarIADefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'GnabarINaP', gnabarINaPDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, paramValues, paramNames, suffixes, varargin{:});
v = iP.Results.v;
celsius = iP.Results.Celsius;
eK = iP.Results.Ek;
eNa = iP.Results.Ena;
eh = iP.Results.Eh;
caOut = iP.Results.CaOut;
caIn = iP.Results.CaIn;
pcabarITSoma = iP.Results.PcabarIT;
pcabarITDend1 = iP.Results.PcabarIT;
pcabarITDend2 = iP.Results.PcabarIT;
shiftmIT = iP.Results.ShiftmIT;
shifthIT = iP.Results.ShifthIT;
slopemIT = iP.Results.SlopemIT;
slopehIT = iP.Results.SlopehIT;
ghbarIhSoma = iP.Results.GhbarIh;
ghbarIhDend1 = iP.Results.GhbarIh;
ghbarIhDend2 = iP.Results.GhbarIh;
shiftmIh = iP.Results.ShiftmIh;
gkbarIKirSoma = iP.Results.GkbarIKir;
gkbarIKirDend1 = iP.Results.GkbarIKir;
gkbarIKirDend2 = iP.Results.GkbarIKir;
gkbarIASoma = iP.Results.GkbarIA;
gkbarIADend1 = iP.Results.GkbarIA;
gkbarIADend2 = iP.Results.GkbarIA;
gnabarINaPSoma = iP.Results.GnabarINaP;
gnabarINaPDend1 = iP.Results.GnabarINaP;
gnabarINaPDend2 = iP.Results.GnabarINaP;
outFolder = iP.Results.OutFolder;

% Construct voltage vector if not provided
if isempty(v)
    v = vMinDefault:vStepDefault:vMaxDefault;
end

%% Add directories to search path for required functions across servers
if exist('/home/Matlab/', 'dir') == 7
    functionsDirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsDirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsDirectory does not exist!');
end
if ~isdeployed
    addpath_custom(fullfile(functionsDirectory, '/Adams_Functions/')); 
                                    % for find_in_strings.m
end

%% Check if needed output directory exist
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory made: %s\n\n', outFolder);
end

% Extract info from arguments
nSets = numel(paramValues);

% Find the parameters that exist in all sets of neuronparams
commonNames = intersect_over_cells(paramNames);

% Make sure there is a parameter to be checked
nCommon = numel(commonNames);
if nCommon < 1
    fprintf('There are no common parameters in the sets provided!\n\n');
    return;
end

% Determine which parameters were changed
nChanged = 0;                       % number of changed parameters
for iCommon = 1:nCommon
    % Get the parameter name of interest
    thisName = commonNames{iCommon};

    % Find the corresponding parameter values in each set
    thisValues = cellfun(@(x, y) y(strcmpi(thisName, x)), ...
                            paramNames, paramValues);

    % Add the parameter and the corresponding values to changedNames and
    %   changedValues if not all the values are the same
    if length(unique(thisValues)) > 1
        nChanged = nChanged + 1;
        changedNames{nChanged} = thisName;
        changedValues{nChanged} = thisValues;
    end
end
if nChanged == 0
    fprintf('No parameters were changed! :P\n\n');
    return;
end

%% Plot graphs according to what parameters were changed
if max(cellfun(@(x) ~isempty(find_in_strings(x, changedNames, ...
                            'SearchMode', 'substrings', ...
                            'IgnoreCase', true)), paramPartsGeo))
                        % if one of the geometry params appears in changedNames
    for iSet = 1:nSets      % for each set of parameters
        % Retrieve updated geometry params, if any
        for iParam = 1:numel(paramNamesGeo)
            % Find the index in paramNames{iSet} for this parameter, if any
            idx = find_in_strings(paramPartsGeo{iParam}, ...
                                        paramNames{iSet}, ...
                                        'SearchMode', 'substrings', ...
                                        'IgnoreCase', true);

            % If this parameter exists exactly once in paramNames{iSet},
            %   update the value using paramValues{iSet}
            if ~isempty(idx) && numel(idx) > 1
                fprintf(['The parameter %s is found ', ...
                        'more than once in Set #%d!!!\n\n'], ...
                        paramNamesGeo{iParam}, iSet);
                return;
            elseif ~isempty(idx)
                eval(sprintf('%s(%d) = %g;\n', ...
                    paramNamesGeo{iParam}, iSet, paramValues{iSet}(idx)));
            end
        end        
     end

    % Plot geometry
    if numel(suffixes) == 2
        % Find all common strings after separation by delimiter '_'
        commonStrs = ...
            intersect_over_cells(cellfun(@(x) strsplit(x, '_'), ...
                                suffixes, 'UniformOutput', false));

        % If there is such a common string, remove it from the second suffix
        if ~isempty(commonStrs) && numel(commonStrs) == 2
            % Note: the first string is alway ''
            suffixesSimplified = ...
                {suffixes{1}, strrep(suffixes{2}, ['_', commonStrs{2}], '')};
        else
            suffixesSimplified = suffixes;
        end

        % Plot invidual plots
        m3ha_compute_and_plot_geometry(outFolder, suffixes{1}, ...
                                   diamSoma(1), LDend1(1), LDend2(1), 'b', 2);
        m3ha_compute_and_plot_geometry(outFolder, suffixes{2}, ...
                                   diamSoma(2), LDend1(2), LDend2(2), 'r', 2); 
    else
        suffixesSimplified = suffixes;
    end
    m3ha_compute_and_plot_geometry2(outFolder, suffixesSimplified, ...
                               diamSoma, LDend1, LDend2, 2);
end

if max(cellfun(@(x) ~isempty(find_in_strings(x, changedNames, ...
                            'SearchMode', 'substrings', ...
                            'IgnoreCase', true)), paramPartsIT))
                        % if one of the IT params appears in changedNames
    for iSet = 1:nSets      % for each set of parameters
        % Retrieve updated IT params, if any
        for iParam = 1:numel(paramNamesIT)  % for each IT parameter
            % Find the index in paramNames{iSet} for this parameter, if any
            idx = find_in_strings(paramPartsIT{iParam}, ...
                                        paramNames{iSet}, ...
                                        'SearchMode', 'substrings', ...
                                        'IgnoreCase', true);

            % If this parameter exists exactly once in paramNames{iSet},
            %   update the value using paramValues{iSet}
            if ~isempty(idx) && numel(idx) > 1
                fprintf(['The parameter %s is found ', ...
                        'more than once in Set #%d!!!\n\n'], ...
                        paramNamesIT{iParam}, iSet);
                return;
            elseif ~isempty(idx)
                eval(sprintf('%s(%d) = %g;\n', ...
                    paramNamesIT{iParam}, iSet, paramValues{iSet}(idx)));
            end
        end
    end
    % Compute and plot IT curves
    m3ha_compute_and_plot_IT2(v, 'OutFolder', outFolder, ...
                        'Suffixes', suffixes, 'Celsius', celsius, ...
                        'Cout', caOut, 'Cin', caIn, ...
                        'Pcabar', pcabarITSoma, ...
                        'Shiftm', shiftmIT, 'Shifth', shifthIT, ...
                        'Slopem', slopemIT, 'Slopeh', slopehIT);
end

if max(cellfun(@(x) ~isempty(find_in_strings(x, changedNames, ...
                            'SearchMode', 'substrings', ...
                            'IgnoreCase', true)), paramPartsIh))
                        % if one of the Ih params appears in changedNames
    for iSet = 1:nSets      % for each set of parameters
        % Retrieve updated Ih params, if any
        for iParam = 1:numel(paramNamesIh)  % for each Ih parameter
            % Find the index in paramNames{iSet} for this parameter, if any
            idx = find_in_strings(paramPartsIh{iParam}, ...
                                        paramNames{iSet}, ...
                                        'SearchMode', 'substrings', ...
                                        'IgnoreCase', true);

            % If this parameter exists exactly once in paramNames{iSet},
            %   update the value using paramValues{iSet}
            if ~isempty(idx) && numel(idx) > 1
                fprintf(['The parameter %s is found ', ...
                        'more than once in Set #%d!!!\n\n'], ...
                        paramNamesIh{iParam}, iSet);
                return;
            elseif ~isempty(idx)
                eval(sprintf('%s(%d) = %g;\n', ...
                    paramNamesIh{iParam}, iSet, paramValues{iSet}(idx)));
            end
        end
    end
    % Compute and plot Ih curves
    m3ha_compute_and_plot_Ih2(v, 'OutFolder', outFolder, ...
                        'Suffixes', suffixes, 'Celsius', celsius, ...
                        'Ghbar', ghbarIhSoma, 'Erev', eh, ...
                        'Shiftm', shiftmIh);
end

if max(cellfun(@(x) ~isempty(find_in_strings(x, changedNames, ...
                            'SearchMode', 'substrings', ...
                            'IgnoreCase', true)), paramPartsIKir))
                        % if one of the IKir params appears in changedNames
    for iSet = 1:nSets      % for each set of parameters
        % Retrieve updated IKir params, if any
        for iParam = 1:numel(paramNamesIKir)  % for each IKir parameter
            % Find the index in paramNames{iSet} for this parameter, if any
            idx = find_in_strings(paramPartsIKir{iParam}, ...
                                        paramNames{iSet}, ...
                                        'SearchMode', 'substrings', ...
                                        'IgnoreCase', true);

            % If this parameter exists exactly once in paramNames{iSet},
            %   update the value using paramValues{iSet}
            if ~isempty(idx) && numel(idx) > 1
                fprintf(['The parameter %s is found ', ...
                        'more than once in Set #%d!!!\n\n'], ...
                        paramNamesIKir{iParam}, iSet);
                return;
            elseif ~isempty(idx)
                eval(sprintf('%s(%d) = %g;\n', ...
                    paramNamesIKir{iParam}, iSet, paramValues{iSet}(idx)));
            end
        end
    end
    % Compute and plot IKir curves
    m3ha_compute_and_plot_IKir2(v, 'OutFolder', outFolder, ...
                            'Suffixes', suffixes, ...
                            'Celsius', celsius, ...
                            'Gkbar', gkbarIKirSoma, 'Ek', eK);
end

if max(cellfun(@(x) ~isempty(find_in_strings(x, changedNames, ...
                            'SearchMode', 'substrings', ...
                            'IgnoreCase', true)), paramPartsINaP))
                        % if one of the IA params appears in changedNames
    for iSet = 1:nSets      % for each set of parameters
        % Retrieve updated IA params, if any
        for iParam = 1:numel(paramNamesIA)  % for each IA parameter
            % Find the index in paramNames{iSet} for this parameter, if any
            idx = find_in_strings(paramPartsIA{iParam}, ...
                                        paramNames{iSet}, ...
                                        'SearchMode', 'substrings', ...
                                        'IgnoreCase', true);

            % If this parameter exists exactly once in paramNames{iSet},
            %   update the value using paramValues{iSet}
            if ~isempty(idx) && numel(idx) > 1
                fprintf(['The parameter %s is found ', ...
                        'more than once in Set #%d!!!\n\n'], ...
                        paramNamesIA{iParam}, iSet);
                return;
            elseif ~isempty(idx)
                eval(sprintf('%s(%d) = %g;\n', ...
                    paramNamesIA{iParam}, iSet, paramValues{iSet}(idx)));
            end
        end
    end
    % Compute and plot IA curves
    m3ha_compute_and_plot_IA2(v, 'OutFolder', outFolder, ...
                        'Suffixes', suffixes, 'Celsius', celsius, ...
                        'Gkbar', gkbarIASoma, 'Ek', eK);
end

if max(cellfun(@(x) ~isempty(find_in_strings(x, changedNames, ...
                            'SearchMode', 'substrings', ...
                            'IgnoreCase', true)), paramPartsINaP))
                        % if one of the INaP params appears in changedNames
    for iSet = 1:nSets      % for each set of parameters
        % Retrieve updated INaP params, if any
        for iParam = 1:numel(paramNamesINaP)  % for each INaP parameter
            % Find the index in paramNames{iSet} for this parameter, if any
            idx = find_in_strings(paramPartsINaP{iParam}, ...
                                        paramNames{iSet}, ...
                                        'SearchMode', 'substrings', ...
                                        'IgnoreCase', true);

            % If this parameter exists exactly once in paramNames{iSet},
            %   update the value using paramValues{iSet}
            if ~isempty(idx) && numel(idx) > 1
                fprintf(['The parameter %s is found ', ...
                        'more than once in Set #%d!!!\n\n'], ...
                        paramNamesINaP{iParam}, iSet);
                return;
            elseif ~isempty(idx)
                eval(sprintf('%s(%d) = %g;\n', ...
                    paramNamesINaP{iParam}, iSet, paramValues{iSet}(idx)));
            end
        end
    end
    % Compute and plot INaP curves
    m3ha_compute_and_plot_INaP2(v, 'OutFolder', outFolder, ...
                        'Suffixes', suffixes, 'Celsius', celsius, ...
                        'GNabar', gnabarINaPSoma, 'ENa', eNa);
end

if max(cellfun(@(x) ~isempty(find_in_strings(x, changedNames, ...
                            'SearchMode', 'substrings', ...
                            'IgnoreCase', true)), paramPartsIV))
                        % if one of the IV params appears in changedNames
    for iSet = 1:nSets      % for each set of parameters
        % Retrieve updated IV params, if any
        for iParam = 1:numel(paramNamesIV)  % for each IV parameter
            % Find the index in paramNames{iSet} for this parameter, if any
            idx = find_in_strings(paramPartsIV{iParam}, ...
                                        paramNames{iSet}, ...
                                        'SearchMode', 'substrings', ...
                                        'IgnoreCase', true);

            % If this parameter exists exactly once in paramNames{iSet},
            %   update the value using paramValues{iSet}
            if ~isempty(idx) && numel(idx) > 1
                fprintf(['The parameter %s is found ', ...
                        'more than once in Set #%d!!!\n\n'], ...
                        paramNamesIV{iParam}, iSet);
                return;
            elseif ~isempty(idx)
                eval(sprintf('%s(%d) = %g;\n', ...
                    paramNamesIV{iParam}, iSet, paramValues{iSet}(idx)));
            end
        end
    end
    % Compute and plot I-V curves
    m3ha_compute_and_plot_all_IV2(v, 'OutFolder', outFolder, ...
                            'Suffixes', cellfun(@(x) [x '_soma'], suffixes, ...
                                                'UniformOutput', false), ...
                            'Celsius', celsius, 'Ek', eK, 'ENa', eNa, ...
                            'Eh', eh, 'CaOut', caOut, 'CaIn', caIn, ...
                            'PcabarIT', pcabarITSoma, ...
                            'ShiftmIT', shiftmIT, 'ShifthIT', shifthIT, ...
                            'SlopemIT', slopemIT, 'SlopehIT', slopehIT, ...
                            'GhbarIh', ghbarIhSoma, ...
                            'ShiftmIh', shiftmIh, ...
                            'GkbarIKir', gkbarIKirSoma, ...
                            'GkbarIA', gkbarIASoma, ...
                            'Gnabar', gnabarINaPSoma);
    m3ha_compute_and_plot_all_IV2(v, 'OutFolder', outFolder, ...
                            'Suffixes', cellfun(@(x) [x '_dend1'], suffixes, ...
                                                'UniformOutput', false), ...
                            'Celsius', celsius, 'Ek', eK, 'ENa', eNa, ...
                            'Eh', eh, 'CaOut', caOut, 'CaIn', caIn, ...
                            'PcabarIT', pcabarITDend1, ...
                            'ShiftmIT', shiftmIT, 'ShifthIT', shifthIT, ...
                            'SlopemIT', slopemIT, 'SlopehIT', slopehIT, ...
                            'GhbarIh', ghbarIhDend1, ...
                            'ShiftmIh', shiftmIh, ...
                            'GkbarIKir', gkbarIKirDend1, ...
                            'GkbarIA', gkbarIADend1, ...
                            'Gnabar', gnabarINaPDend1);
    m3ha_compute_and_plot_all_IV2(v, 'OutFolder', outFolder, ...
                            'Suffixes', cellfun(@(x) [x '_dend2'], suffixes, ...
                                                'UniformOutput', false), ...
                            'Celsius', celsius, 'Ek', eK, 'ENa', eNa, ...
                            'Eh', eh, 'CaOut', caOut, 'CaIn', caIn, ...
                            'PcabarIT', pcabarITDend2, ...
                            'ShiftmIT', shiftmIT, 'ShifthIT', shifthIT, ...
                            'SlopemIT', slopemIT, 'SlopehIT', slopehIT, ...
                            'GhbarIh', ghbarIhDend2, ...
                            'ShiftmIh', shiftmIh, ...
                            'GkbarIKir', gkbarIKirDend2, ...
                            'GkbarIA', gkbarIADend2, ...
                            'Gnabar', gnabarINaPDend2);
end

fprintf('All plots are complete!!! :)\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
                eval(sprintf('%s{%d} = %g;\n', ...
                    paramNamesGeo{iParam}, iSet, paramValues{iSet}(idx)));

paramNamesGeo = {'diamSoma', 'LDend1', ...
                'LDend2', 'distDendPercent'};
paramPartsGeo = {'diamSoma', {'L', 'dend'}, ...
                {'diam', 'dend'}, {'dend', 'perc'}};

m3ha_compute_and_plot_geometry2(outFolder, suffixes, ...
                diamSoma, LDend, diamDendToSoma, distDendPercent, ...
                2);

commonNames = paramNames{1};        % initialize to the first set of paramnames
for iSet = 2:nSets
    % Intersect this set of paramnames with the previous sets of paramnames
    commonNames = intersect(paramNames{iSet}, commonNames);
end

% if nSets < 2
%     fprintf('There is nothing to compare against!\n\n');
%     return;
% end

%}

