function [simParamsTable, simParamsPath] = ...
                m3ha_neuron_create_sim_params (neuronParamsTable, varargin)
%% Generates a table of simulation parameters from table(s) of neuron parameters
% Usage: [simParamsTable, simParamsPath] = ...
%               m3ha_neuron_create_sim_params (neuronParamsTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       simParamsTable  - table of simulation parameters, 
%                           with each row being a simulation
%                       specified as a 2d table
%       simParamsPath   - file path for saved simulation parameters
%                       specified as a character vector
%
% Arguments:
%       neuronParamsTable   
%                   - table(s) of single neuron parameters with 
%                       parameter names as 'RowNames' and with variables:
%                       'Value': value of the parameter
%                       'LowerBound': lower bound of the parameter
%                       'UpperBound': upper bound of the parameter
%                       'JitterPercentage': jitter percentage of the parameter
%                       'IsLog': whether the parameter is 
%                                   to be varied on a log scale
%                   must be a 2d table or a cell array of 2d tables
%       varargin    - 'Prefix': prefix to prepend to file names
%                   must be a character array
%                   default == ''
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SaveParamsFlag': whether to save simulation parameters
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'JitterFlag': whether to introduce noise in parameters
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'CprWindow': current pulse response window in ms
%                   must be a numeric vector with 2 elements
%                   default == [2000, 2360]
%                   - 'IpscrWindow': IPSC response window in ms
%                   must be a numeric vector with 2 elements
%                   default == [2000, 10000]
%                   - 'NSims': number of simulations
%                   must be a positive integer scalar
%                   default == 1
%                   - 'UseHH': whether to use HH channels
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'UseCvode': whether to use CVODE
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false if useHH is true; true otherwise
%                   - 'SecondOrder': whether to use the Crank-Nicholson method
%                               (otherwise the Backward-Euler method is default)
%                       see https://neuron.yale.edu/neuron/static/new_doc/
%                               simctrl/programmatic.html?#secondorder
%                   must be one of the following:
%                           0 - Implicit Backward Euler
%                           1 - Crank-Nicholson 
%                           2 - Crank-Nicholson with ion current correction
%                   default == 0
%                   - 'TauhMode': mode for simulating tauh
%                   must be one of:
%                       0 - original curve
%                       1 - the same as taum
%                       2 - 10 times smaller amplitude
%                       3 - 10 times larger amplitude
%                       consistent with IT.mod
%                   default == 0
%                   - 'BuildMode': TC neuron build mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - insert leak channels only
%                       'active'  - insert both passive and active channels
%                       or a cell array of them TODO
%                   default == 'active'
%                   - 'SimMode': simulation mode
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'passive' - simulate a current pulse response
%                       'active'  - simulate an IPSC response
%                       or a cell array of them TODO
%                   default == 'active'
%                   - 'OutFilePath': path to output file(s)
%                   must be a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == 'auto'
%                   - 'Tstop': simulation end time
%                   must be a numeric vector
%                   default == based on simMode, cprWindow(2) & ipscrWindow(2)
%                   - 'HoldPotential': holding potential for IPSC response (mV)
%                   must be a numeric vector
%                   default == -70 mV
%                   - 'CurrentPulseAmplitude': current pulse amplitude (nA)
%                   must be a numeric vector
%                   default == -0.050 nA
%                   - 'GababAmp': GABA-B IPSC amplitude (uS)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTrise': GABA-B IPSC rising time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTfallFast': GABA-B IPSC falling phase 
%                                           fast component time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababTfallSlow': GABA-B IPSC falling phase 
%                                           slow component time constant (ms)
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'GababWeight': GABA-B IPSC falling phase 
%                                           fast component weight
%                   must be a numeric vector
%                   default == set of all 12 input waveforms
%                   - 'CustomHoldCurrentFlag': whether to use a custom 
%                                               holding current
%                   must be numeric binary vector
%                   default == 0
%                   - 'HoldCurrent': custom holding current 
%                                       for IPSC response (nA)
%                   must be a numeric vector
%                   default == 0 nA
%                   - 'HoldCurrentNoise': custom holding current noise
%                                           for IPSC response (nA)
%                   must be a numeric vector
%                   default == 0 nA
%
% Requires:
%       cd/argfun.m
%       cd/create_simulation_output_filenames.m
%       cd/force_string_end.m
%       cd/isnumericvector.m
%       cd/ispositiveintegerscalar.m
%       cd/m3ha_load_gabab_ipsc_params.m
%       cd/match_row_count.m
%       cd/transpose_table.m
%
% Used by:    
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2018-10-22 Adapted from code in run_neuron_once_4compgabab.m
% 2018-11-12 Now set the default gababAmp for passive fits to be 0, 
%               but retain the other gabab parameters
% 2019-12-04 Changed the default simMode to 'active'
% 2019-12-19 Added 'buildMode'
% 2019-12-27 Added 'useHH'
% 2020-01-02 Added 'tauhMode'
% 2020-01-05 Fixe nSimsDefault
% 2020-01-22 Now uses m3ha_load_gabab_ipsc_params.m
% 2020-04-27 Added 'UseCvode' as an optional argument
% 2020-04-27 Added 'SecondOrder' as an optional argument
% 2019-04-28 Changed timeToStabilize from 2000 to 3000
% TODO: Remove the Cpr parameters and decide on them before
% 

%% Hard-coded parameters
validBuildModes = {'active', 'passive'};
validSimModes = {'active', 'passive'};
simParamsFileName = 'simulation_parameters.csv';
simParamsFromArguments = {'useHH', 'useCvode', 'secondOrder', ...
                            'tauhMode', 'buildMode', 'simMode', ...
                            'outFilePath', 'tstop' ...
                            'holdPotential', 'currentPulseAmplitude', ...
                            'gababAmp', 'gababTrise', ...
                            'gababTfallFast', 'gababTfallSlow', ...
                            'gababWeight', 'customHoldCurrentFlag', ...
                            'holdCurrent', 'holdCurrentNoise'};

% The following must be consistent with both dclampDataExtractor.m & ...
%   singleneuron4compgabab.hoc
cprWinOrig = [0, 360];          % current pulse response window (ms), original
ipscrWinOrig = [0, 8000];       % IPSC response window (ms), original

% The following must be consistent with singleneuron4compgabab.hoc
timeToStabilize = 3000;         % padded time (ms) to make sure initial value 
                                %   of simulations are stabilized

%% Default values for optional arguments
prefixDefault = '';             % prepend nothing to file names by default
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
saveParamsFlagDefault = true;   % save simulation parameters by default
jitterFlagDefault = false;      % no jitter by default
cprWindowDefault = cprWinOrig + timeToStabilize;
ipscrWindowDefault = ipscrWinOrig + timeToStabilize;
nSimsDefault = [];              % set later

%% Default values for simulation parameters
useHHDefault = false;           % don't use HH channels by default
useCvodeDefault = [];           % set later
secondOrderDefault = 0;         % use Backward Euler by default
tauhModeDefault = 0;            % regular tauh by default
buildModeDefault = 'active';    % insert active channels by default
simModeDefault = 'active';      % simulate a IPSC response by default
outFilePathDefault = 'auto';    % set later
tstopDefault = [];              % set later
holdPotentialDefault = -70;     % (mV)
customHoldCurrentFlagDefault = 0; % don't use custom hold current by default
holdCurrentDefault = 0;         % (nA)
holdCurrentNoiseDefault = 0;    % (nA)

%% Default values for current pulse parameters
currentPulseAmplitudeDefault = -0.050;                          % (nA)

%% Default values for GABAB parameters
gababAmpRatioDefault = [100; 200; 400] / 100;
gababAmpDefault = [];           % set later
gababTriseDefault = [];         % set later
gababTfallFastDefault = [];     % set later
gababTfallSlowDefault = [];     % set later
gababWeightDefault = [];        % set later

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
addRequired(iP, 'neuronParamsTable', ...
    @(x) validateattributes(x, {'table', 'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SaveParamsFlag', saveParamsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'JitterFlag', jitterFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'CprWindow', cprWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'IpscrWindow', ipscrWindowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'NSims', nSimsDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'NSims must be empty or a positive integer scalar!'));
addParameter(iP, 'UseHH', useHHDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'UseCvode', useCvodeDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SecondOrder', secondOrderDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'}));
addParameter(iP, 'TauhMode', tauhModeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'}));
addParameter(iP, 'BuildMode', buildModeDefault, ...
    @(x) any(validatestring(x, validBuildModes)));
addParameter(iP, 'SimMode', simModeDefault, ...
    @(x) any(validatestring(x, validSimModes)));
addParameter(iP, 'OutFilePath', outFilePathDefault, ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'Tstop', tstopDefault, ...
    @(x) assert(isnumericvector(x), 'Tstop must be a numeric vector!'));
addParameter(iP, 'HoldPotential', holdPotentialDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'CurrentPulseAmplitude', currentPulseAmplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'GababAmp', gababAmpDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTrise', gababTriseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTfallFast', gababTfallFastDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababTfallSlow', gababTfallSlowDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'GababWeight', gababWeightDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'CustomHoldCurrentFlag', customHoldCurrentFlagDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'binary', 'vector'}));
addParameter(iP, 'HoldCurrent', holdCurrentDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentNoise', holdCurrentNoiseDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, neuronParamsTable, varargin{:});
prefix = iP.Results.Prefix;
outFolder = iP.Results.OutFolder;
saveParamsFlag = iP.Results.SaveParamsFlag;
jitterFlag = iP.Results.JitterFlag;
cprWindow = iP.Results.CprWindow;
ipscrWindow = iP.Results.IpscrWindow;
nSims = iP.Results.NSims;
useHH = iP.Results.UseHH;
useCvode = iP.Results.UseCvode;
secondOrder = iP.Results.SecondOrder;
tauhMode = iP.Results.TauhMode;
buildMode = validatestring(iP.Results.BuildMode, validBuildModes);
simMode = validatestring(iP.Results.SimMode, validSimModes);
outFilePath = iP.Results.OutFilePath;
tstop = iP.Results.Tstop;
holdPotential = iP.Results.HoldPotential;
currentPulseAmplitude = iP.Results.CurrentPulseAmplitude;
gababAmp = iP.Results.GababAmp;
gababTrise = iP.Results.GababTrise;
gababTfallFast = iP.Results.GababTfallFast;
gababTfallSlow = iP.Results.GababTfallSlow;
gababWeight = iP.Results.GababWeight;
customHoldCurrentFlag = iP.Results.CustomHoldCurrentFlag;
holdCurrent = iP.Results.HoldCurrent;
holdCurrentNoise = iP.Results.HoldCurrentNoise;

%% Preparation
% Make sure prefix ends with a '_'
prefix = force_string_end(prefix, '_', 'OnlyIfNonempty', true);

% Count the number of simulation parameters from arguments
nSimParamsFromArguments = numel(simParamsFromArguments);

% Decide on a cprFlag based on simMode
% TODO: case when simMode is already a cell array
switch simMode
    case 'passive'
        cprFlag = true;
    case 'active'
        cprFlag = false;
    otherwise
        error('simMode unrecognized!');
end

% Make buildMode and simMode cell arrays
buildMode = {buildMode};
simMode = {simMode};

% Decide whether to use CVODE
useCvode = set_default_flag(useCvode, ~useHH);

% Decide on the end time of simulations
if isempty(tstop)
    if cprFlag
        % Simulate until the end of the current pulse response
        tstop = cprWindow(2);
    else
        % Simulate until the end of the IPSC response
        tstop = ipscrWindow(2);    
    end
end

% Decide on GABAB parameters
if isempty(gababAmp)
    % Load default GABAB IPSC parameters in uS
    [gababAmpTemplate, gababTriseTemplate, gababTfallFastTemplate, ...
            gababTfallSlowTemplate, gababWeightTemplate] = ...
        m3ha_load_gabab_ipsc_params('AmpScaleFactor', 100, 'AmpUnits', 'uS');

    % Set the amplitude based on whether cprFlag is on
    if cprFlag
        % Set to zero
        gababAmp = 0;
    else
        % Use all possible combinations from the template
        gababAmp = gababAmpRatioDefault * transpose(gababAmpTemplate);
    end

    % Match gababAmp
    [gababTrise, gababTfallFast, gababTfallSlow, gababWeight] = ...
        argfun(@(x) ones(size(gababAmpRatioDefault)) * transpose(x), ...
                gababTriseTemplate, gababTfallFastTemplate, ...
                gababTfallSlowTemplate, gababWeightTemplate);
end

% Force as a column vector
[gababAmp, gababTrise, gababTfallFast, gababTfallSlow, gababWeight] = ...
    argfun(@(x) reshape(x, [], 1), ...
            gababAmp, gababTrise, gababTfallFast, gababTfallSlow, gababWeight);

% Decide on the number of simulations
%   and modify gabab parameters if necessary
if ~isempty(nSims)
    % Count the number of GABAB IPSC inputs
    nGababInputs = length(gababAmp);

    % Make sure the length of gababAmp
    %   match up with the number of simulations
    if mod(nGababInputs, nSims) == 0
        % If nGababInputs is divisible by nSims,
        %   extract the first of every group
        [gababAmp, gababTrise, gababTfallFast, gababTfallSlow, gababWeight] = ...
            argfun(@(x) x((0:(nSims - 1))' * (nGababInputs/nSims) + 1), ...
                    gababAmp, gababTrise, gababTfallFast, ...
                    gababTfallSlow, gababWeight);
    else
        % Use match_row_count.m
        [gababAmp, gababTrise, gababTfallFast, gababTfallSlow, gababWeight] = ...
            argfun(@(x) match_row_count(x, nSims), ...
                    gababAmp, gababTrise, gababTfallFast, ...
                    gababTfallSlow, gababWeight);
    end
else
    % Use the length of either currentPulseAmplitude or gababAmp
    if cprFlag
        nSims = length(currentPulseAmplitude);
    else
        nSims = length(gababAmp);
    end
end

% Generate output file paths
if strcmpi(outFilePath, 'auto')
    outFilePath = ...
        create_simulation_output_filenames(nSims, 'OutFolder', outFolder, ...
                                            'Prefix', prefix);
end

% Check if any other simulation parameter that's not a scalar has the 
%   same number of elements as the number of simulations
for iParam = 1:nSimParamsFromArguments
    % Get this parameter name
    thisParam = simParamsFromArguments{iParam};

    % Get the number of elements
    nOld = eval(sprintf('numel(%s)', thisParam));

    % Match it up to the number of simulations
    if nOld > 1 && nOld ~= nSims
        error('%s should have %d elements!', thisParam, nSims);
    elseif nOld == 1 && nSims > 1
        eval(sprintf('%s = repmat(%s, %d, 1);', thisParam, thisParam, nSims));
    end
end

%% Generate a table of simulation parameters
% Extract neuron parameters and initialize as simParamsTable
simParamsTable = neuronParams2simParams(neuronParamsTable, nSims, jitterFlag);

% Add other simulation parameters to the table
for iParam = 1:nSimParamsFromArguments
    % Get this parameter name
    thisParam = simParamsFromArguments{iParam};

    % Add it to the table
    simParamsTable.(thisParam) = eval(thisParam);
end

%% Save things
% Write simulation parameters to a file if requested
if saveParamsFlag
    % Construct full path
    simParamsPath = fullfile(outFolder, strcat(prefix, simParamsFileName));

    % Write table to file
    writetable(simParamsTable, simParamsPath);
else
    simParamsPath = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simParamsTable = ...
                neuronParams2simParams (neuronParamsTable, nSims, jitterFlag)
%% Extracts NEURON parameters and initialize as simParamsTable

% Check if the number of neuron parameter tables is the 
%   same as the number of simulations
if iscell(neuronParamsTable) && numel(neuronParamsTable) > 1
    if numel(neuronParamsTable) ~= nSims
        error('neuronParamsTable should have %d elements!', nSims);
    end
end

% Extract NEURON parameters and initialize as simParamsTable
if iscell(neuronParamsTable) && numel(neuronParamsTable) > 1
    % Extract just the value column to get simpler tables
    neuronParamValuesTable = ...
        cellfun(@(x) x(:, 'Value'), neuronParamsTable, 'UniformOutput', false);

    % Transpose the tables
    neuronParamValuesTableTransposed = ...
        cellfun(@(x) transpose_table(x, 'RowNames', 'suppress'), ...
                neuronParamValuesTable, 'UniformOutput', false);

    % Vertically concatenate the tables
    simParamsTable = vertcat(neuronParamValuesTableTransposed{:});
else
    % Get the table if it is in a cell array
    if iscell(neuronParamsTable)
        neuronParamsTable = neuronParamsTable{1};
    end

    % Generate up to nSims rows
    if jitterFlag
        % Get the parameter names
        paramNames = neuronParamValuesTable.Properties.RowNames;

        % Count the number of parameters
        nParams = numel(paramNames);

        % Get the original parameter values
        paramValuesOrig = neuronParamsTable.Value;

        % Get the lower bounds
        paramMins = neuronParamsTable.LowerBound;

        % Get the upper bounds
        paramMaxs = neuronParamsTable.UpperBound;

        % Get the jitter percentages
        paramJitters = neuronParamsTable.JitterPercentage;

        % Get whether it is log-scaled
        paramIsLogs = neuronParamsTable.IsLog;

        % Generate up to nSims parameter sets with added jitter 
        paramValuesNew = zeros(nSims, nParams);
        parfor iSim = 1:nSims
            for iParam = 1:nParams
                % Get the current values
                pOrig = paramValuesOrig(iParam);
                pMin = paramMins(iParam);
                pMax = paramMaxs(iParam);
                pJitter = paramJitters(iParam);
                pIsLog = paramIsLogs(iParam);

                % Get a new jitter percentage
                jitter = (pJitter/100) * (-1 + 2*rand);

                % Add jitter to get new value
                if pIsLog
                    pNew = exp(log(pOrig) * (1 + jitter));            
                else
                    pNew =pOrig * (1 + jitter);            
                end

                % Make sure the new value does not exceed bounds
                if pNew < pMin
                    pNew = pMin;
                elseif pNew > pMax
                    pNew = pMax;
                end
                    
                % Store the new value
                paramValuesNew(iSim, iParam) = pNew;
            end
        end

        % Build the table with the parameter names as variables names
        simParamsTable = ...
            array2table(paramValuesNew, 'VariableNames', paramNames);
    else
        % Extract just the value column to get a simpler table
        neuronParamValuesTable = neuronParamsTable(:, 'Value');

        % Transpose the table
        neuronParamValuesTableTransposed = ...
            transpose_table(neuronParamValuesTable, 'RowNames', 'suppress');

        % Just duplicate the rows of the table
        simParamsTable = repmat(neuronParamValuesTableTransposed, nSims, 1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Get all neuron parameter names
neuronParamNames = neuronParamsTable.Properties.RowNames;

% Get all neuron parameter values
neuronParamValues = neuronParamsTable.Value;

% Extract neuron parameters to workspace
cellfun(@(x, y) eval(sprintf('%s = %g;', x, y)), ...
        neuronParamNames, neuronParamValues, 'UniformOutput', false);

neuronParamValuesTableTransposed = ...
    cellfun(@(x, y) transpose_table(x, 'RowNames', y), ...
            neuronParamValuesTable, num2cell(1:nSims), ...
            'UniformOutput', false);

% Unpack info from outparams
currentPulseAmplitude = outparams.currentPulseAmplitude;
gababsyn = outparams.gabab;
cprFlag = outparams.cprFlag;
holdPotential = outparams.holdPotentialCpr;
holdCurrent = outparams.holdCurrentCpr;
holdCurrentNoise = outparams.holdCurrentNoiseCpr;
holdPotential = outparams.holdPotential;
holdCurrent = outparams.holdCurrent;
holdCurrentNoise = outparams.holdCurrentNoise;
tstop = outparams.cprWindow(2);
tstop = outparams.ipscrWindow(2);    

cellfun(@(x) eval('length(%s)', x), simParamsFromArguments)

% Add simulation modes to table
simParamsTable.simMode = repmat({simMode}, nSims, 1);

% Add output file paths to table
simParamsTable.outFilePath = outFilePath;

%                   - 'CprFlag': whether to simulate current pulse responses
%                                   instead of GABA-B IPSC responses
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
cprFlagDefault = true; %false;  % fit active parameters by default
addParameter(iP, 'CprFlag', cprFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
cprFlag = iP.Results.CprFlag;

% Decide on the simulation mode and make it a cell array
if strcmpi(simMode, 'auto')
    if cprFlag
        simMode = 'passive';
    else
        simMode = 'active';
    end
end

eval(sprintf('simParamsTable.%s = %s;', thisParam, thisParam));

if outparams.neuronparams_use(iParam) ~= 0        % if parameter is checked
end

%                outparams.neuronparams(iParam) = outparams.neuronparams(iParam) ^ ...
%                    (1 + (outparams.neuronparams_jit(iParam)/100) * (-1 + 2*rand));
                    paramValuesNew(iSim, iParam) = outparams.neuronparams(iParam) * ...
                        (1 + (outparams.neuronparams_jit(iParam)/100) * (-1 + 2*rand));            

%                   - 'HoldPotentialCpr': holding potential for 
%                                           current pulse response (mV)
%                   must be a numeric vector
%                   default == -70 mV
%                   - 'HoldCurrentCpr': custom holding current 
%                                       for current pulse response (nA)
%                   must be a numeric vector
%                   default == 0 nA
%                   - 'HoldCurrentNoiseCpr': custom holding current noise
%                                               for current pulse response (nA)
%                   must be a numeric vector
%                   default == 0 nA
holdPotentialCprDefault = -70;  % (mV)
holdCurrentCprDefault = 0;      % (nA)
holdCurrentNoiseCprDefault = 0; % (nA)
addParameter(iP, 'HoldPotentialCpr', holdPotentialCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentCpr', holdCurrentCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'HoldCurrentNoiseCpr', holdCurrentNoiseCprDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
holdPotentialCpr = iP.Results.HoldPotentialCpr;
holdCurrentCpr = iP.Results.HoldCurrentCpr;
holdCurrentNoiseCpr = iP.Results.HoldCurrentNoiseCpr;
                        'HoldPotentialCpr', holdPotentialCpr, ...
                        'HoldCurrentCpr', holdCurrentCpr, ...
                        'HoldCurrentNoiseCpr', holdCurrentNoiseCpr, ...

% Decide on the holding potentials and holding currents
if cprFlag
end

% Get the number of elements
nOld = eval(sprintf('numel(%s)', thisParam));

% Match it up to the number of simulations
if nOld > 1 && nOld ~= nSims
    error('%s should have %d elements!', thisParam, nSims);
elseif nOld == 1 && nSims > 1
    eval(sprintf('%s = repmat(%s, %d, 1);', thisParam, thisParam, nSims));
end

gababTrise = 0;
gababTfallFast = 0;
gababTfallSlow = 0;
gababWeight = 0;

% If simulating IPSC responses, make sure the length of gababAmp
%   match up with the number of simulations
if ~cprFlag

%       cd/force_column_vector.m
gababAmp = force_column_vector(gababAmpRatioDefault * gababAmpTemplate');
gababTrise = force_column_vector(gababAmpRatioDefault * gababTriseTemplate');
gababTfallFast = force_column_vector(gababAmpRatioDefault * gababTfallFastTemplate');
gababTfallSlow = force_column_vector(gababAmpRatioDefault * gababTfallSlowTemplate');
gababWeight = force_column_vector(gababAmpRatioDefault * gababWeightTemplate');

elseif mod(nSims, nGababInputs) == 0
    % If nSims is divisible by nGababInputs,
    %   repeat nSims/nGababInputs times
    gababAmp = repmat(gababAmp, nSims/nGababInputs, 1);
    gababTrise = repmat(gababTrise, nSims/nGababInputs, 1);
    gababTfallFast = repmat(gababTfallFast, nSims/nGababInputs, 1);
    gababTfallSlow = repmat(gababTfallSlow, nSims/nGababInputs, 1);
    gababWeight = repmat(gababWeight, nSims/nGababInputs, 1);

gababAmp = gababAmp((1:nSims)' * (nGababInputs/nSims) + 1);
gababTrise = gababTrise((1:nSims)' * (nGababInputs/nSims) + 1);
gababTfallFast = gababTfallFast((1:nSims)' * (nGababInputs/nSims) + 1);
gababTfallSlow = gababTfallSlow((1:nSims)' * (nGababInputs/nSims) + 1);
gababWeight = gababWeight((1:nSims)' * (nGababInputs/nSims) + 1);

gababAmp = match_dimensions(gababAmp, nSims, 1);
gababTrise = match_dimensions(gababTrise, nSims, 1);
gababTfallFast = match_dimensions(gababTfallFast, nSims, 1);
gababTfallSlow = match_dimensions(gababTfallSlow, nSims, 1);
gababWeight = match_dimensions(gababWeight, nSims, 1);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
