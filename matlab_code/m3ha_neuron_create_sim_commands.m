function [simCommands, simCmdsFilePath] = ...
                m3ha_neuron_create_sim_commands (simParamsTable, varargin)
%% Generates simulation commands to be read by NEURON from a table of simulation parameters

% Usage: [simCommands, simCmdsFilePath] = ...
%               m3ha_neuron_create_sim_commands (simParamsTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       simCommands     - simulation commands
%                       specified as a cell array of character vectors
%       simCmdsFilePath - simulation commands file path 
%                       specified as a character vector
% Arguments:    
%       simParamsTable  - table of simulation parameters, 
%                           with each row being a simulation, with variables:
%                               diamSoma
%                               LDend
%                               diamDend
%                               cm
%                               Ra
%                               gpas
%                               epas
%                               corrD
%                               useHH
%                               tauhMode
%                               pcabarITSoma
%                               pcabarITDend1
%                               pcabarITDend2
%                               shiftmIT
%                               shifthIT
%                               slopemIT
%                               slopehIT
%                               ghbarIhSoma
%                               ghbarIhDend1
%                               ghbarIhDend2
%                               ehIh
%                               shiftmIh
%                               gkbarIKirSoma
%                               gkbarIKirDend1
%                               gkbarIKirDend2
%                               gkbarIASoma
%                               gkbarIADend1
%                               gkbarIADend2
%                               gnabarINaPSoma
%                               gnabarINaPDend1
%                               gnabarINaPDend2
%                               buildMode
%                               outFilePath
%                               tstop
%                               holdPotential
%                               currentPulseAmplitude
%                               gababAmp
%                               gababTrise
%                               gababTfallFast
%                               gababTfallSlow
%                               gababWeight
%                               customHoldCurrentFlag
%                               holdCurrent
%                               holdCurrentNoise
%                       must be a 2d table
%       varargin    - 'Prefix': prefix to prepend to file names
%                   must be a character array
%                   default == ''
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'SaveSimCmdsFlag': whether to save simulation commands
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       cd/force_string_end.m
%
% Used by:    
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2018-10-22 Adapted from code in run_neuron_once_4compgabab.m
% 2018-11-01 Made corrD an argument of build()
% 2018-11-23 Fixed placement of adjust_Ih()
% 2018-12-14 Moved T channels to passive parameters
% 2019-10-12 Now allows sim params table to have diamDendToSoma instead 
%               of diamDend
% 2019-11-12 Updated usage of build()
% 2019-11-12 Removed T & h channels from passive fit
% 2019-11-13 Restored usage of build()
% 2019-12-19 Renamed simMode -> buildMode
% 2019-12-27 Made useHH an argument of build()
% 2020-01-02 Made tauhMode an argument of adjust_IT()
% 2020-03-09 Added useHH to sim()
% 2020-04-27 Added useCvode to sim()
% 2020-04-27 Added secondOrder to sim()
% 

%% Hard-coded parameters
simCmdsFileName = 'simulation_commands.txt';

%% Default values for optional arguments
prefixDefault = '';             % prepend nothing to file names by default
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
saveSimCmdsFlagDefault = true;  % save simulation commands by default

%% Default values for simulation parameters


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
addRequired(iP, 'simParamsTable', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'SaveSimCmdsFlag', saveSimCmdsFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, simParamsTable, varargin{:});
prefix = iP.Results.Prefix;
outFolder = iP.Results.OutFolder;
saveSimCmdsFlag = iP.Results.SaveSimCmdsFlag;

% Check relationships between arguments
% TODO

%% Preparation
% Make sure prefix ends with a '_'
prefix = force_string_end(prefix, '_', 'OnlyIfNonempty', true);

% Count the number of simulations
nSims = height(simParamsTable);

% Count the number of parameters
nParams = width(simParamsTable);

% Extract all parameter names
paramNames = simParamsTable.Properties.VariableNames;

% Extract each required variable from the table
diamSoma = simParamsTable.diamSoma;
LDend = simParamsTable.LDend;
if ismember('diamDend', paramNames)
    diamDend = simParamsTable.diamDend;
elseif ismember('diamDendToSoma', paramNames)
    diamDendToSoma = simParamsTable.diamDendToSoma;
    diamDend = diamDendToSoma .* diamSoma;
end
cm = simParamsTable.cm;
Ra = simParamsTable.Ra;
gpas = simParamsTable.gpas;
epas = simParamsTable.epas;
corrD = simParamsTable.corrD;
useHH = simParamsTable.useHH;
useCvode = simParamsTable.useCvode;
secondOrder = simParamsTable.secondOrder;
tauhMode = simParamsTable.tauhMode;
pcabarITSoma = simParamsTable.pcabarITSoma;
pcabarITDend1 = simParamsTable.pcabarITDend1;
pcabarITDend2 = simParamsTable.pcabarITDend2;
shiftmIT = simParamsTable.shiftmIT;
shifthIT = simParamsTable.shifthIT;
slopemIT = simParamsTable.slopemIT;
slopehIT = simParamsTable.slopehIT;
ghbarIhSoma = simParamsTable.ghbarIhSoma;
ghbarIhDend1 = simParamsTable.ghbarIhDend1;
ghbarIhDend2 = simParamsTable.ghbarIhDend2;
ehIh = simParamsTable.ehIh;
shiftmIh = simParamsTable.shiftmIh;
gkbarIKirSoma = simParamsTable.gkbarIKirSoma;
gkbarIKirDend1 = simParamsTable.gkbarIKirDend1;
gkbarIKirDend2 = simParamsTable.gkbarIKirDend2;
gkbarIASoma = simParamsTable.gkbarIASoma;
gkbarIADend1 = simParamsTable.gkbarIADend1;
gkbarIADend2 = simParamsTable.gkbarIADend2;
gnabarINaPSoma = simParamsTable.gnabarINaPSoma;
gnabarINaPDend1 = simParamsTable.gnabarINaPDend1;
gnabarINaPDend2 = simParamsTable.gnabarINaPDend2;
buildMode = simParamsTable.buildMode;
outFilePath = simParamsTable.outFilePath;
tstop = simParamsTable.tstop;
holdPotential = simParamsTable.holdPotential;
currentPulseAmplitude = simParamsTable.currentPulseAmplitude;
gababAmp = simParamsTable.gababAmp;
gababTrise = simParamsTable.gababTrise;
gababTfallFast = simParamsTable.gababTfallFast;
gababTfallSlow = simParamsTable.gababTfallSlow;
gababWeight = simParamsTable.gababWeight;
customHoldCurrentFlag = simParamsTable.customHoldCurrentFlag;
holdCurrent = simParamsTable.holdCurrent;
holdCurrentNoise = simParamsTable.holdCurrentNoise;

%% Generate commands for NEURON
% Build simulation commands to be read by NEURON through the here-document
simCommands = cell(nSims, 1);
parfor iSim = 1:nSims
%for iSim = 1:nSims
    % Start with the build() command in singleneuron4compgabab.hoc
    thisCmds = sprintf('build("%s", %g, %g, %g, %g, %g)\n', ...
                        buildMode{iSim}, diamSoma(iSim), LDend(iSim), ...
                        diamDend(iSim), corrD(iSim), useHH(iSim));

    % Command to adjust global passive parameters
    thisCmds = [thisCmds, sprintf('adjust_globalpas(%g, %g, %g)\n', ...
                                    cm(iSim), Ra(iSim), corrD(iSim))];

    % Command to adjust passive leak channel parameters
    thisCmds = [thisCmds, sprintf('adjust_leak(%g, %g, %g)\n', ...
                                    gpas(iSim), epas(iSim), corrD(iSim))];

    % Add more commands if in active fit mode
    if strcmp(buildMode{iSim}, 'active')
        % Command to adjust H channel parameters
        thisCmds = [thisCmds, ...
                    sprintf('adjust_Ih(%g, %g, %g, %g, %g, %g)\n', ...
                            ghbarIhSoma(iSim), ghbarIhDend1(iSim), ...
                            ghbarIhDend2(iSim), ehIh(iSim), ...
                            shiftmIh(iSim), corrD(iSim))];

        % Command to adjust T-type calcium channel parameters
        thisCmds = [thisCmds, ...
                    sprintf('adjust_IT(%g, %g, %g, %g, %g, %g, %g, %g, %g)\n', ...
                            pcabarITSoma(iSim), pcabarITDend1(iSim), ...
                            pcabarITDend2(iSim), shiftmIT(iSim), ...
                            shifthIT(iSim), slopemIT(iSim), ...
                            slopehIT(iSim), corrD(iSim), tauhMode(iSim))];

        % Command to adjust inward-rectifying potassium channel parameters
        thisCmds = [thisCmds, ...
                    sprintf('adjust_IKir(%g, %g, %g, %g)\n', ...
                            gkbarIKirSoma(iSim), gkbarIKirDend1(iSim), ...
                            gkbarIKirDend2(iSim), corrD(iSim))];

        % Command to adjust A-Type potassium channel parameters
        thisCmds = [thisCmds, ...
                    sprintf('adjust_IA(%g, %g, %g, %g)\n', ...
                            gkbarIASoma(iSim), gkbarIADend1(iSim), ...
                            gkbarIADend2(iSim), corrD(iSim))];

        % Command to adjust persistent sodium channel parameters
        thisCmds = [thisCmds, ...
                    sprintf('adjust_INaP(%g, %g, %g, %g)\n', ...
                            gnabarINaPSoma(iSim), gnabarINaPDend1(iSim), ...
                            gnabarINaPDend2(iSim), corrD(iSim))];
    end

    % The sim() command in singleneuron4compgabab.hoc
    %   NOTE: currentPulseAmplitude must be in nA and gababAmp must be in uS
    thisCmds = [thisCmds, ...
                sprintf(['sim("%s", "%s", %d, %g, %g, ', ...
                        '%g, %g, %g, %g, %g, %d, %g, %g, %g, %g, %g)\n'], ...
                        buildMode{iSim}, outFilePath{iSim}, tstop(iSim), ...
                        holdPotential(iSim), currentPulseAmplitude(iSim), ...
                        gababAmp(iSim), gababTrise(iSim), ...
                        gababTfallFast(iSim), gababTfallSlow(iSim), ...
                        gababWeight(iSim), customHoldCurrentFlag(iSim), ...
                        holdCurrent(iSim), holdCurrentNoise(iSim), ...
                        useHH(iSim), useCvode(iSim), secondOrder(iSim))];

    % Store in cell array
    simCommands{iSim} = thisCmds;
end

%% Save things
% Save simulation commands to a file if requested
if saveSimCmdsFlag
    % Construct full path
    simCmdsFilePath = fullfile(outFolder, strcat(prefix, simCmdsFileName));

    % Open file
    fid = fopen(simCmdsFilePath, 'w');

    % Print simulation commands to file
    for iSim = 1:nSims
        fprintf(fid, '%s\n', simCommands{iSim});
    end

    % Close the file
    fclose(fid);
else
    simCmdsFilePath = '';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fid = fopen(fullfile(outFolder, ...
        [prefix, '_simulation_commands_backup.txt']), 'w');

%   NOTE: GABA-B conductance amplitudes are converted from nS to uS
                    gababAmp(iSim)/1000, gababTrise(iSim), ...

% Note: the following method would preclude a parfor loop from being used
% Extract each variable of the table
for iParam = 1:nParams
    % Get the current parameter name
    thisParamName = paramNames{iParam};

    % Extract the current parameter from the table
    eval(sprintf('%s = simParamsTable.%s;', thisParamName, thisParamName));
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
