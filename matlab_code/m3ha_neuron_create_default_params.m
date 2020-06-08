function defaultTable = m3ha_neuron_create_default_params
%% Creates the default NEURON parameters table
% Usage: defaultTable = m3ha_neuron_create_default_params
% Explanation:
%       TODO
%
% Example(s):
%       defaultTable = m3ha_neuron_create_default_params
%
% Outputs:
%       defaultTable    - the default NEURON parameters table, with columns:
%                           Value
%                           InitValue
%                           LowerBound
%                           UpperBound
%                           Class
%                           IsLog
%                           IsPassive
%                           UseAcrossTrials
%                           UseAcrossCells
%                       specified as a 2D table
%
% Arguments:
%
% Requires:
%
% Used by:
%       cd/m3ha_decide_on_plot_vars.m
%       cd/m3ha_neuron_create_initial_params.m

% File History:
% 2019-12-30 Moved from m3ha_neuron_create_initial_params.m
% 

%% Hard-coded parameters
% Must be consistent with find_passive_params.m
% Default parameter values (most from Destexhe & Neubig 1997)
cmInit = 0.88;          % specific membrane capacitance [uF/cm^2]
RaInit = 173;           % axial resistivity [Ohm-cm]
corrDInit = 1;          % dendritic surface area correction factor
                        %   set to 1 so that curve-fitted parameters
                        %   from find_passive_params.m may be applied directly
                        % Destexhe et al 1998a used 7.954,
                        %   which was estimated by fitting 
                        %   voltage-clamp traces in Destexhe et al 1998a
% Default parameter values from Dextexhe
neuronParamsDefault = [ ...
    38.42; 84.67; 8.50; ...
    cmInit; RaInit; corrDInit; 2e-5; -70; ...
    .2e-3; .2e-3; .2e-3; ...
    1; 1; 1; 1; ...
    2.2e-5; 2.2e-5; 2.2e-5; -28; 0; ...
    2.0e-5; 2.0e-5; 2.0e-5; ...
    5.5e-3; 5.5e-3; 5.5e-3; ...
    5.5e-6; 5.5e-6; 5.5e-6; ...
    ];

% Names for each parameter
neuronParamNames = { ...
    'diamSoma'; 'LDend'; 'diamDend'; ...
    'cm'; 'Ra'; 'corrD'; 'gpas'; 'epas'; ...
    'pcabarITSoma'; 'pcabarITDend1'; 'pcabarITDend2'; ...
    'shiftmIT'; 'shifthIT'; 'slopemIT'; 'slopehIT'; ...
    'ghbarIhSoma'; 'ghbarIhDend1'; 'ghbarIhDend2'; 'ehIh'; 'shiftmIh'; ...
    'gkbarIKirSoma'; 'gkbarIKirDend1'; 'gkbarIKirDend2'; ...
    'gkbarIASoma'; 'gkbarIADend1'; 'gkbarIADend2'; ...
    'gnabarINaPSoma'; 'gnabarINaPDend1'; 'gnabarINaPDend2'; ...
    };

% Lower bounds for each parameter
neuronParamsLowerBound = [ ...
    8; 5; 1; ...
    cmInit; RaInit; corrDInit; 1.0e-6; -95; ...
    1.0e-9; 1.0e-9; 1.0e-9; ...
    -30; -30; 0.1; 0.1; ...
    1.0e-9; 1.0e-9; 1.0e-9; -32; -30; ...
    1.0e-9; 1.0e-9; 1.0e-9; ...
    1.0e-9; 1.0e-9; 1.0e-9; ...
    1.0e-9; 1.0e-9; 1.0e-9; ...
    ];

% Upper bounds for each parameter
neuronParamsUpperBound = [ ...
    300; 400; 100; ...
    cmInit; RaInit; corrDInit; 1.0e-4; -45; ...
    1.0e-1; 1.0e-1; 1.0e-1; ...
    30; 30; 10; 10; ...
    1.0e-2; 1.0e-2; 1.0e-2; -24; 30; ...
    1.0e-2; 1.0e-2; 1.0e-2; ...
    1.0e-1; 1.0e-1; 1.0e-1; ...
    1.0e-2; 1.0e-2; 1.0e-2; ...
    ];

% Parameter class for each parameter
neuronParamsClass = [ ...
    1; 1; 1; ...
    2; 2; 2; 3; 3; ...
    4; 4; 4; ...
    4; 4; 4; 4; ...
    5; 5; 5; 5; 5; ...
    6; 6; 6; ...
    7; 7; 7; ...
    8; 8; 8; ...
    ];

% Whether each parameter will be varied log-scaled
neuronParamsIsLog = logical([ ...
    0; 0; 0; ...
    1; 1; 0; 1; 0; ...
    1; 1; 1; ...
    0; 0; 1; 1; ...
    1; 1; 1; 0; 0; ...
    1; 1; 1; ...
    1; 1; 1; ...
    1; 1; 1; ...
    ]);

% Whether the parameter is considered a 'passive' parameter
%   Note: must be consistent with TC3.tem and m3ha_neuron_create_sim_commands.m
neuronParamsIsPassive = logical([ ...
    1; 1; 1; ...
    1; 1; 1; 1; 1; ...
    0; 0; 0; ...
    0; 0; 0; 0; ...
    0; 0; 0; 0; 0; ...
    0; 0; 0; ...
    0; 0; 0; ...
    0; 0; 0; ...
    ]);

% Whether the parameter is varied when fitting across trials for each cell
%   Note: All passive parameters and all active conductance densities
neuronParamsUseAcrossTrials = logical([ ...
    1; 1; 1; ...
    0; 0; 0; 1; 0; ...
    1; 1; 1; ...
    0; 0; 0; 0; ...
    1; 1; 1; 0; 0; ...
    1; 1; 1; ...
    1; 1; 1; ...
    1; 1; 1; ...
    ]);

% Whether the parameter is varied when fitting across cells
%   Note: All active conductance kinetics
neuronParamsUseAcrossCells = logical([ ...
    0; 0; 0; ...
    0; 0; 0; 0; 0; ...
    0; 0; 0; ...
    1; 1; 1; 1; ...
    0; 0; 0; 1; 1; ...
    0; 0; 0; ...
    0; 0; 0; ...
    0; 0; 0; ...
    ]);

% Column names for the parameters table
columnNames = {'Value'; 'InitValue'; 'LowerBound'; 'UpperBound'; 'Class'; ...
                'IsLog'; 'IsPassive'; 'UseAcrossTrials'; 'UseAcrossCells'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do the job
% Create the default table
defaultTable = table(neuronParamsDefault, neuronParamsDefault, ...
                neuronParamsLowerBound, neuronParamsUpperBound, ...
                neuronParamsClass, neuronParamsIsLog, neuronParamsIsPassive, ...
                neuronParamsUseAcrossTrials, neuronParamsUseAcrossCells, ...
                'RowNames', neuronParamNames, 'VariableNames', columnNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Force as column vectors
[neuronParamsDefault, neuronParamsLowerBound, ...
    neuronParamsUpperBound, neuronParamsClass, ...
    neuronParamsIsLog, neuronParamsIsPassive, ...
    neuronParamsUseAcrossTrials, neuronParamsUseAcrossCells] = ...
    argfun(@force_column_vector, ...
            neuronParamsDefault, neuronParamsLowerBound, ...
            neuronParamsUpperBound, neuronParamsClass, ...
            neuronParamsIsLog, neuronParamsIsPassive, ...
            neuronParamsUseAcrossTrials, neuronParamsUseAcrossCells);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%