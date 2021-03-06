function m3ha_compare_and_plot_across_IC2 (varargin)
%% Plot activation/inactivation and I-V curves across initial conditions
% Usage: m3ha_compare_and_plot_across_IC2 (varargin)
%
% Requires:
%       cd/compare_neuronparams.m
%
% File History:
% 2017-08-09 Created

%% Default output folder
outFolder = '/media/adamX/m3ha/optimizer4gabab/m3ha_compare_neuronparams/';

%% Parameters used in fitting
% Fixed parameters used in simulations
cm_fixed = 0.88;                    % specific membrane capacitance [uF/cm^2] (Destexhe & Neubig 1997)
Ra_fixed = 173;                     % axial resistivity [Ohm-cm] (Destexhe & Neubig 1997)
corrD_fixed = 7.954;                % dendritic surface area correction factor (Destexhe & Neubig 1997)

neuronparamnames = { ...
    'diamSoma', 'LDend1', 'diamDend1ToSoma', ...
    'LDend2', 'diamDend2To1', 'distDendPercent', ...
    'cm', 'Ra', 'corrD', 'gpas', 'epas', ...
    'pcabarITSoma', 'pcabarITDend0', 'pcabarITDend1', 'pcabarITDend2', ...
    'shiftmIT', 'shifthIT', 'slopemIT', 'slopehIT', ...
    'ghbarIhSoma', 'ghbarIhDend0', 'ghbarIhDend1', 'ghbarIhDend2', 'ehIh', 'shiftmIh', ...
    'gkbarIKirSoma', 'gkbarIKirDend0', 'gkbarIKirDend1', 'gkbarIKirDend2', ...
    'gkbarIASoma', 'gkbarIADend0', 'gkbarIADend1', 'gkbarIADend2', ...
    'gnabarINaPSoma', 'gnabarINaPDend0', 'gnabarINaPDend1', 'gnabarINaPDend2', ...
    };
neuronparamsDefaultDestexhe = [ ...
    38.42, 12.49, 0.2676, 84.67, 0.8268, 50, ...
    cm_fixed, Ra_fixed, corrD_fixed, 1e-5, -80, ...
    .2e-3, .2e-3, .2e-3, .2e-3, ...
    1, 1, 1, 1, ...
    2.2e-5, 2.2e-5, 2.2e-5, 2.2e-5, -43, 0, ...
    2.0e-5, 2.0e-5, 2.0e-5, 2.0e-5, ...
    5.5e-3, 5.5e-3, 5.5e-3, 5.5e-3, ...
    5.5e-6, 5.5e-6, 5.5e-6, 5.5e-6, ...
    ];
neuronparamsDefaultChristine = [ ...
    38.42, 12.49, 0.2676, 84.67, 0.8268, 68.6, ...
    0.789, Ra_fixed, corrD_fixed, 8.21e-6, -80.4, ...
    5e-6, 5e-6, 8.91e-6, 3.98e-6, ...
    -13.8, -4.8, 1.4, 1, ...
    1.1e-5, 1.1e-5, 1.1e-5, 1.1e-5, -43, 11.4, ...
    2.0e-5, 2.0e-5, 2.0e-5, 2.0e-5, ...
    5.5e-3, 5.5e-3, 5.5e-3, 5.5e-3, ...
    5.5e-6, 5.5e-6, 5.5e-6, 5.5e-6, ...
    ];
neuronparamsDefaultSingleNeuronFitting5 = [ ...
    36.2403, 120, 0.1, 117.2416, 0.7088, 50, ...
    cm_fixed, Ra_fixed, corrD_fixed, 3.2632e-5, -70.1867, ...
    2.8216e-7, 2.8216e-7, 1.8440e-6, 5.6634e-5, ...
    -13.8, -4.8, 1.4, 1, ...
    3.0206e-7, 3.0206e-7, 2.8128e-6, 1.0226e-6, -43, 11.4, ...
    2.0e-5, 2.0e-5, 2.0e-5, 2.0e-5, ...
    5.5e-3, 5.5e-3, 5.5e-3, 5.5e-3, ...
    5.5e-6, 5.5e-6, 5.5e-6, 5.5e-6, ...
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramValues = {neuronparamsDefaultDestexhe, ...
                neuronparamsDefaultChristine, ...
                neuronparamsDefaultSingleNeuronFitting5};
paramNames = {neuronparamnames, neuronparamnames, neuronparamnames};
suffixes = {'_Destexhe_1998', '_Christine', '_SingleNeuronFitting5'};
m3ha_compare_neuronparams (paramValues, paramNames, suffixes, ...
                        'OutFolder', outFolder);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


