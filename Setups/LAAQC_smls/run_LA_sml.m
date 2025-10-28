%% F0AM-v4.3.0.1 LA Box Simulations - LAAQC-2020/2022 Terpene Study
%
% Purpose:
%   Simulates the diurnal evolution of terpenes under conditions
%   representative of the LAAQC-2020/2022 campaign.
%   Detailed description is provided in
%   "Diurnal Trends Differentiate Anthropogenic and Biogenic Terpenes
%   in the Los Angeles Basin", Section S6.
%   Tasnia et al., in preparation, Geophys. Res. Lett., 2025.
%
% Model setup:
%   Location      : Pasadena, CA
%   Date/time     : 30 June 2022 (UTC–7)
%   Chemical mech.: SAPRC-18
%   Time step     : 1 h (24 hourly steps)
%
% Input data:
%   - Meteorology (T, P, RH) and oxidants (O3, NO, NO2): LAAQC + CITAQS (Parker et al., 2020)
%   - OH: CalNex-2010, scaled x2 (Van Rooy et al., 2021)
%   - PBLH: Zhang et al. (2020), "Summer1" diurnal profile
%   - Emissions:
%       * Species:
%           API: α-pinene, BPI: β-pinene, LIM: limonene, ISP: isoprene
%       * Biogenic (mol km^-2 h^-1):
%           (a) MEGAN-Schlaerth (Schlaerth et al., 2023), "BE_<species>_a"
%           (b) BEIS-Stockwell (Stockwell et al., 2025), "BE_<species>_b"
%       * Anthropogenic (mol km^-2 h^-1), from SUNVEx 2021 inventory
%         (https://csl.noaa.gov/groups/csl7/measurements/2021sunvex/emissions/):
%           (a) SUNVEx-15/10/75, "AE_<species>_a"
%           (b) SUNVEx-25/25/50, "AE_<species>_b"
%
% SAPRC-18 species names:
%   APINE : α-pinene
%   BPINE : β-pinene
%   DLIMO : limonene
%   ISOP  : isoprene
%   See <F0AM_root>/Chem/SAPRC18/ for other species names.
%
% Output:
%   Model results are saved in <F0AM_root>/Runs/<savname>.mat
%
% Usage:
%  1) Load F0AM and add to MATLAB path
%     >>> addpath(genpath('<F0AM_root>'));
%  2) Modify this script as desired
%     (e.g., change savname, emissions, model options)
%  3) Run the script
%     >>> run_LA_sml
%  4) Load and analyze output from <F0AM_root>/Runs/<savname>.mat
%     Species concentrations are available under thier SPARC-18 names
%     (e.g., S.Conc.APINE for α-pinene)

clear

%% Save name
savname = 'LA_output';  % Output run name (under Runs/)

%% Load data
files = {'LAAQC_inputs.mat'};  % Input data files to load

loaded_files = {};
for i = 1:length(files)
    if exist(files{i}, 'file')
        load(files{i});
        loaded_files{end+1} = files{i};
    else
        warning('File %s not found.', files{i});
    end
end
fprintf('Loaded files: %s\n', strjoin(loaded_files, ', '));

%% GEOMETRY, METEOROLOGY, AND EMISSIONS

% Timebase
a_hour = LA_sml.Time;

% Solar geometry inputs
time.year  = 2022;                  % year
time.month = 6;                     % month
time.day   = 30;                    % day
time.hour  = a_hour;                % hour of day
time.min   = zeros(size(a_hour));
time.sec   = zeros(size(a_hour));
time.UTC   = -7;                    % PDT

% Location (Pasadena, CA)
location.latitude  = 34.0549;       % Latitude, degrees (-90 to 90). North of the equator is positive.
location.longitude = -118.2426;     % Longitude, degrees (-180 to 180). East of the meridian is positive.
location.altitude  = 93;            % Altitude, meters above sea level.

% Calculate solar zenith angles
sun = sun_position(time, location);

Met = {...
%   names       %values
    'P'          LA_sml.Pres;           % Pressure, mbar
    'T'          LA_sml.Temp;           % Temperature, K
    'RH'         LA_sml.RH;             % Relative Humidity, %
    'SZA'        sun.zenith;            % solar zenith angle, degrees
    'kdil'       1/(24*60*60);          % dilution constant, /s
    'jcorr'      0.5;                   % optimizes comparison b/w model and observed NO/NO2
    'BLH'        LA_sml.BLH;            % Planetary Boundary Layer height, m

    %% Emissions (mol km^-2 h^-1): biogenic + anthropogenic
    % Bio ('BE_'): MEGAN-Schlaerth ('_a') or BEIS-Stockwell ('_b')
    % Anthro ('AE_'): SUNVEx-15/10/75 ('_a') or SUNVEx-25/25/50 ('_b')
    'EAPI'       double(LA_sml.BE_API_a) + double(LA_sml.AE_API_a);  % alpha-pinene
    'EBPI'       double(LA_sml.BE_BPI_a) + double(LA_sml.AE_BPI_a);  % beta-pinene
    'ELIM'       double(LA_sml.BE_LIM_a) + double(LA_sml.AE_LIM_a);  % limonene
    'EISP'       double(LA_sml.BE_ISP_a) + double(LA_sml.AE_ISP_a);  % isoprene
    };


%% CHEMICAL CONCENTRATIONS

InitConc = {...
    % names           conc(ppb)           HoldMe
    'O3'                LA_sml.O3            1;
    'OH'                LA_sml.OH            1;    
    'NO'                LA_sml.NO            1;
    'NO2'               LA_sml.NO2           1;
    };

%% CHEMISTRY
ChemFiles = {...
    'SAPRC18_k(Met)';
    'SAPRC18_FZS_2_Jvalues(Met)';
    'SAPRC18_AllRxns';
    'EmisRxns'
    };

%% DILUTION CONCENTRATIONS
BkgdConc = {'DEFAULT'       0};

%% OPTIONS
ModelOptions.Verbose        = 3;
ModelOptions.EndPointsOnly  = 1;
ModelOptions.LinkSteps      = 1;
ModelOptions.IntTime        = 3600; % 3600 seconds/hour
ModelOptions.TimeStamp      = LA_sml.Time;
ModelOptions.SavePath       = savname;  % saved under <F0AM_root>/Runs/


%% INPUT REPLICATION AND INTERPOLATION
% For this particular scenario, it might be desirable to modify the inputs in a few ways.
% This sections demonstrates how to do so.

% INTERPOLATION
% Inputs currently have a time resolution of 60 minutes, but this is pretty coarse (the sun can move
% a lot in 60 minutes). The InputInterp function allows you to interpolate all inputs to a finer
% time resolution. NOTES:
%   - If your native data is fast (e.g., 1 Hz), it is generally better practice to bin-average that 
%       data to your desired resolution rather than average down to 60 minutes and then interpolate as done here.
%   - Make sure you adjust ModelOptions.IntTime too!
% To turn this on, set the "0" to "1" below.
if 0
    dt = 1800; %time spacing, seconds
    
    Time_interp = (0:dt:(86400-dt))'/3600; %interpolation timebase, fractional hours (to match LA_sml.Time)
    circularFlag = 1; % time wraps around at midnight
    [Met,InitConc,BkgdConc] = ...
        InputInterp(Met,InitConc,BkgdConc,LA_sml.Time,Time_interp,circularFlag);
    ModelOptions.TimeStamp = Time_interp;
    ModelOptions.IntTime = dt;
end

% REPLICATION
% Sometimes you may want to run the same inputs for multiple times. Typically, this scenario would
% be ground-site observations over one or more days, and you need a "spin-up" for non-measured
% species. The InputReplicate function lets you do this. Note, this only makes sense to use if
% ModelOptions.LinkSteps = 1. This replaces the "ModelOptions.Repeat" functionality in model
% versions prior to F0AMv4.
% Here, we run the same contraints for 3 days.
% The output "repIndex" is used to separate the days with SplitRun later.
nRep = 1; %number of days to repeat
[Met,InitConc,BkgdConc,repIndex] = InputReplicate(Met,InitConc,BkgdConc,nRep);
ModelOptions.TimeStamp = repmat(ModelOptions.TimeStamp,nRep,1);

%% MODEL RUN
% Now we call the model. Note this may take several minutes to run, depending on your system.
% Output will be saved in the "SavePath" above and will also be written to the structure S.

S = F0AM_ModelCore(Met,InitConc,ChemFiles,BkgdConc,ModelOptions);

% clear Met InitConc ChemFiles BkgdConc ModelOptions
