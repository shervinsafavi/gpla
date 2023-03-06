function spikeTrain = gnrt_inhomogeneousPoissonSpikeTrains(FRmodulators, spikeTrainDuraion, SF, varargin)
% spikeTrain = gnrt_homogeneousPoissonSpikeTrains(firingRate, spikeTrainDuraion, SF, nCell)
% Generate an ihomogeneous Poisson spike train with the given modulated firing
% rate
%
% EXAMPLE:
% an example syntax (if it's needed)
%
% ------
% Input:
% 1     FRmodulator: nNeuron * nSample array
%       (no default value)
%       specify how each spiking unit should be modulated 
% 2     spikeTrainDuraion: scalar
%       (no default value)
%       length of your spike train in second
% 3     SF: scalar
%       (no default value)
%       sampling frequency
% (4)   nCell: scalar
%       (default value = 1)
%       number of cells (i.e. generate Poisson spike rrain for a population
% % (5)  method: string
% %       (default value: randBinSpk)
% Output:
% 1     spikeTrain: binary vector
%
% ------
% potential improvments:
% (1) add a method as optional input (see the methods on p30 of Abbott & Dayan
% (2) add an option to pass a function handle for the FR modulator
% (3) check if the size of FRmodulator is matching the rest of the inpuits
% ------
% Code Info:
%   creation: 2017-12-05 by SS (shervin.safavi@tuebingen.mpg.de)
%   modification:
%       $ 201?
% ------
% see also 


%% Handle optional inputs (varargin):
nOpVar = 0; % number of optional variable
nOpVar = nOpVar + 1; optionalVariables.nTr           = []; defaultValues{nOpVar} = size(FRmodulators, 3); 
nOpVar = nOpVar + 1; optionalVariables.nUnit         = []; defaultValues{nOpVar} = size(FRmodulators, 1);

optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);
% nCell = size(FRmodulator, 1);

%%
nBin   = spikeTrainDuraion * SF;
nUnit  = optionalVariables.nUnit;
nTr    = optionalVariables.nTr;

if (size(FRmodulators, 1) > 1 ...
        && optionalVariables.nUnit ~=  size(FRmodulators, 1))
    error('number of cells (nCell) passed by user is vague:')
end



%% Firing rate modulator
FRmodulators_raw = FRmodulators; % clear FRmodulator;
% check if there is any negative firing rate
if sum(FRmodulators_raw(:) < 0) > 0
    warning('You passed negative firing rate(s) in "FRmodulator". All values were shifted to have them all possitive')
    tmpShiftVal = -1 * min(FRmodulators_raw(:)); 
else 
    tmpShiftVal = 0;
end

FRmodulations = FRmodulators_raw + tmpShiftVal;
% r_max = max(FRmodulation(:));

if sum(FRmodulations(:) < 0) > 0
    error('There is something wrong with your FR modulator')
end

r_max = max(FRmodulations, [], 2);
% r_max = squeeze(max(FRmodulations, [], 2));

%% Generate homogeneous Poisson spike train
% nBin = spikeTrainDuraion * SF;
% probilityThrsld = firingRate / SF; 
% % this term corresponds to r_est(t) * dt see p30 of Abbott & Dayan
% spikeTrain = (rand(nCell, nBin) < probilityThrsld);

HP_spikeTrain = gnrt_homogeneousPoissonSpikeTrains(squeeze(r_max), spikeTrainDuraion, SF, ...
    optionalVariables.nTr, optionalVariables.nUnit);

%% Generate inhomogenious poisson process 
% generate inhomogenious poisson process through spike thinning 

if (size(FRmodulations, 1) == 1 && size(FRmodulations, 3) == 1)
%     rejectMat = repmat(FRmodulations, nUnit, 1, nTr) / r_max;
    tmp = repmat(FRmodulations, nUnit, 1, nTr) / r_max;
    rejectMat = tmp > rand(size(tmp));
end

% if we have variability across units, but not across trials
if (size(FRmodulations, 1) ~= 1 && size(FRmodulations, 3) == 1)
%     rejectMat = repmat(FRmodulations, 1, 1, nTr) ./ repmat(r_max, 1, nBin, nTr);
    tmp = repmat(FRmodulations, 1, 1, nTr) ./ repmat(r_max, 1, nBin, nTr);
    rejectMat = tmp > rand(size(tmp));
end

if (size(FRmodulations, 1) ~= 1 && size(FRmodulations, 3) ~= 1)
%     rejectMat = FRmodulations ./ repmat(r_max, 1, nBin);
    tmp = FRmodulations ./ repmat(r_max, 1, nBin);
    rejectMat = tmp > rand(size(tmp));
end

% if size(FRmodulations, 1) > 1
% %     rejectMat = (FRmodulation ./ r_max) > rand(size(HP_spikeTrain));
%     rejectMat = bsxfun(@rdivide, FRmodulations, r_max) > rand(size(HP_spikeTrain));
% else
%     FRmodulation_singleCell = FRmodulations;
%     FRmodulations = repmat(FRmodulation_singleCell, optionalVariables.nCell, 1);
%     rejectMat = bsxfun(@rdivide, FRmodulations, r_max) > rand(size(HP_spikeTrain));
% end

spikeTrain = rejectMat .* HP_spikeTrain;

