function spikeTrain = gnrt_homogeneousPoissonSpikeTrains(firingRates, spikeTrainDuraion, SF, varargin)
% spikeTrain = gnrt_homogeneousPoissonSpikeTrain(firingRate, spikeTrainDuraion, SF, nCell)
% Generate a homogeneous Poisson spike train with the given constant firing
% rate
%
% EXAMPLE:
% an example syntax (if it's needed)
%
% ------
% Input:
% 1     firingRate: scalar
%       (no default value)
% 2     spikeTrainDuraion: scalar
%       (no default value)
%       length of your spike train in second
% 3     SF: scalar
%       (no default value)
%       sampling frequency
% (4)   nCell: scalar
%       (default value = 1)
%       number of cells (i.e. generate Poisson spike rrain for a population
% (5)  method: string
%       (default value: randBinSpk)
% Output:
% 1     spikeTrain: binary vector
%
% ------
% potential improvments:
% (1) add a method as optional input (see the methods on p30 of Abbott&Dayan
% (2) ask user for the output type (sparse, binary, cell array, etc)
% (3) check if the firing rate is horizental transpose it
% ------
% Code Info:
%   creation: 2016-05-16 by SS (shervin.safavi@tuebingen.mpg.de)
%   modification:
%       $ 2018-04-11 multiple spike train w/ different firing (encoded in
%       firingRate)(only for 'based_ISI_wRefractPreiod' method)
% ------
% see also 


%% handle optional inputs (varargin):
nOpVar = 0; % number of optional variable
nOpVar = nOpVar + 1; optionalVariables.nTr           = []; defaultValues{nOpVar} = size(firingRates, 2); 
nOpVar = nOpVar + 1; optionalVariables.nUnit         = []; defaultValues{nOpVar} = size(firingRates, 1); 
nOpVar = nOpVar + 1; optionalVariables.method        = []; defaultValues{nOpVar} = 'based_spikingProbality'; 
nOpVar = nOpVar + 1; optionalVariables.refractPeriod = []; defaultValues{nOpVar} = NaN; 
% nOpVar = nOpVar + 1; optionalVariables.dim           = []; 
% defaultValues{nOpVar} = [size(firingRates, 1) size(firingRates, 2)]; 

optionalVariables = handleVarargin(varargin, optionalVariables, defaultValues);

%% 
nBin   = spikeTrainDuraion * SF;
nUnit  = optionalVariables.nUnit;
nTr    = optionalVariables.nTr;

if (size(firingRates, 1) > 1 ...
        && optionalVariables.nUnit ~=  size(firingRates, 1))
    error('number of units (nUnits) passed by user is vague')
end

% nUnit  = optionalVariables.dim(1);
% nTr    = optionalVariables.dim(2);

%% generate spikes
switch optionalVariables.method
    
    case 'based_spikingProbality'
        % if only of firing rate value is passed, then we have identical
        % time series for all units and trials
        if numel(firingRates) == 1
%             firingRates_all = repmat(firingRates, nUnit, nTr);
            probilityThrslds = repmat(firingRates, nUnit, nBin, nTr) / SF;
        end
        % if only firing rate value for different cell (but not different trials) 
        % is passed, then we have identical
        % time series for all trials
        if size(firingRates, 2) == 1 % the "firingRates" input doesn't have variability in 'trial' direction
%             tmp = nan(optionalVariables.nCell, 1, optionalVariables.nTr);            
%             firingRates_all = repmat(firingRates, 1, nTr);
            probilityThrslds = repmat(firingRates, 1, nBin, nTr) / SF;
        end   
        
        % if we have diff FR for all trials 
        if (sum(size(firingRates) ~= [1 1]) == 2)
            tmp = nan(nUnit, 1, nTr);      
            tmp(:,1,:) = firingRates;
            probilityThrslds = repmat(tmp, 1, nBin) / SF;
        end
        
        % *** if user pass variable firing rates but not diffent cell then
        % inform the how wired it might be 
        
%         probilityThrslds = firingRates_all / SF;
%         if size(firingRate,2) ~= 1, firingRate = firingRate'; end  
        
        % this term corresponds to r_est(t) * dt see p30 of Abbott & Dayan
%         spikeTrain = (rand(optionalVariables.nCell, nBin) < probilityThrsld);
%         spikeTrain = bsxfun(@lt, rand(optionalVariables.nCell, nBin, optionalVariables.nTr), probilityThrsld);
%         spikeTrain = (rand(nUnit, nBin, nTr) ...
%             < repmat(probilityThrslds, 1, nBin));
        spikeTrain = (rand(nUnit, nBin, nTr) < probilityThrslds);
%     case 'based_InterSpikeIntervals'
%         
% 	case 'based_ISI_wRefractPreiod'
%         % the main idea is stolen from here:
%         % https://github.com/dedan/aand/blob/master/prac3/generatePoissonTrains.m
%         for iCell = 1 : optionalVariables.nCell
%             iSpk = 1;
%             spikeTrain_cellArray{iCell}(iSpk) = -1 * log(rand(1)) / firingRates(iCell);
%             while spikeTrain_cellArray{iCell}(iSpk) <= spikeTrainDuraion
%                 iSpk = iSpk + 1;
%                 tmp_isi = (-1 * log(rand(1)) / firingRates(iCell)) + optionalVariables.refractPeriod; 
%                 spikeTrain_cellArray{iCell}(iSpk) = spikeTrain_cellArray{iCell}(iSpk - 1) + tmp_isi;
%             end
%         end
%         
%         % store spikes in a binary matrix
%         % *** (in a sperate loop as the later changes in the funtion might demand it)
%         spikeTrain = zeros(optionalVariables.nCell, nBin);
%         for iCell = 1 : optionalVariables.nCell
%             spikeTrain(iCell, round(spikeTrain_cellArray{iCell} * SF)) = 1;
%         end
        
    otherwise 
end