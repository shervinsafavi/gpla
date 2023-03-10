%% sim_WavdEqu - Simulates the wave equation on a rectagular ...
                                    %medium while introducing arbitrary perturbations.
%
%  SimRes  =  SIM_WAVEEQU(StimVec,varargin)
%
%% Inputs
%
% * StimVec - vector of perturbations (time x 1), if zero no perturbation
% * is added at a given time point, for any integer value between 1 and 7,
% * it creates a perturbation of a specific type at the given time point.
%
%% Optional argument name,value pair
%
% * plot: plot the activity as it is computed   
% * stimdef: new definition of stimulation
% * storesubsampx: spatial subsampling of the stored pattern
% * edp_params: parameters of the PDE
% * storepattern: store the spatial temporal pattern in a 3d vector
%
%% Outputs
%
% * SimRes - spatio-temporal activity generated by the equation (vertical
% * coordinate x horizontal coordinate x time)
%

% Author : Michel Besserve, MPI for Intelligent Systems, MPI for Biological Cybernetics, Tuebingen, GERMANY

function [SimRes,EdpParams, StimVec, SimVarLabels]  =  simNeuralField(StimVec,varargin)
EdpParams = [];
StimStruct = {};
for karg = 1:2:length(varargin)
    switch lower(varargin{karg})
        case {'edp_params'}
            EdpParams = varargin{karg+1};
        case 'stim_struct'
            StimStruct = varargin{karg+1};
        otherwise
            error('unknown argument number %d',karg);

    end
end


StorePattern = 1;
StoreSubSampT = 50;
StoreSubSampX = 4;
% size of spatial domain in metric units
szy = .01; szx = .01;
dx = .0001; % space step;
dt = 0.00002; % time step;


v = 1; %1m/s axonal velocity for non-myelinated axon (local...)
r = .00044; % assume typical range of exc axons is .44mm
b = 1; % assume same input to excitatory

% keep this parameter for now: indicates a damping border size, setting it
% to one allows to insert the stimulaus inside the grid without errors.
DampGrid = 1;
tauE = .02;
tauI = .02;
% note: incresing inhibitory loop increases frequency when one acts on E-I
% synapse, when acting on I-E synapse, the additional effect of lowering
% the rate dyminishes the dynamic gain and thus the frequency too...
%nu = [1.1,.017;0,0];
%nu = 10*[.1,.17;1,0];% kind of baseline rate
%nu = 10*[.1,.1;.1,0]; too much exitation
%nu = 10*[1,1;1,0]; too much inhibiton
%nu = 10*[.8,.5;1,0]; %low pyr activity (8pct max), strong inhibiton
%nu = [8,4;100,0]; % difficult to find optimal point
nu = [4,4;20,0]; % now using S_I with factor 10 : update u(1,2) accordingly
S_E = @(m) 1./(1+exp(-m));
S_I = @(m) 1./(1+exp(-m/10));
fieldeqtype = 'normal';% normal or jirsa (this add a first order derivative term)
memDynType = '1ord';
tauS = 1/185; %synapse response time contant

% list of all parameters
parList = {'dt','dx','DampGrid','nu','tauE','tauI','v','r','b','szy','szx',...
        'StoreSubSampT','StoreSubSampX','S_E','S_I','fieldeqtype','tauS','memDynType'};
if ~isempty(EdpParams)
    EdpParams2 = structpack(parList);
    EdpParams = mergeStructures(EdpParams2,EdpParams);
else
    EdpParams = structpack(parList);
end

% need to unpack again to get the right grids (perhaps make this part
% cleaner), "true" forces overwrite
structunpack(EdpParams,true);
x = 0:dx:(szx-dx);
y = 0:dx:(szy-dx); % space
t = linspace(0,length(StimVec)*dt,length(StimVec));%0:dt:tm; % time


if StorePattern
    SimRes = zeros(length(x(1:StoreSubSampX:end)),length(y(1:StoreSubSampX:end)),...
        length(t(1:StoreSubSampT:end)),7);
end

[X,Y]  =  meshgrid(x,y);

Lx = length(x);
Ly = length(y);

ucurr = zeros(Ly,Lx); % initial value
upast = ucurr; % previose  =  curent  = > velocties  = 0
if isempty(StimStruct)
    StimStruct = generStimStruct(szy);
end
ncurr = ucurr;
vcurr = ucurr;
mEcurr = ucurr;
mIcurr = ucurr;
mEpast = ucurr;
mIpast = ucurr;
lambdaE = ucurr;
lambdaI = ucurr;

%unpack upated sturcture for call to interateEDP
stimState = 0;
ktime = 0;
for tt = t
    ktime = ktime+1;    
    % Stimulation:
    if StimVec(ktime)~= 0
        [ncurr,stimState] = StimulateEdp(0*ncurr,StimStruct{StimVec(ktime)},stimState,EdpParams);
    else
        ncurr = 0*ncurr;
    end
    switch EdpParams.memDynType
        case '2ord'
            [unext,vnext, mEnext, mInext,lambdaE,lambdaI] = ...
                iterateEDP(ucurr,vcurr, upast,mEcurr,mIcurr,ncurr,lambdaE,lambdaI,...
                EdpParams,mEpast,mIpast);
            mEpast = mEcurr;
            mIpast = mIcurr;
        case '1ord'
            [unext,vnext, mEnext, mInext,lambdaE,lambdaI] = ...
                iterateEDP(ucurr,vcurr, upast,mEcurr,mIcurr,ncurr,lambdaE,lambdaI,...
                EdpParams);
        otherwise
            error('unknow mem dyn type')
    end
    upast = ucurr;
    ucurr = unext;
    vcurr = vnext;
    mEcurr = mEnext;
    mIcurr = mInext;
    
     if StorePattern
        if mod(ktime-1,StoreSubSampT) == 0
            SimRes(:,:,ceil((ktime-1)/StoreSubSampT)+1,1) = ucurr(1:StoreSubSampX:end,...
                1:StoreSubSampX:end);
            SimRes(:,:,ceil((ktime-1)/StoreSubSampT)+1,2) = vcurr(1:StoreSubSampX:end,...
                1:StoreSubSampX:end);
            SimRes(:,:,ceil((ktime-1)/StoreSubSampT)+1,3) = mEcurr(1:StoreSubSampX:end,...
                1:StoreSubSampX:end);
            SimRes(:,:,ceil((ktime-1)/StoreSubSampT)+1,4) = mIcurr(1:StoreSubSampX:end,...
                1:StoreSubSampX:end);
            SimRes(:,:,ceil((ktime-1)/StoreSubSampT)+1,5) = lambdaE(1:StoreSubSampX:end,...
                1:StoreSubSampX:end);
            SimRes(:,:,ceil((ktime-1)/StoreSubSampT)+1,6) = lambdaI(1:StoreSubSampX:end,...
                1:StoreSubSampX:end);
            SimRes(:,:,ceil((ktime-1)/StoreSubSampT)+1,7) = ncurr(1:StoreSubSampX:end,...
                1:StoreSubSampX:end);
        end
    end
    
end
SimVarLabels = {'EPSC','IPSP','VE','VI','Erate','I rate','exo Input'};
end



function [u,stimState] = StimulateEdp(u,StimStruct,stimState,EdpParams)
[Zd,Zidx,stimState] = CreateEdpStim(StimStruct,stimState,EdpParams);
u(Zidx{1},Zidx{2}) = u(Zidx{1},Zidx{2})+Zd;
end



function [Zd,Zidx,stimState] = CreateEdpStim(StimStruct,stimState,EdpParams)

structunpack(StimStruct);
dx = EdpParams.dx;
switch lower(stimtype)

    case {'droplet','noise','filt_noise'}
        xd = -2*dsz:dx:2*dsz;
        yd = -2*dsz:dx:2*dsz;
        [Xd,Yd]  =  meshgrid(xd,yd);
        %generates a dropelet either oscillatory or of constant sign
        
        if ~Osc
            Zd = da*exp(-(Xd/dsz).^2-(Yd/dsz).^2);
        else
            Zd = da*exp(-(Xd/dsz).^2-(Yd/dsz).^2).*cos(2*pi*sqrt(Yd.^2+Xd.^2)/dsz);
        end
        
        switch lower(stimtype)
%             case 'bar'
%                 %make the wavelet unidimensional and amplify it a bit
%                  Zd = repmat(ba/da*Zd(ceil(size(Zd,1)/2),:),usz(1),1)*2;
%                 if ~exist('xside')
%                  xside = 1;
%                 end
%                  Zidx{2} = EdpParams.DampGrid+(0*2*dsz+1)*xside+...
%                      (usz(2)-2*EdpParams.DampGrid-0*2*dsz)*(1-xside)+(-ceil(2*dsz):4*dsz-ceil(2*dsz));
%                  Zidx{1} = 1:usz;
            case 'droplet'
                %Zidx = {EdpParams.DampGrid+(y0d-2*dsz:y0d+2*dsz),EdpParams.DampGrid+(x0d-2*dsz:x0d+2*dsz)};
                Zidx = {EdpParams.DampGrid+round((y0d+yd)/dx),EdpParams.DampGrid+round((x0d+xd)/dx)};
            case 'noise'
                %Zidx = {EdpParams.DampGrid+(y0d-2*dsz:y0d+2*dsz),EdpParams.DampGrid+(x0d-2*dsz:x0d+2*dsz)};
                Zidx = {EdpParams.DampGrid+round((y0d+yd)/dx),EdpParams.DampGrid+round((x0d+xd)/dx)};
                Zd = Zd*randn(1);
            case 'filt_noise'
                Zidx = {EdpParams.DampGrid+round((y0d+yd)/dx),EdpParams.DampGrid+round((x0d+xd)/dx)};
                randSamp = randn(1);                
                stimState = stimState+EdpParams.dt/EdpParams.tauE*(-stimState+randSamp);
                Zd = Zd*stimState;

        end
    otherwise
        error('unknown stimulation type %s',StimStruct.type)
end

end

function [unext,vnext, mEnext, mInext,lambdaE,lambdaI] = ...
    iterateEDP(ucurr, vcurr, upast, mEcurr, mIcurr, ncurr, lambdaE, lambdaI,...
    edpPars, mEpast, mIpast)

if nargin<10
    mEpast = [];
    mIpast = [];
end
% list parameters (use notations of Sanz-leon
% v: axonal velocity
% r: spatial effective axonal range
% gamma: temporal damping coefficient, gamma = v/r, 
% b: multiplicative factor for ratio of magnitude of exogenous input to I cells
% relative to E cells
% nui_j: synaptic coupling from population j to population i
parList = {'dt','dx','nu','tauE','tauI','v','r','b','S_E','S_I','fieldeqtype','tauS','memDynType'};
for kpar = 1:length(parList)
    eval([parList{kpar} ' = edpPars.' parList{kpar} ';'])
end
gamma = v/r;

D = [0 1 0; 1 -4 1; 0 1 0]; % 2d laplace operator
R = r/dx;
G = gamma*dt;

% tauE, tauI: E/I membrane time constants

% mE is excitatory membrane potential
% mI inhibitory
% ncurr: exogenous input current (normalize by RL to give potential) to pyr cell
TE = tauE/dt;
TI = tauI/dt;
switch memDynType
    case '1ord'
        mEnext = mEcurr + 1/TE*(-mEcurr+nu(1,1)*ucurr-nu(1,2)*vcurr+ncurr);

        mInext = mIcurr + 1/TI*(-mIcurr+nu(2,1)*ucurr-nu(2,2)*vcurr+b*ncurr);
    case '2ord'
        TS = tauS/dt;
        mEnext = 1/(TS*TE+TS/2+TE/2)*((2*TS*TE-1)*mEcurr-TS*TE*mEpast+(TE+TS)/2*mEpast+...
            nu(1,1)*ucurr-nu(1,2)*vcurr+ncurr);
        mInext = 1/(TS*TI+TS/2+TI/2)*((2*TS*TI-1)*mIcurr-TS*TI*mIpast+(TI+TS)/2*mIpast+...
            nu(2,1)*ucurr-nu(2,2)*vcurr+b*ncurr);

    otherwise
        error('unknown membrane dynamics')
end
lambdaEpast = lambdaE;
lambdaIpast = lambdaI;
lambdaE = S_E(mEnext);
lambdaI = S_I(mInext);


%uconv = conv2(ucurr,D,'same');
%periodic boundary conditions
ucurrExt = [ucurr(end,end),ucurr(end,:),ucurr(end,1);...
    ucurr(:,end),ucurr,ucurr(:,1);...
    ucurr(1,end),ucurr(1,:),ucurr(1,1)];
uconv = conv2(ucurrExt,D,'same');
uconv = uconv(2:end-1,2:end-1);

switch fieldeqtype
    case 'jirsa'
        unext = 1/(G+1) * ( (G-1)*upast + 2*ucurr + G^2*( R^2*uconv  + lambdaE +1/G*(lambdaE-lambdaEpast)-ucurr)); 
    case {'inst','instantaneous'}
        unext = lambdaE + R^2*uconv;
        
    otherwise
        unext = 1/(G+1) * ( (G-1)*upast + 2*ucurr + G^2*( R^2*uconv  + lambdaE -ucurr));
end
vnext = lambdaI;



end


% 
% function lambda = S_E(m)
%     lambda = 1./(1+exp(-m));
% 
% end
% function lambda = S_I(m)
% % use smoother transition to avoid exponential decrease of osc frequency
% % when inhibitory activity increases
%     lambda = 1./(1+exp(-m/10));
% 
% end


