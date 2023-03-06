% *** should have different name 
% *** is not related to the original name 

% load data
% tdp = fullfile(pds.src, 'explorations', '261_shp_fineGrainedHighFreq', 'dat', 'varsForMB.mat');
tdp = fullfile(pds.ldat, 'explorations', '261_shp_fineGrainedHighFreq', 'dat', 'varsForMB.mat');
load(tdp); 

%%
bands = reshape([freqBands{:}],2,[])';
clear dphi
for ifreq = 1:size(spkVec.ca1,2)
    
    ca1_spkVec_exc = spkVec.ca1(cellLabels == 3, ifreq);

    ca1_spkVec_inh = spkVec.ca1(cellLabels == 4, ifreq);
    dphi(ifreq) = (mean(ca1_spkVec_exc)./mean(ca1_spkVec_inh));
end

centPhase = @(x) x+180

%% alpha function for synapses
W = logspace(0,log10(400),100);

tau = .025;
v11 = 0;
v21 = .3;
v22 = 12.0;
v12 = .85;
sigma = .015;

alpha = 1;
delta = tau;
gamma = .002;
  
A = [  -1/tau,         0   ,        0      ,         0         , v11/tau,      -v12/tau  ;...
           0      , -1/delta,         0      ,        0          ,v21/delta, -v22/delta;...
      1/sigma,        0    ,   -1/sigma,           0     ,       0     ,          0;...
        0        ,  1/sigma,            0    ,   -1/sigma ,      0       ,         0;...
        0        ,          0    ,  1/sigma   ,         0      ,  -1/sigma,        0;...
        0        ,          0    ,           0     ,   1/sigma ,       0     ,  -1/sigma];
% outputs; lambdaE, lambdaI, EPSP, iPSP, totLFP
    C = [1,0,0,0,0,0;...
     0,1,0,0,0,0;0,0,0,0,1,0;0,0,0,0,0,1;0,0,0,0,1,1];    
D = [0;0;0;0;1];
 
%DelayT(1) = struct('delay',gamma,'a',Ad,'b',[],'c',[],'d',[]);   % tau1=0.5
sys = ss(A,[1;alpha;0;0;0;0],C,D);


[MAG,PHASE] =bode(sys,W);

%%
for ifreq = 1:size(spkVec.ca1,2)
    
    ca1_spkVec_exc = spkVec.ca1(cellLabels == 3, ifreq);
    ca1_spkVec_inh = spkVec.ca1(cellLabels == 4, ifreq);
    
    phiE(ifreq) = 180/pi*angle(mean(ca1_spkVec_exc));
    phiI(ifreq) = 180/pi*angle(mean(ca1_spkVec_inh));
    
end


simPhaseE = squeeze(PHASE(1,1,:)-PHASE(4,1,:));
simPhaseI = squeeze(PHASE(2,1,:)-PHASE(4,1,:));


clear phaseEVal phaseIVal
for kband = 1:size(bands,1)
    freqInt = W<bands(kband,2) & W>=bands(kband,1);
    phaseEVal(kband) = 180/pi*angle(mean(exp(1i*pi/180*simPhaseE(freqInt))));
    phaseIVal(kband) = 180/pi*angle(mean(exp(1i*pi/180*simPhaseI(freqInt))));
end



unwrapDeg = @(x)180/pi*unwrap(pi/180*x);
centPhase2 = @(x) mod(x,360);

%% plot
semilogx(mean(bands,2),unwrapDeg(centPhase2(phiE)),'b.-',mean(bands,2),unwrapDeg(centPhase2(phiI)),'r.-','MarkerSize',20,'linewidth',2)
hold on
semilogx(mean(bands,2),(centPhase(phaseEVal)),'c+-',mean(bands,2),(centPhase(phaseIVal)),'m+-','MarkerSize',10,'linewidth',2)
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
box off


lh = legend('E GPLA','I GPLA','E MassAlpha filt.','I MassAlpha filt.', 'location','northwest')


displaceFigureStuff(lh, [NaN .02 NaN NaN])
% set(lh, 'location','northwest');
set(lh, 'box','off');



set(gca, ...
    'TickLength', [0.03 0.035], ...
    'TickDir',    'out', ...
    'ytick',      (0:100:400) ...
    );

set(gca, 'fontsize', 8.5)

