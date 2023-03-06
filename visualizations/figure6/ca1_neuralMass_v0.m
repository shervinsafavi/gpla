% load data
tdp = fullfile(pds.src, 'explorations', '261_shp_fineGrainedHighFreq', 'dat', 'varsForMB.mat');
load(tdp); 

%%

bands = reshape([freqBands{:}],2,[])';
clear dphi
for ifreq = 1:size(spkVec.ca1,2)
    
    ca1_spkVec_exc = spkVec.ca1(cellLabels == 3, ifreq);
    ca1_spkVec_inh = spkVec.ca1(cellLabels == 4, ifreq);

    dphi(ifreq) = (mean(ca1_spkVec_exc) ./ mean(ca1_spkVec_inh));
    
end

%% plot GPLA-based results
semilogx(mean(bands,2),180/pi*unwrap(angle(dphi)),'.--', 'markersize',20, 'linewidth', vc.f11.EIphaseShift.lw);
hold all

%% neural
% case that works but criticisable
tau = .02;
delta = .005;
nu21 = .5;
nu22 = .25;
nu12 = .01;


p = tf('s');

% EI ratio up to multiplicative constant
% alpha =0
EIratio_noFFI = 1+delta*p;
% alpha =1
EIratio_FFI = (1+delta*p-nu12)/(1+tau*p+nu21);

W = logspace(0,log10(180),100);
[MAG,PHASE] =bode(EIratio_noFFI,W);
[MAG2,PHASE2] = bode(EIratio_FFI,W); 

semilogx(W,squeeze(PHASE),W,squeeze(PHASE2), 'linewidth', vc.f11.EIphaseShift.lw);
% legend('GPLA outcome','Mass2D without FFI','Mass2d with FFI')
ylabel('E-I phase shift [deg]')
xlabel('Frequency (Hz)')


%% %%%%%%%%%%%%%%%%%alpha function for synapses

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
C = [1,0,0,0,0,0;...
     0,1,0,0,0,0];    

 
%DelayT(1) = struct('delay',gamma,'a',Ad,'b',[],'c',[],'d',[]);   % tau1=0.5
sys = ss(A,[1;alpha;0;0;0;0],C,zeros(2,1));


[MAG,PHASE] =bode(sys,W);

semilogx(W,mod(squeeze(PHASE(1,1,:)-PHASE(2,1,:))+90,360)-90, 'linewidth', vc.f11.EIphaseShift.lw);


phasediff = PHASE(1,1,:)-PHASE(2,1,:);
clear phaseVal
for kband = 1:size(bands,1)
    freqInt = W<bands(kband,2) & W>=bands(kband,1);
    phaseVal(kband) = 180/pi*angle(mean(exp(1i*phasediff(freqInt)*pi/180)));
end


% legend('GPLA outcome','Mass2D without FFI','Mass2d with FFI','Massalpha')

xlim([0, 180])


semilogx(mean(bands,2),phaseVal,'+--','markersize', 9, 'linewidth', vc.f11.EIphaseShift.lw);
lh = legend('GPLA outcome','Mass2D without FFI','Mass2D with FFI', ...
       'MassAlpha','MassAlphaFilt', 'location','northwest')

displaceFigureStuff(lh, [NaN -.1 NaN NaN])
set(lh, 'location','northwest');
set(lh, 'box','off');

box off

set(gca, ...
    'TickLength', [0.03 0.035], ...
    'TickDir',    'out', ...
    'ytick',      (-200:100:300) ...
    );

set(gca, 'fontsize', 8.5)