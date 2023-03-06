%%
% second, same with varying rate instead of trials

nunit = 50;
nlfp = 55;
if nunit>=nlfp, 
    N=nlfp,p=nunit,
else
    N=nunit,p=nlfp,
end
c = N/p;

normPLV = true;

histAxSpec = [0:.1:6];
kappa= .5;
T = 1;
dt = .001;
t = 0:dt:T;
refRate = 5*ones(1,nunit);
refRate(1:2:nunit) = 30;
Ntrial = [1,5,10,20,50,100];
clear NSR*

rateFact = [.1,.2,.5,1,2,5,10];

%%
for krate = 1:length(rateFact)
    kTrial = 1;
    nTrial = Ntrial(kTrial);    
    nSim = 100;


    popN = nunit;
    popIdx = 1:popN;
    Sigma = [.1,.2,.5,1,2,5,10];
    for ksigma = 1:length(Sigma)
        PLV = zeros(nlfp,nunit,nSim);
        Nspk = zeros(nunit,nSim);
        sigma = Sigma(ksigma);
        lfpCorr = zeros(nlfp,nlfp,nSim);
        for ksim = 1:nSim
            for ktrial = 1:nTrial
                spkTimes = rand(length(t),nunit)<(ones(length(t),1)*refRate*dt*rateFact(krate));
                spkTimes(:,popIdx) = rand(length(t),popN)<(exp(kappa*cos(2*pi*t'))*refRate(popIdx)*dt*rateFact(krate));

                % noisy narrow band lfp: filtered noise is a sinusoid with random
                % complex coefficient
                %                lfp = (ones(nlfp,1)*exp(1i*2*pi*t)).*(1+sigma*(randn(nlfp,1)+1i*randn(nlfp,1))*ones(1,length(t)));
                lfp = (ones(nlfp,1)*exp(1i*2*pi*t))+sigma*(randn(nlfp,length(t))+1i*randn(nlfp,length(t)));
                PLV(:,:,ksim) = PLV(:,:,ksim)+lfp*spkTimes;
                Nspk(:,ksim) = Nspk(:,ksim)+sum(spkTimes,1)';
                lfpCorr(:,:,ksim) = lfpCorr(:,:,ksim)+lfp*lfp';
            end
        end
        PLV=PLV/T/nTrial;
        lfpCorr = lfpCorr/nTrial;
        truePLV = ones(nlfp,1)*integral(@(t)exp(1i*2*pi*t).*exp(kappa*cos(2*pi*t)),0,T)*refRate*rateFact(krate)/T;

        svdEst = zeros(nlfp,nunit,nSim);

        for ksim = 1:nSim
            M = (sqrtm(lfpCorr(:,:,ksim)));
            [u,d,v] = svds(inv(M)*PLV(:,:,ksim)*diag(1./sqrt(rate)),2);

            svdEst(:,:,ksim) = M*u(:,1)*d(1,1)*v(:,1)';
            % svdEst(:,:,ksim) = M*u(:,1)*d(1,1)/abs(d(1,1))*sqrt(abs(d(1,1)).^2-abs(d(2,2)).^2)*v(:,1)'*diag(sqrt(rate));
        end

        svdErr = svdEst-truePLV;
        classicErr = PLV- truePLV;

        svdRMS = sum(sum(abs(svdErr).^2,1),2);
        classicRMS = sum(sum(abs(classicErr).^2,1),2);

        NSRclass(ksigma,krate) = sqrt(mean(classicRMS)/mean(sum(sum(abs(truePLV).^2,1),2)));
        NSRsvd(ksigma,krate) = sqrt(mean(svdRMS)/mean(sum(sum(abs(truePLV).^2,1),2)));
    end
end


%%

