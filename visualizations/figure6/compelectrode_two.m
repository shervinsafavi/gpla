function [Electr,Ve] = compelectrode_two(Xsoma,Xdend,Nlabels,varargin)
%This function presents a (simplistic) electrode model for the CA3-CA1
%network models. Here we treat compartments as line sources of current.
%
%We assume the array of neurons to be linear (e.g., CA1) and only the
%length of the compartments is considered. 
%The number of electrode tips is variable, but the inter-tip
%distance is assumed 20um. 
%
%The computation works basically in the spirit of the Ohm's law. 
%
%References:
%   [1] Schomburg et al. (2012), J Neurosci.*
%   [2] Taxidis et al. (2011), Hippocampus.*
%   [3] Cannon et al. (1999), J Comp. Neurol.
%   [4] Gold et al. (2006), J Neurophysiol.
%   [5] Einevoll et al. (2013), Nature
%   [6] Nicholson (1973), IEEE Trans. Biomed. Eng.*
%
%Inputs:
%Xsoma: Somatic currents (both pyramidal and interneurons) (Neurons x Time)
%Xdend: Dendritic currents  (Neurons x Time)
%Nlabels: Neuron labels
%
%
%Outputs:
%Ve: Raw multi-site simulated LFP signal
%
% JF Ramirez-Villegas (c) (May 2014), MPI for Biological Cybernetics
%
% Report bugs/comments to Juan F. Ramirez-Villegas
% juan.ramirez-villegas -at- tuebingen.mpg.de
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
%LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
%IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
%WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH
%THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



if nargin<3
    fprintf('n\Not enough input arguments\n')
    help electrode_two1
    return
end
options.elepos  = 0.5;
options.Tipdist = 20;
options.Idist   = 1;
options.Lsoma   = [80 80];
options.Ldend   = [200 200];
n = size(Xsoma,1);
Maxdist = options.Idist*(n);
centers = linspace(0,Maxdist,n); 
Lsoma1      = 100; %Soma and dendritic compartments are separated by 
Nele        = 32; %number of electrodes
str_tick    = [];
for karg=1:length(varargin)
    switch lower(varargin{karg})
        case 'ypos'
            str_tick = varargin{karg + 1}; 
    end
end
if ~isempty(str_tick)
    cell_ypos = 100 + (str_tick - options.Lsoma(1)).*rand(length(Nlabels),1);
else
    cell_ypos = 100.*ones(length(Nlabels),1);
end

rho = (3.3); %conductivity of the extracellular medium [uS/um] is 0.3 or 
                %3.3 uOhm*um
%Linear probe
Electr1 = linspace(0,options.Tipdist*(Nele(1)-1),Nele(1));
Electr(:,:,2) = [Electr1];
Vsoma = zeros(size(Electr,1),size(Electr,2),1,size(Xsoma,2));
Vdend = zeros(size(Electr,1),size(Electr,2),1,size(Xsoma,2));
Elepos = Maxdist*options.elepos; %Electrode position in the array of neurons (center if L/2)
Electr(:,:,1) = Elepos;

for ktip=1:size(Electr,1)
    for ktip1=1:size(Electr,2)
        for kcell=1:size(Xsoma,1)
            if Nlabels(kcell)==1
                Isyn = Xsoma(kcell,:);
                %compute distance from the electrode to the soma (r)
                r = abs(Electr(ktip,ktip1,1) - centers(kcell)); %from electrode X-axis
                %compute distance from the electrode to soma (h)
                h = ((-Electr(ktip,ktip1,2)) + cell_ypos(kcell)); %from electrode Y-axis / Top of the compartment in zero by default
                %compute distance from the electrode to soma (s)
                s = (-Electr(ktip,ktip1,2) + (options.Lsoma(1) + cell_ypos(kcell)));%s distance is just the difference of h (longitudinal 
                                    %distance from electrode to the top of the compartment)
                                        %and the length of the soma 
                %Compute distances
                d0 = sqrt(h^2 + r^2); d1 = sqrt(s^2 + r^2);
                %compute the actual potential based on Ohm's Law
                Vsoma(ktip,ktip1,1,:) = squeeze((Vsoma(ktip,ktip1,1,:)))' + ...
                        ((rho*Isyn)./(4*pi*options.Lsoma(1))).*log10(abs((d0 - h)./(d1 - s)));
            if ~isempty(Xdend)
                %Compute dendrites
                Isyn = Xdend(kcell,:);
                %compute distance from the electrode to the dendrite (r)
                r = abs(Electr(ktip,ktip1,1) - centers(kcell)); %from electrode X-axis
                %comput distance from the electrode to dendrites (h) (vert)
                h = (-Electr(ktip,ktip1,2) + (Lsoma1 + options.Lsoma(1) + cell_ypos(kcell)));
                %comput distance from the electrode to dendrites (s) (vert)
                s = (-Electr(ktip,ktip1,2) + (Lsoma1 + options.Lsoma(1) + cell_ypos(kcell) + options.Ldend(1)));
                %Compute distances
                d0 = sqrt(h^2 + r^2); d1 = sqrt(s^2 + r^2);
                %compute the actual potential based on Ohm's Law
                Vdend(ktip,ktip1,1,:) = squeeze((Vdend(ktip,ktip1,1,:)))' + ...
                        ((rho*Isyn)./(4*pi*options.Ldend(1))).*log10(abs((d0 - h)/(d1 - s)));
            end
            else %when the label is 2 (from here, obsolete)
                Isyn = Xsoma(kcell,:);
                %compute distance from the electrode to the soma (r)
                r = abs(Electr(ktip,ktip1,1) - centers(kcell)); %from electrode X-axis
                %compute distance from the electrode to soma (h)
                h = (-(Electr(ktip,ktip1,2)) + cell_ypos(kcell)); %from electrode Y-axis / Top of the compartment in zero by default
                %compute distance from the electrode to soma (s)
                s = (-Electr(ktip,ktip1,2) + (cell_ypos(kcell) + options.Lsoma(2)));%s distance is just the difference of h (longitudinal 
                                    %distance from electrode to the top of the compartment)
                                        %and the length of the soma 
                %Compute distances
                d0 = sqrt(h^2 + r^2); d1 = sqrt(s^2 + r^2);
                %compute the actual potential based on Ohm's Law
                Vsoma(ktip,ktip1,1,:) = squeeze((Vsoma(ktip,ktip1,1,:)))' + ...
                    ((rho*Isyn)./(4*pi*options.Lsoma(2))).*log10(abs((d0 - h)./(d1 - s)));
            if ~isempty(Xdend)
                %Compute dendrites
                Isyn = Xdend(kcell,:);
                %compute distance from the electrode to the dendrite (r)
                r = abs(Electr(ktip,ktip1,1) - centers(kcell)); %from electrode X-axis
                %comput distance from the electrode to dendrites (h) 
                h = (-Electr(ktip,ktip1,2) + (Lsoma1 + options.Lsoma(2) + cell_ypos(kcell)));
                %comput distance from the electrode to dendrites (s)
                s = (-Electr(ktip,ktip1,2) + (Lsoma1 + options.Lsoma(2) + cell_ypos(kcell) + options.Ldend(2)));
                %Compute distances
                d0 = sqrt(h^2 + r^2); d1 = sqrt(s^2 + r^2);
                %compute the actual potential based on Ohm's Law
                Vdend(ktip,ktip1,1,:) = squeeze((Vdend(ktip,ktip1,1,:)))' + ...
                    ((rho*Isyn)./(4*pi*options.Ldend(2))).*log10(abs((d0 - h)/(d1 - s)));
            end
            end
        end
    end
end
%Sum contributions of soma, dendrites and all cells considered in the model
Ve = squeeze(sum(Vsoma + Vdend,3));