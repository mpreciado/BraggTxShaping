%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bragg grating simulatinos and results display.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select Bragg grating design file
% 
% The program will generate plots showing the gratings, its sepectral response, and the result from applying a 5ps input gaussian pulse.
% 
% Check parameters at the begining of main.m.
% In particular, bandwidth and grating length can be modified by using the parameter ScalingFactor.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By Miguel Preciado @miguelSciLight
% 
% Please cite:
% Miguel A Preciado, Xuewen Shu, Kate Sugden "Proposal and design of phase-modulated fiber gratings in transmission for pulse shaping", Opt. Lett. 38, 70-72 (2013) 
% 
% More info:
% https://www.researchgate.net/publication/234042470_Proposal_and_design_of_phase-modulated_fiber_gratings_in_transmission_for_pulse_shaping
% https://sites.google.com/view/miguelscilight/resources

clear all;
close all;

%% PARAMETERS
nav=1.452; % average refractive index
LAMBDA0=1550/(2*nav); %Bragg central resonance (1550 nm default)


ScalingFactor=1/5; %This parameter enables the spatial - temporal scaling of gragings, signals, and their spectrums. 
                    %=1  original design and test signals.
                    
                    %>1  scaled gratings/signals. Shorter gratins, shorter
                    %pulses, wider bandwidth
                    
                    %<1  ...the opposite

%%
c=2.9979e8;


[filename, pathname] = uigetfile(  {'*.txt'},'Pick a file with an Bragg grating design');

S=importdata([pathname filename ]);

z=S(:,1)*ScalingFactor;
kappa=S(:,2)/ScalingFactor;
LAMBDA=S(:,3)/ScalingFactor;



%Bragg grating discretized as a sequence of reflectors
Az=mean(diff(z));
ro2=kappa*Az; %Equivalent reflection coefficient for a layer Az

K=2*pi./(1e-9*(LAMBDA+LAMBDA0)); %Spatial frequency
K0=mean(K);
K=K-K0;

phi=K*Az/2;

ro=ro2.*exp(-j*phi); %signo cambiado para corregirlo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation of the grating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ R,T,MA,MB ] = BraggSim(ro);
%%
%%
%PREPARATION OF FUNCTIONS TO PLOT


% Grating period. Trick to avoid problems with initial and final value. Not important.
LAMBDA_PM=LAMBDA;
LAMBDA_PM(1)=LAMBDA_PM(2);
LAMBDA_PM(end)=LAMBDA_PM(end-1);


v=c/(2*nav);
Ts=Az/v;

tiempo=Ts*(1:length(T));
AF=1/Ts;
fs=AF/length(T);
freq=fs*((1:length(T))-length(T)/2);
freq0=v*K0/(2*pi);

lambda=c./(freq+freq0);

rango=abs(lambda-1550e-9)<2e-9/ScalingFactor;
rango2=abs(lambda-1550e-9)<2e-9/ScalingFactor;

aux=(unwrap(angle(T)*2)/2);
aux=aux-mean((diff(aux(rango2))))*(1:length(aux));
aux=aux-mean(aux);


At=20e-12*ScalingFactor;
Np=At/Ts;


%iNPUT PULSE
FWHM=5e-12*ScalingFactor;
sigmat=FWHM/sqrt(4*log(2));
g=exp(-tiempo.^2/(2*sigmat^2));

sigmaf=1/(2*pi*sigmat);
shift_f=0;
G=exp(-((freq-shift_f).^2)/(2*sigmaf^2));

FIG=figure(6789); 
subplot(2,1,1);

[AX,H1,H2] = plotyy(z*100,LAMBDA_PM,z*100,abs(kappa));
set(AX,'FontSize',10, 'FontName','Arial')
set(get(AX(1),'Ylabel'),'String','Grating period variation [nm]','FontSize',12, 'FontName','Arial') 
set(get(AX(2),'Ylabel'),'String','Coupling coefficient [m^{-1}]','FontSize',12, 'FontName','Arial') 
set(get(AX(1),'Xlabel'),'String','z [cm]','FontSize',12, 'FontName','Arial');

set(H1,'LineWidth',2) 

set(H2,'LineStyle','-.')
set(H2,'LineWidth',2)

  set(AX, 'XLim', [0 max(z)*100]);
title('(a)');

subplot(2,1,2);
lambdaSimul=lambda;
hold on

[AX,H1,H2] = plotyy(lambda(rango)*1e9,20*log10(abs(T(rango))),lambda(rango)*1e9,aux(rango));
title('(b)');

hold on;

set(AX,'FontSize',12, 'FontName','Arial')
set(get(AX(1),'Ylabel'),'String','Transmission [dB]','FontSize',12, 'FontName','Arial') 
set(get(AX(2),'Ylabel'),'String','Transmission phase [rad]','FontSize',12, 'FontName','Arial') 
set(get(AX(1),'Xlabel'),'String','Wavelength [nm]','FontSize',12, 'FontName','Arial');
set(H1,'LineWidth',2) 
set(H2,'LineStyle')

set(H2,'LineStyle','-.')
set(H2,'LineWidth',2)

set(FIG,'Position',[0 0 600 500])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

offset=20e-12*ScalingFactor;
fin1=exp(-((tiempo-offset).^2/(2*sigmat^2)));

Fin1=fftshift(fft(fin1));

Fout1=Fin1.*T;
fout1=ifft(fftshift(Fout1));
fout1=fout1;
 
rangot=tiempo>0e-10*ScalingFactor & tiempo < 1.6e-10*ScalingFactor;

FIG=figure(635);
subplot(2,1,1);

[ax,h1,h2]=plotyy(1e12*(tiempo(rangot)),abs(fin1(rangot)).^2,1e12*(tiempo(rangot)),wrapToPi(angle(fin1(rangot))));

set(h1,'linewidth',2)% to change the first line
set(h2,'linewidth',2)% to change the first line

set(ax,'Fontsize',10)

xlabel('Time [ps]');
ylabel('Intensity [a.u.]');

title('Input signal ');


rangot=tiempo>=4.4e-10*ScalingFactor & tiempo <= 5.8e-10*ScalingFactor;
subplot(2,1,2);
[ax,h1,h2]=plotyy(1e12*(tiempo(rangot)),abs(fout1(rangot).^2),1e12*(tiempo(rangot)),wrapToPi(angle(fout1(rangot))));

set(h1,'linewidth',2)% to change the first line
set(h2,'linewidth',2)% to change the first line

set(ax,'Fontsize',10)

xlabel('Time [ps]');
title('Output signal ');
axes(ax(2));
ylabel('Phase [rad]');
%xlim([min(xlim),max(xlim)]);


set(FIG,'Position',[0 0 600 400])



%%
