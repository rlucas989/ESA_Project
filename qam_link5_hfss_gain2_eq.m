clear all
%close all
%
fc=433e6;    %Carrier Frequency
%fc=426e6;
fps=1e-3;    % frame length in sec
M = 16;      % Modulation order
k = log2(M); % Number of bits per symbol
n = 2000;    % Number of symbols per frame
sr = n/fps;  % Symbol rate 
bw = sr*log2(M); % BW 
disp(['Bandwidth = ' num2str(bw/1e6) ' MHz'])
overlap=1.65;  %Sub-band overlap (dB) 
% USF=1;       %SW: USF=30,SF1=1,SF2=1; SW oversampling: USF=30,SF1=1,SF2=30;  
% SF1=36; SF2=1; % R: USF=1,SF1=30,SF2=1; Combine: USF=1,SF1=30,SF2=1; 
% BF=64;
% sps = SF1*BF*USF*SF2;     % Number of samples per symbol (oversampling factor)
% snr=10; % dB
USF=1;   
SF1=36;   %Number of sub-bands
SF2=1;    % R: USF=1,SF1=36,SF2=1; Combine: USF=1,SF1=36,SF2=1; 
BF=10;            % SW oversampling: USF=36,SF1=1,SF2=36;
sps=SF1*BF*SF2;   
fs=sps*bw;  
bw_buffer=20e6;
rng default  % Use default random number generator

dataIn = randi([0 1],n*k,1); % Generate vector of binary data
dataSymbolsIn = bit2int(dataIn,k);
dataMod = qammod(dataSymbolsIn,M); %
tt0=linspace(0,fps,length(dataMod));

% Raised Cosine Filter
RolloffFactor=0.4;
TX_FILT = comm.RaisedCosineTransmitFilter(...
        'RolloffFactor',RolloffFactor, ...
        'FilterSpanInSymbols',6, ...
        'OutputSamplesPerSymbol',sps, ...
        'Gain',sqrt(sps));
txFiltOut0 = TX_FILT(dataMod); 

tx_RF=modulate(real(txFiltOut0),fc,fs,'qam',imag(txFiltOut0));

figure(1);clf
        plotSpectrum(tx_RF,fs,'c',0)

%%
R=3;
switch R
    case 1
        %aa=readmatrix('RealizedGain_Resistive.csv');
        aa=readmatrix('Antenna Params plot 1.csv'); %Gain data from HFSS
        f1=aa(:,1)*1e6;
        gain=aa(:,2);

        freq0 = linspace(0,fs/2,length(tx_RF)/2);
        gain_extrapol=interp1(f1,gain,freq0,'linear','extrap').'; %Extrapolate Gain Data

        figure(25);clf
        plot(f1/1e6,gain,'b','linewidth',2);
        hold on
        plot(freq0/1e6,gain_extrapol,'--r');
        grid on
        xlim([fc-bw*2 fc+bw*2]/1e6)
%% 

        %%%%%%%%%%%
        bb=readmatrix('rE Plot 1.csv');         %Phase Data from HFSS   
        phase_extrapol=interp1(f1,bb(:,4),freq0,'linear','extrap').'; %Interpolate and Extrapolate Data
        idx1=min(find(freq0>380e6));
        idx2=min(find(freq0>480e6));

        phase_extrapol(1:idx1-1)=phase_extrapol(idx1);
        phase_extrapol(idx2+1:end)=phase_extrapol(idx2);
        %%%%%%

        %resp_mag0=10.^(gain_extrapol/10);
        resp_mag0=10.^(gain_extrapol/10).*exp(1j*phase_extrapol*pi/180);
        resp_mag1=vertcat(flipud(resp_mag0),resp_mag0);
        tx_RF_filt1=ifft(ifftshift(fftshift(fft(tx_RF)).*resp_mag1));
        tx_RF_filt=tx_RF_filt1;
        %tx_RF_filt=tx_RF;

        figure(1);clf
        plotSpectrum(tx_RF,fs,'c',0)
        hold on
        plotSpectrum(tx_RF_filt,fs,'b',0)
        xlim([fc-bw_buffer fc+bw_buffer]/1e6)
        legend('Tx signal','After filter/ESA','location','southeast')

%return
%%
    % case 2
    %     aa0=readmatrix('Antenna Params Table 1.csv');
    %     f1=aa0(:,1)*1e6;
    % 
    %     % aa=interpft(aa0(:,2:end),size(aa0,1)*4,1);
    %     % f1=linspace(f1(1),f1(end),size(aa0,1)*4).';
    %     % aa=horzcat(f1,aa);
    %     gain=aa0(:,9);
    % 
    %     freq0 = linspace(0,fs/2,length(tx_RF)/2);
    %     gain_extrapol=interp1(f1,gain,freq0,'linear','extrap').';
    % 
    %     figure(25);clf
    %     plot(f1/1e6,gain,'b','linewidth',2);
    %     hold on
    %     plot(freq0/1e6,gain_extrapol,'--r');
    %     grid on
    %     xlim([fc-bw*2 fc+bw*2]/1e6)
    % 
    %     resp_mag0=10.^(gain_extrapol/10);
    %     resp_mag1=vertcat(flipud(resp_mag0),resp_mag0);
    %     tx_RF_filt1=ifft(ifftshift(fftshift(fft(tx_RF)).*resp_mag1));
    %     tx_RF_filt=tx_RF_filt1;
    %     %tx_RF_filt=tx_RF;
    % 
    %     figure(1);clf
    %     plotSpectrum(tx_RF,fs,'c',0)
    %     hold on
    %     plotSpectrum(tx_RF_filt,fs,'b',0)
    %     xlim([fc-bw_buffer fc+bw_buffer]/1e6)
    %     legend('Tx signal','After filter/ESA','location','southeast')
%%
    case 3
        %aa=readmatrix('Antenna Params Table 1.csv');
        aa0=readmatrix('Antenna Params Plot 2.csv');
        f1=aa0(:,1)*1e6;

        aa=interpft(aa0(:,2:end),size(aa0,1)*4,1);
        f1=linspace(f1(1),f1(end),size(aa0,1)*4).';
        aa=horzcat(f1,aa);

        freq0 = linspace(0,fs/2/USF,length(tx_RF)/2/USF);

        gain0=fliplr(aa(:,2:45));
        for ii=1:size(gain0,2)
            [gainEnv(ii) idx]=max(gain0(:,ii));
            fpeak(ii)=f1(idx);
        end
        gainEnvi=interp1(fpeak,gainEnv,freq0,'linear','extrap').';
        % figure(24);clf
        % plot(freq0/1e6,gainEnvi)
        % %xlim([fc-bw*2 fc+bw*2]/1e6)
        % hold on;grid on
        % plot(f1/1e6,gain0);ylim([-3 5]);

        clear gain_extrapol
        figure(26);clf
        for ii=1:36%size(gain0,2)
            if ii==1
                gain_extrapol(:,1)=interp1(f1,gain0(:,1),freq0,'spline','extrap').';
            else
                gain_extrapol(:,ii)=circshift(gain_extrapol(:,ii-1),gain3dBIdx2-gain3dBIdx1);
                gain_extrapol(:,ii)=gain_extrapol(:,ii)-gainPeak+1*gainEnvi(gainPeakIdx+gain3dBIdx2-gain3dBIdx1);
            end
            [gainPeak gainPeakIdx]=max(gain_extrapol(:,ii));
            gain3dBIdx1=min(find(gain_extrapol(:,ii)>gainPeak-overlap));
            gain3dBIdx2=max(find(gain_extrapol(:,ii)>gainPeak-overlap));

            figure(26);
            plot(freq0/1e6,gain_extrapol(:,ii),'--r','linewidth',1);
            hold on
            plot(f1/1e6,gain0(:,1),'b','linewidth',2);
            grid on;
            ylim([-3 5]);%xlim([fc-2*bw fc+2*bw]/1e6)
        end
        figure(26);
        plot(freq0/1e6,gainEnvi)
        
        idx1=min(find(freq0>300e6));
        idx2=min(find(freq0>500e6));

        gain_extrapol(1:idx1,:)=-1000;
        gain_extrapol(idx2:end,:)=-1000;

        if USF==1
            tx_RF_filt=0;
            clear tx_RF_filt0; figure(2);clf
            for ii=1:size(gain_extrapol,2)
                resp_mag0=10.^(gain_extrapol(:,ii)/10);
                resp_mag1=vertcat(flipud(resp_mag0),resp_mag0);
                tx_RF_filt0(:,ii)=ifft(ifftshift(fftshift(fft(tx_RF)).*resp_mag1));
                
                noise = 500*sqrt(1)* randn(size(tx_RF))/sqrt(length(tx_RF));

                tx_RF_filt=tx_RF_filt+tx_RF_filt0(:,ii)+0*noise;
 
                figure(2);hold on
                plotSpectrum(tx_RF_filt0(:,ii),fs,[],0)
                xlim([fc-bw_buffer fc+bw_buffer]/1e6)
                ylim([-20 120])
            end
        else
            clear tx_RF_filt0;
            for ii=1:size(gain_extrapol,2)
                resp_mag0=10.^(gain_extrapol(:,ii)/10);
                resp_mag1=vertcat(flipud(resp_mag0),resp_mag0);
                sig=tx_RF(ii:USF:end);
                tx_RF_filt0(:,ii)=ifft(ifftshift(fftshift(fft(sig)).*resp_mag1));
            end

            tx_RF_filt00=0;
            kk=1;
            for ii=1:size(tx_RF_filt0,1)
                for jj=1:size(tx_RF_filt0,2)
                    tx_RF_filt00(kk+jj-1)=tx_RF_filt0(ii,jj);
                end
                kk=kk+size(tx_RF_filt0,2);
            end
            tx_RF_filt=tx_RF_filt00.';
        end

        figure(1);clf
        plotSpectrum(tx_RF,fs,'c',0)
        hold on
        plotSpectrum(tx_RF_filt,fs,'b',0)
        xlim([fc-bw_buffer fc+bw_buffer]/1e6)
        legend('Tx signal','After filter/ESA','location','southeast')
end
% % %return
%%
%EbNo = 10;
%snr = convertSNR(EbNo,'ebno', samplespersymbol=sps, bitspersymbol=k);
%rxNoisy = awgn(tx_RF_filt,snr,'measured');
noise = 50*sqrt(1)* randn(size(tx_RF_filt))/sqrt(length(tx_RF_filt));
rxNoisy = tx_RF_filt + 1*noise;

figure(1);
hold on
plotSpectrum(noise,fs,'m',0)
xlim([fc-bw_buffer fc+bw_buffer]/1e6)
legend('Tx signal','After filter/ESA','Noise floor','location','southeast')
ylim([10 130]) 


% idx1=min(find(gainEnvi>=0));
% idx2=max(find(gainEnvi>=0));
% gainEq=zeros(size(gainEnvi));
% gainEq(idx1:idx2)=-gainEnvi(idx1:idx2);
% 
% % figure(77)
% % plot(freq0/1e6,gainEq)
% % xlim([fc-bw fc+bw]/1e6)
% 
% tx_RF_filt=0;
% for ii=1:size(gain_extrapol,2)
%     resp_mag0=10.^(gainEq/10);
%     resp_mag1=vertcat(flipud(resp_mag0),resp_mag0);
%     rxNoisyEq=ifft(ifftshift(fftshift(fft(rxNoisy)).*resp_mag1));
% end
% 
% figure(111);
% plotSpectrum(tx_RF,fs,'b',0)
% hold on
% plotSpectrum(rxNoisyEq,fs,'c',0)
% xlim([fc-bw fc+bw]/1e6)
% ylim([10 130])


[rxI,rxQ]=demod(rxNoisy,fc,fs,'qam');   %Demodulation
rxDemod=rxI+1j*rxQ;

%Raised Cosine Filter Inverse
RX_FILT = comm.RaisedCosineReceiveFilter(...
    'RolloffFactor',RolloffFactor, ...
    'FilterSpanInSymbols',6, ...
    'InputSamplesPerSymbol',sps,'DecimationFactor',sps, ...
    'Gain',1/sqrt(sps));

rxFiltOut = RX_FILT(rxDemod); % Receive filter
rxFiltOut=rxFiltOut(M+1:end);%*0.31;

rx = rxFiltOut/0.062;
eq = comm.LinearEqualizer('Constellation',qammod(0:M-1,M)*1.*exp(1j*0*pi/180),...
    'ReferenceTap',3,'NumTaps',6);
numTrainingSymbols=25;
rxEqOut=eq(rxFiltOut,dataMod(1:numTrainingSymbols));
release(eq)

figure(40);clf
subplot(2,2,1)
scatter(real(dataMod),imag(dataMod),'b.');grid on
xlim([-4 4]);ylim([-4 4]);
subplot(2,2,2)
xlim([-4 4]);ylim([-4 4]);
scatter(real(rx),imag(rx),'c.');grid on
hold on
scatter(real(qammod(0:M-1,M)),imag(qammod(0:M-1,M)),'r+')
%xlim([-4 4]);ylim([-4 4]);
subplot(2,2,4)
scatter(real(rxEqOut),imag(rxEqOut),'g.');grid on
xlim([-4 4]);ylim([-4 4]);

%%
rxFiltOut_plot=rxFiltOut.*exp(1j*4*pi/180); %Rotate Constellation
%rxFiltOut_plot=rxEqOut;

%scale=3.5;
%scale=0.11;
%scale=0.061;

%scale=0.11;
scale=0.062;    %Scale Factor of Reference Constellation Points

evm = comm.EVM;
evm.ReferenceSignalSource = "Estimated from reference constellation";
evm.ReferenceConstellation = qammod([0:M-1],M)*scale;
% %evm.ReferenceConstellation = qammod([0:M-1],M)/USF/1.05.*exp(1j*42*pi/180);
evm.MaximumEVMOutputPort = false;
evm.AveragingDimensions = 2;
e = evm(rxFiltOut_plot);
B = 1/10*ones(10,1);
e_avg = filtfilt(B,1,e);
e_mean=mean(e_avg(1:end));
release(evm)

% figure(3);clf
% plot(e_avg)
% ylim([0 100]);grid on

figure(4);clf
scatter(real(rxFiltOut_plot),imag(rxFiltOut_plot),'c.')
hold on
scatter(real(evm.ReferenceConstellation),imag(evm.ReferenceConstellation),'r+')
grid on;axis equal;
title(['EVM = ' ,num2str(e_mean) ,'%'])
xlabel('In-phase');ylabel('Qaudrature')



