%Initialize
clearvars
clc
close all

rng default  % Use default random number generator


% Basic Parameters
carrierFreq  = 433e6;    %Carrier Frequency

% Symbol Parameters
secondsPerFrame = 1e-3;                         % Seconds per frame
M               = 16;                           % Modulation order
bitsPerSym      = log2(M);                      % Number of bits per symbol
symPerFrame     = 2000;                         % Number of symbols per frame
symbolRate      = symPerFrame/secondsPerFrame;  % Symbol rate per second
bitrate         = bitsPerSym*symbolRate;        % Bit rate per second


bandwidth = bitrate;      % Maximum Bandwidth

%For debug purposes: Output Bandwidth
disp(['Bandwidth = ' num2str(bandwidth/1e6) ' MHz'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USF=1;       %SW: USF=30,SF1=1,SF2=1; SW oversampling: USF=30,SF1=1,SF2=30;    %  
% SF1=36; SF2=1; % R: USF=1,SF1=30,SF2=1; Combine: USF=1,SF1=30,SF2=1;           %
%                                                                                %
% BF=64;                                                                         %
% sps = SF1*BF*USF*SF2;     % Number of samples per symbol (oversampling factor) %
% snr=10; % dB                                                                   % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

overlap=1.65;  %Sub-band overlap (dB) 

USF=1;   
SF1=36;         % Number of sub-bands
SF2=1;          % R: USF=1,SF1=36,SF2=1; Combine: USF=1,SF1=36,SF2=1; 
BF=10;          % SW oversampling: USF=36,SF1=1,SF2=36;
sampPerSym=SF1*BF*SF2;   
sampleRate=sampPerSym*bandwidth;  %Set Sample Rate for transmission
bw_buffer=20e6;


%% Generate Datastream

inputDatastream  = randi([0 1],symPerFrame*bitsPerSym,1);    % Generate vector of random binary bits as our datastream
inputDatastream  = bit2int(inputDatastream,bitsPerSym);      % Convert binary to integers
inputDataSymbols = qammod(inputDatastream,M);                % Convert data ints to symbols

time = linspace(0,secondsPerFrame,length(inputDataSymbols)); % Time Domain t-axis


%% Modulate Datastream

% Raised Cosine Filter
RolloffFactor = 0.4;
TX_FILTER = comm.RaisedCosineTransmitFilter(...
            'RolloffFactor',RolloffFactor, ...
            'FilterSpanInSymbols',6, ...
            'OutputSamplesPerSymbol',sampPerSym, ...
            'Gain',sqrt(sampPerSym));

tx = TX_FILTER(inputDataSymbols); %Complex Baseband RCF Filtered Transmission Sample Stream

tx = modulate(real(tx),carrierFreq,sampleRate,'qam',imag(tx)); % Real-Valued Passband Transmission Stream

% Plot the Transmission Spectrum for debug purposes.
% Should see spikes at the carrier frequency

fig1 = figure(1);
    fig1.Name='PassbandSpectrum';
    plotSpectrum(tx,sampleRate,'c',0);
    title('QAM Passband Spectrum');

%% Things and stuff

R = 1;
switch R
    case 1
        %% Import Gain Data From HFSS
        %  HFSS gain .csv file contains rows of data in two columns
        %  Column 1 is the frequency in MHz
        %  Column 2 is the gain at the corresponding frequency

        antennaGainData = readmatrix('Antenna Params plot 1.csv'); % Gain data from HFSS
        dataFreq        = antennaGainData(:,1)*1e6;     % Frequency points of imported data
        gain            = antennaGainData(:,2);         % Imported Gain Data

        freq = linspace(0,sampleRate/2,length(tx)/2);   % Frequency Domain f-axis

        gainExtrapolated = interp1(dataFreq,gain,freq,'linear','extrap').'; %Extrapolate Gain Data for frequencies outside of range

        fig2 = figure(2);
            gcf.Name='Antenna Gain';
            plot(dataFreq/1e6,gain,'b','linewidth',2);  % Plot imported data
            hold on
            plot(freq/1e6,gainExtrapolated,'--r');      % Plot extrapolated data
            grid on
            xlim([carrierFreq-bandwidth*2 carrierFreq+bandwidth*2]/1e6)
            title('Antenna Gain');
            legend('Imported Gain Data','Extrapolated Gain Data','location','southeast')


        %% Import Phase Data from HFSS
        %  HFSS phase .csv file contains rows of data in four columns
        %  Column 1 is the In-Phase phase (Usually 0)
        %  Column 2 is the Quadrature phase (Usually -90)
        %  Column 3 is frequency in MHz (should be same as in gain file)
        %  Column 4 is the phase data in degrees

        antennaPhaseData    = readmatrix('rE Plot 1.csv');  % Import phase data
        dataFreq            = antennaPhaseData(:,3)*1e6;    % Frequency points of imported data
        phase               = antennaPhaseData(:,4);        % Imported Phase Data

        phaseExtrapolated = interp1(dataFreq,phase,freq,'pchip','extrap').'; %Interpolate and Extrapolate Data
        
        % Extrapolation of phase data generally produces bad results
        % Realistically values should asymptotically approach an offset.
        % Instead pchip produces a cubic-looking function that goes to +/-
        % infinity
        %
        % We can fix this by replacing everything beyond the local min and
        % max of the cubic with that value
        %
        % Comment out the next 4 lines to see what happens if you don't do
        % this

        idx1=find(islocalmax(phaseExtrapolated), 1 );       %Find the local max
        idx2=find(islocalmin(phaseExtrapolated), 1,'last'); %Find the local min
        phaseExtrapolated(1:idx1-1)=phaseExtrapolated(idx1);    %Replace values beyond local max with the local max value
        phaseExtrapolated(idx2+1:end)=phaseExtrapolated(idx2);  %Replace values beyond local min with the local min value

        fig3 = figure(3);
            fig3.Name='Antenna Phase Offset';
            plot(dataFreq/1e6,phase,'b','linewidth',2);  % Plot imported data
            hold on
            plot(freq/1e6,phaseExtrapolated,'--r');      % Plot extrapolated data
            grid on
            xlim([carrierFreq-bandwidth*2 carrierFreq+bandwidth*2]/1e6)
            title('Antenna Phase Offset');
            legend('Imported Phase Data','Extrapolated Phase Data','location','northeast')

        %% Combine Gain and Phase into an Antenna Response to get transmitted signal

        % Combine Gain and Phase into an Antenna Response
        antennaResponse=10.^(gainExtrapolated/10).*exp(1j*phaseExtrapolated*pi/180); 

        %Real-valued signals in the f-domain will be mirrored along the y-axis
        antennaResponse=vertcat(flipud(antennaResponse),antennaResponse);

        % Convolve tx with antennaResponse to get the output after the ESA
        txESA=ifft(ifftshift(fftshift(fft(tx)).*antennaResponse));

        fig4 = figure(4);
            fig4.Name='Tx Spectrum';
            plotSpectrum(tx,sampleRate,'c',0)
            hold on
            plotSpectrum(txESA,sampleRate,'b',0)
            xlim([carrierFreq-bw_buffer carrierFreq+bw_buffer]/1e6)
            title('Tx Spectrum');
            legend('Tx signal','After filter/ESA','location','southeast')
end