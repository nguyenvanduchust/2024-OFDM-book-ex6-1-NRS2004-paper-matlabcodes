%==========================================================================
% "Least Square Channel Estimation Using Special Training Sequences for MIMO
% OFDM Systems in the Presence of Intersymbol Interference"
% Van Duc Nguyen, 15.06.2004, Agder University College, Norway
% Results of SER
%==========================================================================

%clc;
clear all;

NFFT = 64;               % FFT length
G = 0;                   % Guard interval length

M_ary =4;                % Multilevel of M-ary symbol
P_A = sqrt(2);           % Amplitude of pilot symbol
 
D_f = 2;                 % Pilot distance in frequency domain
D_t = 4;                 % pilot distance in time domain
NofZeros = D_f-1;
M = NFFT / (D_f);        % Number of pilot symbol  per OFDM symbol

t_a = 50*10^(-9);       % Sampling duration of HiperLAN/2

%-------------------------------------------------
% Parameters for Monte Carlo channel
%-------------------------------------------------



symbol_duration = NFFT * t_a;   %OFDM symbol duration
number_of_summations = 40;      % Number of summations for Monte-Carlo method




f_dmax = 50.0;                  % Maximum Doppler frequency

load h_decimation.am -ascii;
h11_initial = h_decimation;

h12_initial = h_decimation;


h21_initial = h_decimation;



h22_initial = h_decimation;
                                

N_P  = length(h_decimation);
 





NofOFDMSymbol = 1000;        %Number of data and pilot OFDM symbol





No_Of_OFDM_Data_Symbol = NofOFDMSymbol-ceil(NofOFDMSymbol/D_t);
                            %Number of datay symbols
length_data = (No_Of_OFDM_Data_Symbol) * NFFT;  
                            % The total data length
                            
                            
                        
Number_Relz = 50;
ser_relz = [];
for number_of_relialization= 1: Number_Relz;   


u11 = rand(N_P,number_of_summations); % A random variable
u12 = rand(N_P,number_of_summations); % A random variable
u21 = rand(N_P,number_of_summations); % A random variable
u22 = rand(N_P,number_of_summations); % A random variable
    

%-------------
% Source bites
%-------------
source_data1 = randint(length_data,2);
source_data2 = randint(length_data,2);
%--------------------
% bit to symbol coder
%--------------------

symbols1 = bi2de(source_data1);  
symbols2 = bi2de(source_data2);  
%----------------------------
% QPSK modulator in base band
%----------------------------------------------
QASK_Symbol1 = dmodce(symbols1,1,1,'qask',M_ary);
QASK_Symbol2 = dmodce(symbols2,1,1,'qask',M_ary);

%-----------------------------------------------
% Preparing data pattern
%
%-----------------------------------------------

Data_Pattern1 = []; % Transmitted Signal before IFFT
m = 0;
for i=0:No_Of_OFDM_Data_Symbol-1;
    QASK_tem = [];
    for n=1:NFFT;
          QASK_tem = [QASK_tem,QASK_Symbol1(i*NFFT+n)];
    end;
    Data_Pattern1 = [Data_Pattern1;QASK_tem];
    
    clear QASK_tem;

end;

Data_Pattern2 = []; % Transmitted Signal before IFFT
m = 0;
for i=0:No_Of_OFDM_Data_Symbol-1;
    QASK_tem = [];
    for n=1:NFFT;
          QASK_tem = [QASK_tem,QASK_Symbol2(i*NFFT+n)];
    end;
    Data_Pattern2 = [Data_Pattern2;QASK_tem];
    
    clear QASK_tem;

end;


%------------------------
% Preparing pilot pattern for Antenna 1
%------------------------

PP_A1 = []; 


for m = 0:M-1; 
    
    PP_A1 = [PP_A1,P_A*exp(j*D_f*pi*(m)^2/NFFT)];  
    for l = 1:D_f -1;
    PP_A1=[PP_A1,zeros(1,NofZeros)]; 
    end;
    
end;


%------------------------
% preparing pilot pattern for Antenna 2
%------------------------

PP_A2 = [];

for m = 0:M-1;
  
    PP_A2 = [PP_A2,P_A*exp(j*D_f*pi*(m+M/2)^2/NFFT)];

    for l = 1:D_f -1;
    PP_A2=[PP_A2,zeros(1,NofZeros)]; 
    end;
   
   
end;

%----------------------------------------
% FFT matrix
%----------------------------------------

F = [];

for k=0:NFFT-1
    W_tem = [];
    for n = 0:NFFT-1;
        W_tem = [W_tem,exp(-j*2*pi*n*k/NFFT)];
    end;
    F = [F;W_tem];
end;

%---------------------------------------
% Least square filter coefficients
%---------------------------------------


PP = [diag(PP_A1)*F(:,1:N_P),diag(PP_A2)*F(:,1:N_P)];
Q = inv(PP'*PP);

%******************************************************************
% Transmitted Signal of Antenna 1 (Insert pilot sequence into data
% sequence)
%------------------------------------------------------------------

TS1_BeforeIFFT = Insert_PilotSymbol(PP_A1,Data_Pattern1,D_t,NofOFDMSymbol,NFFT);


%******************************************************************
% Transmitted Signal of Antenna 2 (Insert pilot sequence into data
% sequence)
%------------------------------------------------------------------

TS2_BeforeIFFT = Insert_PilotSymbol(PP_A2,Data_Pattern2,D_t,NofOFDMSymbol,NFFT);




%------------------------------------------------------------------------
% Transmitted signal of trasmitt antena 1 to receive antenna 1 (channel: h11)
% and transmitted signal of trasmitt antena 1 to receive antenna 2 (channel: h12)
%-----------------------------------------------------------------------
ser_without_isic = [];

snr_min =0;
snr_max =50;
step = 5;
for snr = snr_min:step:snr_max; 
    



rs11_frame = [];
rs12_frame = [];
initial_time=0;                 % Initial time

for i=0:NofOFDMSymbol-1;
   OFDM_signal_tem = OFDM_Modulator(TS1_BeforeIFFT(i+1,:),NFFT,G);
   % OFDM signal from the first transmitt antenna is created
   
   [h11, t] = MCM_channel_model(u11, initial_time, number_of_summations, symbol_duration, ...,
       f_dmax, h11_initial);

   [h12, t] = MCM_channel_model(u12, initial_time, number_of_summations, symbol_duration, ...,
       f_dmax, h12_initial);
   

   
   initial_time = t;
   
   rs11 = conv(OFDM_signal_tem, h11);
   rs12 = conv(OFDM_signal_tem, h12);
   % The received signal over multhipath channel is created
   rs11 = awgn(rs11,snr,'measured','dB');
   rs12 = awgn(rs12,snr,'measured','dB');
   
   rs11_frame = [rs11_frame; rs11];
   rs12_frame = [rs12_frame; rs12];
   clear OFDM_signal_tem;
end;


%------------------------------------------------------------------------
% Transmitted signal of antenna 2 to receive antenna 1 (channel: h21)
% and transmitted signal of antenna 2 to receive antenna 2 (channel: h22)
%-----------------------------------------------------------------------


rs21_frame = [];
rs22_frame = [];

initial_time=0;                 % Initial time
for i=0:NofOFDMSymbol-1;
   OFDM_signal_tem = OFDM_Modulator(TS2_BeforeIFFT(i+1,:),NFFT,G);
   % OFDM signal from the second antenna is created
   
   [h21, t] = MCM_channel_model(u21, initial_time, number_of_summations, symbol_duration, ...,
       f_dmax, h21_initial);
   
   [h22, t] = MCM_channel_model(u22, initial_time, number_of_summations, symbol_duration, ...,
       f_dmax, h22_initial);
   

   
   initial_time = t;
   
   rs21 = conv(OFDM_signal_tem, h21);
   rs22 = conv(OFDM_signal_tem, h22);
   % The received signal over multhipath channel is created
   rs21 = awgn(rs21,snr,'measured','dB');
   rs22 = awgn(rs22,snr,'measured','dB');
   rs21_frame = [rs21_frame; rs21];
   rs22_frame = [rs22_frame; rs22];
   clear OFDM_signal_tem;
end;


%-----------------
% Recever 1: OFDM demodulator, channel estimation
%-----------------

estimated_h11_21 = [];
Received_PP= []; % Prepare a matrix for reveived pilot symbols
Receiver_Data = []; %  Prepare a matrix for reveived data symbols

SignalPostFFT1 = [];
SignalPostFFT2 = [];

d1 = []; % Received signal of the first receive antenna
d2 = []; % Received signal of the second receive antenna
data_symbol_1 = [];
data_symbol_2 = [];

for i=1:NofOFDMSymbol;
   if (N_P > G+1) & (i>1)
          % if it is not the first symbol and the length of CIR is longer than
          % the gaurd interval length, then the ISI term must be taken into 
          % account
          previous_symbol11 = rs11_frame(i-1,:);
                           % previous OFDM symbol of the first transmitt
                           % antenna to the first receive antenna
          previous_symbol21 = rs21_frame(i-1,:);
                           % previous OFDM symbol of the second transmitt antenna
                           % to the first receive antenna
          previous_symbol12 = rs12_frame(i-1,:);
                           % previous OFDM symbol of the first transmitt
                           % antenna to the second receive antenna
          previous_symbol22 = rs22_frame(i-1,:);
                           % previous OFDM symbol of the second transmitt antenna
                           % to the second receive antenna                 
                           
                           
          ISI_term11 = previous_symbol11(NFFT+2*G+1:NFFT+G+N_P-1); 
                           % the position from NFFT+2G+1: NFFT+G+N_P-1 is ISI term
          ISI_11 = [ISI_term11,zeros(1,length(previous_symbol11)-length(ISI_term11))];                  
          
          
          ISI_term21 = previous_symbol21(NFFT+2*G+1:NFFT+G+N_P-1); 
          
          ISI_21 = [ISI_term21,zeros(1,length(previous_symbol21)-length(ISI_term21))];  
          
          
          ISI_term12 = previous_symbol12(NFFT+2*G+1:NFFT+G+N_P-1); 
                          % the position from NFFT+2G+1: NFFT+G+N_P-1 is ISI term
          ISI_12 = [ISI_term12,zeros(1,length(previous_symbol12)-length(ISI_term12))];                  
          
          
          ISI_term22 = previous_symbol22(NFFT+2*G+1:NFFT+G+N_P-1); 
          
          ISI_22 = [ISI_term22,zeros(1,length(previous_symbol22)-length(ISI_term22))];
          
          %rs1_i = rs11_frame(i,:) + rs21_frame(i,:);  
                          
          %rs2_i = rs12_frame(i,:) + rs22_frame(i,:);
          
          rs1_i = rs11_frame(i,:) + rs21_frame(i,:) + ISI_11 +  ISI_21;  
          rs2_i = rs12_frame(i,:) + rs22_frame(i,:) + ISI_12 +  ISI_22;
                        % the ISI term is added to the current OFDM symbol
      else
          
          rs1_i = rs11_frame(i,:) + rs21_frame(i,:);
          rs2_i = rs12_frame(i,:) + rs22_frame(i,:);

   end;  
   
   
   if (mod(i-1,D_t)==0)
        %Demodulated_Pilot1 = PilotSymbolExtractor(rs1_i,NFFT,NFFT,G,D_f); 
        Demodulated_Pilot1 = OFDM_Demodulator(rs1_i,NFFT,NFFT,G);
        
        SignalPostFFT1 = [SignalPostFFT1; Demodulated_Pilot1];
        %Demodulated_Pilot2 = PilotSymbolExtractor(rs2_i,NFFT,NFFT,G,D_f);
        Demodulated_Pilot2 = OFDM_Demodulator(rs2_i,NFFT,NFFT,G);
        
        SignalPostFFT2 = [SignalPostFFT2; Demodulated_Pilot2];
        Demodulated_P1 = [];
        Demodulated_P2 = [];
        for i = 1:NFFT;
            Demodulated_P1 = [Demodulated_P1; Demodulated_Pilot1(i)];
            Demodulated_P2 = [Demodulated_P2; Demodulated_Pilot2(i)];
        end;
      
        
        estimated_h11_21_i = Q * PP' * Demodulated_P1;

        
        he11_21_i = estimated_h11_21_i;
        
        he11_i = he11_21_i(1:length(estimated_h11_21_i)/2);
        H11_i = fft([he11_i;zeros(NFFT-N_P,1)]);
        he21_i = he11_21_i(length(he11_21_i)/2+1:length(he11_21_i));
        H21_i = fft([he21_i;zeros(NFFT-N_P,1)]);
        
        estimated_h12_22_i = Q * PP' * Demodulated_P2;
        
        he12_22_i = estimated_h12_22_i;
        he12_i = he12_22_i(1:length(he12_22_i)/2);
        H12_i = fft([he12_i;zeros(NFFT-N_P,1)]);
        
        he22_i = he12_22_i(length(he12_22_i)/2+1:length(he12_22_i));
        H22_i = fft([he22_i;zeros(NFFT-N_P,1)]);        
        
        
        
   else
        %-----------------------------------------------------------
        % OFDM demodulator
        %-----------------------------------------------------------
        Demodulated_signal1_i = OFDM_Demodulator(rs1_i,NFFT,NFFT,G);  
        SignalPostFFT1 = [SignalPostFFT1; Demodulated_signal1_i];
        Demodulated_signal2_i = OFDM_Demodulator(rs2_i,NFFT,NFFT,G); 
        SignalPostFFT2 = [SignalPostFFT2; Demodulated_signal2_i];
       
        d1_i = [];
        d2_i = [];
        for k = 1:NFFT;
            H_k = [H11_i(k),H21_i(k); H12_i(k),H22_i(k)];
            y = [Demodulated_signal1_i(k);Demodulated_signal2_i(k)];
            x = inv(H_k) * y;
            d1_i = [d1_i,x(1)];
            d2_i = [d2_i,x(2)];
        end;
        
        %--------------------------
        % Demodulated signal 
        %--------------------------
        d1 = [d1; d1_i];
        
        d2 = [d2; d2_i];
        
        demodulated_symbol_1i = ddemodce(d1_i,1,1,'qask',M_ary);
        demodulated_symbol_2i = ddemodce(d2_i,1,1,'qask',M_ary);
        data_symbol_1 = [data_symbol_1, demodulated_symbol_1i];
        
        data_symbol_2 = [data_symbol_2, demodulated_symbol_2i];
        

        
    end;

end;



    data_symbol_1 = data_symbol_1';
    data_symbol_2 = data_symbol_2';
    [number1_without_isic, ratio1_without_isic] = symerr(symbols1,data_symbol_1);
    

    
    ser_without_isic = [ser_without_isic, ratio1_without_isic];

    
end;



ser_relz = [ser_relz;ser_without_isic];
end;

ser = sum(ser_relz)/Number_Relz;



snr = snr_min:step:snr_max;




semilogy(snr, ser,'b*');

ylabel('SER');
xlabel('SNR');


data = [snr; ser];
save ser_conventional_method.am data -ascii;
