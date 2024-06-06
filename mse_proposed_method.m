%==========================================================================
% "Least Square Channel Estimation Using Special Training Sequences for MIMO
% OFDM Systems in the Presence of Intersymbol Interference"
% Van Duc Nguyen, 15.06.2004, Agder University College, Norway
% Results of MSE of estimated channel
%==========================================================================

clear all;
%-------------------------------------------------
% Parameters for OFDM system
%-------------------------------------------------

NFFT = 64;              % FFT length
G = 0;                  % Guard interval length

M_ary =4;               % Multilevel of M-ary symbol
P_A = sqrt(2);          % Amplitude of pilot symbol
 
D_f = 2;                % Pilot distance in frequency domain
D_t = 4;                % pilot distance in time domain
NofZeros = D_f-1;
M = NFFT / (D_f);       % number of pilot symbol  per OFDM symbol
t_a = 50*10^(-9);       % Sampling duration of HiperLAN/2

%-------------------------------------------------
% Parameters for Monte Carlo channel
%-------------------------------------------------


N_P  = 4;
symbol_duration = NFFT * t_a;   %OFDM symbol duration
number_of_summations = 40;      % Number of summations for Monte-Carlo method



f_dmax = 50.0;                  % Maximum Doppler frequency

load h_decimation.am -ascii;

h11_initial = h_decimation;

h12_initial = h_decimation;


h21_initial = h_decimation;



h22_initial = h_decimation;
                                

N_P  = length(h_decimation);



 


%--------------------------------------
% Preparing pilot pattern for Antenna 1
%--------------------------------------

PP_A1 = []; 


for m = 0:M-1; 
    
    PP_A1 = [PP_A1,P_A*exp(j*D_f*pi*(m)^2/NFFT)];
    %randm = rand(1,1); 
    %PP_A1 = [PP_A1,P_A*exp(j*2*pi*randm)];
    for l = 1:D_f -1;
    PP_A1=[PP_A1,zeros(1,NofZeros)]; 
    end;
    
end;


%------------------------
% Preparing pilot pattern for Antenna 2
%------------------------

PP_A2 = [];

for m = 0:M-1;
  
    PP_A2 = [PP_A2,P_A*exp(j*D_f*pi*(m+M/2)^2/NFFT)];
    %randm = rand(1,1); 
    %PP_A2 = [PP_A2,P_A*exp(j*2*pi*randm)];

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

R = Q * PP';



NofOFDMSymbol = 1000;    % Number of OFDM symbole to be processed
No_Of_OFDM_Data_Symbol = NofOFDMSymbol-ceil(NofOFDMSymbol/D_t);
                        % Number of data symbol
length_data = NofOFDMSymbol * NFFT;  
                        % The total data length
                        
                        
Number_Relz = 100;
mse_relz = [];
for number_of_relialization= 1: Number_Relz;         
    
u1 = rand(N_P,number_of_summations); % A random variable
u2 = rand(N_P,number_of_summations); % A random variable

%--------------------------
% Generation of source bits
%--------------------------
source_data1 = randint(length_data,2);
source_data2 = randint(length_data,2);

%--------------------
% bit to symbol coder
%--------------------

symbols1 = bi2de(source_data1);  
symbols2 = bi2de(source_data2);  
%----------------------------
% QPSK modulator in base band
%----------------------------
QASK_Symbol1 = dmodce(symbols1,1,1,'qask',M_ary);
QASK_Symbol2 = dmodce(symbols2,1,1,'qask',M_ary);


%----------------------------
% Preparation of data pattern
%----------------------------

Data_Pattern1 = []; % Transmitted Signal before IFFT from the first antenna
m = 0;
for i=0:No_Of_OFDM_Data_Symbol-1;
    QASK_tem = [];
    for n=1:NFFT;
          QASK_tem = [QASK_tem,QASK_Symbol1(i*NFFT+n)];
    end;
    Data_Pattern1 = [Data_Pattern1;QASK_tem];
    
    clear QASK_tem;

end;

Data_Pattern2 = []; % Transmitted Signal before IFFT from the second antenna
m = 0;
for i=0:No_Of_OFDM_Data_Symbol-1;
    QASK_tem = [];
    for n=1:NFFT;
          QASK_tem = [QASK_tem,QASK_Symbol2(i*NFFT+n)];
    end;
    Data_Pattern2 = [Data_Pattern2;QASK_tem];
    clear QASK_tem;
  
end;



%******************************************************************
%
% Transmitted Signal of Antenna 1
%------------------------------------------------------------------

TS1_BeforeIFFT = Insert_PilotSymbol(PP_A1,Data_Pattern1,D_t,NofOFDMSymbol,NFFT);


%******************************************************************
%
% Transmitted Signal of Antenna 2
%------------------------------------------------------------------

TS2_BeforeIFFT = Insert_PilotSymbol(PP_A2,Data_Pattern2,D_t,NofOFDMSymbol,NFFT);




mse = []

snr_min = 0;
snr_max = 50;
step = 5;

for snr = snr_min:step:snr_max;




initial_time=0;                 % Initial time



    
    

%--------------------------------------------------------------------------
% Transmitted signal of trasmitt antena 1 to receive antenna 1 (channel: h11)
%--------------------------------------------------------------------------


rs1_t_frame = [];
rs2_t_frame = [];
h11_frame = [];
h21_frame = [];


for i=0:NofOFDMSymbol-1;
   OFDM_signal_tem1 = OFDM_Modulator(TS1_BeforeIFFT(i+1,:),NFFT,G);
   % OFDM signal from the first antenna is created
   [h11, t] = MCM_channel_model(u1, initial_time, number_of_summations, symbol_duration, ...,
       f_dmax, h11_initial);
   h11_frame = [h11_frame; h11];
   
   rs1_t = conv(OFDM_signal_tem1, h11);
   % The received signal over multhipath channel is created
   rs1_t = awgn(rs1_t,snr,'measured','dB');
   %rs1_t = awgn(rs1_t,snr);
   % The received signal over multhipath channel with additive noise is created
   rs1_t_frame = [rs1_t_frame; rs1_t];
   clear OFDM_signal_tem1;
   

%------------------------------------------------------------------------
% Transmitted signal of antenna 2 to receive antenna 1 (channel: h21)
%-----------------------------------------------------------------------


   OFDM_signal_tem2 = OFDM_Modulator(TS2_BeforeIFFT(i+1,:),NFFT,G);
   % OFDM signal from the second antenna is created
   
   [h21, t] = MCM_channel_model(u2, initial_time, number_of_summations, symbol_duration, ...,
       f_dmax, h21_initial);
   h21_frame = [h21_frame; h21];
   initial_time = t;
   
   rs2_t = conv(OFDM_signal_tem2, h21);
   % The received signal over multhipath channel is created
   rs2_t = awgn(rs2_t,snr,'measured','dB');
   %rs2_t = awgn(rs2_t,snr);
   % The received signal over multhipath channel with additive noise is created
   rs2_t_frame = [rs2_t_frame; rs2_t];
   clear OFDM_signal_tem2;
end;


%------------------------------------------------
% Recever 1: OFDM demodulator, channel estimation
%------------------------------------------------
Estimated_CTF = [];
mse_v = [];
estimated_cir = [];
estimated_h11 = [];
estimated_h21 = [];

Received_PP= [];    % Prepare a matrix for reveived pilot symbols
Receiver_Data = []; %  Prepare a matrix for reveived data symbols

rs = [];
rs_f_frame = [];
for i=1:NofOFDMSymbol;
   if (N_P > G+1) & (i>1)
          % If it is not the first symbol and the length of CIR is longer than
          % the gaurd interval length, then the ISI term must be taken into 
          % account
          previous_symbol1 = rs1_t_frame(i-1,:);
                           % previous OFDM symbol of the first antenna
          previous_symbol2 = rs2_t_frame(i-1,:);
                           % previous OFDM symbol of the second antenna 
                           
          %----------------------------------------------------
          % Extract the ISI term from the previous symbol 
          %----------------------------------------------------
          
          ISI_term1 = previous_symbol1(NFFT+2*G+1:NFFT+G+N_P-1); 
                           % The position from NFFT+2G+1: NFFT+G+N_P-1 is ISI term
          ISI_extended1 = [ISI_term1,zeros(1,length(previous_symbol1)-length(ISI_term1))];                  
          
          
          ISI_term2 = previous_symbol2(NFFT+2*G+1:NFFT+G+N_P-1); 
          
          ISI_extended2 = [ISI_term2,zeros(1,length(previous_symbol2)-length(ISI_term2))];  
          
          rs_t = rs1_t_frame(i,:) + rs2_t_frame(i,:) + ISI_extended1 +  ISI_extended2;          
                           % The ISI term is added to the current OFDM symbol
          
 
          
      else
          
          rs_t = rs1_t_frame(i,:) + rs2_t_frame(i,:);
          
   end;  
   %--------------------------------
   % Extract the pilot symbols
   %--------------------------------
   
   if (mod(i-1,D_t)==0)
        Demodulated_Pilot = PilotSymbolExtractor(rs_t,NFFT,NFFT,G,D_f); 
        %Demodulated_Pilot = OFDM_Demodulator(rs_t,NFFT,NFFT,G);

        
        Demodulated_P = [];
        for i = 1:NFFT;
            Demodulated_P = [Demodulated_P; Demodulated_Pilot(i)];
        end;
      
       
        
        estimated_cir_i_tem = R * Demodulated_P;
        estimated_cir_i = [];
        for i = 1:length(estimated_cir_i_tem);
            estimated_cir_i = [estimated_cir_i, estimated_cir_i_tem(i)];
        end;
        
        estimated_h11_i = [estimated_cir_i(1:N_P)];
        estimated_h21_i = [estimated_cir_i(N_P+1:length(estimated_cir_i))];
           
        estimated_h11 = [estimated_h11;  estimated_h11_i];
        estimated_h21 = [estimated_h21;  estimated_h21_i];


   else
        Demodulated_signal = OFDM_Demodulator(rs_t,NFFT,NFFT,G);  
        % OFDM demodulator
        Receiver_Data = [Receiver_Data; Demodulated_signal];
        
   end;

end;
 mse_v1 = sum((abs(h11_frame(1:D_t:NofOFDMSymbol,:)-estimated_h11).^2)')/N_P;
 
 mse_awgn1 = sum(mse_v1)/length(mse_v1);
 
 mse_v2 = sum((abs(h21_frame(1:D_t:NofOFDMSymbol,:)-estimated_h21).^2)');
 
 mse_awgn2 = sum(mse_v2)/length(mse_v2);
 mse_awgn = (mse_awgn1 + mse_awgn2)/2;
 
mse = [mse, mse_awgn];
end;

mse_relz = [mse_relz;mse];
end;

mse = sum(mse_relz)/Number_Relz;

snr = snr_min:step:snr_max;

plot(snr,10*log10(mse),'bo');
data = [snr; mse];
save mse_proposed_method.am data -ascii;

