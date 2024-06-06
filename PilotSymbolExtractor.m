
%---------------------------------------------------------------------------
%  OFDM Demodulator and GI Extractor for Pilot Symbols of Kim and Stüber
%  NFFT: FFT length
%  chnr: number of subcarrier
%  G: guard length
%  N_P: channel impulse response length
%---------------------------------------------------------------------------

function [y] = PilotSymbolExtractor(data,chnr,NFFT,G,D_f);


x_remove_guard_interval = [data(G+1:NFFT+G)]; % remover the guard interval

x_remove_first_part_of_cir = [x_remove_guard_interval(NFFT/D_f+1:NFFT)]; % remover the guard interval

x_t = [x_remove_first_part_of_cir, x_remove_first_part_of_cir];
x_f = fft(x_t);

y = x_f(1:chnr); %Zero removing



	
