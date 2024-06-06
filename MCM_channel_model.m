
%=============================================
% Monte Carlo method for a time-variant channel modelling 
% Van Duc Nguyen, Agder Unversity College
% 12.06.04, Grimstad, Norway
%=============================================
function [h,t_next] = MCM_channel_model(u, initial_time, number_of_summations, symbol_duration, f_dmax, channel_coefficients);


t = initial_time;



Channel_Length = length(channel_coefficients);





h_vector = [];


for k=1:Channel_Length;
    u_k = u(k,:); % A random variable
    phi = 2 * pi * u_k; % Phase coefficients are created
    f_d = f_dmax * sin(2*pi*u_k); % Doppler frequency after Monte Carlo method is created
    h_tem= channel_coefficients(k)* 1/(sqrt(number_of_summations)) * sum(exp(j*phi).*exp(j*2*pi*f_d*t));
    h_vector = [h_vector,  h_tem];
end;

h = h_vector;
t_next = initial_time + symbol_duration; %Coherent time for the next symbol

    