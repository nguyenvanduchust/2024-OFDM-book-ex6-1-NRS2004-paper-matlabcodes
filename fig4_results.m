clear;

load mse_conventional_method.am  -ascii;

c_x = mse_conventional_method(1,:);

c_y = mse_conventional_method(2,:);



load  mse_proposed_method.am -ascii;

p_x = mse_proposed_method(1,:);

p_y = mse_proposed_method(2,:);




plot(c_x, 10*log10(c_y),'k*-');
hold on;
plot(p_x, 10*log10(p_y),'ko-');

hold off;

legend('\fontsize{12}LS channel estimator', 'Proposed channel estimator')
xlabel('SNR in dB', 'FontSize',12);
ylabel('MSE in dB', 'FontSize',12);

