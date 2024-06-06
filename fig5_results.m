clear;
load ser_conventional_method.am -ascii;

ser_c_x = ser_conventional_method(1,:);
ser_c_y = ser_conventional_method(2,:);

load ser_proposed_method.am -ascii;

ser_p_x = ser_proposed_method(1,:);
ser_p_y = ser_proposed_method(2,:);



load ser_sufficientGL.am -ascii;
ser_sufficientGL_x = ser_sufficientGL(1,:);

ser_sufficientGL_y = ser_sufficientGL(2,:);


semilogy(ser_c_x, ser_c_y ,'ro-');
hold on;

semilogy(ser_p_x, ser_p_y ,'b*-');




semilogy(ser_sufficientGL_x,ser_sufficientGL_y,'kx-');

hold off;
xlabel('SNR in dB','FontSize',12);
ylabel('SER','FontSize',12);
legend('\fontsize{12}Conventional LS estimator, without GI','\fontsize{12}Proposed method, without GI', '\fontsize{12}Conventional LS estimator, sufficient GI length');