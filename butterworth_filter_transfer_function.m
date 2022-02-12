%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth Order Butterworth Filter Design
% Transfer Function Discretization
% Embedded Scientific Computing 4450:410
% The University of Akron
% Wilson Woods
% 12.8.2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;

% variable definitions
N = 10001;
omega = logspace(-4,pi,N);
z = exp(1i*omega);              % z-domain independent variable
fc = 880;                       % -3dB frequency 880 Hz         

% test values for delta t
dt = [1/2000 1/5000 1/8000 1/10000 1/20000 1/30000 1/40000 1/50000];
len = length(dt);
wc = 2*pi*fc;                   % -3dB radian frequency
k = 4;

found = 0;

% construct transfer function terms
p1 = wc*exp(1i*(5*pi/8));
p2 = wc*exp(1i*(7*pi/8));
p3 = wc*exp(1i*(9*pi/8));
p4 = wc*exp(1i*(11*pi/8));

% s-domain poles
a0 = (p1*p2*p3*p4);
a1 = (p1*p3*p4+p2*p3*p4+p1*p2*p3+p1*p2*p4);
a2 = (p3*p4+p1*p3+p1*p4+p2*p3+p2*p4+p1*p2);
a3 = (p1+p2+p3+p4);

for g = 1:1:len

    s = (log(z)/dt(g));
    
    % exact (s-domain) transfer function H(s)
    Hd_exact = (k*wc.^4)./(s.^4-s.^3*a3+s.^2*a2-s*a1+a0);

    c4 =  16- 8*dt(g)*a3 + 4*a2*dt(g).^2 - 2*dt(g).^3*a1 + dt(g).^4*a0;
    
    c3 = -64+16*a3*dt(g)-4*a1*dt(g).^3+4*a0*dt(g).^4;
    
    c2 = 96-8*a2*dt(g).^2+6*a0*dt(g).^4;
    
    c1 = -64 -16*a3*dt(g)+4*a1*dt(g).^3+4*a0*dt(g).^4;
    
    c0 = 16+8*a3*dt(g)+4*a2*dt(g).^2+2*a1*dt(g).^3+a0*dt(g).^4;
    
    % tustin (z-domain) approximation H(z)
    Hd_T_num = (k*wc.^4*dt(g).^4*(z.^4+4*z.^3+6*z.^2+4*z+1));
    
    Hd_T_denom = z.^4*c4 + z.^3*c3 + z.^2*c2 + z*c1 + c0;
    
    Hd_T = Hd_T_num./Hd_T_denom;
    
    % calculate tustin approximation error
    for a = 1:1:N
        error(a) = abs(Hd_T(a)-Hd_exact(a))/max(abs(Hd_exact));
    end
   
    % plot approximation using various sampling rates
    figure(1)
    subplot(2,1,1)
    title('Sampling Rate (dt) Selection');
    txt1 = ['Hd(Exact) dt = ', num2str(dt(g))];
    txt2 = ['Hd(Tustin) dt = ', num2str(dt(g))];
    txt3 = ['% Error dt = ', num2str(dt(g))];
    
    plot(omega, abs(Hd_exact), 'LineWidth' ,2, 'DisplayName', txt1)
    hold on
    plot(omega, abs(Hd_T), '--', 'DisplayName', txt2)
    xlabel('\Omega (rad/s)')
    ylabel('magnitude (dB)')
    legend('Location','northeastoutside')

    subplot(2,1,2)
    title('Tustin Approximation Error for Sampling Rates');
    if max(error) < 0.02
        plot (omega, error, 'DisplayName', txt3) 
        hold on
        xlabel('\Omega (rad/s)')
        ylabel('Normalized Error')
        legend('Location','southeastoutside')
    end

    % find sampling rate that yields < 1% error
    if max(error) < 0.01
        while found == 0

            % plot selected approximation
            figure (2)
            subplot(2,1,1)
            plot(omega, abs(Hd_exact), 'k', 'LineWidth' , 2 , 'DisplayName', txt1)
            hold on
            plot(omega, abs(Hd_T), '--y', 'LineWidth' , 2 ,'DisplayName', txt2)
            title('Exact and Tustin Comparison for Selected Sampling Rate');
            xlabel('\Omega (rad/s)')
            ylabel('magnitude (V)')
            legend('Location','northeastoutside')
    
            subplot(2,1,2)
            
            plot (omega, error, 'DisplayName', txt3) 
            title('Tustin Approximation Error for Selected Sampling Rate');
            xlabel('\Omega (rad/s)')
            ylabel('Normalized Error')
            legend('Location','southeastoutside')
            
            found = 1;
        end
    end
end

%% from above, sampling rate is chosen to be 1/30000

% determining Hc(iw) and Hd(exp(iwdt))
% Hc(iw):           continuous-time function of i * omega
% Hd(exp(iwdt)):    discrete-time function of i * omega * (1 / 30000)

s_c = 1i*wc;
dt_chosen = 1/30000;
z_wc = exp(1i*wc*dt_chosen);

% continuous-time transfer function
Hc_wc = abs((k*wc.^4)./(s_c.^4-s_c.^3*a3+s_c.^2*a2-s_c*a1+a0));

c4_wc =  16- 8*dt_chosen*a3 + 4*a2*dt_chosen.^2 - 2*dt_chosen.^3*a1 + dt_chosen.^4*a0;

c3_wc = -64+16*a3*dt_chosen-4*a1*dt_chosen.^3+4*a0*dt_chosen.^4;

c2_wc = 96-8*a2*dt_chosen.^2+6*a0*dt_chosen.^4;

c1_wc = -64 -16*a3*dt_chosen+4*a1*dt_chosen.^3+4*a0*dt_chosen.^4;

c0_wc = 16+8*a3*dt_chosen+4*a2*dt_chosen.^2+2*a1*dt_chosen.^3+a0*dt_chosen.^4;

Hd_wc_num = (k*wc.^4*dt_chosen.^4*(z_wc.^4+4*z_wc.^3+6*z_wc.^2+4*z_wc+1));

Hd_wc_denom = z_wc.^4*c4_wc + z_wc.^3*c3_wc + z_wc.^2*c2_wc + z_wc*c1_wc + c0_wc;

% discrete-time (tustin approximation) transfer function
Hd_wc = abs(Hd_wc_num./Hd_wc_denom);

% plots at chosen dt
s_dt = (log(z)/dt_chosen);

Hc_dt = (k*wc.^4)./(s_dt.^4-s_dt.^3*a3+s_dt.^2*a2-s_dt*a1+a0);

Hd_dt_num = (k*wc.^4*dt_chosen.^4*(z.^4+4*z.^3+6*z.^2+4*z+1));

Hd_dt_denom = z.^4*c4_wc + z.^3*c3_wc + z.^2*c2_wc + z*c1_wc + c0_wc;

Hd_dt = Hd_dt_num./Hd_dt_denom;


% calculate DC gain, required to equal 4
s_c2 = 1i*(wc*1.275);
z_wc2 = exp(1i*(wc*1.275)*dt_chosen);

Hc_wc2 = abs((k*wc.^4)./(s_c2.^4-s_c2.^3*a3+s_c2.^2*a2-s_c2*a1+a0));

Hd_wc2_num = (k*wc.^4*dt_chosen.^4*(z_wc2.^4+4*z_wc2.^3+6*z_wc2.^2+4*z_wc2+1));

Hd_wc2_denom = z_wc2.^4*c4_wc + z_wc2.^3*c3_wc + z_wc2.^2*c2_wc + z_wc2*c1_wc + c0_wc;

Hd_wc2 = abs(Hd_wc2_num./Hd_wc2_denom);

% Legend Text
txt4 = ['Hc(i*wc) = ' num2str(Hc_wc)];
txt5 = ['Hd(exp(i*wc*dt)) = ' num2str(Hd_wc)];
txt6 = ['Hc(s) @ dt = ', num2str(dt_chosen)];
txt7 = ['Hd(z) @ dt = ', num2str(dt_chosen)];
txt8 = ['-3dB of Hc(i*wc) = ' num2str(20*log10(Hc_wc))];
txt9 = ['-3dB of Hd(exp(i*wc*dt)) = ' num2str(20*log10(Hd_wc))];
txt10 = ['max dB of Hc(i*wc) = ' num2str(20*log10(max(abs(Hc_dt))))];
txt11 = ['max dB of Hd(exp(i*wc*dt)) = ' num2str(20*log10(max(abs(Hd_dt))))];
txt12 = ['Hc(i*wc) = ' num2str(Hc_wc2)];
txt13 = ['Hd(exp(i*wc*dt)) = ' num2str(Hd_wc2)];

% Figure
figure (3)
subplot(2,1,1)

plot(omega/dt_chosen, abs(Hc_dt), 'k', 'LineWidth' , 2 , 'DisplayName', txt6)
hold on
plot(omega/dt_chosen, abs(Hd_dt), '--y', 'LineWidth' , 2 ,'DisplayName', txt7)
plot(wc, Hc_wc, 'r*', 'MarkerSize', 8, 'DisplayName', txt4)
plot(wc, Hd_wc, 'go', 'MarkerSize', 8, 'DisplayName', txt5)
plot(wc*1.275, Hc_wc2, 'c*', 'MarkerSize', 8, 'DisplayName', txt12)
plot(wc*1.275, Hd_wc2, 'Color', [0.5 0 0.8] , 'marker' , 'o', 'MarkerSize', 8, 'DisplayName', txt13)
title('Exact and Tustin Transfer Functions Max and -3dB Points');
xlim([0 pi/3/dt_chosen])
ax = gca;
ax.XAxis.Exponent = 0;
xtickformat('%.0f')
xlabel('\omega (Hz)')
ylabel('magnitude (V)')
legend('Location','northeastoutside')

subplot(2,1,2)
semilogx(omega/dt_chosen, 20*log10(abs(Hc_dt)), 'k', 'LineWidth' , 2 , 'DisplayName', txt6)
hold on
semilogx(omega/dt_chosen, 20*log10(abs(Hd_dt)), '--y', 'LineWidth' , 2 ,'DisplayName', txt7)
semilogx(wc, 20*log10(Hc_wc), 'r*', 'MarkerSize', 8, 'DisplayName', txt8)
semilogx(wc, 20*log10(Hd_wc), 'go', 'MarkerSize', 8, 'DisplayName', txt9)
semilogx(3, 20*log10(max(abs(Hc_dt))), 'c*', 'MarkerSize', 8, 'DisplayName', txt10)
semilogx(3, 20*log10(max(abs(Hd_dt))), 'Color', [0.5 0 0.8], 'Marker', 'o', 'MarkerSize', 8, 'DisplayName', txt11)
ylim([-40 20])
xlabel('\omega (Hz)')
ylabel('magnitude (V)')
legend('Location','northeastoutside')
