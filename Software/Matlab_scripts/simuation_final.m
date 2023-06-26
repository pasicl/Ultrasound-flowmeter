%% Simulation of wave superposition and height change
clearvars
load('angles.mat');
f = 1E6; %freq in Hz
alfa = pi/3;
H_0 = 0.05;
H_c = 0.01;
attenuation_coefficient = 0.0022;
k_w = 0.001*1E-6*10^(attenuation_coefficient/10);%attenuation in water%m^-1*Hz^-1 
%k=100*1E-6; 
t = 0:1E-8:1E-4;
c_w = 1500;
Z1 = 10; Z2 = 1.5; Z3 = 10;
c_c = 2750;
MM=5000;
phi_M_1=zeros(1,MM);
phi_M_2=zeros(1,MM);

N = 20;

D = H_0/tan(alfa);

fi_vec_dns = zeros(N,1);
fi_vec_ups = zeros(N,1);

path_j = zeros(N,1);
s_mat_1 = zeros(N,length(t));s_mat_2 = s_mat_1;
angle_vec = angles(1:N);

v = 2e-6 * [0:MM-1];
H_v = H_0 + v/10;
att_vec = zeros(N,1);
Rpvc = ((Z2-Z3)/(Z2+Z3))^2;
Rg = ((Z2-Z1)/(Z2+Z1))^2;
for ii=1:MM
    
H = H_v(ii);
%H = H_0;
c_dns = c_w + v(ii);
c_ups = c_w - v(ii);
    for jj = 1:N
        
        path_j(jj) = H*((jj*2)-1)/sin(angle_vec(jj));
        
        fi_vec_dns(jj) = pi + path_j(jj)/c_dns*(2*pi*f);
        fi_vec_ups(jj) = pi + path_j(jj)/c_ups*(2*pi*f);
        att = Rpvc^(jj-1)*Rg^(jj-1)*exp(-k_w*f*path_j(jj));
        
        s_mat_1(jj,:) = cos(2*pi*f*t+fi_vec_dns(jj))*att;
        s_mat_2(jj,:) = cos(2*pi*f*t+fi_vec_ups(jj))*att;
        
        att_vec(jj) = att;
    end
    
s_sum_1 = sum(s_mat_1);
s_sum_2 = sum(s_mat_2);


ref_signal = cos(2*pi*f*t');
ref_signal_90 = cos(2*pi*f*t'+pi/2);

re_1 = s_sum_1*ref_signal;
im_1 = s_sum_1*ref_signal_90;

phi_i_1 =atan(im_1/re_1);
    if im_1 > 0 & re_1 < 0
        phi_i_1 = pi - abs(phi_i_1);
    elseif im_1 < 0 & re_1 < 0
        phi_i_1 = phi_i_1 + pi;
    elseif im_1 < 0 & re_1 > 0
        phi_i_1 = 2*pi - abs(phi_i_1);
    end
phi_M_1(ii) = phi_i_1;

re_2 = s_sum_2*ref_signal;
im_2 = s_sum_2*ref_signal_90;

phi_i_2 =atan(im_2/re_2);
    if im_2 > 0 & re_2 < 0
        phi_i_2 = pi - abs(phi_i_2);
    elseif im_2 < 0 & re_2 < 0
        phi_i_2 = phi_i_2 + pi;
    elseif im_2 < 0 & re_2 > 0
        phi_i_2 = 2*pi - abs(phi_i_2);
    end
phi_M_2(ii) = phi_i_2;
%[~,index] = max(s_sum(1:100));
%phi_M(ii) = index;
end
phi_M_1 (246:end) = phi_M_1(246:end)+2*pi;
phi_M_2 (246:end) = phi_M_2(246:end)+2*pi;

f1 = figure;plot(v,phi_M_1);hold on;plot(v,phi_M_2);legend("Pulse upstream", "Pulse downstream");xlabel("Fluid velocity (cm/s)");ylabel("Phase (rad)");grid on;
title("Evolution of phase of pulses in both directions when increasing velocity and height")
%%
figure;yyaxis left;plot([0:N-1],angles(1:N)*180/pi,'-o');ylabel("Output angle \alpha with respect to medium interface (º)");
yyaxis right;plot([0:N-1],att_vec,'-sq');xlabel("Path number (n)")
ylabel('Normalized amplitude');grid on
%% correct the phase repetition
dif = phi_M_1-phi_M_2;
figure
%plot(v,phi_M_1);
figure;plot(dif);
for ii =1:MM
    if dif(ii) < pi
        dif(ii) = dif(ii)+2*pi;
    end
    if dif(ii) > pi
        dif(ii) = dif(ii) - 2 * pi;
    end
end
f2 = figure;plot(v,dif);xlabel("Fluid velocity (cm/s)");ylabel("Phase (rad)");grid on;title("Phase difference")
%%
f3 = figure;
subplot(1,2,1)
plot(100*v,phi_M_1);hold on;plot(100*v,phi_M_2);legend("Pulse upstream", "Pulse downstream");xlabel("Fluid velocity (cm/s)");ylabel("Phase (rad)");grid on;title("Absolute phase in both directions")
subplot(1,2,2)
plot(100*v,-dif);xlabel("Fluid velocity (cm/s)");ylabel("Phase difference (rad)");grid on;title("Phase difference")
sgtitle("Simulation of wave superposition and cover movement")