clear; clc; close all;

%% Monte Carlo setting
Monte_Carlo = 60;
rng(12,'twister');                
init_dispersion = 0.3*randn(Monte_Carlo,1);  


RR = 0.8;                       
zL = 0; zU = 1;
t0 = 0; tf = 3;

d = 1;                          
h1 = 0.08;
h2 = 0.12;

K1 = -1.1076;
K2 = -1.1076;

t_sample = 0.001;
z_sample = 0.05;
N_t = round((tf-t0+d)/t_sample) + 1;
N_z = round((zU-zL)/z_sample) + 1;
M_t = round(d/t_sample);
ttt = (0:N_t-1)*t_sample;
zzz = (0:N_z-1)*z_sample;


alpha_bar = 0.8;                  
eta_bar   = 0.5;                
kappa_bar = 0.4;                  

zoh1 = round(h1/t_sample);
zoh2 = round(h2/t_sample);

san1 = 4;
san2 = 6;
sensor_1 = round(san1*N_z/10);
sensor_2 = round(san2*N_z/10);
subdomain_mid = round(5*N_z/10);

J = 0;

final_Error1 = zeros(Monte_Carlo,N_t);


for j = 1:Monte_Carlo

  
    xbar1   = zeros(N_t,N_z);    
    xhbar1  = zeros(N_t,N_z);   
    ebar1   = zeros(N_t,N_z);    
    ybar1   = zeros(N_t,N_z);    
    yhbar1  = zeros(N_t,N_z);     
    yybar1  = zeros(N_t,N_z);    
    yyhbar1 = zeros(N_t,N_z);     

    u1      = zeros(1,N_t);      
    u2      = zeros(1,N_t);      
    u1_nom  = zeros(1,N_t);       
    u2_nom  = zeros(1,N_t);

  
    amp_hat = 2 + init_dispersion(j);

    for it = 1:M_t
        for iz = 1:N_z
            c = zzz(iz);

            xbar1(it,iz)  = -2*sin(pi*c);
            xhbar1(it,iz) = amp_hat*sin(pi*c);
            ebar1(it,iz)  = xhbar1(it,iz) - xbar1(it,iz);

            ybar1(it,iz)   = -2*pi*cos(pi*c);
            yhbar1(it,iz)  =  amp_hat*pi*cos(pi*c);
            yybar1(it,iz)  =  2*pi^2*sin(pi*c);
            yyhbar1(it,iz) = -amp_hat*pi^2*sin(pi*c);
        end
    end

  
    sample_index = zeros(1,N_t);  
    sample_type  = zeros(1,N_t); 
    eta1         = zeros(1,N_t);
    kappa        = zeros(1,N_t);

    k_sample = M_t;
    while k_sample <= N_t
        if rand <= alpha_bar
            h_step = zoh1;
            h_type = 1;
        else
            h_step = zoh2;
            h_type = 2;
        end

        idx_block = k_sample:min(k_sample+h_step-1,N_t);
        sample_index(idx_block) = k_sample;
        sample_type(idx_block)  = h_type;


        eta_k   = double(rand <= eta_bar);
        kappa_k = double(rand <= kappa_bar);
        eta1(idx_block)  = eta_k;
        kappa(idx_block) = kappa_k;

        k_sample = k_sample + h_step;
    end

    sample_index(sample_index==0) = M_t;
    sample_type(sample_type==0)   = 1;

  
    for it = M_t:(N_t-1)
    
        a_now = max(ttt(it)-d,0);

    
        tau_now = 0.5 + 0.5*sin(a_now);
        id_tau = max(1, it-round(tau_now/t_sample));


        id_samp = sample_index(it);
        ys1 = ebar1(id_samp,sensor_1);
        ys2 = ebar1(id_samp,sensor_2);

     
        F1_ys1 = tanh(0.2*ys1);
        F2_ys1 = tanh(0.1*ys1);
        F1_ys2 = tanh(0.2*ys2);
        F2_ys2 = tanh(0.1*ys2);

   
        u1_nom(it) = K1*ys1;
        u2_nom(it) = K2*ys2;
        u1(it) = K1*(eta1(it)*ys1 - (1-eta1(it))*(kappa(it)*F1_ys1 + (1-kappa(it))*F2_ys1));
        u2(it) = K2*(eta1(it)*ys2 - (1-eta1(it))*(kappa(it)*F1_ys2 + (1-kappa(it))*F2_ys2));


        xbar1(it+1,1) = 0;       xbar1(it+1,N_z) = 0;
        xhbar1(it+1,1) = 0;      xhbar1(it+1,N_z) = 0;
        ebar1(it+1,1) = 0;       ebar1(it+1,N_z) = 0;

        for iz = 2:(N_z-1)

            fff1  = (abs(xbar1(it,iz)+1)  - abs(xbar1(it,iz)-1))/2;
            fffh1 = (abs(xhbar1(it,iz)+1) - abs(xhbar1(it,iz)-1))/2;
            ffd1  = (abs(xbar1(id_tau,iz)+1)  - abs(xbar1(id_tau,iz)-1))/2;
            ffhd1 = (abs(xhbar1(id_tau,iz)+1) - abs(xhbar1(id_tau,iz)-1))/2;

            if iz <= subdomain_mid
                us = u1(it);
            else
                us = u2(it);
            end

            xbar1(it+1,iz) = xbar1(it,iz) + t_sample*( ...
                RR*yybar1(it,iz) - 0.6*xbar1(it,iz) + 0.1*fff1 + 0.4*ffd1 + J);

            xhbar1(it+1,iz) = xhbar1(it,iz) + t_sample*( ...
                RR*yyhbar1(it,iz) - 0.6*xhbar1(it,iz) + 0.1*fffh1 + 0.4*ffhd1 + us + J);

            ebar1(it+1,iz) = xhbar1(it+1,iz) - xbar1(it+1,iz);
        end


        for iz = 2:N_z
            ybar1(it+1,iz)  = (xbar1(it+1,iz)  - xbar1(it+1,iz-1))/z_sample;
            yhbar1(it+1,iz) = (xhbar1(it+1,iz) - xhbar1(it+1,iz-1))/z_sample;
        end
        ybar1(it+1,1)  = (xbar1(it+1,2)  - xbar1(it+1,1))/z_sample;
        yhbar1(it+1,1) = (xhbar1(it+1,2) - xhbar1(it+1,1))/z_sample;

        for iz = 2:(N_z-1)
            yybar1(it+1,iz)  = (ybar1(it+1,iz+1)  - ybar1(it+1,iz))/z_sample;
            yyhbar1(it+1,iz) = (yhbar1(it+1,iz+1) - yhbar1(it+1,iz))/z_sample;
        end
        yybar1(it+1,1)   = (ybar1(it+1,2)   - ybar1(it+1,1))/z_sample;
        yybar1(it+1,N_z) = (ybar1(it+1,N_z) - ybar1(it+1,N_z-1))/z_sample;
        yyhbar1(it+1,1)   = (yhbar1(it+1,2)   - yhbar1(it+1,1))/z_sample;
        yyhbar1(it+1,N_z) = (yhbar1(it+1,N_z) - yhbar1(it+1,N_z-1))/z_sample;
    end


    yy1 = sqrt(sum(ebar1.^2,2)'*z_sample);
    final_Error1(j,:) = yy1;
end


figure(1)
set(1,'position',[131 400 344 179])
hold on
for j = 1:Monte_Carlo
    plot(ttt-d, final_Error1(j,:), 'LineWidth', 0.5)
end
xlim([0 tf]);
ylim([0 4]);
xlabel('time (sec.)', 'fontname','Times New Roman', 'FontSize',8);
set(gca,'fontname','Times New Roman', 'FontSize',8);
grid off;
box on;

