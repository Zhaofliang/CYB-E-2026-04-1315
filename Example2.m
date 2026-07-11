clear all;
clc;
rng(11,'twister');


RR = 1;                 
zL = 0; zU = 1;
t0 = 0; tf = 2;
d  = 1;               

h1 = 0.05;
h2 = 0.10;
alpha_bar = 0.75;       
eta_bar   = 0.60;        
kappa_bar = 0.30;        

J = 5;
K1 = [-1.1350  0.1882];
K2 = [-1.1478  0.1896];

t_sample = 0.002;
z_sample = 0.08;
N_t = round((tf-t0+d)/t_sample)+1;
N_z = round((zU-zL)/z_sample)+1;
M_t = round(d/t_sample);
ttt = (0:N_t-1)*t_sample;
t_plot = ttt - ttt(M_t);


xbar1  = zeros(N_t,N_z);
xbar2  = zeros(N_t,N_z);
xhbar1 = zeros(N_t,N_z);
xhbar2 = zeros(N_t,N_z);
ebar1  = zeros(N_t,N_z);
ebar2  = zeros(N_t,N_z);

ybar1  = zeros(N_t,N_z);
ybar2  = zeros(N_t,N_z);
yhbar1 = zeros(N_t,N_z);
yhbar2 = zeros(N_t,N_z);
yybar1  = zeros(N_t,N_z);
yybar2  = zeros(N_t,N_z);
yyhbar1 = zeros(N_t,N_z);
yyhbar2 = zeros(N_t,N_z);

u1 = zeros(1,N_t);
u2 = zeros(1,N_t);
u1_nom = zeros(1,N_t);
u2_nom = zeros(1,N_t);
eta1 = zeros(1,N_t);
kappa = zeros(1,N_t);
sample_index = zeros(1,N_t);
sample_type  = zeros(1,N_t);

zzz = zeros(1,N_z);
for iz = 1:N_z
    zzz(iz) = (iz-1)*z_sample;
end


for it = 1:M_t
    for iz = 1:N_z
        c = zzz(iz);
        xbar1(it,iz)  = -5*sin(pi*c);
        xbar2(it,iz)  = -1*sin(pi*c);
        xhbar1(it,iz) =  5*sin(pi*c);
        xhbar2(it,iz) =  1*sin(pi*c);

        ebar1(it,iz) = xhbar1(it,iz)-xbar1(it,iz);
        ebar2(it,iz) = xhbar2(it,iz)-xbar2(it,iz);

        ybar1(it,iz)  = -5*pi*cos(pi*c);
        ybar2(it,iz)  = -1*pi*cos(pi*c);
        yhbar1(it,iz) =  5*pi*cos(pi*c);
        yhbar2(it,iz) =  1*pi*cos(pi*c);

        yybar1(it,iz)  =  5*pi^2*sin(pi*c);
        yybar2(it,iz)  =  1*pi^2*sin(pi*c);
        yyhbar1(it,iz) = -5*pi^2*sin(pi*c);
        yyhbar2(it,iz) = -1*pi^2*sin(pi*c);
    end
end


xbar1(1:M_t,[1,N_z])  = 0;  xbar2(1:M_t,[1,N_z])  = 0;
xhbar1(1:M_t,[1,N_z]) = 0;  xhbar2(1:M_t,[1,N_z]) = 0;
ebar1(1:M_t,[1,N_z])  = 0;  ebar2(1:M_t,[1,N_z])  = 0;


zoh1 = round(h1/t_sample);
zoh2 = round(h2/t_sample);

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
sample_type(sample_type==0) = 1;

san1 = 4;
san2 = 6;
idx_s1 = round(san1*N_z/10);
idx_s2 = round(san2*N_z/10);

for it = M_t:(N_t-1)
    a_now = max((it-M_t)*t_sample,0);
    tau_now = exp(a_now)/(exp(a_now)+1);       
    id_tau = max(1,it-round(tau_now/t_sample));


    id_samp = sample_index(it);

    ys1 = [ebar1(id_samp,idx_s1); ebar2(id_samp,idx_s1)];
    ys2 = [ebar1(id_samp,idx_s2); ebar2(id_samp,idx_s2)];

    F1_ys1 = [tanh(0.2*ys1(1)); tanh(0.2*ys1(2))];
    F2_ys1 = [tanh(0.3*ys1(1)); tanh(0.1*ys1(2))];

    F1_ys2 = [tanh(0.2*ys2(1)); tanh(0.2*ys2(2))];
    F2_ys2 = [tanh(0.3*ys2(1)); tanh(0.1*ys2(2))];

    u1_nom(it) = K1*ys1;
    u2_nom(it) = K2*ys2;

    u1(it) = K1*(eta1(it)*ys1 ...
        - (1-eta1(it))*(kappa(it)*F1_ys1 + (1-kappa(it))*F2_ys1));

    u2(it) = K2*(eta1(it)*ys2 ...
        - (1-eta1(it))*(kappa(it)*F1_ys2 + (1-kappa(it))*F2_ys2));


    for iz = 2:(N_z-1)
        c = zzz(iz);

        f1  = tanh(xbar1(it,iz));
        f2  = tanh(xbar2(it,iz));
        fh1 = tanh(xhbar1(it,iz));
        fh2 = tanh(xhbar2(it,iz));

        fd1  = tanh(xbar1(id_tau,iz));
        fd2  = tanh(xbar2(id_tau,iz));
        fhd1 = tanh(xhbar1(id_tau,iz));
        fhd2 = tanh(xhbar2(id_tau,iz));

        if c < 0.4
            us = u1(it);
        else
            us = u2(it);
        end

        xbar1(it+1,iz) = xbar1(it,iz) + t_sample*( ...
            RR*yybar1(it,iz) - xbar1(it,iz) ...
            + 0.66*f1 - 0.06*f2 - 0.10*fd1 - 0.30*fd2 + J);

        xbar2(it+1,iz) = xbar2(it,iz) + t_sample*( ...
            RR*yybar2(it,iz) - xbar2(it,iz) ...
            - 0.12*f1 + 0.54*f2 - 0.05*fd1 - 0.35*fd2);

        xhbar1(it+1,iz) = xhbar1(it,iz) + t_sample*( ...
            RR*yyhbar1(it,iz) - xhbar1(it,iz) ...
            + 0.66*fh1 - 0.06*fh2 - 0.10*fhd1 - 0.30*fhd2 + us + J);

        xhbar2(it+1,iz) = xhbar2(it,iz) + t_sample*( ...
            RR*yyhbar2(it,iz) - xhbar2(it,iz) ...
            - 0.12*fh1 + 0.54*fh2 - 0.05*fhd1 - 0.35*fhd2);

        ebar1(it+1,iz) = xhbar1(it+1,iz)-xbar1(it+1,iz);
        ebar2(it+1,iz) = xhbar2(it+1,iz)-xbar2(it+1,iz);
    end


    xbar1(it+1,1)=0;  xbar1(it+1,N_z)=0;
    xbar2(it+1,1)=0;  xbar2(it+1,N_z)=0;
    xhbar1(it+1,1)=0; xhbar1(it+1,N_z)=0;
    xhbar2(it+1,1)=0; xhbar2(it+1,N_z)=0;
    ebar1(it+1,1)=0;  ebar1(it+1,N_z)=0;
    ebar2(it+1,1)=0;  ebar2(it+1,N_z)=0;


    for iz = 2:N_z
        ybar1(it+1,iz)  = (xbar1(it+1,iz)-xbar1(it+1,iz-1))/z_sample;
        ybar2(it+1,iz)  = (xbar2(it+1,iz)-xbar2(it+1,iz-1))/z_sample;
        yhbar1(it+1,iz) = (xhbar1(it+1,iz)-xhbar1(it+1,iz-1))/z_sample;
        yhbar2(it+1,iz) = (xhbar2(it+1,iz)-xhbar2(it+1,iz-1))/z_sample;
    end
    ybar1(it+1,1)  = (xbar1(it+1,2)-xbar1(it+1,1))/z_sample;
    ybar2(it+1,1)  = (xbar2(it+1,2)-xbar2(it+1,1))/z_sample;
    yhbar1(it+1,1) = (xhbar1(it+1,2)-xhbar1(it+1,1))/z_sample;
    yhbar2(it+1,1) = (xhbar2(it+1,2)-xhbar2(it+1,1))/z_sample;

    for iz = 2:(N_z-1)
        yybar1(it+1,iz)  = (ybar1(it+1,iz+1)-ybar1(it+1,iz))/z_sample;
        yybar2(it+1,iz)  = (ybar2(it+1,iz+1)-ybar2(it+1,iz))/z_sample;
        yyhbar1(it+1,iz) = (yhbar1(it+1,iz+1)-yhbar1(it+1,iz))/z_sample;
        yyhbar2(it+1,iz) = (yhbar2(it+1,iz+1)-yhbar2(it+1,iz))/z_sample;
    end
    yybar1(it+1,1)  = (ybar1(it+1,2)-ybar1(it+1,1))/z_sample;
    yybar2(it+1,1)  = (ybar2(it+1,2)-ybar2(it+1,1))/z_sample;
    yyhbar1(it+1,1) = (yhbar1(it+1,2)-yhbar1(it+1,1))/z_sample;
    yyhbar2(it+1,1) = (yhbar2(it+1,2)-yhbar2(it+1,1))/z_sample;
    yybar1(it+1,N_z)  = (ybar1(it+1,N_z)-ybar1(it+1,N_z-1))/z_sample;
    yybar2(it+1,N_z)  = (ybar2(it+1,N_z)-ybar2(it+1,N_z-1))/z_sample;
    yyhbar1(it+1,N_z) = (yhbar1(it+1,N_z)-yhbar1(it+1,N_z-1))/z_sample;
    yyhbar2(it+1,N_z) = (yhbar2(it+1,N_z)-yhbar2(it+1,N_z-1))/z_sample;
end


yy1 = zeros(1,N_t);
yy2 = zeros(1,N_t);
for it = 1:N_t
    yy1(it) = sqrt(sum(ebar1(it,:).^2)*z_sample);
    yy2(it) = sqrt(sum(ebar2(it,:).^2)*z_sample);
end



idx_u = M_t:(N_t-1);
t_u = t_plot(idx_u);
t_u = t_u - t_u(1);

figure(5)
set(5,'position',[131 400 344 179])
stairs(t_u,u1_nom(idx_u),'b')
xlim([0 tf]);
xlabel('time (sec.)', 'fontname','Times New Roman', 'FontSize',8);
set(gca,'fontname','Times New Roman', 'FontSize',8);
set(gca,'YGrid','off');

figure(6)
set(6,'position',[131 400 344 179])
stairs(t_u,u2_nom(idx_u),'b')
xlim([0 tf]);
xlabel('time (sec.)', 'fontname','Times New Roman', 'FontSize',8);
grid off;
set(gca,'YGrid','off');
set(gca,'fontname','Times New Roman', 'FontSize',8);

figure(7)
plot(t_plot, yy1,'k-',t_plot, yy2,'k:')
xlim([0 tf]);
xlabel('time (sec.)', 'fontname','Times New Roman', 'FontSize',8);
set(gca,'fontname','Times New Roman', 'FontSize',8);
h=legend('$$\|\textbf{e}_1(\cdot,a)\|_2$$','$$\|\textbf{e}_2(\cdot,a)\|_2$$','1');
set(h,'Interpreter','latex');


h_interval = h1*ones(1,N_t);
h_interval(sample_type==2) = h2;

idx_h_all = M_t:N_t;
t_h_all = ttt(idx_h_all) - ttt(M_t);

idx_keep = (t_h_all < tf);
t_h_plot = t_h_all(idx_keep);
h_plot   = h_interval(idx_h_all(idx_keep));


t_h_plot(end+1) = tf;
h_plot(end+1)   = h_plot(end);

figure(8)
set(8,'position',[131 400 344 179])
stairs(t_h_plot,h_plot,'b','LineWidth',1.5)
xlim([0 tf]);
ylim([0 1.2*max(h1,h2)]);
xlabel('time (sec.)', 'fontname','Times New Roman', 'FontSize',8);
ylabel('$$h_k$$','Interpreter','latex','fontname','Times New Roman','FontSize',8);
set(gca,'fontname','Times New Roman', 'FontSize',8);
grid on;
set(gca,'YGrid','on','XGrid','on');

