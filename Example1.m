clear all;
clc;
%%
RR=0.8;

zL=0;zU=1;
t0=0;tf=3;

d=1;
h1=0.08;
h2=0.12;

K1 =-1.1076;
K2 =-1.1076;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_sample=0.001;
z_sample=0.05;
N_t=round((tf-t0+d)/t_sample)+1;
N_z=round((zU-zL)/z_sample)+1;
M_t=round(d/t_sample);
time = t0 : t_sample : tf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for it=1:M_t
for iz=1:N_z
    zzz(iz)=(iz-1)*z_sample;
    xbar1(it,iz)=-2*sin(pi*(iz-1)*z_sample); 
 
    xhbar1(it,iz)=2*sin(pi*(iz-1)*z_sample); 

    ebar1(it,iz)=xhbar1(it,iz)-xbar1(it,iz); 

                                
   ybar1(it,iz)=-(pi)*2*cos(pi*(iz-1)*z_sample); 
                 
   yhbar1(it,iz)=(pi)*2*cos(pi*(iz-1)*z_sample); 


   yybar1(it,iz)=(pi^2)*2*sin(pi*(iz-1)*z_sample);

   yyhbar1(it,iz)=-(pi^2)*2*sin(pi*(iz-1)*z_sample);

end
end

% ================================================================
rng(12,'twister');       
eta_bar   = 0.5;        
kappa_bar = 0.4;       
alpha_bar = 0.8;        

eta1      = zeros(1,N_t);   
kappa     = zeros(1,N_t);  
u1_nom    = zeros(1,N_t);   
u2_nom    = zeros(1,N_t);
for it=1:(M_t-1)
for iz=1:N_z
    zzz(iz)=(iz-1)*z_sample;
    ttt(it)=(it-1)*t_sample;

    ybar1(it,1)=0;
    ybar1(it,N_z)=0;
    ybar2(it,1)=0;
    ybar2(it,N_z)=0; 
   u1(it)=0;   
   u2(it)=0;    
end
end
for it=(M_t):N_t
    ttt(it)=(it-1)*t_sample;
    ybar1(it,1)=0;
    ybar1(it,N_z)=0;
    ybar2(it,1)=0;
    ybar2(it,N_z)=0; 
     u1(it)=0;   
     u2(it)=0;    
end
for it=1:N_t
      sum1(it)=0; 
     yy1(it)=0;  
end


ys1_trace     = zeros(1,N_t);   
ys2_trace     = zeros(1,N_t);   
F1_ys1_trace  = zeros(1,N_t);   
F2_ys1_trace  = zeros(1,N_t);  
Fmix1_trace   = zeros(1,N_t);  
yrecv1_trace  = zeros(1,N_t);   
yrecv2_trace  = zeros(1,N_t);  
attack_flag   = zeros(1,N_t);  
attack_source = zeros(1,N_t);   

%%
J=0;
 zoh1=round(h1/t_sample);
 zoh2=round(h2/t_sample);


valid_sequence = false;
max_attempts = 500;
attempt = 0;

while ~valid_sequence && attempt < max_attempts
    attempt = attempt + 1;

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

    idx_show = M_t:N_t;
    g1_count = sum(eta1(idx_show)==0 & kappa(idx_show)==1);
    g2_count = sum(eta1(idx_show)==0 & kappa(idx_show)==0);
    no_count = sum(eta1(idx_show)==1);

    valid_sequence = (g1_count >= round(0.12/t_sample)) && ...
                     (g2_count >= round(0.12/t_sample)) && ...
                     (no_count >= round(0.20/t_sample));
end

%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it=M_t:(N_t-1)
  for iz=1:N_z
dd(it)=exp(sin(ttt(it)))/(exp(sin(ttt(it)))+1);
  end
%%%%%%%%%%%%%%%%%
  for iz=1:N_z
%%  
 san1=4;
 san2=6;
 
 a_now = max(ttt(it)-d,0);
 tau_now = 0.5 + 0.5*sin(a_now);
 id_tau = max(1, it-round(tau_now/t_sample));
 fff1  = (abs(xbar1(it,iz)+1)-abs(xbar1(it,iz)-1))/2;
 fffh1 = (abs(xhbar1(it,iz)+1)-abs(xhbar1(it,iz)-1))/2;
 ffd1  = (abs(xbar1(id_tau,iz)+1)-abs(xbar1(id_tau,iz)-1))/2;
 ffhd1 = (abs(xhbar1(id_tau,iz)+1)-abs(xhbar1(id_tau,iz)-1))/2;
 id_samp = sample_index(it);
 ys1 = ebar1(id_samp,round(san1*N_z/10));
 ys2 = ebar1(id_samp,round(san2*N_z/10));
 F1_ys1 = tanh(0.2*ys1);
 F2_ys1 = tanh(0.1*ys1);
 F1_ys2 = tanh(0.2*ys2);
 F2_ys2 = tanh(0.1*ys2);
 u1_nom(it) = K1*ys1;
 u2_nom(it) = K2*ys2;
 u1(it) = K1*(eta1(it)*ys1 - (1-eta1(it))*(kappa(it)*F1_ys1 + (1-kappa(it))*F2_ys1));
 u2(it) = K2*(eta1(it)*ys2 - (1-eta1(it))*(kappa(it)*F1_ys2 + (1-kappa(it))*F2_ys2));


 ys1_trace(it)     = ys1;
 ys2_trace(it)     = ys2;
 F1_ys1_trace(it)  = F1_ys1;
 F2_ys1_trace(it)  = F2_ys1;
 Fmix1_trace(it)   = kappa(it)*F1_ys1 + (1-kappa(it))*F2_ys1;
 yrecv1_trace(it)  = eta1(it)*ys1 - (1-eta1(it))*Fmix1_trace(it);
 yrecv2_trace(it)  = eta1(it)*ys2 - (1-eta1(it))*(kappa(it)*F1_ys2 + (1-kappa(it))*F2_ys2);
 attack_flag(it)   = 1 - eta1(it);
 if attack_flag(it) == 0
     attack_source(it) = 0;
 elseif kappa(it) == 1
     attack_source(it) = 1;
 else
     attack_source(it) = 2;
 end

 
           if iz>=2 && iz<=round(5*N_z/10) 

             
      xbar1(it+1,iz) = xbar1(it,iz)+t_sample*(RR*yybar1(it,iz)-(0.6)*xbar1(it,iz)+0.1*fff1+0.4*ffd1+J);
   
     
      xhbar1(it+1,iz) = xhbar1(it,iz)+t_sample*(RR*yyhbar1(it,iz)-(0.6)*xhbar1(it,iz)+0.1*fffh1+0.4*ffhd1+u1(it)+J);

      
      ebar1(it+1,iz) = xhbar1(it+1,iz)-xbar1(it+1,iz);


            end
           
      
        if iz>=round(5*N_z/10)+1 &&iz<=N_z-1 
                
             
      xbar1(it+1,iz) = xbar1(it,iz)+t_sample*(RR*yybar1(it,iz)-(0.6)*xbar1(it,iz)+0.1*fff1+0.4*ffd1+J);
     
     
      xhbar1(it+1,iz) = xhbar1(it,iz)+t_sample*(RR*yyhbar1(it,iz)-(0.6)*xhbar1(it,iz)+0.1*fffh1+0.4*ffhd1+u2(it)+J);
      
      
      ebar1(it+1,iz) = xhbar1(it+1,iz)-xbar1(it+1,iz);
     
         end
  end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iz=2:N_z
     if iz>=2&&iz<=round(5*N_z/10)  
           ybar1(it+1,iz)=(xbar1(it+1,iz)-xbar1(it+1,iz-1))/z_sample;

           
           yhbar1(it+1,iz)=(xhbar1(it+1,iz)-xhbar1(it+1,iz-1))/z_sample;%%

     end
   if iz>=floor(5*N_z/10)+1&&iz<=N_z
         ybar1(it+1,iz)=(xbar1(it+1,iz)-xbar1(it+1,iz-1))/z_sample;

         
         yhbar1(it+1,iz)=(xhbar1(it+1,iz)-xhbar1(it+1,iz-1))/z_sample;%%

   end
      
    
end
    ybar1(it+1,1)=(xbar1(it+1,2)-xbar1(it+1,1))/z_sample;
    ybar1(it+1,N_z)=(xbar1(it+1,N_z)-xbar1(it+1,N_z-1))/z_sample;
    
    yhbar1(it+1,1)=(xhbar1(it+1,2)-xhbar1(it+1,1))/z_sample;%%
    yhbar1(it+1,N_z)=(xhbar1(it+1,N_z)-xhbar1(it+1,N_z-1))/z_sample;%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
for iz=2:(N_z-1)
    if iz>=2&&iz<=round(5*N_z/10)  
      yybar1(it+1,iz)=(ybar1(it+1,iz+1)-ybar1(it+1,iz))/z_sample;

      
      yyhbar1(it+1,iz)=(yhbar1(it+1,iz+1)-yhbar1(it+1,iz))/z_sample;

    end
       if iz>=floor(5*N_z/10)+1&&iz<=N_z
      yybar1(it+1,iz)=(ybar1(it+1,iz+1)-ybar1(it+1,iz))/z_sample;
  
      
       yyhbar1(it+1,iz)=(yhbar1(it+1,iz+1)-yhbar1(it+1,iz))/z_sample;

       end
         
end
    yybar1(it+1,1)=(ybar1(it+1,2)-ybar1(it+1,1))/z_sample;
    yybar1(it+1,N_z)=(ybar1(it+1,N_z)-ybar1(it+1,N_z-1))/z_sample;
    yyhbar1(it+1,1)=(yhbar1(it+1,2)-yhbar1(it+1,1))/z_sample;
    yyhbar1(it+1,N_z)=(yhbar1(it+1,N_z)-yhbar1(it+1,N_z-1))/z_sample;
    
 end
    
 for it=1:N_t
  for iz=1:N_z
   
    sum1(it)=sum1(it)+ebar1(it,iz)*ebar1(it,iz)*z_sample; 
   
  end
  
     yy1(it)=sqrt(sum1(it));
   
 end

 


figure(4)
set(4,'position',[131 400 344 179])
plot(ttt-d,yy1,'b')
xlim([0 tf]);
xlabel('time (sec.)', 'fontname','Times New Roman', 'FontSize',8);
grid off;
set(gca, 'YGrid', 'off');
set(gca,'fontname','Times New Roman', 'FontSize',8);


h_interval = h1*ones(1,N_t);
h_interval(sample_type==2) = h2;
idx_h = M_t:N_t;
t_h = ttt(idx_h) - ttt(M_t);

figure(5)
set(5,'position',[131 400 344 179])
stairs(t_h,h_interval(idx_h),'b','LineWidth',1.5)
xlim([0 tf]);
ylim([0 1.2*max(h1,h2)]);
xlabel('time (sec.)', 'fontname','Times New Roman', 'FontSize',8);
ylabel('$$h_k$$','Interpreter','latex','fontname','Times New Roman','FontSize',8);
set(gca,'fontname','Times New Roman', 'FontSize',8);
grid on;
set(gca,'YGrid','on','XGrid','on');




figure(6)
set(6,'position',[150 80 430 430])

subplot(2,1,1)
hold on;

t_normal = ttt(M_t:N_t) - d;
normal_meas = ys1_trace(M_t:N_t);

x_step = zeros(1,2*length(t_normal)-1);
y_step = zeros(1,2*length(normal_meas)-1);

x_step(1:2:end) = t_normal;
x_step(2:2:end) = t_normal(2:end);

y_step(1:2:end) = normal_meas;
y_step(2:2:end) = normal_meas(1:end-1);

plot(x_step, y_step, '-', 'Color',[0 1 0], 'LineWidth',1.8);

xlim([0 tf]);


ylim([min(normal_meas)-0.1, max(normal_meas)+0.2]);

ylabel('Normal measurement', ...
       'fontname','Times New Roman', ...
       'FontSize',8);

h_clean = line(nan,nan, ...
    'Color',[0 1 0], ...
    'LineStyle','-', ...
    'LineWidth',1.8);

legend(h_clean,'normal measurement','Location','best');

set(gca,'fontname','Times New Roman','FontSize',8);
grid off;
box on;

text(0.50,-0.22,'(a)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'fontname','Times New Roman', ...
    'FontSize',9);

subplot(2,1,2)
hold on;

t_attack = ttt - d;
recv_meas = yrecv1_trace;
mode_trace = attack_source;  

for it = M_t:(N_t-1)
    x1 = t_attack(it);
    x2 = t_attack(it+1);
    y1 = recv_meas(it);
    y2 = recv_meas(it+1);

    if mode_trace(it) == 0
        line_color = [0 1 0];   
    elseif mode_trace(it) == 1
        line_color = [1 0 0];     
    else
        line_color = [0 0 1];    
    end
    line([x1 x2], [y1 y1], 'Color', line_color, 'LineStyle', '-', 'LineWidth', 1.5);
    if abs(y2-y1) > 1e-12
        line([x2 x2], [y1 y2], 'Color', [0 1 0], 'LineStyle', '-', 'LineWidth', 1.5);
    end
end
h_no = line(nan,nan,'Color',[0 1 0],'LineStyle','-','LineWidth',1.5);
h_g1 = line(nan,nan,'Color',[1 0 0],'LineStyle','-','LineWidth',1.5);
h_g2 = line(nan,nan,'Color',[0 0 1],'LineStyle','-','LineWidth',1.5);
xlim([0 tf]);
xlabel('time (sec.)', 'fontname','Times New Roman', 'FontSize',8);
ylabel('Received measurement', 'fontname','Times New Roman', 'FontSize',8);
h_leg1 = legend([h_no h_g1 h_g2], ...
       'normal  measurement', ...
       '$\mathcal{G}_1(\cdot)$-attacked measurement', ...
       '$\mathcal{G}_2(\cdot)$-attacked measurement', ...
       'Location','best');
set(h_leg1,'Interpreter','latex','FontSize',8);
set(gca,'fontname','Times New Roman', 'FontSize',8);
grid off; box on;
text(0.50,-0.28,'(b)','Units','normalized','HorizontalAlignment','center', ...
     'fontname','Times New Roman','FontSize',9)
