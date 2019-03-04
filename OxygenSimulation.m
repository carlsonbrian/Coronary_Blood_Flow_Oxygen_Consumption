
% initial masses
x0(1) = 1*V_pa(1);
x0(2) = 1*V_11(1);
x0(3) = 0.25*V_21(1); % Rest 0.17, ExwUniReg 0.39, ExwVarReg 0.25
x0(4) = 1*V_12(1);
x0(5) = 0.25*V_22(1); % Rest 0.20, ExwUniReg 0.25, ExwVarReg 0.25
x0(6) = 1.0*V_13(1);
x0(7) = 0.05*V_23(1); % Rest 0.26, ExwUniReg 0.00, ExwVarReg 0.05
x0(8) = 0.2*V_pv(1);
x0(9) = 0.25;         % Rest 0.17, ExwUniReg 0.39, ExwVarReg 0.25
x0(10) = 0.25;        % Rest 0.20, ExwUniReg 0.25, ExwVarReg 0.25
x0(11) = 0.05;        % Rest 0.26, ExwUniReg 0.00, ExwVarReg 0.05

% simulation
T = t;
[tOxy,xOxy] = ode15s(@dCdT_oxygen,[0 9.99],x0,[],T,V_pa,V_11, ...
    V_12,V_21,V_22,V_13,V_23,V_pv,Q_pa,Q_11,Q_m1,Q_12,Q_21, ...
    Q_m2,Q_22,Q_13,Q_m3,Q_23,Q_pv);

% state variables
m_pa = xOxy(:,1);
m_11 = xOxy(:,2);
m_21 = xOxy(:,3);
m_12 = xOxy(:,4);
m_22 = xOxy(:,5);
m_13 = xOxy(:,6);
m_23 = xOxy(:,7);
m_pv = xOxy(:,8);
c_t1 = xOxy(:,9);
c_t2 = xOxy(:,10);
c_t3 = xOxy(:,11);

% volumes
v_pa = interp1(T,V_pa,tOxy);
v_11 = interp1(T,V_11,tOxy);
v_12 = interp1(T,V_12,tOxy);
v_21 = interp1(T,V_21,tOxy);
v_22 = interp1(T,V_22,tOxy);
v_13 = interp1(T,V_13,tOxy);
v_23 = interp1(T,V_23,tOxy);
v_pv = interp1(T,V_pv,tOxy) + 0.1; % add offset to v_pv to keep it stays positive

% concentrations
c_pa = m_pa./v_pa;
c_11 = m_11./v_11;
c_21 = m_21./v_21;
c_12 = m_12./v_12;
c_22 = m_22./v_22;
c_13 = m_13./v_13;
c_23 = m_23./v_23;
c_pv = m_pv./v_pv;

figure(6); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
plot(tOxy,c_13,'--r','linewidth',1.5); 
plot(tOxy,c_12,'--g','linewidth',1.5); 
plot(tOxy,c_11,'--b','linewidth',1.5); 
plot(tOxy,c_23,'-r','linewidth',1.5); 
plot(tOxy,c_22,'-g','linewidth',1.5); 
plot(tOxy,c_21,'-b','linewidth',1.5); 
l = legend('endo in','mid in','epi in','endo out','mid out','epi out'); 
set(l,'fontsize',12,'location','east');
set(gca,'fontsize',14); box on;
ylabel('Oxygen Saturation','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([5 10 -0.1 1.1]); grid

figure(7); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
plot(tOxy,c_23,'-r','linewidth',1.5); 
plot(tOxy,c_22,'-g','linewidth',1.5); 
plot(tOxy,c_21,'-b','linewidth',1.5); 
l = legend('endo out','mid out','epi out'); 
set(l,'fontsize',12,'location','southeast');
set(gca,'fontsize',14); box on;
ylabel('Oxygen Saturation','interpreter','latex','fontsize',16);
xlabel('time (sec)','interpreter','latex','fontsize',16);
axis([0 10 -0.1 0.35]); grid
%%
figure(8); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
plot(tOxy,c_23,'-r','linewidth',3); 
plot(tOxy,c_22,'-g','linewidth',3); 
plot(tOxy,c_21,'-b','linewidth',3); 
l = legend('endo out','mid out','epi out'); 
set(l,'fontsize',16,'location','southeast');
set(gca,'fontsize',20); box on;
ylabel('Oxygen Saturation','interpreter','latex','fontsize',24);
xlabel('time (sec)','interpreter','latex','fontsize',24);
axis([5 10 -0.2 0.55]); grid


% figure(7); clf; axes('position',[0.15 0.15 0.75 0.75]);  hold on;
% plot(tOxy,c_23,'color','r','linewidth',1.5); 
% plot(tOxy,c_22,'color','g','linewidth',1.5); 
% plot(tOxy,c_21,'color','b','linewidth',1.5); 
% l = legend('endo','mid','epi'); 
% set(l,'fontsize',12,'location','northeast');
% set(gca,'fontsize',14); box on;
% ylabel('Norm Oxygen Conc Out','interpreter','latex','fontsize',16);
% xlabel('time (sec)','interpreter','latex','fontsize',16);
% axis([5 10 0 1.1]); grid

% venous oxygen content
mean(interp1(tOxy,c_pv,4:0.01:9.99))
