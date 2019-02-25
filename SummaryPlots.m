%           REST    EXER.
F_Exp     = [70.8  102.3];

%           REST
F_Model   = [71.8  72.3  102    100]; 
R_EpiEndo = [1.14  0.80  0.51   0.85];
    
figure(10); clf; axes('position',[0.15 0.20 0.30 0.75]);  hold on;
ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',16);
plot([0 1],[F_Exp(1) F_Exp(1)],'k--','linewidth',1.5)
plot([0.2],F_Model(1),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',[1 1 1]);
set(gca,'fontsize',14); box on;
axis([0 1 0 120]);
set(gca,'xtick',[0.2 ],'xticklabel',{'Rest'})

axes('position',[0.65 0.20 0.30 0.75]);  hold on;
ylabel('Endo/Epi Flow Ratio','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); box on;
plot([0.2],R_EpiEndo(1),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',[1 1 1]);
axis([0 1 0 1.2]);
set(gca,'xtick',[0.2 ],'xticklabel',{'Rest'})


%
figure(11); clf; axes('position',[0.15 0.20 0.30 0.75]);  hold on;
ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',16);
plot([0 1],[F_Exp(2) F_Exp(2)],'k--','linewidth',1.5)
plot([0.2],F_Model(1),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',[1 1 1]);
plot([0.8],F_Model(2),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',0.5*[1 1 1]);
set(gca,'fontsize',14); box on;
axis([0 1 0 120]);
set(gca,'xtick',[0.2 0.8],'xticklabel',{'Rest','Exercise'})

axes('position',[0.65 0.20 0.30 0.75]);  hold on;
ylabel('Endo/Epi Flow Ratio','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); box on;
plot([0.2],R_EpiEndo(1),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',[1 1 1]);
plot([0.8],R_EpiEndo(2),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',0.5*[1 1 1]);
axis([0 1 0 1.2]);
set(gca,'xtick',[0.2 0.8],'xticklabel',{'Rest','Exercise'})

%
figure(12); clf; axes('position',[0.15 0.20 0.30 0.75]);  hold on;
ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',16);
plot([0 1],[F_Exp(2) F_Exp(2)],'k--','linewidth',1.5)
plot([0.2],F_Model(1),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',[1 1 1]);
plot([0.8],F_Model(3),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',0.5*[1 1 1]);
set(gca,'fontsize',14); box on;
axis([0 1 0 120]);
set(gca,'xtick',[0.2 0.8],'xticklabel',{'Rest','Exercise'})

axes('position',[0.65 0.20 0.30 0.75]);  hold on;
ylabel('Endo/Epi Flow Ratio','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); box on;
plot([0.2],R_EpiEndo(1),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',[1 1 1]);
plot([0.8],R_EpiEndo(3),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',0.5*[1 1 1]);
axis([0 1 0 1.2]);
set(gca,'xtick',[0.2 0.8],'xticklabel',{'Rest','Exercise'})

%
figure(13); clf; axes('position',[0.15 0.20 0.30 0.75]);  hold on;
ylabel('Myocardial flow (mL/min)','interpreter','latex','fontsize',16);
plot([0 1],[F_Exp(2) F_Exp(2)],'k--','linewidth',1.5)
plot([0.2],F_Model(1),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',[1 1 1]);
plot([0.8],F_Model(4),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',0.5*[1 1 1]);
set(gca,'fontsize',14); box on;
axis([0 1 0 120]);
set(gca,'xtick',[0.2 0.8],'xticklabel',{'Rest','Exercise'})

axes('position',[0.65 0.20 0.30 0.75]);  hold on;
ylabel('Endo/Epi Flow Ratio','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); box on;
plot([0.2],R_EpiEndo(1),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',[1 1 1]);
plot([0.8],R_EpiEndo(4),'ko','linewidth',1.5,'Markersize',12,'Markerfacecolor',0.5*[1 1 1]);
axis([0 1 0 1.2]);
set(gca,'xtick',[0.2 0.8],'xticklabel',{'Rest','Exercise'})