%% ME - 221 Dynamical systems, Spring 2026
% Solution to MATLAB PROBLEM SET 1
% Student Name: Haruki Morita
% SCIPER: [421834]
% Submission date: 26_03_2026

%% Initiate the script
clc; clear; close all

%% Define the parameters
l_0 = 0.11; % m 
M = 2.5; % kg 
k = 95; % N/m 
c = 45; % Ns/m 
g = 10; % m/s^2 

%% Q2 Simulation of dynamic response
% Harmonic force parameters
A = 17; % Force amplitude [N] 
w1 = 1; % rad/s 
w2 = 3; % rad/s 
w3 = 21; % rad/s 
w = [w1 w2 w3];	% rad/s

% Initial conditions
y0 = l_0 - (M*g)/k; % m 
v0 = 0; % m/s 
IC = [y0; v0];

% Define the time vector
tt = linspace(0, 10, 1000); % s 

% Simulate the three cases
opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
[~,z_1] = ode45(@(t,z)  A1_P1_ode_2026(t,z,M,k,c,A,w(1),l_0,g), tt, IC, opts);
[~,z_2] = ode45(@(t,z)  A1_P1_ode_2026(t,z,M,k,c,A,w(2),l_0,g), tt, IC, opts);
[~,z_3] = ode45(@(t,z)  A1_P1_ode_2026(t,z,M,k,c,A,w(3),l_0,g), tt, IC, opts);

% Data visualization 
figure(1);
yyaxis left
    plot(tt,z_1(:,1),'linewidth',2);
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$x$(m)','interpreter','Latex')
    set(gca,'FontSize',11)
yyaxis right
    plot(tt,z_1(:,2),'linewidth',2);
    ylabel('$\dot{x}$(m/s)','interpreter','Latex')
    set(gca,'FontSize',11)
set(gcf,'color','w','units','centimeters','outerposition',[5 5 15 10])
exportgraphics(gcf,'Figure_Q2_1.png','Resolution',300);    

figure(2);
yyaxis left
    plot(tt,z_2(:,1),'linewidth',2);
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$x$(m)','interpreter','Latex')
    set(gca,'FontSize',11)
yyaxis right
    plot(tt,z_2(:,2),'linewidth',2);
    ylabel('$\dot{x}$(m/s)','interpreter','Latex')
    set(gca,'FontSize',11)
set(gcf,'color','w','units','centimeters','outerposition',[5 5 15 10])
exportgraphics(gcf,'Figure_Q2_2.png','Resolution',300);  

figure(3);
yyaxis left
    plot(tt,z_3(:,1),'linewidth',2);
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$x$(m)','interpreter','Latex')
    set(gca,'FontSize',11)
yyaxis right
    plot(tt,z_3(:,2),'linewidth',2);
    ylabel('$\dot{x}$(m/s)','interpreter','Latex')
    set(gca,'FontSize',11)
set(gcf,'color','w','units','centimeters','outerposition',[5 5 15 10])
exportgraphics(gcf,'Figure_Q2_3.png','Resolution',300);  

%% Q4 Plot the rail model
L0 = 0.32; % m 
l_01 = 0.05; % m 
l_02 = 0.23; % m 
l_03 = 0.95; % m 
l_0V = [l_01 l_02 l_03];  
alpha = 0.23; % m 
beta = 0.42; % 1/m 

xx = linspace(-4*pi, 4*pi, 1000); % m 
yy = alpha.*(cos(beta.*xx)-1); % m 

for i = 1:3
    l0 = l_0V(i);
    V = (1/2).*k.*((-1).*l0+L0+alpha.*((-1)+cos(beta.*xx))).^2+alpha.*g.*M.*(( ...
          -1)+cos(beta.*xx));
    dVf = @(x) alpha.*beta.*((-1).*g.*M+k.*(alpha+l0+(-1).*L0+(-1).*alpha.*cos( ...
          beta.*x))).*sin(beta.*x);
    ddVf = @(x) alpha.*beta.^2.*((-1).*g.*M.*cos(beta.*x)+k.*cos(beta.*x).*( ...
          alpha+l0+(-1).*L0+(-1).*alpha.*cos(beta.*x))+alpha.*k.* ...
          sin(beta.*x).^2);

    %Finding the stable equilibrium points
    S = csapi(xx,dVf(xx));      
    z = fnzeros(S);             
    temp = ~(z(1,:)-z(2,:));    
    roots = z(1,temp).';        
    roots(ddVf(roots)<=0) = [];  

    % Print to command window
    disp("Equilibrium points for l_0 = "+l0+"m:")
    disp(roots)
  
    figure;clf
    yyaxis left
    plot(xx,yy,'linewidth',1.5);hold all 
    for jj = 1:length(roots)
        xline(roots(jj),'--k')
    end
        xlabel('$x$(m)','interpreter','Latex')
        ylabel('$y$(m)','interpreter','Latex')
        set(gca,'fontsize',8)
        xlim([xx(1) xx(end)])
    yyaxis right
    plot(xx,V,'linewidth',1.5);hold all 
        xlabel('$x$(m)','interpreter','Latex')
        ylabel('$V$(J)','interpreter','Latex')
        set(gca,'fontsize',8)
        set(gcf,'color','w','units','centimeters','outerposition',[5 5 12 8])
        xlim([xx(1) xx(end)])
    
    exportgraphics(gcf,strcat('Figure_Q4_',num2str(i),'.png'),'Resolution',300);  

end

%% Q5 Simulation
% Harmonic force parameters
l0 = 0.23; % m 

V = (1/2).*k.*((-1).*l0+L0+alpha.*((-1)+cos(beta.*xx))).^2+alpha.*g.*M.*(( ...
      -1)+cos(beta.*xx));
dVf = @(x) alpha.*beta.*((-1).*g.*M+k.*(alpha+l0+(-1).*L0+(-1).*alpha.*cos( ...
      beta.*x))).*sin(beta.*x);
ddVf = @(x) alpha.*beta.^2.*((-1).*g.*M.*cos(beta.*x)+k.*cos(beta.*x).*( ...
      alpha+l0+(-1).*L0+(-1).*alpha.*cos(beta.*x))+alpha.*k.* ...
      sin(beta.*x).^2);

% Finding the stable equilibrium points
S = csapi(xx,dVf(xx));      
z = fnzeros(S);             
temp = ~(z(1,:)-z(2,:));    
roots = z(1,temp).';        
roots(ddVf(roots)<=0) = [];  

% Print to command window
disp("Equilibrium points for l_0 = "+l0+"m:")
disp(roots)

% Initial conditions for both simulations    
ICM = [-0.05, 0; 7.5, 1.6]; 

%% Initial conditions (1)
i=1;
IC = ICM(i,:)';
tt = linspace(0, 400, 2000);  % s 
opts = odeset('RelTol',1e-4,'AbsTol',1e-6); 

% Simulate
[~,z_1] = ode45(@(t,z) A1_P2_ode_2026(t,z,M,k,c,L0,l0,alpha,beta), tt, IC, opts); 

% Plot position and velocity of particle over time.
figure(7); clf
subplot(1,2,1)
yyaxis left
    plot(tt,z_1(:,1),'linewidth',1.5); 
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$x$(m)','interpreter','Latex')
    for jj = 1:length(roots)
        yline(roots(jj),'--k')
    end
    set(gca,'FontSize',8)
    ylim([0.9*min(z_1(:,1)) 1.1*max(z_1(:,1))])
yyaxis right
    plot(tt,z_1(:,2),'linewidth',1.5); 
    ylim([1.1*min(z_1(:,2)) 1.1*max(z_1(:,2))])
    ylabel('$\dot{x}$(m/s)','interpreter','Latex')
    set(gca,'FontSize',8)
set(gcf,'color','w','units','centimeters','outerposition',[5 5 12*2 8])

%Plot particle position on rail and potential energy in space.
figure(7);
subplot(1,2,2)
yyaxis left
plot(xx,yy,'linewidth',1.5);hold all
COL = cool(size(z_1,1));
for jj = 1:length(roots)
    xline(roots(jj),'--k')
end
for kk = 1:size(z_1,1)/50:size(z_1,1)
    plot(z_1(kk,1),alpha*(cos(beta*z_1(kk,1))-1),'o','color',COL(kk,:),'MarkerFaceColor',COL(kk,:));
end
    xlim([xx(1) xx(end)])
    xlabel('${x}$(m)','interpreter','Latex')
    ylabel('$y$(m)','interpreter','Latex')
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','outerposition',[5 5 12*2 8])
yyaxis right
    plot(xx,V,'linewidth',2);hold all 
    xlabel('$x$(m)','interpreter','Latex')
    ylabel('$V$(J)','interpreter','Latex')
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','outerposition',[5 5 12*2 8])
    xlim([xx(1) xx(end)]) 

exportgraphics(gcf,strcat('Figure_Q5_1.png'),'Resolution',300); 
        
%% Initial conditions (2)
i=2;
IC = ICM(i,:)';
opts = odeset('RelTol',1e-4,'AbsTol',1e-6); 

% Simulate
[~,z_1] = ode45(@(t,z) A1_P2_ode_2026(t,z,M,k,c,L0,l0,alpha,beta), tt, IC, opts); 

% Plot position and velocity of particle over time.
figure(8); clf
subplot(1,2,1)
yyaxis left
    plot(tt,z_1(:,1),'linewidth',1.5); 
    xlabel('$t$(s)','interpreter','Latex')
    ylabel('$x$(m)','interpreter','Latex')
    for jj = 1:length(roots)
        yline(roots(jj),'--k')
    end
    set(gca,'FontSize',8)
    ylim([0.9*min(z_1(:,1)) 1.1*max(z_1(:,1))])
    set(gca,'FontSize',8)
yyaxis right
    plot(tt,z_1(:,2),'linewidth',1.5); 
    ylim([1.1*min(z_1(:,2)) 1.1*max(z_1(:,2))])
    ylabel('$\dot{x}$(m/s)','interpreter','Latex')
    set(gca,'FontSize',8)
set(gcf,'color','w','units','centimeters','outerposition',[5 5 12*2 8])

%Plot particle position on rail and potential energy in space.
figure(8);
subplot(1,2,2)

yyaxis left
plot(xx,yy,'linewidth',1.5);hold all 
COL = cool(size(z_1,1));
for jj = 1:length(roots)
    xline(roots(jj),'--k')
end
for kk = 1:size(z_1,1)/50:size(z_1,1)
    plot(z_1(kk,1),alpha*(cos(beta*z_1(kk,1))-1),'o','color',COL(kk,:),'MarkerFaceColor',COL(kk,:));
end
    xlim([xx(1) xx(end)])
    xlabel('${x}$(m)','interpreter','Latex')
    ylabel('$y$(m)','interpreter','Latex')
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','outerposition',[5 5 12*2 8])

yyaxis right
    plot(xx,V,'linewidth',2);hold all 
    xlabel('$x$(m)','interpreter','Latex')
    ylabel('$V$(J)','interpreter','Latex')
    set(gca,'fontsize',8)
    set(gcf,'color','w','units','centimeters','outerposition',[5 5 12*2 8])
    xlim([xx(1) xx(end)]) 
exportgraphics(gcf,strcat('Figure_Q5_2.png'),'Resolution',300);