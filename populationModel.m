% Ariana Pineda, CAAM 210, SPRING 2022, Population Models
% populationModel.m
% this script helps models populations of different animals using logistic
% models and predator-prey models
% Last modified: March 7, 2022

function populationModel
%this driver plots logistic models and predator-prey models of rabbits and foxes

%Ra: plot rabbit population vs time

%initialize parameters
K = 100;
r = 1;
R0 = 20;
delta = 0.01;
T = 15;
%solve ODE using logistic model
[t,R] = logisticRabbits(R0,K,T,r,delta);

%plot rabbit population vs time
figure(1)
hold on
grid on
plot(t,R, 'LineWidth',1.8)
title('Logistic Model of Rabbits R(0)=20')
xlabel('Time, t')
ylabel('Rabbit Population, R(t)')
hold off

%Rb: plot rabbit pop vs time

%initialize parameters
K = 100;
r = 1;
R0 = 150;
delta = 0.01;
T = 15;

%solve ODE using logistic model
[t,R] = logisticRabbits(R0,K,T,r,delta);

%plot solution
figure(2)
plot(t,R, 'LineWidth', 1.8)
title('Logistic Model of Rabbits R(0)=150')
xlabel('Time, t')
ylabel('Rabbit Population, R(t)')
hold off

% given any positive initial population, the rabbit population will plateau
% at the carrying capacity


%RFa: Predator-Prey Model with a Euler Step Size of 0.01

%initialize parameters
k1 = 3;
k2 = 0.003;
k3 = 0.0006;
k4 = 0.5;
R0 = 1000;
F0 = 500;
T = 15;
delta = 0.01;

figure(3)
tiledlayout(2,2)

%solve ODE using Predator-Prey model
[t,R,F] = PredatorPreyModel(delta, T, R0, F0, k1, k2, k3, k4);
nexttile
hold on
grid on

%plot solutions
plot(t, R,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.8)
plot(t,F,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.8)
title('Predator-Prey Model with a Euler Step Size of 0.01')
xlabel('Time, t')
ylabel('Populations, R and F')
legend('Rabbit Population, R', 'Fox Population, F')
hold off

%RFa: Predator-Prey Model with a Euler Step Size of 0.001

%initialize parameters
k1 = 3;
k2 = 0.003;
k3 = 0.0006;
k4 = 0.5;
R0 = 1000;
F0 = 500;
T = 15;
delta = 0.001;

%solve ODE using Predator-Prey model
[t,R,F] = PredatorPreyModel(delta, T, R0, F0, k1, k2, k3, k4);
nexttile
hold on
grid on

%plot solution
plot(t, R,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.8)
plot(t,F,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.8)
legend('Rabbit Population, R', 'Fox Population, F')
title('Predator-Prey Model with a Euler Step Size of 0.001')
xlabel('Time, t')
ylabel('Populations, R and F')
hold off

%RFb

%solve ODE system using ode45
dPdt=@(t,P) [k1*P(1)-k2*P(1)*P(2);
    k3*P(1)*P(2)-k4*P(2)];
[t,Solutions]=ODESolver(dPdt,15,1000,500);


%RFc
%the populations of the two grow together for some time, then become
% inversely related as the population of the fox dominates the rabbit population. The population of the rabbits
% then dominate the fox population. Both populations seem to almost be
% growing and declining together. The solutions seem to closely resemble
% that using ode45. This tells me that the populations are not mutually
% beneficial. The introduction of the foxes limits the population of the
% rabbits.

%calculate dPdt
dPdt=@(t,P) [k1*P(1)-k2*P(1)*P(2);
    k3*P(1)*P(2)-k4*P(2)];

%solve ODE using ode45
[t3,Solutions] = ODESolver(dPdt,15,1000,500);
nexttile([1 2])
hold on
grid on

%plot solutions
plot(t3,Solutions(:,1),'Color','#4DBEEE', 'LineWidth', 1.8)
plot(t3,Solutions(:,2), 'Color', '#EDB120', 'LineWidth', 1.8)
xlabel('Time, t')
ylabel('Populations, R and F')
legend('Rabbit Population, R', 'Fox Population, F')
title('Predator-Prey Model ODE45')
hold off


% %RFd plot 
% I noticed that the population of rabbits and foxes are almost inversely related.
% As the population of rabbit`s increases, the population of foxes
% decreases. As in figure from RFc, the populations of the two grow
% together for some time, then become inversely related as the population
% of the fox dominates the rabbit population. The population of the rabbits
% then dominate the fox population. Both populations seem to almost be
% growing and declining together.

figure(4)

%initialize parameters
k1 = 3;
k2 = 3*exp(-3);
k3 = 6*exp(-4);
k4 = 0.5;
dPdt=@(t,P) [k1*P(1)-k2*P(1)*P(2);
    k3*P(1)*P(2)-k4*P(2)];

%solve ODE using ode45
[t,Solutions]=ODESolver(dPdt,15,1000,500);
Sol = Solutions;

hold on
%iterate over time to plot each population per time
for i=1:length(t)
    plot(Sol(:,1),Sol(:,2),'r',Sol(i,1),Sol(i,2),'ro');
    grid on

end
%titles and labels
title('Rabbit vs Fox Populations');
xlabel('Rabbit Population, R(t)','fontsize',15)
ylabel('Fox Population, F(t)','fontsize',15)
hold off

%RFe

%initialize parameters
k1 = 3;
k2 = 0.003;
k3 = 0.0006;
k4 = 0.5;
dPdt=@(t,P) [k1*P(1)-k2*P(1)*P(2);
    k3*P(1)*P(2)-k4*P(2)];

%solve ODE
[t,Solutions]=ODESolver(dPdt,15,1000,500);
Sol=Solutions;
figure(5)
PopVid = VideoWriter('C:\Users\arianapineda\Desktop\Populationvid.avi');
open(PopVid)

%iterate over time to plot each population per time
for i=1:length(t)
    plot(Sol(:,1),Sol(:,2),'r',Sol(i,1),Sol(i,2),'ro');
    grid on
    title(['Rabbit vs Fox Populations at t=',num2str(t(i))])
    xlabel('Rabbit Population, S=R(t)','fontsize',15)
    ylabel('Fox Population, F(t)','fontsize',15)
    writeVideo(PopVid, getframe(gcf));
end
%end video
close(PopVid);

end

function [t,R] = logisticRabbits(R0,K,T,r,delta)
%this function solves the ODE for logistic population growth using Euler's
%method
%inputs: R0: initial rabbit population, K: carrying capacity, T: time span,
% r: positive constant related to fecundity, delta: step size
%outputs: t:time, R:vector representing rabbit population at time t

t = 0:delta:T;
R = zeros(1,length(t));
R(1) = R0;
%iterate over time to calculate rabbit population for each time
for i=1:length(t)-1
    dR=(0.45*R(i)*(1-R(i)/K))*delta;
    R(i+1)=R(i)+dR;
end
end


function [t,R,F] = PredatorPreyModel(delta, T, R0, F0, k1, k2, k3, k4)
%this function solves the ODE using Predator-Prey population model
%inputs: delta: step size, T: time span, R0: initial rabbit population, F0:
% initial fox population, k1: constant, k2: constant, k3: constant, k4: constant
%outputs: t:time, R:vector representing rabbit population at time tk F:
%vector representing fox population at time t

%initialize variables
t=0:delta:T;
R=zeros(1,length(t));
F=zeros(1,length(t));
R(1)=R0;
F(1)=F0;

%iterate over time to calculate rabbit and fox population for each time
for i=1:length(t)-1
    dR=(k1*R(i)-k2*R(i)*F(i))*delta;
    dF=(k3*R(i)*F(i)-k4*F(i))*delta;
    R(i+1)=R(i)+dR;
    F(i+1)=F(i)+dF;
end
end


function [t,Solutions]=ODESolver(dPdt,Tfinal,ini1,ini2)
%this function solves an ODE using ode45
%Inputs: dPdt: derivative of population, Tfinal: final time, ini1: initial
%population of rabbits, ini2: initial populationof foxes
%Outputs: t:time, Solutions: a vector containing rabbit and fox populations
%at time t

options=odeset('RelTol',0.000001);
[t,Solutions]=ode45(dPdt,[0 Tfinal],[ini1 ini2],options);
end
