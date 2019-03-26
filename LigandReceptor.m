clear all
close all

%ODE solver

totaltime = 200; %time in seconds
dt = 0.01;

timeaxis = [0:dt:(totaltime-dt)]';

kon = 0.1667; %1/(uM.s)
koff = 0.001; %1/s


L = zeros(totaltime/dt,1);
R = zeros(totaltime/dt,1);
LR = zeros(totaltime/dt,1);

% set intitial conditions

L(1,1) = 1;
R(1,1) = 1;
LR(1,1) = 0;

%start simulation

for loop = 2:totaltime/dt
    
    L(loop) = L(loop-1) + dt*(-kon*L(loop-1,1)*R(loop-1,1) + koff*LR(loop-1,1));
    
    R(loop) = R(loop-1) + dt*(-kon*L(loop-1,1)*R(loop-1,1) + koff*LR(loop-1,1));
    
    LR(loop) = LR(loop-1) + dt*(kon*L(loop-1,1)*R(loop-1,1) - koff*LR(loop-1,1));
    
end

axismax = max( [L ; R ; LR] );

subplot(1,3,1)
plot(timeaxis,L);
xlabel('Time (s)');ylabel('L (\mu M)')
ylim([0 , axismax])
axis square

subplot(1,3,2)
plot(timeaxis,R);
xlabel('Time (s)');ylabel('R (\mu M)')
ylim([0 , axismax])
axis square

subplot(1,3,3)
plot(timeaxis,LR);
xlabel('Time (s)');ylabel('LR (\mu M)')
ylim([0 , axismax])
axis square