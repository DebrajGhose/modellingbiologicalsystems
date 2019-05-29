% Biophysics workshop solving ODEs
clear 
%close all

%% L + R -> LR. LR -> L + R. Ligand-receptor interaction.

totaltime = 50; %total simulation time  in seconds
dt = 0.1; %time step seconds
nsteps = totaltime/dt; %simulation steps

%parameters

kon = 0.008; %L + R -> LR
koff = 0.005; %LR -> L + R


%set up matrices
L = zeros(nsteps,1);
R = zeros(nsteps,1);
LR = zeros(nsteps,1);
timeaxis = (0:nsteps-1)*dt; %timeaxis for plotting purposes

%Intial conditions


LR(1) = 0;
L(1) = 100; %in uMolar
R(1) = 70;


for tt = 2:nsteps
    
    %set up and update equations
    
    L0 = L(tt-1); R0 = R(tt-1); LR0 = LR(tt-1);
    dLdt = @(L) -kon*L*R0 + koff*LR0; %L + R -> LR
    dRdt = @(R) -kon*L0*R + koff*LR0; %L + R -> LR
    dLRdt = @(LR) kon*L0*R0 - koff*LR; %LR -> L + R
    
    %use numerical integration
    
    L(tt) = L(tt-1) + Runge_Kut(dLdt,dt,L(tt-1));
    R(tt) = R(tt-1) + Runge_Kut(dRdt,dt,R(tt-1));
    LR(tt) = LR(tt-1) + Runge_Kut(dLRdt,dt,LR(tt-1));
    
end

%figure

hold on
plot(timeaxis,L);
plot(timeaxis,R);
plot(timeaxis,LR);
xlabel('Time (s)');ylabel('Concentration (\mu M)')
legend('L','R','LR');
axis square



%% Functions

%% Euler's method

function [df] = Eulers_met(dfdt,dt,X)

% dfdt = function for derivative
% dt = timestep
% X = thing being changed

df = dfdt(X)*dt;

end

function [df] = Runge_Kut(dfdt,dt,X) %the y is to keep track of what is changing right now

k1 = dfdt(X)*dt;
k2 = dfdt(X+0.5*k1)*dt;
k3 = dfdt(X+0.5*k2)*dt;
k4 = dfdt(X+k3)*dt;

df = (k1 + 2*k2 + 2*k3 + k4)/6;

end



