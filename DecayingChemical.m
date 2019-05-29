% Biophysics workshop solving ODEs
clear 
close all

%% A -> B. Nuclear decay.

totaltime = 10; %total simulation time  in seconds
dt = 0.1; %time step seconds
nsteps = totaltime/dt; %simulation steps

%parameters

koff = 0.8; %A -> B

%set up matrices

AE = zeros(nsteps,1); %store solution you get by Euler's method here
AR = zeros(nsteps,1); %store solution you get by RK method here

timeaxis = (0:nsteps-1)*dt; %timeaxis for plotting purposes

%Intial conditions
A0 = 100; %in uMolar
AE(1) = A0; %in uMolar
AR(1) = A0; %in uMolar

for ii = 2:nsteps
    
    %set up and update equations
    
    dAdt = @(A) -koff*A; %A -> B
    
        AE(ii) = AE(ii-1) + Eulers_met(dAdt,dt,AE(ii-1));
        AR(ii) = AR(ii-1) + Runge_Kut(dAdt,dt,AR(ii-1));
        
    
end

% Plot how ligand and receptor evolution with time.

% figure

hold on
plot(timeaxis,AE);
plot(timeaxis,AR);
xlabel('Time (s)');ylabel('L (\mu M)')
axis square

%plot analytical solution

plot(timeaxis,A0*exp(-koff*timeaxis),'--k')

legend('Euler','RK','Analytical')

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



