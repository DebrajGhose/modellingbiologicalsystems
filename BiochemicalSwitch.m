% Biophysics workshop solving ODEs
%% A -> B. B -> A. Protein phosphorylation.

clear
close all



totaltime = 200; %total simulation time  in seconds
dt = 0.01; %time step seconds
nsteps = totaltime/dt; %simulation steps

%parameters

simtype = 3; %1 simple transduction. 2 linear feedback. 3 ultrasenstive feedback.

switch simtype
    
    case 1 %simple transduction
        
        k1 = 0.08; %A -> B
        k2 = 0.05; %B -> A
        S = 1; %Stimulus.
        I = 1;
        kf = 0; %Feedback strength.
        Km = 1;
        n = 1; %Hill Coefficient
        
    case 2 %linear feedback
        
        k1 = 0.08; %A -> B
        k2 = 0.05; %B -> A
        S = 0; %Stimulus.
        I = 1;
        kf = 5e18; %Feedback strength.
        Km = 10e18;
        n = 1; %Hill Coefficient
        
    case 3 %ultrasensitive feedback
        
        k1 = 0.08; %A -> B
        k2 = 0.05; %B -> A
        S = 0; %Stimulus.
        I = 1;
        kf = 0.5; %Feedback strength.
        Km = 1;
        n = 3; %Hill Coefficient
        
end



%first plot rate balance plot 

subplot(1,2,1)

rB = [0:0.01:1];
rA = 1 - rB;
FR = (k1*S + kf.*rB.^n./(rB.^n + Km.^n)).*rA;
BR = k2*I*rB;

plot(rB,FR,'g');
hold on
plot(rB,BR,'r');

axis square
legend('Forward reaction','Backward reaction')
xlabel('B'); ylabel('Rate (Conc/s)');


%% do simulations

for simnum = 1:10
    
    %set up matrices
    A = zeros(nsteps,1);
    B = zeros(nsteps,1);
    timeaxis = (0:nsteps-1)*dt; %timeaxis for plotting purposes
    
    %Intial conditions
    T = 1; %Total protein
    
    A(1) = 0.1*simnum;%rand(); 
    B(1) = T - A(1);
    
    for tt = 2:nsteps
        
        % provide stimulus
        
        if 0 %provide stimulus
        
        if tt*dt>100 && tt*dt<150, S = 1; 
        else, S = 0;
        end
        
        end
        
        %set up and update equations
        
        A0 = A(tt-1); B0 = B(tt-1);
        dAdt = @(A) -(k1*S + kf*B0^n/(B0^n + Km^n))*A + k2*I*B0; 
        dBdt = @(B) (k1*S + kf*B^n/(B^n + Km^n))*A0 - k2*I*B; 
        
        %use numerical integration
        
        A(tt) = A(tt-1) + Runge_Kut(dAdt,dt,A(tt-1));
        B(tt) = B(tt-1) + Runge_Kut(dBdt,dt,B(tt-1));
        
    end
    
    subplot(1,2,2)
    hold on
    %plot(timeaxis,A,'r');
    plot(timeaxis,B,'b');
    xlabel('Time (s)');ylabel('B (Conc.)')
    axis square
    
    %{
    subplot(1,2,2)
    hold on
    plot(A,B,'b')
    plot(A(end),B(end),'xr')
    axis square
    %}
    
end



%% Functions

%% Euler's method

function [df] = Eulers_met(dfdt,dt,X)

% dfdt = function for derivative
% dt = timestep
% X = thing being changed

df = dfdt(X)*dt;

end

function [df] = Runge_Kut(dfdt,dt,X) %the y is to keep track of what is changing right now

K1 = dfdt(X)*dt;
K2 = dfdt(X+0.5*K1)*dt;
K3 = dfdt(X+0.5*K2)*dt;
K4 = dfdt(X+K3)*dt;

df = (K1 + 2*K2 + 2*K3 + K4)/6;

end



