clear all
close all

%% diffusing particle(s)

numberofparticles = 10000;
totaltime = 10; %seconds
Dconst = 1; %um2/s
dt = 0.1;
timeaxis = [0:dt:(totaltime-dt)]';

particles = zeros( totaltime/dt , numberofparticles , 2 ); %row is time, column is particle, depth is x/y axes

%simulate diffusion
distmoved = sqrt(4*Dconst*dt);

for loop = 2:totaltime/dt
    
    direct = rand( 1 , numberofparticles , 1  )*2*pi ;
    
    [x , y ] = pol2cart(  direct , distmoved );
    
    particles(loop, : , :) = particles(loop-1, : , :) + cat(3,x,y);
    
end

subplot(1,3,1)

plot( particles(:,1:10,1) , particles(:,1:10,2)  );

axis equal
axis square

subplot(1,3,2)

msds = mean(( particles(:,:,1).^2 + particles(:,:,2).^2 ),2);

plot(timeaxis,msds);
xlabel('Time (s)')
ylabel('MSD (\mu m^2)')

%do a fit

para = polyfit(timeaxis,msds,1);
hold on
plot(timeaxis, para(1)*timeaxis + para(2));

legend('Data','Fit')

text( 0.1,0.7 ,[ 'D = ' num2str( para(1)/4  )] ,'Units','normalized');

axis square

subplot(1,3,3)

inspecttime = 6; %inspect particle distribution at this time

particlesalongaxis = particles(inspecttime/dt,:,1);

%histogram(particlesalongaxis);
histfit(particlesalongaxis(:));

disp(['Mean distance traveled=', num2str(mean(particlesalongaxis(:)))]);
disp(['Standard deviation=', num2str(std(particlesalongaxis(:)))]);

%pd = fitdist(particlesalongaxis(:),'Normal')

sqrt(2*Dconst*inspecttime)

axis square
