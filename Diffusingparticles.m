clear all
close all

%% diffusing particle(s)

numberofparticles = 1000;
totaltime = 10; %seconds
Dconst = 1; %um2/s
dt = 0.1;
timeaxis = [0:dt:(totaltime-dt)]';

particles = zeros( totaltime/dt , numberofparticles , 2 );

%simulate diffusion
distmoved = sqrt(4*Dconst*dt);

for loop = 2:totaltime/dt
    
    direct = rand( 1 , numberofparticles , 1  )*2*pi ;
    
    [x , y ] = pol2cart(  direct , distmoved );
    
    
    particles(loop, : , :) = particles(loop-1, : , :) + cat(3,x,y);
    
end

subplot(1,2,1)

%plot( particles(:,:,1) , particles(:,:,2)  );

subplot(1,2,2)

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

