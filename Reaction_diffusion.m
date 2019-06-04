%% 2d reaction diffusion equations
clear
close all

% simulation parameters

N = 100;
dt = 1;
dx = 0.01;
tottaltime = 30000;
Du = 2*10^-5;
Dv = 10^-5;
F = 0.05;
k = 0.063;

%oscillating spots [0.021,0.05]
%dividing spots [0.021,0.06]
%stripes [0.0]

reactionstep = 1; % run with reaction reaction of not

%% generate adjacency matrix

adjmat = sparse(N^2,N^2);


for ii = 1:N
    for jj = 1:N
        
        linindex = sub2ind([N,N],ii,jj);
        
        up = mod(ii-1-1,N)+1; upi = sub2ind([N,N],up,jj); adjmat(upi,linindex) = 1; adjmat(linindex,upi) = 1;
        down = mod(ii+1-1,N)+1; downi = sub2ind([N,N],down,jj); adjmat(downi,linindex) = 1; adjmat(linindex,downi) = 1;
        left = mod(jj-1-1,N)+1; lefti = sub2ind([N,N],ii,left); adjmat(lefti,linindex) = 1; adjmat(linindex,lefti) = 1;
        right = mod(jj+1-1,N)+1; righti = sub2ind([N,N],ii,right); adjmat(righti,linindex) = 1; adjmat(linindex,righti) = 1;
        
        adjmat(linindex,linindex) = -4;
        
    end
end

%% initialize U and V

if reactionstep == 1
    
    U = ones(N,N);
    U = U(:)';
    
    V = zeros(N,N);
    V = V(:)';
    
else
    
    U = zeros(N,N); 
    U = U(:)';
    
    V = zeros(N,N); V(round(N/2),round(N/2)) = 1000; %introduce impulse function
    V = V(:)';
    
end

%% diffuse particles

for tt = 1:tottaltime/dt
   
    
    U = U + dt*reactionstep*(-U.*V.^2 + F*(1-U)) + Du*dt/dx^2*U*adjmat; 
    V = V + dt*reactionstep*(U.*V.^2 - (F + k)*V) + Dv*dt/dx^2*V*adjmat;
    
    if tt == 100 % add perturbation
        
        if reactionstep == 1 %do this only if you are running reactions
            
            %circular perturbation
            
            %{
            xs = repmat([1:N]-N/2',N,1);
            ys = xs';
            [ths , rs ] = cart2pol(xs,ys);
            pert = zeros(N,N);
            pert(rs<20) = 1;
            pert = logical( pert(:));
            U( pert ) = 0.5;
            V( pert ) = 0.25;
            U = U + rand(size(U))*0.5;
            V = V + rand(size(V))*0.5;
            %}
            
            
            boxsize = round(10);
            
            rowstuff = repmat([1 : boxsize]' , boxsize , 1 );
            colstuff = repmat( [1 : boxsize] , boxsize , 1  ); colstuff = colstuff(:);
            
            U( sub2ind([N,N] , rowstuff, colstuff ) ) = 0.5;
            V( sub2ind([N,N] , rowstuff, colstuff ) ) = 0.25;
            
            U = U + rand(size(U))*0.1;
            V = V + rand(size(V))*0.1;
            
            
        end
        
    end
   
   if mod(tt,10)==0
       
       if reactionstep == 1
           
           Vdisp = reshape(V',N,N);
           imagesc(Vdisp)
           %colorbar
           %caxis([0 1])
           axis square
           drawnow
           
       else
          
           subplot(1,2,1)
           Vdisp = reshape(V',N,N);
           imagesc(Vdisp)
           colorbar
           caxis([0 10])
           axis square
           drawnow
           
           subplot(1,2,2)
           xdist = [1:N]*dx;
           plot( xdist , exp( -(xdist-N/2*dx).^2/(4*Dv*tt*dt)), 'r' )
           hold on
           plot(xdist,Vdisp(:,round(N/2))/max(Vdisp(:)) , '--k');
           legend('Theoretical Conc.' , 'Simulated Conc.'  )
           ylabel('NOrmalized Conc.'); xlabel('Distance (a.u.)')
           hold off
           %ylim([0 10])
           axis square
           
       end
       
       
   end
   
end




