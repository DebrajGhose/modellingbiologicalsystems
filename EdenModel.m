%% Eden model
clear all
close all

N = 1000; %size of domain -- make this an even number
agar = zeros(N,N); %
agar(N/2,N/2) = 1; % 1 - red. 0 - none. 2 - white.
mutat = 0.01;
boundary = [];

%create initial boundary around seed

boundary = adjustb(boundary,agar,[N/2 N/2]);

for ii = 1:10000
   
    [agar,boundary] = growme(agar,boundary,mutat);
    
    if mod(ii,100)==0
        
        imagesc(agar)
        %plot(boundary(:,1),boundary(:,2),'.')
        axis square
        drawnow
        
    end
end

function [agar,boundary] = growme(agar,boundary,mutat)

%choose random cell in boundary to fill

randnum = randi( [ 1 , size(boundary,1) ] , 1 , 1 ) ; %pick a location to grow from
cell2fill = boundary( randnum , : );
boundary( randnum , : ) = [];

if rand() < mutat %small chance of mutating stuff
    agar(cell2fill(1),cell2fill(2)) = 2;
else
    agar(cell2fill(1),cell2fill(2)) = 1;
end % small chance of causing a mutation

%readjust boundary

boundary = adjustb(boundary,agar,[cell2fill(1),cell2fill(2)]);

end


function [boundary] = adjustb(boundary,agar,poi)

for ii = [1 0 -1]
    
    for jj = [1 0 -1]
        
        if ~( ii == 0 && jj == 0 )
            
            if agar(poi(1)+ ii,poi(2)+ jj) < 1 %check for empty spaces to extend boundary into
                
                if isempty(boundary)  %if empty then ignore ismember command to avoid errors
                    
                    
                    boundary = [boundary ; [poi(1) + ii , poi(2) + jj ]   ];
                    
                else
                    
                    
                    if ~ismember( [poi(1)+ ii , poi(2)+ jj]  , boundary , 'rows'  ) %make sure boundary value dose not already exist
                        
                        boundary = [boundary ; [poi(1) + ii , poi(2) + jj ]   ];
                        
                        
                    end
                    
                end
                
            end
        end
    end 
end

end
