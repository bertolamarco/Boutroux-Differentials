% A Chebotarov continuum with poles E. To have a single continuum we
% require the number of stagnation points, L, to be zero. In general the
% number of connected components we obtain is L+1, as it follows from
% elementary Morse theory applied to the function $\Im \int^z \sqrt{Q}$. 
close all;
phi=1; %t_0=1, all higher ``times'' set to zero. 
L=0;
e = [0,1,-1+1i,2i,3i];
Precision = 1e-10; % The threshold of Fun
PRINT=false; %change to true if you want the picture to be saved.
Draw = true; %turn to false if you don't want to see any picture (but why would you do that?)
epsilon = 0.001; %The step used in plotting;
Length  = 4;% a proxy for the length of the trajectories that we plot at the end (not quite a length of any sort, though).
% 
% 

[S,D, e,tR]= G_functions(phi,L,e,Precision, PRINT,  Draw, epsilon, Length );

ax1 = findall(0,'type','axes');
ax2 = copyobj(ax1, figure(2));

%Now the same with 2 continua;
L=1;
[S,D, e,tR]= G_functions(phi,L,e,Precision, PRINT,  Draw, epsilon, Length );


ax1 = findall(0,'type','axes');
ax2 = copyobj(ax1, figure(3));

%Now the same with 2 continua; it likely will give a different solution (if
%not, rerun)
L=1;
[S,D, e,tR]= G_functions(phi,L,e,Precision, PRINT,  Draw, epsilon, Length );


ax1 = findall(0,'type','axes');
ax2 = copyobj(ax1, figure(4));
% 
%An example with external field 
e = [0,1,2i,3i, 2];
phi = [3i+1,0]; %The external field is Im((3i+1) z), with $t_0=1$.
L=1;
[S,D, e,tR]= G_functions(phi,L,e,Precision, PRINT,  Draw, epsilon, Length );


ax1 = findall(0,'type','axes');
ax2 = copyobj(ax1, figure(5));
