function [S,D, e,tR]= G_functions(phi,L,e,Precision, PRINT,  Draw, epsilon, Length )
%[S,D, e,tR]= G_functions(phi,L,e,Precision, PRINT,  Draw, epsilon, Length )
%
% Finds the quadratic differential $Q= S^2 Delta/E$ where E = prod (z-e),
% deg S = L, deg Delta = M and \sqrt{Q} \simeq phi' + t_0/z + ..., where
% sqrt{Q} has purely imaginary periods. 
% The quadratic differential is Q = S^2 D/E where E is the poly with roots
% at e and deg S =L, deg D =M, deg E = N, N  = 2L+M+2;
% The argument phi=[t_R, ..., t_1,t_0]; t_0 is the residue of
% \sqrt Q. 
% e=[list of points]
% L= degree of L. If the passed L is too large it is defaulted to L=0, with
% a complaint. 
%
%PRINT=true/false; if true (default false), it will save the picture in the
%same directory and also a small .tex wrapper with the same filename
%
%Precision= threshold of the functional  below which we stop the
%computation (Default = 1e-8)
%
%Draw, epsilon, Length are parameters for the external function that draws
%the critical trajectories. They don't affect the data output. 
%
%OUTPUT:  S= the polynomial S;D = the polynomial Delta; e= the same as the
% one passed; tR= phi(1) (leading coefficient). If needed one can change
% the code to return the polynomials S, D(elta) themselves. 
%
%DEFAULTS: if phi is not defined or empty it defaults to phi=1 (i.e. t_0=1;
%Chebotarov case); if e is not defined or empty, it defaults to N=randi(10)+1 random points;
%if L is not defined it defaults to L=0; DRAW=true, epsilon = 0.004; Length=20

% Written for Matlab R2024a but possibly working for older versions as
% well. 
% 
%This package is the accompanying code to https://www.arxiv.org/abs/2408.15234
%
% August 2024
% (c) Marco Bertola
% email: Marco.Bertola@concordia.ca


%%% Input validation %%%
if ~exist('phi','var') || size(phi,1)==0
    phi = 1;
end
if ~exist('PRINT','var')
    PRINT = 0;
end
t0=phi(end);
tR=phi(1);
R = size(phi,2)-1;
if ~exist('e','var') || size(e,1)==0 
    N=randi(10)+1;
    MM= randn(N) + 1i*randn(N);
    e = eig(MM).';
    %e = exp(2i*pi*(2/3:2/3:2));
    %e = exp(2i*pi*(2/3:1/3:2));
     %e=exp(1i*pi*(1/2:1/2:2));
    %e= [1,2,1i,3i, 1+1i];
    
end
if ~exist('Draw','var')
    Draw= true;
end
if ~exist('L','var')
    L=0;
    %No stagnation points, the standard Chebotarov problem
end
if ~exist ('epsilon','var')
    epsilon =0.004;
end
if ~exist ('Length','var')
    Length =20;
end
N = size(e,2);
M = N-2*L-2+2*R;
g = (M+N)/2-1;
if g<0  || M< 0
    disp('Error: there are too many stagnation points');
    L = 0 ;
    M = N-2*L-2+2*R;
    %g = (M+N)/2-1;
end
%%% End Input Validation %%%


% M is The number of simple zeros (generically) of the quadratic differential
% L is the number of stagnation points. 

S = phi(1)*poly(randn(1,L) + 1i*randn(1,L));
D = poly(randn(1,M)  + 1i*randn(1,M));

% The quadratic differential is Q = D/E dz^2 and we start with a D such
% that the roots are in the convex hull; by a theorem , the derivative does
% that.
% 
% We will choose the cuts ....
droot = roots(D).';
sroot = roots(S).';


ITER = 1;
MAXITER = 1000;
Fun =3;
Funold=Fun;
if ~exist('Precision','var')
    Precision = 1e-8;
end

if Draw
    figure(1);
    clf;
    ZZ = cat(2,e,droot, sroot);
    dd=plot(droot,'.k','MarkerSize',10);hold on;
    ee=plot(real(e),imag(e),'.r','MarkerSize',20);
    ss=plot(real(sroot),imag(sroot), '.g', 'MarkerSize',12);
    axis ([min(real(ZZ))-1,max(real(ZZ))+1,min(imag(ZZ))-1,max(imag(ZZ))+1]); hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bounce=0;
Sign =1;
t_ells = polyder(phi);
t_ells = cat(2,t_ells,t0);
t_ells = t_ells(2:end);
while( (Fun>Precision && ITER<MAXITER))
    M = size(droot,2);
    %g = (N+M)/2-1;
    [Inty, P]  = Cycles(e,droot, sroot,phi(1));
    % Inty(1:R) contains the quantities denoted $T_\ell[Q]$ in the paper.
    % $\ell=0,...,R-1$;
    % Inty(R+1:end) contains the quantities denoted $P_j[Q]$ (the periods
    % of the square root of the quadratic differential).
    AA = P(1:R,:);
    BB = P(R+1:end,:);
    PP = [ imag(AA), real(AA); real(AA), -imag(AA); real(BB), -imag(BB)];
    % PP is the matrix of real/imag parts of the Holomorphic periods as
    % well as differentials z^{g-1+k}/sqrt(DE), k=1,... 2*L)
    V1 = Inty(1:R) - t_ells.';
    V2 = Inty(R+1:end);
    V = [imag(V1);real(V1); real(V2)];
    dotNum_re_im = PP^(-1)*V;
    dotNum = dotNum_re_im(1:end/2) + 1i*dotNum_re_im(end/2+1:end);
    dotNum =-Sign*cat(2,0, dotNum.'); 
    %dotNum is the numerator in 
    %\dot \sqrt{Q} = (dot S D+ 1/2 \dot D S)/sqrt(D E);
    MM = mysylvester(S,D);
    T1 =dotNum(2:end)*MM^(-1);
    dotS = T1(end-L+1:end);
    dotD = 2*T1(1:M);
    dotS = cat(2,0,dotS);
    dotD = cat(2,0,dotD);

    c1 = sum(abs(dotD))/sum(abs(D)) + sum(abs(dotS))/sum(abs(S));
    dt = (1/2 + c1/(1+c1)/4);

    D = D + dt*dotD;
    S = S + dt*dotS;
    droot = roots(D).';
    sroot = roots(S).';

    Fun =sqrt( sum(abs(V1).^2) + sum(real(V2).^2));
    %The functional we use is the square root of the one in the paper. 
    if Funold<Fun
        Bounce=Bounce+1;
    else
        Bounce=0;
    end
    if Bounce>10
        Sign=-Sign;
        Bounce=0;
    end
    Funold=Fun;
    ITER = ITER+1;

    if Draw
        delete(dd);
        delete(ee);
        delete(ss);
        dd=plot(real(cat(2,droot)),imag(cat(2,droot)) ,'.k','MarkerSize',10);hold on;
        ss=plot(real(cat(2,sroot)),imag(cat(2,sroot)) ,'.g','MarkerSize',10);
        ee=plot(real(e),imag(e),'.r','MarkerSize',10);
        title(['Fun=',num2str(Fun), '. Iteration=', num2str(ITER)]);
        drawnow;
    else
        disp(['L^2 norm of imaginary periods=',num2str(Fun),'. Iteration=',num2str(ITER)]);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% END MAIN LOOP %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(ITER==MAXITER)
    disp(['Max iter! Tolerance=',num2str(Fun)]);
    if Fun>1
        stop;
    end
end
if Draw
    figure(1); clf;
    ZZ = cat(2,sroot, droot);
    muld = ones(size(droot));
    muls = ones(size(sroot));
    mul=cat(2,muls*2,muld);
    CriticalTrajectories_Rational(-1,ZZ,mul,...
        e,ones(size(e)),epsilon,Length);
    win = [min(real(cat(2,ZZ,e)))-2, max(real(cat(2,ZZ,e)))+2, min(imag(cat(2,ZZ,e)))-1, max(imag(cat(2,ZZ,e)))+1];
    plot(real(e),imag(e), 'r.','MarkerSize',24);
    plot(real(droot),imag(droot), 'k.','MarkerSize',18);
    
    axis square;
    axis equal; 
    axis(win);

    grid on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exit preparations and printout %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = N-2*L-2+2*R;
g = (M+N)/2-1;
disp(['The final genus is g=',num2str(g)]);
disp(['We converged in ',num2str(ITER),' iterations and we have Fun<=',num2str(Fun)]);

if PRINT==true
    if R>0
        Caption=['In this case $\\Phi = \\sum_{\\ell=1}^',num2str(R),...
            ' t_\\ell z^\\ell$ with $[t_',num2str(R),',\\dots,t_1]=[',...
            regexprep(num2str(phi(1:end-1)),'\s+',',')...
            ,']$. Here  $t_0=',num2str(t0),...
            '$. The set $E=\\{',...
            regexprep(num2str(e),'\s+',','),...
            '\\}$ and $L=',num2str(L),'$.'];
    else
        Caption=['Chebotarov continuum for the set $E=\\{',...
            regexprep(num2str(e),'\s+',','),'\\}$ with $L=',num2str(L),'$.'];
    end


    symbols = ['a':'z' 'A':'Z' '0':'9'];
    stLength = 10;
    nums = randi(numel(symbols),[1 stLength]);
    st = symbols (nums);

    %name = ['/Users/bertola/Library/CloudStorage/OneDrive-ConcordiaUniversity-Canada/Noether/Research in progress/Numerics_Chebotarov/Example_',st];

    name = ['Example_',st];

    TexInclude = [ ...
        '\\begin{minipage}{0.5\\textwidth} \n',...
        '\\includegraphics[width=1\\textwidth]{',name,'.pdf} \n',...
        '\\captionof{figure}{\\small ',Caption,'} \n',...
        '\\end{minipage} '];


    title('');
    exportgraphics(gcf,cat(2,name,'.pdf'),'ContentType','vector');
    FileTex = fopen(cat(2,name,'.tex'),'w');
    fprintf(FileTex, TexInclude);
    fclose(FileTex);
    %
end
end

%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Inty,PowerPeriods] = Cycles (e,droot,sroot,tR)
    %
    % Computes the integral of y dz = sqrt(Q)  along the cycles as well as the holomorphic ones
    
    N=size(e,2);
    L=size(sroot,2);
    D = poly(droot);
    M= size(droot,2);
    R = (M+2*L+2-N)/2;
    S = tR*poly(sroot);
    Bpts = cat(2,e,droot);
    %%% DEBUG
    % DEBUG=0;
    % if DEBUG
    %     Bpts = (1:N+M);
    %     Bpts(2:2:end) = Bpts(2:2:end)+2i;
    % end
    %%%
    
    Cuts = BranchCuts(Bpts);
    %Cuts = [Bpts(1:2:end).', Bpts(2:2:end).'];
    % Rad is the radical of the Riemann surface; y = T/R dz
    
    ttt1 = conv(S, D);
    T = @(z) polyval (ttt1,z);
    H = (N+M)/2;
    t = Cuts(:,2)-Cuts(:,1);
    t = t./abs(t); %Tangent vectors to cuts
    tt = prod(-conj(t.^(-1)));
    % Rad = @(z)  prod(...
    %      -conj(t.^(-1))*ones(size(z)).*...
    %      sqrt(-conj(t)*ones(size(z)).*( ones(H,1)*z-Cuts(:,1)*ones(size(z)) ) )...
    %     .*sqrt(-conj(t)*ones(size(z)).*( ones(H,1)*z-Cuts(:,2)*ones(size(z)) ) )...
    %     ) ;
    
    Rad = @(z) tt* prod(...
        sqrt(-conj(t)*ones(size(z)).*( ones(H,1)*z-Cuts(:,1)*ones(size(z)) ) )...
        .*sqrt(-conj(t)*ones(size(z)).*( ones(H,1)*z-Cuts(:,2)*ones(size(z)) ) )...
        ) ;
    if H==1 %it means we are in genus 0 and the code above does not work
        Rad= @(z) tt*sqrt(-conj(t)*(z-Cuts(1,1))).*sqrt(-conj(t)*(z-Cuts(1,2)));
    end
    y= @(z)  T(z)./Rad(z);
    r=0.001;
    g= (N+M)/2-1;
    
    % clf;
    % figure(1);
    %
    % for j=1:g+1
    %     plot(Cuts(j,:),'r');hold on;
    % end
    
    Path=cell(2*g+R,1);
    % each path has a  component on the first and on the second sheet (possibly none of the latter)
    %Don't put any endpoint on a cut!
    KK=40;
    for j=1:R
        Path{j}   =KK*[ 1+1i, -1+1i, -1-1i, 1-1i, 1+1i];
    end
    for j=1:g
        %Path{j}   =[Cuts(j,1)-r*1i*t(j)*exp(1i*pi*( (0:-1/5:-1))), Cuts(j,2)-r*1i*t(j)*exp(1i*pi*( (1:-1/5:0))), Cuts(j,1)- r*1i*t(j)];
        Path{j+R}   =[Cuts(j,1)+ r*1i*t(j), ...
            Cuts(j,2)+r*1i*t(j)*exp(1i*pi*( (2:-1/5:0))), Cuts(j,1)+ r*1i*t(j)*exp(1i*pi*( (2:-1/5:0)))];
        %A-cycles
    end
    
    t = Cuts(2:end,1) - Cuts(1:end-1,2);
    t = t./abs(t);
    for j=1:g
        Path{g+R+j}   =[ Cuts(j,2)+ r*t(j),Cuts(j+1,1)-r*t(j)*exp(1i*pi*( (0:-2/7:-2))),...
            Cuts(j,2)+r*t(j)*exp(1i*pi*( (2:-2/7:0)))];
        % Path{g+j}   =[ Cuts(j,2)+ r*1i*t(j),Cuts(j+1,1)+ r*1i*t(j)*exp(1i*pi*( (2:-1/5:0))),...
        %      Cuts(j,2)+r*1i*t(j)*exp(1i*pi*( (2:-1/5:0)))];
        %pseudo B-cycles; not respecting the intersection number;
    end
    
    
    % %%%%%%%%%% DEBUG
    % if DEBUG
    %     clf;
    %     plot(real(Bpts), imag(Bpts),'r.','MarkerSize',12);
    %     hold on
    %     % tt = 0:0.01:4;
    %     % figure(2);
    %     % plot(tt, real(Rad(tt+0i)));
    %     % grid on;
    %     % figure(1);
    %
    %     for j=1:g+1
    %         CC = [Cuts(j,1), Cuts(j,2)];
    %         plot(real(CC), imag(CC));
    %     end
    %
    %     for j=1:g
    %         plot(Path{j},'r');
    %     end
    %     for j=g+1:2*g
    %         plot(Path{j},'b');
    %     end
    %     axis equal;
    %     grid on;
    % end
    % %%%%%%%%%%
    PowerPeriods= zeros(2*g+R,g+R);
    Inty=zeros(2*g+R,1); %
    for k=1:g+R
        F = @(z) z.^(g+R-k)./Rad(z);
        % first the residues at infinity...
        for j=1:R
            PowerPeriods(j, k) = Hyperelliptic_Integral( @(z)F(z)./z.^(R-j), Path{j}, Cuts)/2i/pi;
        end
        % then the periods
        for j=R+1:2*g+R
            PowerPeriods(j, k) =  Hyperelliptic_Integral (@(z) F(z),Path{j}, Cuts);
        end
    end
    for j=1:R
        Inty(j) =  Hyperelliptic_Integral (@(z) y(z)./z.^(R-j) ,Path{j}, Cuts)/2i/pi;
    end
    for j=R+1:2*g+R
        Inty(j) =  Hyperelliptic_Integral (@(z) y(z),Path{j}, Cuts);
    end
    % %The B--periods as of now are  not correct because we choose the A-period
    % %as the cycle around z_{2j-1}, z_{2j} and B around z_{2j}, z_{2j+1}. Thus
    % %the B_j period as of now is B_j - B_{j+1}.
    % %To rectify the situation  we take the current B and multiply appropriately
    % %the A periods.
    % A = PowerPeriods(end-2*g+1:end-g,end-g+1:end);
    % G = ones(g);
    % G = triu(G);
    % B = G*PowerPeriods(end-g+1:end, end-g+1:end) ;
    % tau = B*A^(-1);
    % if sum(sum(abs(tau-tau.')))>1e-3
    %     disp("oufh: instability detected in the tau matrix, it is not exactly symmetric. Too large genus?");
    %     % tau should be symmetric; if it is not it is a bad day...
    % end
end
    
function Cuts= BranchCuts( Q )
    % Given a set of points $Q$ in the plane, finds non intersecting segments
    % between pairs. We take the convex hull, remove all the sides, repeat.
    Q0=Q;
    Cuts = zeros(size(Q0,2)/2,2);
    count=1;
    while(size(Q0,2))>2
        K = convhull ([real(Q0.'), imag(Q0.')]);
        ct= floor((size(K,1)-1)/2);
        Cuts(count:count+ct-1,:) = [Q0(K(1:2:2*ct)).', Q0(K(2:2:2*ct)).'];
        Q0(K(1:2*ct))=[];
        count=count+ct;
    end
    if size(Q0,2)==2
        Cuts(count,:) = [Q0(1), Q0(2)];
    end
end

function Syl= mysylvester(S,D)
 % Given the two polynomials S, D it returns the Sylvester matrix. No
 % sanity check is done. YOLO! NOTE: if S,D are constants, it returns a 1x1
 % matrix and not an empty matrix. But why would you compute the Sylvester
 % matrix of two constants?

 s = size(S,2);
 d = size(D,2);
 
 %Syl = zeros(s+d-2,s+d-2);%
 if s==1
     Syl=toeplitz([S(1);zeros(d-2,1)], [S,zeros(1,d-2)]);
 else
     if d==1
         Syl=toeplitz([D(1);zeros(s-2,1)], [D,zeros(1,s-2)]);
     else
         Syl=[toeplitz([S(1);zeros(d-2,1)], [S,zeros(1,d-2)]);
             toeplitz([D(1);zeros(s-2,1)], [D,zeros(1,s-2)])];
     end
 end

end