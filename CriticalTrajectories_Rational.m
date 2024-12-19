function [CritTraj,Adjacency] = CriticalTrajectories_Rational(varargin)
%% CriticalTrajectories_Rational(Const, Zeros, mu, Poles,nu,epsilon,L, varargin)
%
% August 2024
% (c) Marco Bertola
% email: Marco.Bertola@concordia.ca
%% Usage 
%% Crittraj (Const, Z, mu, Poles, nu, epsilon,Length,interactive,Points)
%% Plots the critical trajectories of the quadratic differential
%% R(x) dx^2 where R(x) is a rational function.
%% Usage Crittraj(Const, Zeros, mu, Poles, nu);
%% Zeros is a vector of complex numbers containing the position of the zeroes;\n\
%% mu is a vector of integers containing the multiplicities of the zeroes
%%  Poles, nu are  similar (can be empty).
%% Const is the overall multiplicative constant that multiplies 
%% the rational function with those zeroes and poles.
%% epsilon is the step for the numerical integration
%% Length is the length of the arc drawn
%% interactive is a flag (0 or 1); if 1 then it prompts for points where to draw horizontal trajectories.
%% The vector Points contains points through which we want to draw
%% additional horizontal trajectories.
%% The parameter Thickness is used to scale the thickness of the lines by
%% this common factor.
%% The function returns the adjacency matrix of the critical points.

if ~exist('interactive', 'var')
    interactive=0;
end
p = inputParser;
% if ~exist('Const','var')
%     %Const=-1;
% end

clf ;
CritTraj=cell(1,100); % We assume there won't be more than 100 critical trajectories! Change as needed.
countcrit=0;

DefaultColorCrit= 'b';
un=get(gca,'Units');
set(gca,'Units','pixels');
tmp=get(gca,'Position');
set(gca,'Units',un);
DefaultThick = tmp(4)/280;

addOptional(p,'Const',-1);
addOptional(p,'Zeros',[1i,-1i]);
addOptional(p,'mu',[2,2]);
addOptional(p,'Poles',[-1,1]);
addOptional(p,'nu',[1,1]);
addOptional(p,'epsilon',0.04);
addOptional(p,'L',4);
addParameter(p,"Interactive", 0);
addParameter(p,"ColorCritical",DefaultColorCrit);
addParameter(p,"Points",[]);
addParameter(p,"Draw",'true');
addParameter(p,"Thick",DefaultThick);
parse(p, varargin{:});
Const=p.Results.Const;
Zeros= p.Results.Zeros;
Poles= p.Results.Poles;
mu=p.Results.mu;
nu=p.Results.nu;
Points = p.Results.Points;
L=p.Results.L;
ColorCrit=p.Results.ColorCritical;
Draw=p.Results.Draw;
Thick=p.Results.Thick;
epsilon=p.Results.epsilon;
interactive=p.Results.Interactive;

%global Den Num  Zeros_Poles ;
CurrAxes = gca;
Length=L;


%%%%%%
%This is just really a front-end for the external function QTrajectory,
%which actually computes the horizonal trajectroy of a quadratic
%differential, but does none of the plotting or interface.

%Adjust multiplicities: if two roots are closer than epsilon to each other
%we merge them adding the corresponding multiplicities.


tmp0=length(Zeros);
if tmp0>0
    tmp = ones(tmp0,1)*Zeros- Zeros.'*ones(1, tmp0);
    tmp=abs(tmp) + eye(tmp0);
    mindist= min(min(tmp));
    Colors= zeros(length(Zeros),3);
    for k=1:length(Zeros)
        if rem(mu(k),2)==1
            Colors(k,:)=[1 0 0];
        else
            Colors(k,:)=[0 1 0];
        end
    end
    while(mindist<2*epsilon)
        [t1,t2]=min(tmp);
        [~,tt2]=min(t1);
        %Zeros(t2(tt2)) and Zeros(tt2) are closer than epsilon; we need to merge
        Zeros(tt2)=[];
        mu(t2(tt2)) = mu(t2(tt2))+mu(tt2);
        if rem(mu(tt2),2)==1 || rem (mu(t2(tt2)),2)==1
            Colors(t2(tt2),:)=[ 0 0 1];
        end
        Colors(tt2,:)=[];
        mu(tt2)=[];
        tmp0=tmp0-1;
        tmp = ones(tmp0,1)*Zeros- Zeros.'*ones(1, tmp0) + eye(tmp0);
        tmp=abs(tmp);
        mindist= min(min(tmp));
    end
end
%We should do the same with poles but it is a rare instance.. oh well
SpecialPoints = cat(2,Zeros,Poles);
ZZ= zeros(1,sum(mu));
for i=1:length(mu)
  for j=1:mu(i)
    ZZ(sum(mu(1:i-1))+j) = Zeros(i);
  end
end
PP = zeros(1,sum(nu));
for i=1:length(nu)
  for j=1:nu(i)
    PP(sum(nu(1:i-1))+j) = Poles(i);
  end
end
Num = Const * poly(ZZ);
Den = poly  (PP);
Q= @(z) polyval(Num,z)./polyval(Den,z);

%clf;
%legend("off");
%grid on;
%Now we plot also the critical trajectories issuing from the
%simple poles   
for k=1:length(Poles) %We also plot the critical trajectories from simple poles of the quadratic differential
  if (nu(k)==1) 
    tmp=polyder(Den);
    tmp1 = (polyval(Num,Poles(k))/polyval(tmp,Poles(k)));
    tmp2 = tmp1;
    tmp2 = 1/(tmp2/abs(tmp2));
    StartAngle = tmp2;
    Phase = StartAngle;
    pt = Poles(k) + epsilon/9*Phase;
    %tmp3 = sqrt(polyval(Num,pt)/polyval(Den,pt));
    %tmp3= tmp3/abs(tmp3);    
    %Y = Trajectory(Num, Den, Length,pt ,epsilon/10,1-flag);
    Y = QTrajectory(Q,SpecialPoints,Length,pt,epsilon);
    countcrit=countcrit+1;
    CritTraj{countcrit} = Y;
    if Draw 
        plot(CurrAxes, real(Y(1:end-1)),imag(Y(1:end-1)),ColorCrit,'LineWidth',Thick);
        hold on;
    end
  end
end


plot(real(Poles), imag(Poles),'r.','MarkerSize',12);
Adjacency = zeros(length(Zeros)); 
%this is the adjacency matrix telling us which critical points are connected to which other.

for k=1:length(Zeros)
    h=1e-3;
    C=Q(Zeros(k)+h)/h^(mu(k));
    StartAngle=exp(1i*angle(C)/2);

    for j=0:(mu(k)+2)-1
        Phase = StartAngle*exp(2i*pi/(mu(k)+2)*(j));
        pt = Zeros(k) + epsilon/50*Phase;
        %Y = Trajectory(Num, Den, Length,pt ,epsilon,flag);
        Y = QTrajectory(Q,SpecialPoints,Length,pt,epsilon);
        countcrit=countcrit+1;
        CritTraj{countcrit} = Y;
        %What remains to do now is to find if this trajectory connects two
        %critical points...The current point where the trajectory started is
        %Zeros(k), so we need to make a list with this point put very far so that it is not counted.

        TMP = Zeros;
        TMP(k) = 1e5;
        Distances = abs (TMP - Y(end)); % the distances of each critical point from the last point of the trajectory
        [Delta,ell] = min(Distances);

        if (Delta<5*epsilon)
            % we are significantly close to the point Zeros(ell) so that Zeros(k) and Zeros(ell) are connected
            Adjacency(k,ell)=1;
            % Thickness =Thick;
            Color='b';
        else
            Color = 'b';
            % Thickness = Thick;
        end

        %We finally plot the trajectory in the appropriate color: thick red for
        %cuts, thick blue for (critica) shorelines, thin blue for
        if (mod(mu(k),2)==0)
            if sum(Colors(k,:)==[0 1 0])==3
                if Draw 
                    plot (CurrAxes, real(Y), imag(Y),'g', 'LineWidth',Thick/2);
                end
            else
                if Draw 
                    plot(CurrAxes, real(Y), imag(Y),ColorCrit, 'LineWidth',Thick);
                end
            end

        else
            if Draw 
                plot(CurrAxes,real(Y), imag(Y),Color,'LineWidth',Thick);
            end
        end
        if Draw
            hold on;
            plot(real(Zeros), imag(Zeros),'b.','markersize',10);
        end

    end

end


if (interactive==1)
    while(pt ~= 0)
        pt = input('Enter a regular point or enter 0 to exit: ');
        if(pt~=0)
            %Y = Trajectory(Num, Den, Length,pt ,epsilon,1);
            Y = QTrajectory(Q,SpecialPoints,Length,pt,epsilon);
            plot(CurrAxes,real(Y), imag(Y));
            %Y = Trajectory(Num, Den, Length,pt ,epsilon,0);
            Y = QTrajectory(Q,SpecialPoints,Length,pt,epsilon);
            plot(CurrAxes, real(Y), imag(Y));
        end
    end
end
 if (interactive ==0) && exist('Points','var')
   for k=1:length (Points)
     %Y = Trajectory(Num, Den, Length,Points(k) ,epsilon,1);
     Y = QTrajectory(Q,SpecialPoints,Length,Points(k),epsilon);
     plot(CurrAxes,real(Y), imag(Y), 'b', 'LineWidth',2*Thick);
     %Y = Trajectory(Num, Den, Length,Points(k) ,epsilon,0);
     Y = QTrajectory(Q,SpecialPoints,Length,Points(k),epsilon);
     plot(CurrAxes,real(Y), imag(Y), 'b', 'LineWidth',2*Thick);
   end
 end
 %Adjacency;   
 %axis( [X0-Size X0+Size Y0-Size Y0+Size],'equal');
 OddZeros=Zeros(rem(mu,2)==1);
 EvenZeros=Zeros(rem(mu,2)==0);

 if Draw 
     plot(real(OddZeros), imag(OddZeros), 'r.', 'markersize',10);
     plot(real(EvenZeros), imag(EvenZeros), 'g.', 'markersize',10);
     if (isempty(Poles)>0)
         plot(CurrAxes,real(Poles),imag(Poles),'.r','MarkerSize',10);hold on;
     end
 end

 CritTraj=CritTraj(1:countcrit); 

end


% %%%%%%%%%%%%%%%%%%%%%%%Functions%%%%%%%%%%%%%%%%%%%%%%%
% function X=Trajectory(Num, Den, Length,Startpoint, epsilon,flag)
% %% Computes the horizontal trajectory of the quadratic differential Q = Num/Den dz^2
% %% through Startpoint, of approximate length Length, with step epsilon. 
% %No = floor(Length/epsilon)+1;
% L = Length/epsilon;
% N=floor(L);
% %N should be the number of steps so that the final length of the trajectory
% %is what we requested in Length, given our adaptive step.
% 
% 
% X=zeros(1,N);
% X(1)= Startpoint;
% x=X(1);
% y = sqrt( polyval(Num,x)/polyval(Den,x));
% Direction = conj(y)/abs(y);
% ZEROS= cat(2,roots(Num), roots(Den));
% %Wrong name; anyways, it contains all zeros and poles
% %remlist= find();
% 
% % remlist = [];
% % for jt = 1:length(ZEROS)
% %     %All critical points closer than 5(3)epsilon from the starting one are
% %     %excluded from the proximity test (because presumably have been already
% %     %counted...)
% %     if (abs(ZEROS(jt) - Startpoint)<3*epsilon)
% %        remlist=[remlist,jt];
% %     end
% % end
% tmpZEROS=ZEROS(abs(ZEROS- Startpoint)>=3*epsilon).'; % We will check proximity only with these.
% 
% 
% if (flag==1) 
%     Direction=-Direction;
% end
% for kt=1:N
%     dt=epsilon;
% 
%     K1 =1i* Direction;
%     ynew =  sqrt( polyval(Num,x+dt*K1/2)/polyval(Den,x+dt*K1/2));
%     DirectionNew = conj(ynew)/abs(ynew);
%     if (abs(DirectionNew+Direction)<abs(DirectionNew-Direction)) %Continuity check on the determination of the root.
%         Direction1=-DirectionNew;
%     else
%         Direction1=DirectionNew;
%     end
% 
%     K2 =1i* Direction1;
%     ynew =  sqrt( polyval(Num,x+dt*K2/2)/polyval(Den,x+dt*K2/2));
%     DirectionNew = conj(ynew)/abs(ynew);
%     if (abs(DirectionNew+Direction)<abs(DirectionNew-Direction))
%         Direction2=-DirectionNew;
%     else
%         Direction2=DirectionNew;
%     end
% 
%     K3 =1i* Direction2;
%     ynew =  sqrt( polyval(Num,x+dt*K3)/polyval(Den,x+dt*K3));
%     DirectionNew = conj(ynew)/abs(ynew);
%     if (abs(DirectionNew+Direction)<abs(DirectionNew-Direction))
%         Direction3=-DirectionNew;
%     else
%         Direction3=DirectionNew;
%     end
%     K4=1i*Direction3;
% 
%     Direction= 1/6i*(K1+2*K2+2*K3+K4);
%     x= x+dt*1i*Direction;
%     Direction= Direction/abs(Direction);
%     % Runge Kutta 4
%     % End Runge Kutta
% 
%     X(kt) = x;
% end
% 
% %Now we find out if our trajectory came close to some of the other points,
% %in which case we truncate.
% ProximityRadius = 2*epsilon;
% tmp=length(tmpZEROS);
% tmp2 = abs((X.')*ones(1,tmp)  - ones(N,1)*tmpZEROS); %Contains the distances of trajectory X from each of the zeros/poles
% %We now find the first occurrence when X gets closer than ProximityRadius
% %to any of the zeros/poles.
% 
% [row,~]=find (tmp2<ProximityRadius); %Find indices where the trajectory is near a special point
% numpoints = min(kt-1,min(row)); %Find the first instance. The tail after that we truncate below.
% X(numpoints+1:end)=[];
% end
