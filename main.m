% The construction of the Taylor state coded below is presented in
% detail in A.J. Cerfon and M. O'Neil, Exact axisymmetric Taylor states 
%for shaped plasmas, Physics of Plasmas 21, 064501 (2014)


% Copyright (C) 2014: Antoine Cerfon
% Contact: cerfon@cims.nyu.edu
% 
% This program is free software; you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation; either version 2 of the License, or 
% (at your option) any later version.  This program is distributed in 
% the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE.  See the GNU General Public License for more 
% details. You should have received a copy of the GNU General Public 
% License along with this program; 
% if not, see <http://www.gnu.org/licenses/>.


clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%               Define parameters for equilibrium of interest
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = 0.9;%Inverse aspect ratio
delta = 0.15;%Triangularity
kappa = 1.2;%Elongation
alpha = asin(delta);
N1 = -(1+alpha)^2/(epsilon*kappa^2);%Curvature at outer equatorial point
N2 = (1-alpha)^2/(epsilon*kappa^2);%Curvature at inner equatorial point
N3 = kappa/(epsilon*(cos(alpha))^2);%Curvature at top
lambda=2.4/epsilon;%Initial guess for lambda/mu
k = 1.6/epsilon;%Initial guess for kappa
xsep = 1+0.5*epsilon;%X location of the separatrix
ysep = 1.25*epsilon*kappa;%Y location of the separatrix
p=1; %Boolean: 0 if one wants to use a prespecified vector D0 for the initial guess 
     %         1 if one wants to compute D0 by solving the linear system with fixed k and lambda 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Compute initial guess D0 from boundary conditions
%                       (done if p==1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
if p==1
A = [psi1(1+epsilon,0,lambda) psi2(1+epsilon,0,k,lambda) psi3(1+epsilon,0,k,lambda) psi4(1+epsilon,0,lambda) psi5(1+epsilon,0,lambda) psi6(1+epsilon,0,lambda) psi7(1+epsilon,0,lambda) psi8(1+epsilon,0,k,lambda) psi9(1+epsilon,0,k,lambda) psi10(1+epsilon,0,lambda)
psi1(1-epsilon,0,lambda) psi2(1-epsilon,0,k,lambda) psi3(1-epsilon,0,k,lambda) psi4(1-epsilon,0,lambda) psi5(1-epsilon,0,lambda) psi6(1-epsilon,0,lambda) psi7(1-epsilon,0,lambda) psi8(1-epsilon,0,k,lambda) psi9(1-epsilon,0,k,lambda) psi10(1-epsilon,0,lambda)
psi1(1-epsilon*delta,-kappa*epsilon,lambda) psi2(1-epsilon*delta,-kappa*epsilon,k,lambda) psi3(1-epsilon*delta,-kappa*epsilon,k,lambda) psi4(1-epsilon*delta,-kappa*epsilon,lambda) psi5(1-epsilon*delta,-kappa*epsilon,lambda) psi6(1-epsilon*delta,-kappa*epsilon,lambda) psi7(1-epsilon*delta,-kappa*epsilon,lambda) psi8(1-epsilon*delta,-kappa*epsilon,k,lambda) psi9(1-epsilon*delta,-kappa*epsilon,k,lambda) psi10(1-epsilon*delta,-kappa*epsilon,lambda)
psi1(xsep,ysep,lambda) psi2(xsep,ysep,k,lambda) psi3(xsep,ysep,k,lambda) psi4(xsep,ysep,lambda) psi5(xsep,ysep,lambda) psi6(xsep,ysep,lambda) psi7(xsep,ysep,lambda) psi8(xsep,ysep,k,lambda) psi9(xsep,ysep,k,lambda) psi10(xsep,ysep,lambda)
psi1_x(1-epsilon*delta,-kappa*epsilon,lambda) psi2_x(1-epsilon*delta,-kappa*epsilon,k,lambda) psi3_x(1-epsilon*delta,-kappa*epsilon,k,lambda) psi4_x(1-epsilon*delta,-kappa*epsilon,lambda) psi5_x(1-epsilon*delta,-kappa*epsilon,lambda) psi6_x(1-epsilon*delta,-kappa*epsilon,lambda) psi7_x(1-epsilon*delta,-kappa*epsilon,lambda) psi8_x(1-epsilon*delta,-kappa*epsilon,k,lambda) psi9_x(1-epsilon*delta,-kappa*epsilon,k,lambda) psi10_x(1-epsilon*delta,-kappa*epsilon,lambda)
psi1_y(1+epsilon,0,lambda) psi2_y(1+epsilon,0,k,lambda) psi3_y(1+epsilon,0,k,lambda) psi4_y(1+epsilon,0,lambda) psi5_y(1+epsilon,0,lambda) psi6_y(1+epsilon,0,lambda) psi7_y(1+epsilon,0,lambda) psi8_y(1+epsilon,0,k,lambda) psi9_y(1+epsilon,0,k,lambda) psi10_y(1+epsilon,0,lambda)
psi1_y(1-epsilon,0,lambda) psi2_y(1-epsilon,0,k,lambda) psi3_y(1-epsilon,0,k,lambda) psi4_y(1-epsilon,0,lambda) psi5_y(1-epsilon,0,lambda) psi6_y(1-epsilon,0,lambda) psi7_y(1-epsilon,0,lambda) psi8_y(1-epsilon,0,k,lambda) psi9_y(1-epsilon,0,k,lambda) psi10_y(1-epsilon,0,lambda)
psi1_x(xsep,ysep,lambda) psi2_x(xsep,ysep,k,lambda) psi3_x(xsep,ysep,k,lambda) psi4_x(xsep,ysep,lambda) psi5_x(xsep,ysep,lambda) psi6_x(xsep,ysep,lambda) psi7_x(xsep,ysep,lambda) psi8_x(xsep,ysep,k,lambda) psi9_x(xsep,ysep,k,lambda) psi10_x(xsep,ysep,lambda)
psi1_y(xsep,ysep,lambda) psi2_y(xsep,ysep,k,lambda) psi3_y(xsep,ysep,k,lambda) psi4_y(xsep,ysep,lambda) psi5_y(xsep,ysep,lambda) psi6_y(xsep,ysep,lambda) psi7_y(xsep,ysep,lambda) psi8_y(xsep,ysep,k,lambda) psi9_y(xsep,ysep,k,lambda) psi10_y(xsep,ysep,lambda)
psi1_yy(1+epsilon,0,lambda)+N1*psi1_x(1+epsilon,0,lambda) psi2_yy(1+epsilon,0,k,lambda)+N1*psi2_x(1+epsilon,0,k,lambda) psi3_yy(1+epsilon,0,k,lambda)+N1*psi3_x(1+epsilon,0,k,lambda) psi4_yy(1+epsilon,0,lambda)+N1*psi4_x(1+epsilon,0,lambda) psi5_yy(1+epsilon,0,lambda)+N1*psi5_x(1+epsilon,0,lambda) psi6_yy(1+epsilon,0,lambda)+N1*psi6_x(1+epsilon,0,lambda) psi7_yy(1+epsilon,0,lambda)+N1*psi7_x(1+epsilon,0,lambda) psi8_yy(1+epsilon,0,k,lambda)+N1*psi8_x(1+epsilon,0,k,lambda) psi9_yy(1+epsilon,0,k,lambda)+N1*psi9_x(1+epsilon,0,k,lambda) psi10_yy(1+epsilon,0,lambda)+N1*psi10_x(1+epsilon,0,lambda)
psi1_yy(1-epsilon,0,lambda)+N2*psi1_x(1-epsilon,0,lambda) psi2_yy(1-epsilon,0,k,lambda)+N2*psi2_x(1-epsilon,0,k,lambda) psi3_yy(1-epsilon,0,k,lambda)+N2*psi3_x(1-epsilon,0,k,lambda) psi4_yy(1-epsilon,0,lambda)+N2*psi4_x(1-epsilon,0,lambda) psi5_yy(1-epsilon,0,lambda)+N2*psi5_x(1-epsilon,0,lambda) psi6_yy(1-epsilon,0,lambda)+N2*psi6_x(1-epsilon,0,lambda) psi7_yy(1-epsilon,0,lambda)+N2*psi7_x(1-epsilon,0,lambda) psi8_yy(1-epsilon,0,k,lambda)+N2*psi8_x(1-epsilon,0,k,lambda) psi9_yy(1-epsilon,0,k,lambda)+N2*psi9_x(1-epsilon,0,k,lambda) psi10_yy(1-epsilon,0,lambda)+N2*psi10_x(1-epsilon,0,lambda)
psi1_xx(1-epsilon*delta,-epsilon*kappa,lambda)+N3*psi1_y(1-epsilon*delta,-epsilon*kappa,lambda) psi2_xx(1-epsilon*delta,-epsilon*kappa,k,lambda)+N3*psi2_y(1-epsilon*delta,-epsilon*kappa,k,lambda) psi3_xx(1-epsilon*delta,-epsilon*kappa,k,lambda)+N3*psi3_y(1-epsilon*delta,-epsilon*kappa,k,lambda) psi4_xx(1-epsilon*delta,-epsilon*kappa,lambda)+N3*psi4_y(1-epsilon*delta,-epsilon*kappa,lambda) psi5_xx(1-epsilon*delta,-epsilon*kappa,lambda)+N3*psi5_y(1-epsilon*delta,-epsilon*kappa,lambda) psi6_xx(1-epsilon*delta,-epsilon*kappa,lambda)+N3*psi6_y(1-epsilon*delta,-epsilon*kappa,lambda) psi7_xx(1-epsilon*delta,-epsilon*kappa,lambda)+N3*psi7_y(1-epsilon*delta,-epsilon*kappa,lambda) psi8_xx(1-epsilon*delta,-epsilon*kappa,k,lambda)+N3*psi8_y(1-epsilon*delta,-epsilon*kappa,k,lambda) psi9_xx(1-epsilon*delta,-epsilon*kappa,k,lambda)+N3*psi9_y(1-epsilon*delta,-epsilon*kappa,k,lambda) psi10_xx(1-epsilon*delta,-epsilon*kappa,lambda)+N3*psi10_y(1-epsilon*delta,-epsilon*kappa,lambda)];


B = -[psi0(1+epsilon,0,lambda)
    psi0(1-epsilon,0,lambda)
    psi0(1-epsilon*delta,-kappa*epsilon,lambda)
    psi0(xsep,ysep,lambda)
    psi0_x(1-epsilon*delta,-kappa*epsilon,lambda)
    psi0_y(1+epsilon,0,lambda)
    psi0_y(1-epsilon,0,lambda)
    psi0_x(xsep,ysep,lambda)
    psi0_y(xsep,ysep,lambda)
    psi0_yy(1+epsilon,0,lambda)+N1*psi0_x(1+epsilon,0,lambda)
    psi0_yy(1-epsilon,0,lambda)+N2*psi0_x(1-epsilon,0,lambda)
    psi0_xx(1-epsilon*delta,-epsilon*kappa,lambda)+N3*psi0_y(1-epsilon*delta,-epsilon*kappa,lambda)];

C0 = A\B;

D0 = [C0(1) C0(2) C0(3) C0(4) C0(5) k lambda C0(6) C0(7) C0(8) C0(9) C0(10)];
else
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Specify initial guess D0 by hand
%                       (done if p==0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D0 = [7.2803    1.6935   -6.6120   -6.4564    5.7232 1.8465    2.6635   -2.5153    3.7255    1.2570 -1.8505   -0.0806];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Nonlinear solve for actual unknowns D
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options=optimset('Display','iter');
D = fsolve(@(C)mysystem(C,epsilon,kappa,delta,N1,N2,N3,xsep,ysep),D0,options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Construct flux function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(-(1+epsilon)-0.3:.005:1+epsilon+0.3,-kappa*epsilon-0.2:.005:ysep+0.2);
Z = eq_psi(abs(X),Y,D);
I = find(Z<0);
Z(I)=NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Plot flux contours
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
contour(X, Y, Z, 10,'LineWidth',1.5,'color','black');
hold on
pcolor(X, Y, Z);
shading interp
axis equal
xlim([-(1+epsilon+0.3) 1+epsilon+0.3])
ylim([-kappa*epsilon-0.2 ysep+0.2])
colorbar
grid on
xlabel('R','FontSize',14,'FontName','Times')
ylabel('Z','FontSize',14,'FontName','Times')
set(gca,'FontSize',14,'FontName','Times')