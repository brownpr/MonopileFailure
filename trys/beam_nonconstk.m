clf; clear all; close all; clc

%PARAMETERS
x_depth = 50;       %Depth of pile
tt = 0.06;          %Thickness of pile
ODP = 9;            %Outer pile diameter
Ep = 210*10^6;      %Modulus of steel (pile)
Fs = 18.4*10^3;     %Force at surface
Ms = 0;             %Moment at surface
Fb = 0;             %Force at base
Mb = 0;             %Moment at base
Px = 0;             %Applied lateral load


%SECTIONS
N=200;             %Desired sections, must be multiple of 10
NN=N+5;             %Number of nodes
NE=N+4;             %Number of elements
h = x_depth/NE;     %Height of elements


%EQUATIONS
%y1dot = (1/h)*(-0.5*y(ii-1) + 0.5*y(ii+1));                                 %Rotation
%y2dot = (1/(h^2))*(y(ii-1) - 2*y(ii) + y(ii+1));                            %Moment
%y3dot = (1/(h^3))*(-0.5*y(ii-2) + y(ii-1) - y(ii+1) +0.5*y(ii+2));          %Shear
%y4dot = (1/(h^4))*(y(ii-2) - 4*y(ii-1) + 6*y(ii) - 4*y(ii+1) + y(ii+2));    %Soil Reaction

Ip = (pi/4)*(((ODP/2)^4)-((ODP/2)-tt)^4);                                    %Second moment of Inertia of Pile

%p = -k*y(NN);                                                               %
%EBE = Ep*Ip*y4dot + p;                                                      %Simplified Euler Bernouli equation.


%SOLUTION VECTOR
SOL_EBE = zeros(NN,1);        %Empty vector of solutions to the problem.
SOL_EBE(NN-3,1) = Fs;         %Stores Force at surface in solution vector
SOL_EBE(NN-2,1) = Ms;         %Stores Moment at surface in solution vector
SOL_EBE(NN-1,1) = Fb;         %Stores Force at base in solution vector
SOL_EBE(NN,1) = Mb;           %Stores Moment at base in solution vector


%SETTING UP VARYING K

K = zeros(N,1);                 %Empty array for k values
K(1,1) = 300*10^3;              %Initial spring stiffness
dy = 4/length(K);               %Change in y for calculating k values
for ik = 2:length(K)
    K(ik,1) = K(1,1)*(tanh(ik*dy)-tanh((ik-1)*dy))/tanh(dy);
    
end
plot(K)
%SETTING UP EQUATION MATRIX
EQ_M = zeros(NN);       %Empty equation matrix


%from boundary conditions
EQ_M(NN-3,1:5) = (1/(h^3))*(Ep*Ip)*[0.5, 2, 0, -4, 5];                  %Shear at surface using BC, NN=3
EQ_M(NN-2,2:4) = (1/(h^2))*(Ep*Ip)*[1, -2, 1];                          %Moment at surface using BC, NN=3
EQ_M(NN-1,NN-4:NN) = (1/(h^3))*(Ep*Ip)*[-0.5, 2, 0, -4, 0.5];           %Shear at base using BC, NN=NN-3
EQ_M(NN,NN-3:NN-1) = (1/(h^2))*(Ep*Ip)*[1, -2, 1];                      %Moment at surface using BC, NN = NN-3


%SOLVER FOR VARYING K

y = cell(N,1);           	%Empty cell of y values of length N
for itr = 1:N
    k = K(itr,1);
    jj = 1;                 %Counter for columns 
    for ii = 3:NN-2
        %Filling equation matrix
        EQ_M(jj,ii-2:ii+2) = (1/(h^4))*(Ep*Ip)*[1, -4, 6, -4, 1] + k*[0, 0, 1, 0, 0] + Px*(1/(h^2))*(Ep*Ip)*[0, 1, -2, 1, 0];
    
        jj = jj+1;       	%Increase counter
    end
    
    
    y{itr,1} = EQ_M\SOL_EBE;
    
end

%PLOT
x = linspace(0,(-1)*x_depth,NN);
%Plots the first iteration and every N/10 iteration after that.
c = itr/(N/10);
jj_fig = 1;                     %Figure counter
for itr_fig = 1:N/10:length(y)
	figure(jj_fig)
	plot(y{itr_fig,1},x);
	xlabel('Y displacement')
	ylabel('Depth')
    K_title =  num2str(round((K(itr_fig,1)/1000),3));
    title(['Displacement of pile, k = ' K_title 'kN/m^2'])
	jj_fig = jj_fig +1;
end
