clear all
close all
clc

%SECTIONS
N=20;               %Desired sections
NN=N+5;             %Number of nodes
NE=N+4;             %Number of elements


%PARAMETERS

x_depth = 50;                               %Depth of pile
tt = 0.06;                                  %Thickness of pile
ODP = 9;                                    %Outer pile diameter
k = 300*10^3;                               %Spring stiffness
Ep = 210*10^6;                              %Modulus of steel (pile)
Fs_min = 0;                                 %minimum lateral load at surface
Fs_max = 18.4*10^3;                         %maximum lateral load at surface
Fs = linspace(Fs_min , Fs_max , N);         %Force at surface
Ms = 0;                                     %Moment at surface
Fb = 0;                                     %Force at base
Mb = 0;                                     %Moment at base
Px = 0;                                     %Applied lateral load
Ip = (pi/4)*(((ODP/2)^4)-((ODP/2)-tt)^4);   %Second moment of Inertia of Pile
h = x_depth/NE;                             %Height of elements
depth = linspace(0,x_depth,N);              %Depth of the monopile
P_cell = cell(N,1);
P_cell(1,1) = zeros(N,1);

B = 0.1;                                    %soil constant
phi = 35;                                   %soil constant for sand typically between 30 and 60 (degrees)
beta = 45 + phi/2;
alpha = phi/2;
K0 = 0.4;
Kp = (tan(45 + phi/2))^2;
Ka = (tan(45 - phi/2))^2;
C1 = tan(beta)*(Kp*tan(alpha) + K0*(tan(phi)*sin(beta)*(1/cos(alpha) + 1) - tan(alpha)));
C2 = Kp - Ka;
C3 = (Kp^2)*(Kp + K0*tan(phi)) - Ka;
gamma = 10^3;		%specific weight of soil (N/m^3)


%EQUATIONS
%y1dot = (1/h)*(-0.5*y(ii-1) + 0.5*y(ii+1));                                 %Rotation
%y2dot = (1/(h^2))*(y(ii-1) - 2*y(ii) + y(ii+1));                            %Moment
%y3dot = (1/(h^3))*(-0.5*y(ii-2) + y(ii-1) - y(ii+1) +0.5*y(ii+2));          %Shear
%y4dot = (1/(h^4))*(y(ii-2) - 4*y(ii-1) + 6*y(ii) - 4*y(ii+1) + y(ii+2));    %Soil Reaction


%SOLUTION CELL
SOL_EBE = cell(N,1);		%Empty cell for solutions

%EQUATION CELL
EQ_C = cell(N,1);
SOL_EBE{1,1} = zeros(NN,1);          %initial sol are 0

%y CELL
y_cell = cell(N,1);
y_cell{1,1} = zeros(N,1);

for i = 2:N

    %Computing new solution matrix
	fs = Fs(i);
	SOL_init = zeros(NN,1);        	%Empty vector of solutions to the problem.
	SOL_init(NN-3,1) = fs;         	%Stores Force at surface in solution vector
	SOL_init(NN-2,1) = Ms;         	%Stores Moment at surface in solution vector
	SOL_init(NN-1,1) = Fb;         	%Stores Force at base in solution vector
	SOL_init(NN,1) = Mb;           	%Stores Moment at base in solution vector
	SOL_EBE{i,1} = SOL_init;

	P_vals = P_cell(i-1,1);
	
    
	k_star = zeros(N,1);
	for j = 1:N

		z = depth(j);

		Pus = (C1*z +C2*ODP)*gamma*z;
		Pud = C3*ODP*gamma*z;
		Pu = min(Pus,Pud);

		A = (3 - 0.8*(z/ODP));
		if A < 0.9
			A = 0.9;
		end
		
		P = P_vals(j);
		H = (1/B)*((A*Pu - P)^2)/A*Pu;
		k_star(j,1) = k*H/(k+H);
	end


	%SETTING UP EQUATION MATRIX
	EQ_M = zeros(NN);       %Empty equation matrix

	jj = 1;                 	%Counter for rows
	for ii = 3:NN-2

		
    		%Setting up matrix with 4th order solutions
    		EQ_M(jj,ii-2:ii+2) = (1/(h^4))*(Ep*Ip)*[1, -4, 6, -4, 1] + k_star(ii-2,1)*[0, 0, 1, 0, 0] + Px*(1/(h^2))*(Ep*Ip)*[0, 1, -2, 1, 0];
    
    		jj = jj+1;                                  %Increase counter
	end

	EQ_M(NN-3,1:5) = (1/(h^3))*(Ep*Ip)*[0.5, 2, 0, -4, 5];                  %Shear at surface using BC, NN=3
	EQ_M(NN-2,2:4) = (1/(h^2))*(Ep*Ip)*[1, -2, 1];                          %Moment at surface using BC, NN=3
	EQ_M(NN-1,NN-4:NN) = (1/(h^3))*(Ep*Ip)*[-0.5, 2, 0, -4, 0.5];           %Shear at base using BC, NN=NN-3
	EQ_M(NN,NN-3:NN-1) = (1/(h^2))*(Ep*Ip)*[1, -2, 1];                      %Moment at surface using BC, NN = NN-3
	
	EQ_C(i,1) = EQ_M;
	

	%SOLVER
	y_m = zeros(NN,1);                                                        %Empty vector of y values of length NN
	y_m = EQ_C(i,1)\SOL_EBE(i,1);
	y_cell(i,1) = y_m;

	y_delta = y_cell(i,1) - y_cell(i-1,1);
	P_delta = k_star.*y_delta;
	P_cell(i,1) = P_cell(i-1,1) + P_delta;
	
end

a = 10; %number of desired graphs
hold on
for uu = 1:N/a:N;
	y_vals = zeros(N,1);
	x_vals = zeros(N,1);

	for u = 1:N
		plot_y_axis = P_cell(u,1);
		y_vals(u,1)= plot_y_axis(uu,1);
	
		plot_x_axis = y_cell(u,1);
		x_vals(u,1) = plot_x_axis(uu,1);
	end
	plot(x_vals,y_vals)

end

hold off