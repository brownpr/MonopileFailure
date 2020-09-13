clear all;clc

%SECTIONS
N=500;                              %Desired sections
NN=N+5;                             %Number of nodes
NE=N+4;                             %Number of elements


%PARAMETERS

ODP = 9;                            %Outer pile diameter [m] 
k = 30*10^6;                        %Spring stiffness [N/m^2 or Pa]
Ep = 190*10^9;                      %Modulus of steel (pile) [N/m^2 or Pa]
Fs_min = 0;                         %minimum lateral load at surface
Fs_max = 18.4*10^6;                 %maximum lateral load at surface
NLoad = 5*N;                       	%Separate scalling for Nload than N
Fs = linspace(Fs_min , Fs_max , NLoad);	%Force at surface
Ms = 0;                             %Moment at surface
Fb = 0;                             %Force at base
Mb = 0;                             %Moment at base
Px = 0;                             %Applied lateral load
Ip = (pi/4)*(((ODP/2)^4)-((ODP/2)-tt)^4);                                       %Second moment of Inertia of Pile
h = x_depth/N;                      %Height of elements
depth = linspace(0,x_depth,N+1);    %Depth of the monopile
P_cell = cell(NLoad,1);
%P_cell{1} = zeros(N+1,1);

B = 0.1;                            %soil constant
phi = 35;                           %soil constant for sand typically between 30 and 60 (degrees)
beta = 45 + phi/2;
alpha = phi/2;
K0 = 0.4;
Kp = (tand(45 + phi/2))^2;
Ka = (tand(45 - phi/2))^2;
C1 = tand(beta)*(Kp*tand(alpha) + K0*(tand(phi)*sind(beta)*(1/cosd(alpha) + 1) - tand(alpha)));
C2 = Kp - Ka;
C3 = (Kp^2)*(Kp + K0*tand(phi)) - Ka;
gamma = 10^4;                       %specific weight of soil (N/m^3)

a = 10;                             %number of desired plots, will also plot inititial values

%EQUATIONS
%y1dot = (1/h)*(-0.5*y(ii-1) + 0.5*y(ii+1));                                    %Rotation
%y2dot = (1/(h^2))*(y(ii-1) - 2*y(ii) + y(ii+1));                               %Moment
%y3dot = (1/(h^3))*(-0.5*y(ii-2) + y(ii-1) - y(ii+1) +0.5*y(ii+2));          	%Shear
%y4dot = (1/(h^4))*(y(ii-2) - 4*y(ii-1) + 6*y(ii) - 4*y(ii+1) + y(ii+2));       %Soil Reaction


%SOLUTION CELL
SOL_EBE = cell(NLoad,1);                %Empty cell for solutions

%EQUATION CELL
EQ_C = cell(NLoad,1);

%y CELL
y_cell = cell(NLoad,1);

%k_star cell
k_cell = cell(NLoad,1);

%ITTERATION VARYING FORCE
for i = 1:NLoad

	fs = Fs(i);

	SOL_init = zeros(NN,1);        	%Empty vector of solutions to the problem.
	SOL_init(NN-3,1) = fs;         	%Stores Force at surface in solution vector
	SOL_init(NN-2,1) = Ms;         	%Stores Moment at surface in solution vector
	SOL_init(NN-1,1) = Fb;         	%Stores Force at base in solution vector
	SOL_init(NN,1) = Mb;           	%Stores Moment at base in solution vector
	SOL_EBE{i} = SOL_init;
    
    if i == 1
        P_vals = zeros(N+1, 1);
    else
        P_vals = P_cell{i-1};
    end

    %CALCULATION K_STAR FOR VARYING DEPTH
	k_star = zeros(N,1);
	for j = 1:N+1

		z = depth(j);

		Pus = (C1*z +C2*ODP)*gamma*z;
		Pud = C3*ODP*gamma*z;
		Pu = min(Pus,Pud);
        
        %Pu cannot = 0, instead set as a small number.  
        if Pu == 0 
            Pu = 1;
        end
        
		A = (3 - 0.8*(z/ODP));
		if A < 0.9
 			A = 0.9;
        end 
        
		P = P_vals(j);
		H = (1/B)*((A*Pu - P)^2)/(A*Pu);
		k_star(j,1) = (k*H)/(k+H);
        
        if P >= Pu
        	k_star(j,1) = ones(N,1).*0.00001;
        end
        
        
	end
    k_cell{i} = k_star;                                                     %Saving k_star in k_cell
    
	%SETTING UP EQUATION MATRIX
	EQ_M = zeros(NN);                                                       %Empty equation matrix

	jj = 1;                                                                 %Counter for columns 
    for ii = 3:NN-2

		k_act = k_star(ii-2,1);
        %Setting up matrix with 4th order solutions
        EQ_M(jj,ii-2:ii+2) = (1/(h^4))*(Ep*Ip)*[1, -4, 6, -4, 1] + k_act*[0, 0, 1, 0, 0] + Px*(1/(h^2))*(Ep*Ip)*[0, 1, -2, 1, 0];
    
    	jj = jj+1;                                                          %Increase counter
	end

	EQ_M(NN-3,1:5) = (1/(h^3))*(Ep*Ip)*[0.5, 2, 0, -4, 5];                  %Shear at surface using BC, NN=3
	EQ_M(NN-2,2:4) = (1/(h^2))*(Ep*Ip)*[1, -2, 1];                          %Moment at surface using BC, NN=3
	EQ_M(NN-1,NN-4:NN) = (1/(h^3))*(Ep*Ip)*[-0.5, 2, 0, -4, 0.5];           %Shear at base using BC, NN=NN-3
	EQ_M(NN,NN-3:NN-1) = (1/(h^2))*(Ep*Ip)*[1, -2, 1];                      %Moment at surface using BC, NN = NN-3
	
	EQ_C{i} = EQ_M;
	
	%SOLVER                                         %Empty vector of y values of length NN
	equ = EQ_C{i};
    sols = SOL_EBE{i};
    y_m = equ\sols;
	y_cell{i} = y_m;
    
    %Initial y values are 0 (0 displacement) else they are equal to
    %previous y values
    if i == 1
        y_cell_temp = zeros(NN,1);
    else
        y_cell_temp = y_cell{i-1};
    end
    
    
	y_delta = y_cell{i} - y_cell_temp;
    
    pimp = y_delta(3:end -2);
	P_delta = k_star.*pimp;                         %y_delta;
	
    if i == 1
        P_temp = zeros(N+1,1);
    else
        P_temp = P_cell{i-1};
    end
    P_cell{i} = P_temp + P_delta;
end

%DIFFERENTIATING FOR ROTATION, R
R_cell = cell(NLoad,1);
for rot = 1:length(y_cell)
    R_cell{rot} = diff(y_cell{rot}).*(Ep*Ip);
end


%DIFFERENTIATING FOR BENDING MOMENT, M
M_cell = cell(NLoad,1);
for mom = 1:length(R_cell)
    M_cell{mom} = diff(R_cell{mom});
end


%DIFFERENTIATING FOR SHEAR, S
S_cell = cell(NLoad,1);
for she = 1:length(M_cell)
    S_cell{she} = diff(M_cell{she});
end


%DIFFERENTIATING FOR SOIL REACTION, P
SR_P_cell = cell(NLoad,1);
for sr = 1:length(S_cell)
    SR_P_cell{sr} = diff(S_cell{sr});
end

%PLOTTING
figure
clf;
subplot(1,5,1);
hold on
x = linspace(0,(-1)*x_depth,N+1);
for vv = 1:NLoad/a:NLoad
    y_temp = y_cell{vv};
    Force = round(Fs(vv)/(10^6),4,'significant');
    plot(y_temp(3:NN-2),x,'DisplayName',['Force: ' num2str(Force) ' MN']);

end
y_last = y_cell{end};
Force = round(Fs_max/(10^6),4,'significant');
plot(y_last(3:NN-2),x,'DisplayName',['Force: ' num2str(Force) ' MN']);
title('Y displacement');
xlabel('Y displacement [m]')
ylabel('Depth [m]')
hold off
legend show
legend('location','best') 

subplot(1,5,2);
hold on
for vv = 1:NLoad/a:NLoad
    R_temp = R_cell{vv};
    plot(R_temp(2:end-2),x);

end
R_last = R_cell{end};
plot(R_last(3:NN-2),x);

xlabel('Rotation [deg]')
ylabel('Depth [m]')
title('Rotation, R');
hold off

subplot(1,5,3);
hold on
for vv = 1:NLoad/a:NLoad
    M_temp = M_cell{vv};
    plot(M_temp(2:end-1),x);
    
end
M_last = M_cell{end};
plot(M_last(2:end-1),x);

xlabel('Moment [Nm]')
ylabel('Depth [m]')
title('Moment, M');
hold off


subplot(1,5,4);
hold on
for vv = 1:NLoad/a:NLoad
    S_temp = S_cell{vv};
    plot(S_temp(1:end-1),x);
    
    
end
S_last = S_cell{end};
plot(S_last(1:end-1),x);

xlabel('Shear [N]')
ylabel('Depth [m]')
title('Shear, S');
hold off

subplot(1,5,5);
hold on
for vv = 1:NLoad/a:NLoad
    SR_P_temp = SR_P_cell{vv};
    plot(SR_P_temp(1:end),x);
    
end
SR_P_last = SR_P_cell{end};
plot(SR_P_last(1:end),x);

xlabel('Soil Reaction [N/m]')
ylabel('Depth [m]')
title('Soil Reaction, P');
hold off







% Plotting p-y curves. Not so effective 
% figure(2)
% hold on
% for uu = 1:N/a:N
% 	y_vals = zeros(N,1);
% 	x_vals = zeros(N,1);
% 
% 	for u = 1:N
% 		x_axis = y_cell{u};
%         x_vals(u) = x_axis(uu);
% 		
%         y_axis = P_cell{u};
%         y_vals(u) = y_axis(uu);
% 	end
% 	plot(x_vals,y_vals,'DisplayName',['Itr:' num2str(uu)])
% end
% title('P-Y Curves');
% hold off
% legend show
% legend('location','se') 