function [y_vals,k_vals,P_vals,R_vals,M_vals,S_vals,SR_P_vals] = callBeamState(Fs,N,P_init,y_init)

global x_depth tt ODP k Ep Ms Fb Mb Px B phi gamma NLoad NN

%PARAMETERS
Ip = (pi/4)*(((ODP/2)^4)-((ODP/2)-tt)^4);   %Second moment of Inertia of Pile
h = x_depth/N;                              %Height of elements
depth = linspace(0,x_depth,N+1);            %Depth vector of the monopile


beta = 45 + phi/2;
alpha = phi/2;
K0 = 0.4;
Kp = (tand(45 + phi/2))^2;
Ka = (tand(45 - phi/2))^2;
C1 = tand(beta)*(Kp*tand(alpha) + K0*(tand(phi)*sind(beta)*(1/cosd(alpha) + 1) - tand(alpha)));
C2 = Kp - Ka;
C3 = (Kp^2)*(Kp + K0*tand(phi)) - Ka;


%CREATING CELLS
k_vals = zeros(N,NLoad);                    %k_star values cell
P_vals = zeros(N+1,NLoad);                  %P value cell
y_vals = zeros(NN,NLoad);                   %Displacement values cell
R_vals = zeros(NN-1,NLoad);                 %Rotation values cell
M_vals = zeros(NN-2,NLoad);                 %Moment values cell
S_vals = zeros(NN-3,NLoad);                 %Shear values cell
SR_P_vals = zeros(NN-4,NLoad);          	%Soil reaction values cell


%CREATING SOLUTION AND EQUATION CELL
SOL_C = cell(NLoad,1);                      %Solution cell
EQ_C = cell(NLoad,1);                       %Equation cell


%ITTERATION VARYING FORCE
for i = 1:NLoad

    
    %RESETTING EQUATION AND SOLUTION TO ZERO
	EQ_M = zeros(NN);                       %Equation matrix
    SOL_init = zeros(NN,1);                 %Empty vector of solutions to the problem.
    
    
    %SETTING UP SOLUTIONS
	fs = Fs(i);
	SOL_init(NN-3,1) = fs;         	%Stores Force at surface in solution vector
	SOL_init(NN-2,1) = Ms;         	%Stores Moment at surface in solution vector
	SOL_init(NN-1,1) = Fb;         	%Stores Force at base in solution vector
	SOL_init(NN,1) = Mb;           	%Stores Moment at base in solution vector
	SOL_C{i} = SOL_init;
    
    if i == 1
        P_init_vals = P_init;
    else
        P_init_vals = P_vals(1:end,i-1);
    end

    %CALCULATION K_STAR FOR VARYING DEPTH
	k_star = zeros(N,NLoad);
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
        
		P = P_init_vals(j);
		H = (1/B)*((A*Pu - P)^2)/(A*Pu);
		k_star(j,i) = (k*H)/(k+H);
        
        if P >= Pu
        	k_star(j,i) = 0.00001;
        end
        
        
    end

    
    %SETTING UP EQUATION MATRIX
	jj = 1;                                                                 %Counter for columns 
    for ii = 3:NN-2

		k_act = k_star(ii-2,i);
        %Setting up matrix with 4th order solutions
        EQ_M(jj,ii-2:ii+2) = (1/(h^4))*(Ep*Ip)*[1, -4, 6, -4, 1] + k_act*[0, 0, 1, 0, 0] + Px*(1/(h^2))*(Ep*Ip)*[0, 1, -2, 1, 0];
    	jj = jj+1;                                                          %Increase counter
	end

	EQ_M(NN-3,1:5) = (1/(h^3))*(Ep*Ip).*[-0.5, 1, 0, -1, 0.5];               %Shear at surface using BC, NN=3
	EQ_M(NN-2,1:5) = (1/(h^2))*(Ep*Ip).*[0, 1, -2, 1, 0];                    %Moment at surface using BC, NN=3
	EQ_M(NN-1,NN-4:NN) = (1/(h^3))*(Ep*Ip).*[-0.5, 1, 0, -1, 0.5];           %Shear at base using BC, NN=NN-3
	EQ_M(NN,NN-4:NN) = (1/(h^2))*(Ep*Ip).*[0, 1, -2, 1, 0];                  %Moment at base using BC, NN = NN-3
	
	EQ_C{i} = EQ_M;
	
	%SOLVER                                         
	equ = EQ_C{i};
    sols = SOL_C{i};
    y_m = equ\sols;
	y_vals(1:end,i) = y_m;
    
    %Initial y values are 0 (0 displacement) else they are equal to
    %previous y values
    if i == 1
        y_vals_temp = y_init;
    else
        y_vals_temp = y_vals(1:end,i-1);
    end
    
	y_delta = y_vals(1:end,i) - y_vals_temp;        %y_delta
    
    y_d = y_delta(3:end -2);
	P_delta = k_star(1:end,i).*y_d;              	%P_delta
	
    if i == 1
        P_temp = P_init;
    else
        P_temp = P_vals(1:end,i-1);
    end
    P_vals(1:end,i) = P_temp + P_delta;

end

%DIFFERENTIATING FOR ROTATION, R
for rot = 1:length(y_vals)
    R_vals(1:end,rot) = diff(y_vals(1:end,rot)).*(Ep*Ip);
end


%DIFFERENTIATING FOR BENDING MOMENT, M
for mom = 1:length(R_vals)
    M_vals(1:end,mom) = diff(R_vals(1:end,mom));
end


%DIFFERENTIATING FOR SHEAR, S
for she = 1:length(M_vals)
    S_vals(1:end,she) = diff(M_vals(1:end,she));
end


%DIFFERENTIATING FOR SOIL REACTION, P
for sr = 1:length(S_vals)
    SR_P_vals(1:end,sr) = diff(S_vals(1:end,sr));
end

end
