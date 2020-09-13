function [y_vals,k_star,P_vals,R_vals,M_vals,S_vals,SR_P_vals] = callBeamState(Fs,N,P_init,y_init, cycNum)

global depth Ip ODP k Ep Ms Fb Mb Px B gamma NLoad NN C1 C2 C3 h

%CREATING CELLS
k_star = zeros(NN,NLoad);               	%k_star values cell
P_vals = zeros(NN,NLoad);                   %P value cell
y_vals = zeros(NN,NLoad);                   %Displacement values cell
R_vals = zeros(N+1,NLoad);                  %Rotation values cell
M_vals = zeros(N+1,NLoad);                  %Moment values cell
S_vals = zeros(N+1,NLoad);                  %Shear values cell
SR_P_vals = zeros(N+1,NLoad);               %Soil reaction values cell

fail = 0;                                   %fail condition P>Pu
LFF = -1;                                   %?

%ITTERATION VARYING FORCE
for i = 1:NLoad

    %RESETTING EQUATION AND SOLUTION TO ZERO
	EQ_M = zeros(NN);                       %Equation matrix
    SOL_init = zeros(NN,1);                 %Empty vector of solutions to the problem.
    
    fs = Fs(i);

    if i == 1
        P_temp_vals = P_init;
        y_vals_temp = y_init;
    else
        P_temp_vals = P_vals(1:end,i-1);
        y_vals_temp = y_vals(1:end,i-1);
    end

    %CALCULATION K_STAR FOR VARYING DEPTH
	for j = 1:NN

		z = depth(j);
        
        %Pu cannot = 0, instead set as a small number.
        %Calculating kk fot varying depth
        if z < 0.0001 
            Pu = 0.001;
            kk=k*0.001;
        else 
            Pus = (C1*z +C2*ODP)*gamma*z;
            Pud = C3*ODP*gamma*z;
            Pu = min(Pus,Pud);
            kk = k*z;
        end
        
		A = (3 - 0.8*(z/ODP));
		if A < 0.9
 			A = 0.9;
        end 
        
        P = P_temp_vals(j);
        if P < 0 
            Pu = (-1)*Pu;
        end
        
        if abs(P) >= abs(Pu*A)
        	P_temp_vals(j)=A*Pu;
            %Only print out fail once per iteration. 
%             if fail == 0 
%                 warning(['During cycle: ' num2str(cycNum) ' and force load: ' sprintf('%.4g', fs) 'N, P value has exceeded Pu_max at depth: ' num2str(z)]) 
%                 fail = 1;
%             end
        end
        
		H = (1/B)*((A*Pu*LFF - P)^2)/(abs(A*Pu)); 
		k_star(j,i) = (kk*H)/(kk+H);
    end
    
    %SETTING UP EQUATION MATRIX
	jj = 1;                                                                 %Counter for columns 
    for ii = 3:NN-2
		k_act = k_star(ii,i);
        %Setting up matrix with 4th order solutions
        EQ_M(jj,ii-2:ii+2) = (1/(h^4))*(Ep*Ip)*[1, -4, 6, -4, 1] + k_act*[0, 0, 1, 0, 0] + Px*(1/(h^2))*(Ep*Ip)*[0, 1, -2, 1, 0];
    	SOL_init(jj)= (-1)*(P_temp_vals(ii,1)-k_act*y_vals(ii,i));
        jj = jj+1;                                                          %Increase counter
    end

    %SETTING UP SOLUTIONS
	SOL_init(NN-3,1) = fs;         	%Stores Force at surface in solution vector
	SOL_init(NN-2,1) = Ms;         	%Stores Moment at surface in solution vector
	SOL_init(NN-1,1) = Fb;         	%Stores Force at base in solution vector
	SOL_init(NN,1) = Mb;           	%Stores Moment at base in solution vector
	
	EQ_M(NN-3,1:5) = (1/(h^3))*(Ep*Ip).*[-0.5, 1, 0, -1, 0.5];               %Shear at surface using BC, NN=3
	EQ_M(NN-2,1:5) = (1/(h^2))*(Ep*Ip).*[0, 1, -2, 1, 0];                    %Moment at surface using BC, NN=3
	EQ_M(NN-1,NN-4:NN) = (1/(h^3))*(Ep*Ip).*[-0.5, 1, 0, -1, 0.5];           %Shear at base using BC, NN=NN-3
	EQ_M(NN,NN-4:NN) = (1/(h^2))*(Ep*Ip).*[0, 1, -2, 1, 0];                  %Moment at base using BC, NN = NN-3
	
	%SOLVER                                         
	equ = EQ_M;
    sols = SOL_init;
    y_m = equ\sols;
	y_vals(1:end,i) = y_m;
    
	y_delta = y_vals(1:end,i) - y_vals_temp;        %y_delta
	P_delta = k_star(1:end,i).*y_delta;            	%P_delta
    P_vals(1:end,i) = P_temp_vals + P_delta;        %Calculating P
    
    %CALCULATING ROTATION, MOMENT, SHEAR AND P
    for cnt = 1:N+1
        cnt2 = cnt+2;
        R_vals(cnt,i) = (1/h)*(-0.5*y_vals(cnt2-1,i)+0.5*y_vals(cnt2+1,i));
        M_vals(cnt,i) = ((Ep*Ip)/(h^2))*(y_vals(cnt2-1,i)-2*y_vals(cnt2,i)+y_vals(cnt2+1,i));
        S_vals(cnt,i) = ((Ep*Ip)/(h^3))*(-0.5*y_vals(cnt2-2,i)+y_vals(cnt2-1,i)-y_vals(cnt2+1,i)+0.5*y_vals(cnt2+2,i));
        %SR_P_vals(cnt,i) = Ep*Ip/h^4*(y_vals(cnt2-2)-4*y_vals(cnt2-1)+6*y_vals(cnt2)-4*y_vals(cnt2+1)+y_vals(cnt2+2));
        SR_P_vals(cnt,i) = P_vals(cnt,i);
    end
end
end