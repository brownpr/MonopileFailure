clear all;close all; clc;
global x_depth tt ODP k Ep Ms Fb Mb Px B phi gamma NLoad NN


a = 10;                             %number of desired plots, will also plot inititial values

%SECTIONS
N = 50;                            %Number of desired pile sections 
NN=N+5;                             %Number of nodes
NE=N+4;                             %Number of elements
NLoad = 5*N;                      	%Separate scalling for Nload than N

%PILE PARAMETERS
x_depth = 50;                       %Depth of pile [m]
tt = 0.06;                          %Thickness of pile [m]
ODP = 9;                            %Outer pile diameter [m] 
Ep = 190*10^9;                      %Modulus of pile, set as steel [N/m^2 or Pa]

%BOUNDARY CONDITION PARAMETERS
Fs_min = 0;                         %Min Force at surface [N]
Fs_max = 18.4*(10^6);               %Max Force at surface [N]
Ms = 0;                             %Moment at surface [Nm]
Fb = 0;                             %Force at base [N]
Mb = 0;                             %Moment at base [Nm]
Px = 0;                             %Applied lateral load [N]

%SOIL PARAMETERS
k = 30*10^6;                        %Elastic spring stiffness [N/m^2 or Pa]
B = 0.1;                            %soil constant
phi = 35;                           %soil constant for sand typically between 30 and 60 [deg]
gamma = 10^4;                       %specific weight of soil [N/m^3]


%CYCLIC STUFF
cycle = [Fs_min, Fs_max, Fs_min, -1*Fs_max, Fs_min];
%y_vals = cell(NLoad, length(cycle)-1);

P_init = zeros(N+1, 1);
y_init = zeros(NN,1);


%COMPUTE CYCLIC
for i = 1:length(cycle)-1
    Fs_1 = cycle(i);
    Fs_2 = cycle(i+1);
    Fs = linspace(Fs_1, Fs_2, NLoad);
    [y_vals,k_vals,P_vals,R_vals,M_vals,S_vals,SR_P_vals] = callBeamState(Fs,N, P_init, y_init);

    y_cell(1:NLoad,i) = y_vals;
    k_cell(1:NLoad,i) = k_vals;
    P_cell(1:NLoad,i) = P_vals;
    R_cell(1:NLoad,i) = R_vals;
    M_cell(1:NLoad,i) = M_vals;
    S_cell(1:NLoad,i) = S_vals;
    SR_P_cell(1:NLoad,i) = SR_P_vals;
    
    %Using initial positions and P for beam for initial iterations.
    P_init = P_cell{end,i};
    y_init = y_cell{end,i};
end




%PLOTTING
for ii = 1:4
    
    figure
    clf;
    subplot(1,5,1);
    hold on
    x = linspace(0,(-1)*x_depth,N+1);
    force = linspace(cycle(ii),cycle(ii+1),NLoad);
    for vv = 1:NLoad/a:NLoad
        y_temp = y_cell{vv,ii};
        Force = round(force(vv)/(10^6),4,'significant');
        plot(y_temp(3:NN-2),x,'DisplayName',['Force: ' num2str(Force) ' MN']);
        
    end
    y_last = y_cell{end,ii};
    Force = round(force(end)/(10^6),4,'significant');
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
        R_temp = R_cell{vv,ii};
        plot(R_temp(2:end-2),x);
        
    end
    R_last = R_cell{end,ii};
    plot(R_last(3:NN-2),x);
    
    xlabel('Rotation [deg]')
    ylabel('Depth [m]')
    title('Rotation, R');
    hold off
    
    subplot(1,5,3);
    hold on
    for vv = 1:NLoad/a:NLoad
        M_temp = M_cell{vv,ii};
        plot(M_temp(2:end-1),x);
        
    end
    M_last = M_cell{end,ii};
    plot(M_last(2:end-1),x);
    
    xlabel('Moment [Nm]')
    ylabel('Depth [m]')
    title('Moment, M');
    hold off
    
    
    subplot(1,5,4);
    hold on
    for vv = 1:NLoad/a:NLoad
        S_temp = S_cell{vv,ii};
        plot(S_temp(1:end-1),x);
        
        
    end
    S_last = S_cell{end,ii};
    plot(S_last(1:end-1),x);
    
    xlabel('Shear [N]')
    ylabel('Depth [m]')
    title('Shear, S');
    hold off
    
    subplot(1,5,5);
    hold on
    for vv = 1:NLoad/a:NLoad
        SR_P_temp = SR_P_cell{vv,ii};
        plot(SR_P_temp(1:end),x);
        
    end
    SR_P_last = SR_P_cell{end,ii};
    plot(SR_P_last(1:end),x);
    
    xlabel('Soil Reaction [N/m]')
    ylabel('Depth [m]')
    title('Soil Reaction, P');
    hold off

end