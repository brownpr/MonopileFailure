tic
clear all; clc;
global depth Ip ODP k Ep Ms Fb Mb Px B phi gamma NLoad NN C1 C2 C3 h  

plts = 5;                          %number of desired plots, will also plot inititial values
graphs = 5;                      	%number of desired graphs (cycles), will also plot first cycle.

%SECTIONS AND CYCLE VALUES
N=100;                          	%Number of desired pile sections 
NN=N+5;                             %Number of nodes
NE=N+4;                             %Number of elements
NLoad = 100;                      	%Separate scalling for Nload than N
NCycles = 200;                    	%Number of cycles

%PILE PARAMETERS
x_depth = 50;                       %Depth of pile [m]
tt = 0.06;                          %Thickness of pile [m]
ODP = 9;                            %Outer pile diameter [m] 
Ep = 190*10^9;                      %Modulus of pile, set as steel [N/m^2 or Pa]
Ip = (pi/4)*(((ODP/2)^4)-((ODP/2)-tt)^4);   %Second moment of Inertia of Pile
h = x_depth/N;                              %Height of elements
depth = linspace(-2*h,x_depth+2*h,NN);      %Depth vector of the monopile

%BOUNDARY CONDITION PARAMETERS
Fs_min = 0;                         %Min Force at surface [N]
Fs_max = 20*(10^6);                %Max Force at surface [N]
Ms = 0;                             %Moment at surface [Nm]
Fb = 0;                             %Force at base [N]
Mb = 0;                             %Moment at base [Nm]
Px = 0;                             %Applied lateral load [N]

%SOIL PARAMETERS
k = 30*10^6;                        %Elastic spring stiffness [N/m^2 or Pa]
B = 0.1;                            %soil constant
phi = 35;                           %soil constant for sand typically between 30 and 60 [deg]
gamma = 10^4;                       %specific weight of soil [N/m^3]
beta = 45 + phi/2;
alpha = phi/2;
K0 = 0.4;
Kp = (tand(45 + phi/2))^2;
Ka = (tand(45 - phi/2))^2;
C1 = tand(beta)*(Kp*tand(alpha) + K0*(tand(phi)*sind(beta)*(1/cosd(alpha) + 1) - tand(alpha)));
C2 = Kp - Ka;
C3 = (Kp^2)*(Kp + K0*tand(phi)) - Ka;

%CYCLIC STUFF
cycle = [Fs_min, Fs_max, Fs_min];    %Cycling vector 
%y_vals = cell(NLoad, length(cycle)-1);

%INITIALIZING CELLS
y_cell = cell(NCycles,2);
k_cell = cell(NCycles,2);
P_cell = cell(NCycles,2);
R_cell = cell(NCycles,2);
M_cell = cell(NCycles,2);
S_cell = cell(NCycles,2);
SR_P_cell = cell(NCycles,2);

%Initial conditions of system are zero. Else will be assined as last
%iterated values. 
P_init = zeros(NN,1);
y_init = zeros(NN,1);


for j = 1:NCycles
    cycNum = j; %Cycle number, used for error messages. 
    for i = [1,2]
        
        %Loading or unloading
%         if i == 2
%             LFF = 1;
%         else
%             LFF = -1;
%         end
        
        Fs_1 = cycle(i);
        Fs_2 = cycle(i+1);
        Fs = linspace(Fs_1, Fs_2, NLoad);
        [y_vals,k_vals,P_vals,R_vals,M_vals,S_vals,SR_P_vals] = callBeamState(Fs, N, P_init, y_init, cycNum);
        
        %Storing values in cells. 
        y_cell{j,i} = y_vals;
        k_cell{j,i} = k_vals;
        P_cell{j,i} = P_vals;
        R_cell{j,i} = R_vals;
        M_cell{j,i} = M_vals;
        S_cell{j,i} = S_vals;
        SR_P_cell{j,i} = SR_P_vals;
        
        %Collecting initial values for next cycle.
        P_init_vals = P_cell{j,i};
        y_init_vals = y_cell{j,i};
        P_init = P_init_vals(:,end);
        y_init = y_init_vals(:,end);
    end
end

'Cycle time'
t = toc

%PLOTS

%Max vals plot
figure
y_axis = linspace(0,(-1)*x_depth,N+1);
colours = ['k','b','r','m','c'];
cnt = 1;
for i = 1:NCycles/graphs:NCycles+1
    if cnt < 6
        col = colours(cnt);
        cnt = cnt +1;
    else 
        col = [rand, rand, rand];
    end
        
    if i == 1
        i_ac = 1;
    else
        i_ac = i - 1;
    end
    for j = [1,2]
        if j == 1
            dir = 'loading';
            lin = '-';
        elseif j == 2
            dir ='unloading';
            lin = '--';
        end
        for val = 1:5
            subplot(1,5,val);
            hold on
            if val == 1
                y_plot_vals = y_cell{i_ac,j};
                x_axis = y_plot_vals(3:end-2,end); 
                tit = 'Displacement';
                lab = 'Y Displacement of pile, [m]';
            elseif val == 2
                y_plot_vals = R_cell{i_ac,j};
                x_axis = y_plot_vals(1:end,end); 
                tit = 'Rotation';
                lab = 'Rotation of pile, [deg]';
            elseif val == 3
                y_plot_vals = M_cell{i_ac,j};
                x_axis = y_plot_vals(1:end,end); 
                tit = 'Moment';
                lab = 'Moment across pile, [Nm]';
            elseif val == 4
                y_plot_vals = S_cell{i_ac,j};
                x_axis = y_plot_vals(1:end,end); 
                tit = 'Shear';
                lab = 'Shear across pile, [Pa]';
            elseif val == 5
                y_plot_vals = SR_P_cell{i_ac,j};
                x_axis = y_plot_vals(1:end,end); 
                tit = 'Soil Reaction';
                lab = 'Soil Reaction across pile, [N/m^2]';
            end
            plot(x_axis,y_axis,'Color',col,'linestyle',lin,'DisplayName',['Cycle:' num2str(i_ac) ' ' dir ], 'LineWidth', 1.5);
            xlabel(lab,'FontSize',10,'FontWeight','bold');
            ylabel('Depth [m]','FontSize',10,'FontWeight','bold');
            title(tit)
            %xline(0,':k'); 
            hold off
        end
    end
end
legend show
legend('location','best') 

%p-y plot for different cycle at depth dd (loading only)

% dd = 10;
% hold on
% for i = 1:NCycles %/plts:NCycles+1
%     figure
%     if cnt < 6
%         col = colours(cnt);
%         cnt = cnt +1;
%     else
%         col = [rand, rand, rand];
%     end
% 
%     hold on
%     p_temp = P_cell{i,1};
%     y_temp = y_cell{i,1};
%     y_axis = p_temp(dd,:);
%     x_axis = y_temp(dd,:);
%     plot(x_axis,y_axis,'Color',col);
%     hold off
% end
% hold off

% %displacement plot
% figure
% hold on
% x = linspace(0,(-1)*x_depth,N+1);
% for i = 1:NCycles/graphs:NCycles+1
%     col = [rand, rand, rand];
%     if i == 1
%         i = 1;
%     else
%         i = i - 1;
%     end
%     for j = 1:4
%         y_plot_vals = y_cell{i,j};
%         for ii = 1:NLoad/plts:NLoad+1
%             if ii == 1
%                 ii = 1;
%             else
%                 ii = ii - 1;
%             end
%             y_axis = y_plot_vals(3:end-2,ii);
%             plot(y_axis,x,'Color',col)
%         end
%     end
% 
% end
% hold off

% %Animated displacement plot 
% figure
% hold on
% x = linspace(0,(-1)*x_depth,N+1);
% 
% for i = 1:NCycles/graphs:NCycles+1
%     col = [rand, rand, rand];
%     if i == 1
%         i = 1;
%     else
%         i = i - 1;
%     end
%     curve = animatedline('Color',col);
%     for j = 1:4
%         y_plot_vals = y_cell{i,j};
%         for ii = 1:NLoad/plts:NLoad+1
%             if ii == 1
%                 ii = 1;
%             else
%                 ii = ii - 1;
%             end
%             y_axis = y_plot_vals(3:end-2,ii);
%             for uu = 1:length(y_axis)
%                 addpoints(curve,y_axis(uu),x(uu));
%                 drawnow
%             end
%         end
%     end
% end
% hold off

% % %Animated displacement plot and save as vid, 
% % NOT WORKING as too many loops
% 
% figure
% hold on
% x = linspace(0,(-1)*x_depth,N+1);
% cnt = 1;
% for i = 1:NCycles/graphs:NCycles+1
%     cnt = 1
%     col = [rand, rand, rand];
%     if i == 1
%         i = 1;
%     else
%         i = i - 1;
%     end
%     curve = animatedline('Color',col);
%     for j = 1:4
%         y_plot_vals = y_cell{i,j};
%         for ii = 1:NLoad/plts:NLoad+1
%             if ii == 1
%                 ii = 1;
%             else
%                 ii = ii - 1;
%             end
%             y_axis = y_plot_vals(3:end-2,ii);
%             for uu = 1:length(y_axis)
%                 addpoints(curve,y_axis(uu),x(uu));
%                 drawnow
%                 pause(0.001); %time betweeen each point plot
%                 
%                 frm = [(cnt-1)*(length(y_axis)*4*plts)+cnt:(cnt)*(length(y_axis)*4*plts)+cnt];
%                 for vv = 1:length(frm)
%                     u = frm(vv);
%                     Frames(u) = getframe(gcf);
%                 end
%             end
%         end
%     end
% end
% hold off
% video = VideoWriter('CyclicLoading.avi','Uncompressed AVI');
% video.FrameRate=60;
% open(video)
% writeVideo(video, Frames);
% close(video)

'Plot time'
plotime= toc - t