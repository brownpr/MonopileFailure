tic
clear all;close all; clc;
global x_depth tt ODP k Ep Ms Fb Mb Px B phi gamma NLoad NN


plts = 10;                          %number of desired plots, will also plot inititial values
graphs = 10;                      	%number of desired graphs.

%SECTIONS
N = 50;                             %Number of desired pile sections 
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
NCycles = N;                        %Number of cycles
cycle = [Fs_min, Fs_max, Fs_min, -1*Fs_max, Fs_min];    %Cycling vector 
%y_vals = cell(NLoad, length(cycle)-1);

P_init = zeros(N+1, 1);
y_init = zeros(NN,1);

for j = 1:NCycles
    for i = 1:length(cycle)-1
        Fs_1 = cycle(i);
        Fs_2 = cycle(i+1);
        Fs = linspace(Fs_1, Fs_2, NLoad);
        [y_vals,k_vals,P_vals,R_vals,M_vals,S_vals,SR_P_vals] = callBeamState(Fs,N, P_init, y_init);
        
        y_cell{j,i} = y_vals;
        k_cell{j,i}  = k_vals;
        P_cell{j,i}  = P_vals;
        R_cell{j,i}  = R_vals;
        M_cell{j,i}  = M_vals;
        S_cell{j,i}  = S_vals;
        SR_P_cell{j,i}  = SR_P_vals;
        
        %Using initial positions and P for beam for initial iterations.
        P_init_vals = P_cell{j,i};
        y_init_vals = y_cell{j,i};
        
        P_init = P_init_vals(1:end,end);
        y_init = y_init_vals(1:end,end);
    end
end

'Cycle time'
t = toc

%PLOTS

%Max displacement plot
figure
hold on
y_axis = linspace(0,(-1)*x_depth,N+1);
for i = 1:NCycles/graphs:NCycles+1
    col = [rand, rand, rand];
    if i == 1
        i_ac = 1;
    else
        i_ac = i - 1;
    end
    for j = [1,3]
        y_plot_vals = y_cell{i_ac,j};
        x_axis = y_plot_vals(3:end-2,end);
        plot(x_axis,y_axis,'Color',col,'DisplayName',['Itr:' num2str(i)]);
    end
end
xline(0,'--k','DisplayName','y = 0'); 
legend show
legend('location','best') 
hold off


%displacement plot
figure
hold on
x = linspace(0,(-1)*x_depth,N+1);
for i = 1:NCycles/graphs:NCycles+1
    col = [rand, rand, rand];
    if i == 1
        i = 1;
    else
        i = i - 1;
    end
    for j = 1:4
        y_plot_vals = y_cell{i,j};
        for ii = 1:NLoad/plts:NLoad+1
            if ii == 1
                ii = 1;
            else
                ii = ii - 1;
            end
            y_axis = y_plot_vals(3:end-2,ii);
            plot(y_axis,x,'Color',col)
        end
    end

end
hold off

 
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