function [x ,y] = call_beamstate()

%Variables and constants 
Ep = 200*10^6;
Ip = 20;
Epy = 50*10^6;
x_depth = 50;
y_0 = [1,0,0,0];
xspan=[0,x_depth];


%ODE45 solver
[x,y] = ode45(@beamstate, xspan, y_0);
    %ODE
    function ydot = beamstate(x,y)
        ydot = zeros(4,1);
        ydot(1) = y(2);
        ydot(2) = y(3);
        ydot(3) = y(4);
        ydot(4) = (1/(Ep*Ip))*(-1*Epy*y(1));
    end

p = Epy*y;
plot(p,y);
xlabel('y'); ylabel('p');
legend('y1','y2','y3','y4','location','southeast')
end



