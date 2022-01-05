clear;clc

D = 1; %Diffusion coefficient
L = 1; %Length of rod
m = 2000; %Number of points in time
n = 30; %Number of points in space

tmin = 0; tmax = 0.5;
xmin = 0; xmax = L;

x = linspace(xmin, xmax, n+1);
dx = x(2)-x(1);


t = linspace(tmin, tmax, m+1);
dt = t(2)-t(1);

if dt > D*dx^2/2
    disp('WARNING: Forward Euler Method is not stable.')
end

alpha = D * dt/dx^2;


% Forward Euler 2D
U = zeros(length(x),length(x));

% Initial and boundary conditions
U(1,:) = 0; U(:,1) = 0; U(:,length(x)) = 1; U(length(x),:) = 1; U(17:23,17:23)=0;

U2{1} = U;

for k = 1:length(t)
    for i = 2:length(x)-1
        for j = 2:length(x)-1
            U(i,j) = U(i,j) + alpha*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)-4*U(i,j));
        end
    end
    U2{k+1} = U;
end

figure(1)
for k = 1:length(t)
    M = U2{k};
    surf(x,x,M)
    zlim([0,1])
    pause(0.0001)
end
title('2D Diffusion Animation')
xlabel('x')
ylabel('y')
zlabel('U')

figure(2)
hold on
 
k = floor(0.01*length(t)/10);
subplot(2,2,1)
M = U2{k};
surf(x,x,M)
zlim([0,1])
xlabel('y', 'Interpreter', 'Latex')
ylabel('x', 'Interpreter', 'Latex')
zlabel('U(x,t)', 'Interpreter', 'Latex')
title(['Forward Euler at $t = $',  num2str(t(k))], ...
    'Interpreter', 'Latex')
 
k = floor(0.1*length(t)/10);
subplot(2,2,2)
M = U2{k};
surf(x,x,M)
zlim([0,1])
xlabel('y', 'Interpreter', 'Latex')
ylabel('x', 'Interpreter', 'Latex')
zlabel('U(x,t)', 'Interpreter', 'Latex')
title(['Forward Euler at $t = $',  num2str(t(k))], ...
    'Interpreter', 'Latex')
 
k = floor(0.6*length(t)/10);
subplot(2,2,3)
M = U2{k};
surf(x,x,M)
zlim([0,1])
xlabel('y', 'Interpreter', 'Latex')
ylabel('x', 'Interpreter', 'Latex')
zlabel('U(x,t)', 'Interpreter', 'Latex')
title(['Forward Euler at $t = $',  num2str(t(k))], ...
    'Interpreter', 'Latex')
 
k = floor(9*length(t)/10);
subplot(2,2,4)
M = U2{k};
surf(x,x,M)
zlim([0,1])
xlabel('t', 'Interpreter', 'Latex')
ylabel('x', 'Interpreter', 'Latex')
zlabel('U(x,t)', 'Interpreter', 'Latex')
title(['Forward Euler at $t = $',  num2str(t(k))], ...
    'Interpreter', 'Latex')





%Analytical
Analytical = zeros(n+1,n+1);

%U3{1} = U;


%for k = 1:length(t)
%    for i = 1:length(x)
%        for j = 1:length(x)
%            Analytical(i,j) = Uanalytical(x(i),x(j),t(k),L,300); 
%        end
%    end
%    U3{k+1} = Analytical;
%end


%Plot analytical solution
%figure(3)
%surf(x,x,M)




%Analytical function (2D)
function [result] = Uanalytical(x, y, t, L, nPoints)
temp = 0;
for l = 1:nPoints
    for k = 1:nPoints
        temp = temp + 8/(k*l*pi^2*L^2)*((-1)^(l+1)*(1-(-1)^k)+(1-(-1)^l)*(-1)^(k+1))...
            * sin(l*pi*x/L) * sin(k*pi*y/L) ...
            * exp(-(l^2+k^2)*pi^2/L^2*t);
    end
end
result = temp;
end



