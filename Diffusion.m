clear;clc

D = 1; %Diffusion coefficient
L = 1; %Length of rod
m = 100; %Number of points in time
n = 10; %Number of points in space

tmin = 0; tmax = 0.15;
xmin = 0; xmax = L;

x = linspace(xmin, xmax, n+1);
dx = x(2)-x(1);


t = linspace(tmin, tmax, m+1);
dt = t(2)-t(1);

if dt > D*dx^2/2
    disp('WARNING: Forward Euler Method is not stable.')
end

alpha = D * dt/dx^2;

% Explicit Forward Euler Algorithm (stable if dt < dx^2/(2D)).

U = zeros(length(x),length(t));

% Initial and boundary conditions
U(:,1) = 0; U(1,:) = 0; U(length(x),:) = 1;

for j = 1:m
    
    for i = 2:n
        
        U(i,j+1) = (1-2*alpha)*U(i,j) + alpha*(U(i+1,j) + U(i-1,j));

    end
    
end

Uf = U;

%{
figure(1)
surf(t, x, Uf)
xlabel('t')
ylabel('x')
zlabel('U')
title('Forward Euler')
%}

% Implicit Backward Euler Method

A = zeros(length(x)-2, length(x)-2);
b = zeros(n+1,1);

U = zeros(length(x),length(t));
U(:,1) = 0; U(1,:) = 0; U(length(x),:) = 1;

for i = 2:length(x)-3
    A(i,i-1) = -alpha;
    A(i,i) = 1+2*alpha;
    A(i,i+1) = -alpha;
end
A(1,1) = 1+2*alpha; A(1,2) = -alpha; A(length(x)-2, length(x)-3) = -alpha;
A(length(x)-2, length(x)-2) = 1+2*alpha;


Umid = zeros(length(x)-2, length(t));
b = zeros(length(x)-2, length(t)); 
for i = 1:length(t)
    b(1,i) = alpha*U(1,i);
    b(length(x)-2,i) = alpha*U(length(x),i);
end

for i = 1:m
    Umid(:,i+1) = inv(A)*Umid(:,i) + inv(A)*b(:,i);
end 

for i = 2:length(x)-1
    for j = 1:length(t)
        U(i,j) = Umid(i-1,j);
    end
end


Ub = U;

%{
figure(2)
surf(t, x, Ub)
%shading interp 
%   view(2);
xlabel('t')
ylabel('x')
zlabel('U')
title('Backward Euler')
%}


%Crank-Nicholson Method

U = zeros(length(x),length(t));
U(:,1) = 0; U(1,:) = 0; U(length(x),:) = 1;
    
A = zeros(length(x)-2, length(x)-2);
B = zeros(length(x)-2, length(x)-2);

for i = 2:length(x)-3
    A(i,i-1) = -alpha;
    A(i,i) = 2*(1+alpha);
    A(i,i+1) = -alpha;
end
A(1,1) = 2*(1+alpha); A(1,2) = -alpha; 
A(length(x)-2, length(x)-3) = -alpha;
A(length(x)-2, length(x)-2) = 2*(1+alpha);

for i = 2:length(x)-3
    B(i,i-1) = alpha;
    B(i,i) = 2*(1-alpha);
    B(i,i+1) = alpha;
end
B(1,1) = 2*(1-alpha); B(1,2) = alpha; B(length(x)-2, length(x)-3) = alpha;
B(length(x)-2, length(x)-2) = 2*(1-alpha);
    

Umid = zeros(length(x)-2, length(t));
b = zeros(length(x)-2, length(t)); 
for i = 1:length(t)
    b(1,i) = alpha*U(1,i);
    b(length(x)-2,i) = alpha*U(length(x),i);
end

for i = 1:m
    Umid(:,i+1) = inv(A)*B*Umid(:,i) + 2*inv(A)*b(:,i);
end 

for i = 2:length(x)-1
    for j = 1:length(t)
        U(i,j) = Umid(i-1,j);
    end
end

Ucn = U;


%{
figure(3)
surf(t, x, Ucn)
%shading interp 
%   view(2);
xlabel('t')
ylabel('x')
zlabel('U')
title('Crank-Nicholson')
%}


%Analytical

Analytical = zeros(n+1,m+1);

for i = 1:n+1
    for j = 1:m+1
        Analytical(i,j) = Uanalytical(x(i),t(j),L,5000);
    end
end


%{
figure(4)
surf(t,x,Analytical)
xlabel('t')
ylabel('x')
zlabel('U')
title('Analytical Solution')
zlim([0,1])
%}


%All methods on same plot
figure(1)
hold on

subplot(2,2,1)
surf(t, x, Uf)
xlabel('t', 'Interpreter', 'Latex')
ylabel('x', 'Interpreter', 'Latex')
zlabel('U(x,t)', 'Interpreter', 'Latex')
title('Forward Euler')

subplot(2,2,2)
surf(t, x, Ub)
title('Backward Euler')
xlabel('t', 'Interpreter', 'Latex')
ylabel('x', 'Interpreter', 'Latex')
zlabel('U(x,t)', 'Interpreter', 'Latex')

subplot(2,2,3)
surf(t, x, Ucn)
title('Crank-Nicholson')
xlabel('t', 'Interpreter', 'Latex')
ylabel('x', 'Interpreter', 'Latex')
zlabel('U(x,t)', 'Interpreter', 'Latex')

subplot(2,2,4)
surf(t,x,Analytical)
title('Analytical')
xlabel('t', 'Interpreter', 'Latex')
ylabel('x', 'Interpreter', 'Latex')
zlabel('U(x,t)', 'Interpreter', 'Latex')
hold off

%Compare analytical and numerical solutions

%Early time (still curved)
earlyTime = floor(length(t)/10);

U_AnalyticalEarly = U(:,earlyTime);
U_fEulerEarly = Uf(:,earlyTime);
U_bEulerEarly = Ub(:,earlyTime);
U_cnEarly = Ucn(:,earlyTime);


%Late time (flat steady state)
lateTime = floor(9*length(t)/10);

U_AnalyticalLate = U(:,lateTime);
U_fEulerLate = Uf(:,lateTime);
U_bEulerLate = Ub(:,lateTime);
U_cnLate = Ucn(:,lateTime);

%Compare them all on the same plot
figure(2)
subplot(2,1,1)
hold on 
plot(x,U_AnalyticalEarly, 'r')
plot(x,U_fEulerEarly, 'k--')
plot(x,U_bEulerEarly, 'g--')
plot(x,U_cnEarly, 'b--')
xlabel('x', 'Interpreter', 'Latex')
ylabel('U(x,t)', 'Interpreter', 'Latex')
title('Early Time Solution')
legend({'Analytical', 'Forward Euler', 'Backward Euler', ...
    'Crank-Nicholson'}, 'Location', 'northwest')

subplot(2,1,2)
hold on
plot(x,U_AnalyticalLate, 'r')
plot(x,U_fEulerLate, 'k--')
plot(x,U_bEulerLate, 'g--')
plot(x,U_cnLate, 'b--')
xlabel('x', 'Interpreter', 'Latex')
ylabel('U(x,t)', 'Interpreter', 'Latex')
title('Late Time Solution')
legend({'Analytical', 'Forward Euler', 'Backward Euler', ...
    'Crank-Nicholson'}, 'Location', 'northwest')


%Analytical function (1D)
function [result] = Uanalytical(x, t, L, nPoints)
temp = 0;
for l = 1:nPoints
    temp = temp + -2*(-1)^(l+1)/(pi*l*L) ...
        * sin(l*pi*x/L) * exp(-(l*pi/L)^2*t);
end
result = temp + x/L;
end

