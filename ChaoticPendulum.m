clear;clc

%Simulation Parameters
N = 1000000; Delta = 1/1000;
tau_min = 0; tau_max = Delta*N;

%Variables
nuArray = [0.5]; 
l = 1; 
w_true = 2/3;  
A_trueArray = [1.2]; 
m = 1; 
g_Earth = 1; 
w_0 = g_Earth/l;

%Dimensionless variables
w = w_true/w_0; 

%Couple differential equations in phi and theta.
f = @(phi) phi;
g = @(phi, theta, tau, Q, A) -1/Q*phi - sin(theta) + A*cos(w*tau);

%Analytic solutions for small angles
theta = @(tau, A, Q, w) A*((1-w^2)*cos(w*tau) + w/Q*sin(w*tau)) / ((1-w^2)^2 + w^2/Q^2);

%Arrays for the equations of motion
tau = linspace(tau_min, tau_max, N);
phiArray = zeros(1,N); thetaArray = zeros(1,N);

nPoincarePoints = floor(tau_max*w / (2*pi));
PoincareTimes = zeros(length(nPoincarePoints), 1);
PoincareStep = floor((2*pi/w)/Delta);
index = 1;
for i = 1:nPoincarePoints
    PoincareTimes(i) = index;
    index = index + PoincareStep;
end

%Initial Conditions for theta (If there's more than one it will plot them
%all)
theta_0 = [0.2, 0.2000001];

%Plots various values of A
for h = 1:length(A_trueArray)

    A = A_trueArray(h)/(m*g_Earth); 

    %Plots various values of nu
    for r = 1:length(nuArray)

        Q = m*g_Earth/(nuArray(r)*w_0); 

        %Inital conditions
        for j = 1:length(theta_0)

            phiArray(1) = 0; thetaArray = theta_0(j);

            %Implement RK4
            for i = 1:N-1

                k1 = Delta * f(phiArray(i));
                l1 = Delta * g(phiArray(i), thetaArray(i), tau(i), Q, A);

                k2 = Delta * f(phiArray(i) + l1/2);
                l2 = Delta * g(phiArray(i) + l1/2, thetaArray(i) + k1/2, tau(i) + Delta*tau(i)/2, Q, A);

                k3 = Delta * f(phiArray(i) + l2/2);
                l3 = Delta * g(phiArray(i) + l2/2, thetaArray(i) + k2/2, tau(i) + Delta*tau(i)/2, Q, A);

                k4 = Delta * f(phiArray(i) + l3);
                l4 = Delta * g(phiArray(i) + l3, thetaArray(i) + k3, tau(i) + Delta*tau(i), Q, A);

                phiArray(i+1) = 1/6 * (l1+2*l2+2*l3+l4) + phiArray(i);
                thetaArray(i+1) = 1/6 * (k1+2*k2+2*k3+k4) + thetaArray(i);

            end

            figure(1)
            hold on
            plot(tau, thetaArray)
            xlabel('$\tau$', 'Interpreter', 'LaTeX')
            ylabel('$\theta(\tau)$', 'Interpreter', 'LaTeX')
            title('Pendulum Angle, $\Delta=1/1000$', 'Interpreter', 'LaTeX')
            hold off

            figure(2)
            hold on
            plot(phiArray, thetaArray)
            xlabel('$\dot{\theta}(\tau)$', 'Interpreter', 'LaTeX')
            ylabel('$\theta(\tau)$', 'Interpreter', 'LaTeX')
            title('Phase Space, $\Delta=1/1000$', 'Interpreter', 'LaTeX')
            hold off
            
            
            %Poincare Maps
            
            PoincareTheta = zeros(nPoincarePoints, 1);
            PoincarePhi = zeros(nPoincarePoints, 1);
            for i = 1:nPoincarePoints
                PoincareTheta(i) = thetaArray(PoincareTimes(i));
                PoincarePhi(i) = phiArray(PoincareTimes(i));
            end
            
            figure(3)
            hold on
            plot(PoincarePhi, PoincareTheta, 'o')
            xlabel('$\dot{\theta}(\tau)$', 'Interpreter', 'LaTeX')
            ylabel('$\theta(\tau)$', 'Interpreter', 'LaTeX')
            title('Poincare Map, A=1.2')
            hold off
            
            %Energy (Unforced)
            E = zeros(length(tau),1);
            for i = 1:length(tau)
                E(i) = 1 + 0.5*phiArray(i)^2 - cos(thetaArray(i));
            end
            figure(4)
            hold on 
            plot(tau, E)
            xlabel('$\tau$', 'Interpreter', 'LaTeX')
            ylabel('$E(\tau)$', 'Interpreter', 'LaTeX')
            title('Energy')
            hold off
            
            Analytic = zeros(length(tau),1);
            for i = 1:length(tau)
                Analytic(i) = theta(tau(i), A, Q, w);
            end
            
            %Comparison to analytic solution. Only holds for certain
            %parameters.
            figure(5)
            hold on
            plot(tau, thetaArray)
            plot(tau, Analytic, 'r--')
            title('Comparison to Analytic Solution')
            xlabel('$\tau$', 'Interpreter', 'LaTeX')
            ylabel('$\theta(\tau)$', 'Interpreter', 'LaTeX')
            hold off
        end
    end
end
    
    
    