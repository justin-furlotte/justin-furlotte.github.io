clear;clc

%Simulation parameters
L = 100; M = 100; J = 1; h = 0; kT = [1]; %linspace(0,10,100); 
N = L*M; %Number of spins
sweepsPerSimulation = 50000; %Iterate a large amount more than the number of spins
numberOfSimulations = 1;

%Define lots of empty arrays/matrices
E_static = zeros(numberOfSimulations,length(kT)); 
M_static = zeros(numberOfSimulations,length(kT));
Chi = zeros(1,length(kT));
C = zeros(1,length(kT));

E_transient = zeros(numberOfSimulations,sweepsPerSimulation); 
M_transient = zeros(numberOfSimulations,sweepsPerSimulation);
C_transient = zeros(1,sweepsPerSimulation); 
Chi_transient = zeros(1,sweepsPerSimulation);

acceptedConfigs = zeros(1,sweepsPerSimulation); counter = 0;


for l = 1:length(kT) %Perform Monte Carlo at each temperature
    
    beta = 1/kT(l);
    
    for p = 1:numberOfSimulations %Perform certain amount of simulations

        %Begin Monte Carlo
        x = 1:L; y = 1:M; [X,Y] = meshgrid(x,y); %Mesh grid for plotting later

        %Startw with empty lattice
        Snew = zeros(L,M);

        %Fill the lattice with random spins
        for i = 1:L
            for j = 1:M
                Snew(i,j) = RS;
            end
        end

        for k = 1:sweepsPerSimulation
            %Begin with an empty L by M lattice
            S = Snew;

            %Select a random spin from the lattice
            n = randi(L,1); m = randi(M,1);
            s = S(n,m);

            Eb = Espin(S,n,m,L,M,J,h); %Energy in this configuration of spins

            %Change the spin s_i |--> -s_i
            Strial = S;
            Strial(n,m) = -s;

            %Calculate the energy from the trial configuration and 
            %the change in energy from the old configuration
            Et = Espin(Strial,n,m,L,M,J,h); %trial energy
            dE = Et - Eb; %Change in energy due to this flip

            %See if the new lattice is acceptable of not
            Pacc = exp(-beta * dE); r = rand;
            if (dE < 0) || (r < Pacc)
                Snew = Strial; %Keep the new lattice, keep dE
                counter = counter+1;
                acceptedConfigs(k) = counter;
            else
                Snew = S; 
                dE = 0; %Keep the old lattice, dE must be zero
                acceptedConfigs(k) = counter;
            end

            %Animation
            
            if rem(k, 1000)==0
                figure(1)
                hold on
                pcolor(X,Y,Snew);
                drawnow;
                title('Magentization of Spin Lattice')
            end
            

            %Only measure the observables after they have had a chance
            %to equilibrate (i.e. after the desired number of cycles)
            if k == sweepsPerSimulation
                E_static(p,l) = H(Snew,L,M,J,h); %Per spin
                M_static(p,l) = Mmean(Snew); %Per spin
            end

        %Transient part (Plot this vs. number of cycles)
        E_transient(p,k) = H(Snew,L,M,J,h);
        M_transient(p,k) = Mmean(Snew); %Mmean(Snew);
        
        end        
    end
    %Calculate observables versus number of cycles
    Energy = mean(E_static); 
    Magnetization = mean(M_static);
    C(l) = var(E_static(:,l)) .* beta^2; %Really C/k
    Chi(l) = var(M_static(:,l)) .* beta;
end

%For plotting C and Chi versus cycles
for k = 1:sweepsPerSimulation
    C_transient(k) = (mean(E_transient(:,k).^2) - mean(E_transient(:,k)).^2) * beta^2;
    Chi_transient(k) = (mean(M_transient(:,k).^2) - mean(M_transient(:,k)).^2) * beta;
end




%Plot number of accepted configurations vs cycles
%{
figure(1)
plot(1:sweepsPerSimulation, acceptedConfigs)
title('Number of Accepted Configurations, $T=1$', 'Interpreter', 'LaTeX')
xlabel('Cycles')
ylabel('Configurations Accepted')
%}




%Plot probability of configuration
%{
E_stable = E_transient(:,4000:end);
figure(2)
histogram(E_stable./N, 'Normalization', 'probability')
xlabel('$E$', 'Interpreter', 'LaTeX')
ylabel('$P(E)$', 'Interpreter', 'LaTeX')
title(['Configuration Probability at $k_BT =$ ', num2str(kT)], ...
    'Interpreter', 'LaTeX')

figure(3)
plot(4000:6000, var(E_stable)/N)
xlabel('Cycles', 'Interpreter', 'LaTeX')
ylabel('$\sigma_E^2$', 'Interpreter', 'LaTeX')
title(['Energy Variance at $k_BT =$ ', num2str(kT)], ...
    'Interpreter', 'LaTeX')
%}




%Plot observables vs cycles
%{
figure(4)
subplot(2,2,1)
plot(1:sweepsPerSimulation, mean(E_transient)./N, 'b')
title('Energy')
xlabel('Cycles', 'Interpreter', 'LaTeX')
ylabel('$E/N$', 'Interpreter', 'LaTeX')

subplot(2,2,2)
plot(1:sweepsPerSimulation, mean(M_transient), 'b')
title('Magnetization')
xlabel('Cycles', 'Interpreter', 'LaTeX')
ylabel('$\langle |M| \rangle / N$', 'Interpreter', 'LaTeX')

subplot(2,2,3)
plot(1:sweepsPerSimulation, C_transient./N, 'b')
title('Heat Capaticity')
xlabel('Cycles', 'Interpreter', 'LaTeX')
ylabel('$C_V/(Nk_B)$', 'Interpreter', 'LaTeX')

subplot(2,2,4)
plot(1:sweepsPerSimulation, Chi_transient, 'b')
title('Magnetic Susceptibility')
xlabel('Cycles', 'Interpreter', 'LaTeX')
ylabel('$\chi/N$', 'Interpreter', 'LaTeX')
%}





%Plot observables vs temperature
figure(5)
subplot(2,2,1)
plot(kT, Energy./N, 'b')
title('Energy')
xlabel('$k_BT$', 'Interpreter', 'LaTeX')
ylabel('$E$', 'Interpreter', 'LaTeX')

subplot(2,2,2)
plot(kT, Magnetization, 'b')
title('Magnetization')
xlabel('$k_BT$', 'Interpreter', 'LaTeX')
ylabel('$\langle |M| \rangle / N$', 'Interpreter', 'LaTeX')

subplot(2,2,3)
plot(kT, C./N, 'b')
title('Heat Capaticity')
xlabel('$k_BT$', 'Interpreter', 'LaTeX')
ylabel('$C_V/k_B$', 'Interpreter', 'LaTeX')
ylim([0,2])

subplot(2,2,4)
plot(kT, Chi, 'b')
title('Magnetic Susceptibility')
xlabel('$k_BT$', 'Interpreter', 'LaTeX')
ylabel('$\chi$', 'Interpreter', 'LaTeX')






%Plot Numerical vs Analytical 
%{
E_true = zeros(1,length(kT));
M_true = zeros(1,length(kT));
Cv_true = zeros(1,length(kT));
Chi_true = zeros(1,length(kT));
for i = 1:length(kT)
    beta = 1/kT(i);
    E_true(i) = E_analytic(beta,J,h,N);
    M_true(i) = M_analytic(beta,J,h,N);
    Cv_true(i) = Cv_analytic(beta,J,h,N);
    Chi_true(i) = Chi_analytic(beta,J,h,N);
end

figure(4)
subplot(2,2,1)
hold on
plot(kT, Energy./N, 'b')
plot(kT, E_true, 'r')
hold off
title('Energy')
xlabel('$k_BT$', 'Interpreter', 'LaTeX')
ylabel('$E$', 'Interpreter', 'LaTeX')

subplot(2,2,2)
hold on
plot(kT, Magnetization, 'b')
plot(kT, M_true, 'r')
hold off
title('Magnetization')
xlabel('$k_BT$', 'Interpreter', 'LaTeX')
ylabel('$\langle |M| \rangle / N$', 'Interpreter', 'LaTeX')

subplot(2,2,3)
hold on
plot(kT, C./N, 'b')
plot(kT, Cv_true, 'r')
hold off
title('Heat Capaticity')
xlabel('$k_BT$', 'Interpreter', 'LaTeX')
ylabel('$C_V/k_B$', 'Interpreter', 'LaTeX')
ylim([0,0.5])

subplot(2,2,4)
hold on
plot(kT, Chi, 'b')
plot(kT, Chi_true, 'r')
hold off
title('Magnetic Susceptibility')
xlabel('$k_BT$', 'Interpreter', 'LaTeX')
ylabel('$\chi$', 'Interpreter', 'LaTeX')
ylim([0,0.05])
%}




%Analytic expression for the energy of the lattice (2x2 case)
function [result] = E_analytic(beta,J,h,N)

result = -(1/Z(beta,J,h)) * ((8*J+4*h)*exp(beta*(8*J+4*h)) + ...
    8*h*exp(2*beta*h) - 16*J*exp(-8*beta*J) - 8*h*exp(-2*beta*h) ...
    + (8*J - 4*h)*exp(beta*(8*J-4*h))) / N;

end

%Analytic expression for the mean magnetization of the lattice (2x2 case)
function [result] = M_analytic(beta,J,h,N)

result = 1/(beta * Z(beta,J,h)) * (4*beta*exp(beta*(8*J+4*h)) + ...
    8*beta*exp(2*beta*h) + 8*beta*exp(-2*beta*h) + ...
    4*beta*exp(beta*(8*J-4*h))) / N;

end

function [result] = Cv_analytic(beta,J,h,N)

result = (((8*J+4*h)^2*exp(beta*(8*J+4*h)) + ...
    16*h^2*exp(2*beta*h) - 128*J^2*exp(-8*beta*J) - 16*h^2*exp(-2*beta*h) ...
    + (8*J - 4*h)^2*exp(beta*(8*J-4*h))) * Z(beta,J,h) - ...
    ((8*J+4*h)*exp(beta*(8*J+4*h)) + ...
    8*h*exp(2*beta*h) - 16*J*exp(-8*beta*J) - 8*h*exp(-2*beta*h) ...
    + (8*J - 4*h)*exp(beta*(8*J-4*h)))^2) / Z(beta,J,h)^2 * beta^2/N;

end

function [result] = Chi_analytic(beta,J,h,N)

result = ((16*beta^2*exp(beta*(8*J+4*h)) + ...
    16*beta^2*exp(2*beta*h) + 16*beta^2*exp(-2*beta*h) + ...
    16*beta^2*exp(beta*(8*J-4*h)))*Z(beta,J,h) - ...
    ((4*beta*exp(beta*(8*J+4*h)) + ...
    8*beta*exp(2*beta*h) + 8*beta*exp(-2*beta*h) + ...
    4*beta*exp(beta*(8*J-4*h))))^2 ) / (beta * N^2 * Z(beta,J,h)^2);

end

%Analytic partition function (2x2 case)
function [result] = Z(beta,J,h)

result = exp(beta*(8*J+4*h)) + 4*exp(2*beta*h) + 4 + 2*exp(-8*beta*J) ...
    + 4*exp(-2*beta*h) + exp(beta*(8*J-4*h));

end

%Hamiltonian, J * sum_{i,j}(s_i*s_j) - h * sum_i(s_i)
function [result] = H(S,L,M,J,h)

%Sum over "nearest neighbours" at the spin S(n,m). Recall that
%   a = n-1, b = n+1, c = m-1, d = m+1.
neighbours = @(n,m) S(n,m)*S(a(n,L),m) + S(n,m)*S(b(n,L),m) + ...
    S(n,m)*S(n,c(m,M)) + S(n,m)*S(n,d(m,M));

%Sum over all nearest neighbors for each spin in the lattice 
%and sum over all spins. The factor of 1/2 accounts for the double counting.

P = 0; G = 0;
for i = 1:L
    for j = 1:M
        P = P + neighbours(i,j);
        G = G + S(i,j);
    end
end
        
result = -1/2 * J*P - h*G;

end

%Energy of a single spin S(n,m)
function [result] = Espin(S,n,m,L,M,J,h)

%Sum over "nearest neighbours" at the spin S(n,m). Recall that
%   a = n-1, b = n+1, c = m-1, d = m+1.
%The factor of 1/2 accounts for the double counting.
neighbours = (S(n,m)*S(a(n,L),m) + S(n,m)*S(b(n,L),m) + ...
    S(n,m)*S(n,c(m,M)) + S(n,m)*S(n,d(m,M)));
        
G = S(n,m);

result = -J*neighbours - h*G;

end

%Mean magnetization 
function [result] = Mmean(S)

result = abs(mean(mean(S)));

end

%The following four functions allow me to implement the periodic boundary
%conditions
function [result] = a(n,L)
if n == 1
    result = L;
else
    result = n-1;
end
end
function [result] = b(n,L)
if n == L
    result = 1;
else
    result = n+1;
end
end
function [result] = c(m,M)
if m == 1
    result = M;
else
    result = m-1;
end
end
function [result] = d(m,M)
if m == M
    result = 1;
else
    result = m+1;
end
end

%Random spins, either -1 or 1
function [result] = RS
x = rand;
if x<0.5
    result = -1;
else
    result = 1;
end
end