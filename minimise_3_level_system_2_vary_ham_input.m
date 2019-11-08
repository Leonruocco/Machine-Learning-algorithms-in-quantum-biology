%%% driving function %%%
function minimise_3_level_system_2_vary_ham_input
%defining some constants
global h epsilon0 clight hbar N ti tf tspan options a
h=4.135667517e-15;
hbar=h/(2*pi);
clight=299792458;
epsilon0=h*clight; %energy scaling
%%% for comparison %%%
options=odeset('RelTol',1e-4);
N=1000;
ti=0;
tf=0.00000000005;
tspan=linspace(ti,tf,N);
a=31;

options=optimset('TolFun',1e-4);
options=optimset('TolX',1e-3);

global psi8_full

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data to supervise ML algorithm 
fid = fopen('prob_psi_8_full.txt');
if (fid < 0)
    error('could not open file "prob_psi_8_full.txt"');
end
prob_psi_8_full = fscanf(fid, '%g %g', [2 inf]);
fclose(fid);
psi8_full=reshape(prob_psi_8_full,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ham_vec=zeros(a,4);
fvalvec=zeros(1,a);

for n=2:a
    tic
    
    % define parameter space to explore in optimisiation routine
    ham_input=n*[10,10,2,2]; %elements 1,2,3,4 are E1,E2,omega13,omega23 respectively
    
    % call optimisation routine
    [opt_ham,fval,exitflag,output]=fminsearch(@ham_iter,ham_input,options);
%     disp(fval)
%     disp(exitflag)
%     disp(output)
    
    global ham gamma;
    
    gamma=[1 0 0; 0 1 0; 0 0 1];
    gamma=0; %no dephasing
    
    %build hamiltonian matrix from optimisation parameters
    %copy and pasted from plenio pdf
    
    row1=[opt_ham(1)   opt_ham(3)   0];
    row2=[opt_ham(3)   -1i*(5/1.88)   opt_ham(4)];
    row3=[0            opt_ham(4)   opt_ham(2)];
    
    ham=[row1;row2;row3];
    disp(ham)
    ham=epsilon0*100*ham/hbar; %scale
    % ham=epsilon0*100*ham/hbar;
    % ham=hbar*ham/(100*epsilon0);
    %size of energy increment
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve schrodinger equation using test Hamiltonian
    % to produce final output???
    %wavefunction to be acted upon by hamiltonian
    rho_0=zeros(9,1);
    rho_0(1,1)=1;
    
    [t,y]=ode45(@fmo7,tspan,rho_0,options);
    
    for k=1:N
        for x=1:3
            %arithmetic progression to access
            %every 4th element of vector
            rhodot(k,x)=y(k,((x-1)*4)+1);
        end
    end
    
    norm=sum(rhodot,2);
    prob_psi_8=1-norm;
    ham_vec(n,:)=opt_ham;
    fvalvec(1,n)=fval;
    
    toc
end

t=1e12*t;
plot_vec=[prob_psi_8, psi8_full];

figure(1)
clf
plot(t,plot_vec)
xlabel('Time (ps)')
ylabel('|\psi_{8}|^{2}')
legend('3-level system','7-level system','location','northwest');

xlswrite('ham_vec.xls', ham_vec);
dlmwrite('fmo_coupled_unscaled_mintec_vary_gamma.dlm',fvalvec);
% open the file with write permission
fid = fopen('fmo_coupled_unscaled_mintec_vary_gamma_fvalvec.txt', 'w');
if (fid < 0)
    error('could not open file "fmo_coupled_unscaled_mintec_vary_gamma.txt"');
end
fprintf(fid, '%12.8f %12.8f\r\n', fvalvec);
fclose(fid);

end

%%% comparison function %%%
function avgdiff=ham_iter(ham_input)

global psi8_full epsilon0 hbar N tspan options

global ham gamma;

gamma=[1 0 0; 0 1 0; 0 0 1];
gamma=0; %no dephasing
 
%build hamiltonian matrix
%copy and pasted from plenio pdf

row1=[ham_input(1)   ham_input(3)   0];
row2=[ham_input(3)   -1i*5/1.88   ham_input(4)];
row3=[0            ham_input(4)   ham_input(2)];

ham=[row1;row2;row3];
ham=epsilon0*100*ham/hbar;
        
        %wavefunction to be acted upon by hamiltonian
        rho_0=zeros(9,1);
        rho_0(1,1)=1;
        
        %solve calling function below
        
        [t,y]=ode45(@fmo7,tspan,rho_0,options);

        for k=1:N
            for x=1:3
                %arithmetic progression to access
                %every 4th element of vector
                rhodot(k,x)=y(k,((x-1)*4)+1);
            end
        end
        
        norm=sum(rhodot,2);
        prob_psi_8=1-norm;


%%% compare %%%

diff=abs(psi8_full-prob_psi_8);
num=numel(prob_psi_8);
avgdiff=sum(diff)/num;
% disp(avgdiff)


end

%%% math function %%%
function rhodot=fmo7(t,y)

global ham gamma;

rho=reshape(y,3,3);    

for n=1:3
    rhodiag(n,n)=rho(n,n);
end

rho_t=-1i*(ham*rho-rho*ham')-gamma*rho+gamma*rhodiag; %master equation

rhodot=reshape(rho_t,9,1);

end