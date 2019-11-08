%%% driving function %%%
function minimise_3_level_system_2
%defining some constants
global h epsilon0 clight hbar N ti tf tspan options 
h=4.135667517e-15;
hbar=h/(2*pi);
clight=299792458;
epsilon0=h*clight; %energy scaling
%%% for comparison %%%
options=odeset('RelTol',1e-2);
N=1000;
ti=0;
tf=0.0000000001;
tspan=linspace(ti,tf,N);

global psi8_full

fid = fopen('prob_psi_8_full.txt');
if (fid < 0)
    error('could not open file "prob_psi_8_full.txt"');
end
prob_psi_8_full = fscanf(fid, '%g %g', [2 inf]);
fclose(fid);
psi8_full=reshape(prob_psi_8_full,[],1);
tic
% minimisation 
ham_vec=10*[10,10,1,1]; %elements 1,2,3,4 are E1,E2,omega13,omega23 respectively

options=optimset('TolFun',1e-2);
options=optimset('TolX',1e-1);
[opt_ham,fval,exitflag,output]=fminsearch(@ham_iter,ham_vec,options);

toc
disp(fval)
global ham gamma;

gamma=[10 0 0; 0 10 0; 0 0 10];
gamma=0; %no dephasing
 
%build hamiltonian matrix

row1=[opt_ham(1)   opt_ham(3)   0];
row2=[opt_ham(3)   -1i*(5/1.88)   opt_ham(4)];
row3=[0            opt_ham(4)   opt_ham(2)];

ham=[row1;row2;row3];
disp(ham)
ham=epsilon0*100*ham/hbar;
% ham=epsilon0*100*ham/hbar;
% ham=hbar*ham/(100*epsilon0);
%size of energy increment

        
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
t=1e12*t;
plot_vec=[prob_psi_8, psi8_full];
% fprintf(t,plot_vec);
% % 
% fileID = fopen('plot_vec','w');
% fprintf(fileID,'%f',t,plot_vec);
% fclose(fileID);

fid = fopen('plot_vec.txt_2', 'w');
if (fid < 0)
    error('could not open file "plot_vec.txt_2"');
end
fprintf(fid, '%12.8f %12.8f\r\n', prob_psi_8);

% fid = fopen('plot_vec.txt_2', 'w');
% if (fid < 0)
%     error('could not open file "plot_vec.txt_2"');
% end
% fprintf(fid, '%12.8f %12.8f\r\n', [prob_psi_8, psi8_full]);


% disp(prob_psi_8)
% disp(psi8_full)
fclose(fid);
% disp(min_avgdiff)
% disp(fval)

figure(1)
clf
loglog(t,plot_vec)
xlabel('Time (ps)')
ylabel('|\psi_{8}|^{2}')
legend('3-level system','7-level system','location','northwest');

end

%%% comparison function %%%
function avgdiff=ham_iter(ham_vec)

global psi8_full epsilon0 hbar N tspan options min_avgdiff

global ham gamma;

gamma=[1 0 0; 0 1 0; 0 0 1];
gamma=0; %no dephasing
 
%build hamiltonian matrix

row1=[ham_vec(1)   ham_vec(3)   0];
row2=[ham_vec(3)   -1i*5/1.88   ham_vec(4)];
row3=[0            ham_vec(4)   ham_vec(2)];

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
min_avgdiff=100;
   
%     disp('atb')
    
    fid = fopen('diff_vec.txt', 'w');
    if (fid < 0)
        error('could not open file');
    end
    fprintf(fid, '%12.8f %12.8f\r\n', diff);



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
