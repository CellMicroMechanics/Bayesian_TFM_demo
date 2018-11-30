% This program calculates cellular traction forces using Bayesian regularization
% The program solves the linear relation u = Xf + noise_u, where X is a matrix and u, f 
% and noise_u are vetors.
%
% The code also requires the free Matlab package "Regularization tools" by P.C.
% Hansen Math. Model. Comput.(1998) that can be found at
% https://www.mathworks.com/matlabcentral/fileexchange/52-regtools.
% To use this tool, download the package and add it to the Matlab path
% ,e.g, by entering the command "path(path,'Hansen')".
%
% Yunfei Huang and Benedikt Sabass, Forschungszentrum Jülich, 2018
% Reference: Huang et al. https://arxiv.org/abs/1810.05848.

function  F  = reconstruction_BL2(X1, u1, noise_u,beta1)
%%% > X is the matrix connecting discretized traction to discretized
%%% displacement.
%%% > u is the vector of gel displacements caused by cells.
%%% > noise_u is the noise vector. Here you can input displacement
%%% measurements that are done far away from cells.
%%% >beta is an optional input.

% The inverse variance of the noise can either be prescribed as input
% or it can be calculated from noise_u
if nargin < 4
   beta1 = 1/var(noise_u); 
end

global beta X u U s V XX aa C
beta=beta1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% standardize the input data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
sd = std(X1);
X= (X1-repmat(mean(X1,1),size(X1,1),1))./repmat(std(X1),size(X1,1),1);
u  = u1-mean(u1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% singular value decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('SVD');
[U,s,V] = csvd(X);
XX=X'*X;
aa=size(X); 
c=ones(aa(2),1);
C=diag(c);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% search for optimal regularization parameter lambda_2=alpha/beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We perform a Golden section search method to determine alpha. Note that alpha must lie in the 
% interval between alpha1 and alpha2. 
alpha1 = 200;      %  initial left alpha
alpha2 = 40000;    %  initial right alpha
alpha_opt = fminbnd(@logevidence,alpha1,alpha2);


%%%% Plot the evidence
plot_alpha = [alpha_opt*0.03:alpha_opt*0.02:alpha_opt*2.5]; 
a = size(plot_alpha);
lambda_p = plot_alpha./beta;

   for i = 1:a(2)
     evidence(i) = -logevidence(plot_alpha(i));
   end
 
figure(1);
plot(lambda_p, evidence,'o');
xlabel('\lambda_2');
ylabel('Log(Evidence)');
hold on;
lambda_2 = alpha_opt/beta;
evidence_one = -logevidence(alpha_opt);
plot(lambda_2,evidence_one,'r*','MarkerSize',10);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% L2 regularization and rescaling the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 [F] = tikhonov(U,s,V,u,lambda_2);
  F = F./sd';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% calculating log evidence function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function evidence_value= logevidence(alpha)
% Equation of log(Evidence) can be found in ref. Huang et al. Sci. Rep. (2018).

global beta X u U s V XX aa C
flambda = alpha/beta;
[F] = tikhonov(U,s,V,u,flambda); 
   
% calculate log(det(A))
A = alpha*C+beta*XX;
L = chol(A);
logdetA = 2*sum(log(diag(L)));
   
% calculate log(Evidence)
evidence_value = -(-0.5*alpha*F'*F - 0.5*beta*(X*F-u)'*(X*F-u) ...
                 -0.5*logdetA + 0.5*aa(1)*log(beta) + 0.5*aa(2)*log(alpha)...
                 -0.5*aa(1)*log(2*pi));
end