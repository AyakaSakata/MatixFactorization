% AMP for dictionary learning  
clear variables;
close all;

% dimension of signals
N = input('N = ');
% compression ratio
alpha = input('alpha = ');
% density of non-zero components
rho = input('rho = ');
% number of training samples/signal dimensions
gamma = input('gamma = ');

% number of measurements
M = ceil(alpha*N);
% number of training samples
P = ceil(gamma*N);
% number of non-zero components
K = round(N*rho);
% variance of the original non-zero signal
sigma_x2 = 1;
% noise strength
tau = 1e-10;

% damping factor
dmp = 0.1;
% number of samples
N_SAMPLE = 1;
% maximum number of BP iteration
STEP_MAX = 5000;

sqN = sqrt(N);
sqM = sqrt(M);

% initialization of quantities
m_F_mem = zeros(STEP_MAX,N_SAMPLE);
m_X_mem = m_F_mem;
q_F_mem = m_F_mem;
q_X_mem = m_F_mem;
MSE_y_mem = m_F_mem;
permind_mem = zeros(N,STEP_MAX,N_SAMPLE);

% Controling random number generator 
%rng(3);

for n_s = 1: N_SAMPLE

    % Problem set up 
    F0 = randn(M,N);
    F0_norm = sqrt(diag(F0'*F0));
    F0 = F0./(ones(M,1)*F0_norm')*sqM;

    X0 = sqrt(sigma_x2)*randn(N,P);
    non_active_set = zeros(K,P);
    for l = 1: P
        non_active_set = randperm(N,N-K);
        X0(non_active_set,l) = 0;
    end

    Y = F0*X0/sqN;

%%% Start BP 

    % initial condition
    F = F_initial(Y,N);
    X = pinv(F)*sqN*Y;
   
    %%% BP starting from random initial condition does not work!
    %F = randn(M,N);
    %F_norm = sqrt(diag(F*F'));
    %F = F./(F_norm*ones(1,N))*sqM;   
    %X = rho*sigma_x2*randn(N,P);
    %%%
    
    s_x = rho*ones(N,P);
    s_F = ones(M,N);
   
    [permindx,parity,overlap] = perm(F0,F,N);
    
    F_perm = F(:,permindx)*diag(parity);
    X_perm = diag(parity)*X(permindx,:);

    g_out = zeros(M,P);
    Qh_F = ones(N,1);
    Qh_F_mat = ones(M,1)*Qh_F';
    
    m_F_mem(1,n_s) = sum(diag(F0'*F_perm))/(N*M);
    m_X_mem(1,n_s) = sum(diag(X0*X_perm'))/(N*P);
    MSE_y_mem(1,n_s) = trace((F*X/sqN-Y)'*(F*X/sqN-Y))/(M*P);
    permind_mem(:,1,n_s) = permindx; 
    
    for bp_step = 1: STEP_MAX 

        X_old = X;
        F_old = F;
        g_out_old = g_out;
        
        %%% Macroscopic quantities
        q_x = mean(X(:).^2);
        Q_x = q_x+sum(s_x(:))/(N*P);
        
        q_F = mean(F(:).^2);
        Q_F = q_F+sum(s_F(:))/(M*N);
        
        %%% Row-dependent overlap
        q_x_vec = mean(X.^2,2);
        Q_x_vec = q_x_vec+mean(s_x,2);
        
        %%% Column-dependent overlap
        q_F_vec = (mean(F.^2))';
        Q_F_vec = q_F_vec+(mean(s_F))';
        
        % Variable -> Factor
        w = F*X/sqN-g_out_old*(q_F*(Q_x-q_x)+q_x*(Q_F-q_F));
        s_w = Q_F*Q_x-q_F*q_x;
        g_out = (Y-w)/(s_w+tau);
        g_out_p = -1/(s_w+tau);
        q_hat = trace(g_out'*g_out)/(M*P);
        chi_hat = -g_out_p;
        %chi_hat = q_hat;  %%% Nishimori condition
        
        % Factor -> Variable
        Sigma_x = (1./(alpha*(Q_F_vec*chi_hat-(Q_F_vec-q_F_vec)*q_hat)))*ones(1,P);
        Sigma_F = ones(M,1)*(1./(gamma*(Q_x_vec*chi_hat-(Q_x_vec-q_x_vec)*q_hat)))';
 
        mu_x = Sigma_x.*(F'*g_out/sqN+alpha*X*(q_F*chi_hat-(Q_F-q_F)*q_hat));
        mu_F = Sigma_F.*(g_out*X'/sqN+gamma*F*(q_x*chi_hat-(Q_x-q_x)*q_hat));

        rhoh_n = (1-rho)./sqrt(Sigma_x).*exp(-0.5*(mu_x.^2)./Sigma_x);
        rhoh_p = rho./sqrt(Sigma_x+sigma_x2).*exp(-0.5*(mu_x.^2)./(Sigma_x+sigma_x2));
        rhoh = rhoh_p./(rhoh_n+rhoh_p);
        
        % Matrix Variable -> Hyperparameter
        muF2_mean = ones(M,1)*mean((mu_F./Sigma_F).^2);
        Qh_F_mat = (1+sqrt(1+4*muF2_mean))/2-(1./Sigma_F);

        % Update variables
        X_new = rhoh.*mu_x*sigma_x2./(Sigma_x+sigma_x2);
        s_x_new = Sigma_x.*rhoh*sigma_x2./(Sigma_x+sigma_x2)+...
            ((rhoh_p.*rhoh_n)./(rhoh_p+rhoh_n).^2).*(mu_x*sigma_x2./(Sigma_x+sigma_x2)).^2;
        F_new = mu_F./(Qh_F_mat.*Sigma_F+1);
        s_F_new = Sigma_F./(Qh_F_mat.*Sigma_F+1);

        % damping
        X = X + dmp*(X_new-X);
        s_x = s_x + dmp*(s_x_new-s_x);
        F = F + dmp*(F_new-F);
        s_F = s_F + dmp*(s_F_new-s_F);       

        % fix the parity and the permutation indices
        [permindx,parity,overlap] = perm(F0,F,N);
            
        F_perm = F(:,permindx)*diag(parity);
        X_perm = diag(parity)*X(permindx,:);
            
        % time evolution of overlap between estimates and the true of F
        m_F_mem(bp_step+1,n_s) = sum(diag(F0'*F_perm))/(N*M);
        % time evolution of overlap between estimates and the true of X
        m_X_mem(bp_step+1,n_s) = sum(diag(X0*X_perm'))/(N*P);
        
        % time evolution of sg order parameter of F
        q_F_mem(bp_step+1,n_s) = mean(F(:).^2);
        % time evolution of sg order parameter of X
        q_X_mem(bp_step+1,n_s) = mean(X(:).^2);
        
        % time evolution of mse 
        MSE_y_mem(bp_step+1,n_s) = trace((F*X/sqN-Y)'*(F*X/sqN-Y))/(M*P);
        
        % time evolution of permutation index
        permind_mem(:,bp_step+1,n_s) = permindx; 

    end
    
end
