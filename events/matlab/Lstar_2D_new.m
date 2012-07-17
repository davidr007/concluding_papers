function [f,grad] = Lstar_2D(x)
% This function takes a bunch of locations and returns the value of L* and 
% its gradient at x (the location of N events with some minor
% modification - see below). 
%
% INPUTS: 
% x         [double - vector] A column vector containing the coordinates 
%           for the events of interest. Note that if there are N events
%           then length(x) = N-5 and has form 
%             x = [x2, x3,y3,x4,y4,...,xN,yN]'
%           Note that x1 = y1 = 0 sets the first event to be the
%           origin and y2 =0 sets the second event to be on the
%           arbitrarily oriented x-axis. The two contstraints remove 
%           non-uniqueness associated with translation and rotation 
%           of the event set, respectively. 
%
% OUTPUTS:
% f         [scalar] value of L* at x
% grad      [vector] this is the gradient of L* at x
% EXTRAS        **  Note that this function will search for a file called
%                   CWI_stat in the local directory and load it. This file 
%                   must contain a matrix with same name CWI_stat with 8 
%                   columns as given by 
%                       col1 = eventid1
%                       col2 = eventid2
%                       col3 = mean
%                       col4 = std
%                       col5 = median
%                       col6 = mu_N for non-zero Gaussian by fitting
%                       col7 = sigma_N for non-zero Gaussian by fitting
%                       col8 = fdom
%               **  Note that the order of the events in CWI_stat is very
%               important. CWI_stat(1,1) gives event1 (the master event).
%               CWI_stat(1,2) gives the second event which is fixed on the
%               local x axis. Then events 3,4,5,..,N are given by
%               CWI_stat(2:N,2). 
%
% David Robinson
% 16 June 2008


vs = 3300;

x = x(:);


% setup
xvec = [0; x(1); x(2:2:end)];
yvec = [0; 0; x(3:2:end)];
n = length(xvec);

load CWI_stat.mat
MUN = zeros(n,n);
SIGMAN = zeros(n,n);
FDOM = zeros(n,n);
U = zeros(n,n);   %An upper triangular matrix
count = 1;
for i = 1:n-1
    MUN(i,i+1:n) = CWI_stat(count:count+n-i-1,6);
    SIGMAN(i,i+1:n) = CWI_stat(count:count+n-i-1,7);
    FDOM(i,i+1:n) = CWI_stat(count:count+n-i-1,8);
    U(i,i+1:n) = ones(size([i+1:n]));
    count = count +n-i;
end

% Build an upper triangular matrix of the correct size.
% U = ones(n);
% for i = 1:n
%     U(i,1:i) = zeros(size([1:i]));
% end

XVEC_P = repmat(xvec,1,length(xvec));
XVEC_P = XVEC_P.*U;
XVEC_Q = repmat(xvec',length(xvec),1);
XVEC_Q = XVEC_Q.*U;
YVEC_P = repmat(yvec,1,length(yvec));
YVEC_P = YVEC_P.*U;
YVEC_Q = repmat(yvec',length(yvec),1);
YVEC_Q = YVEC_Q.*U;
DELTA_t = FDOM/vs.*sqrt((XVEC_P - XVEC_Q).^2 + (YVEC_P - YVEC_Q).^2);
%DELTA_t = DELTA_t.*U;

deriv_deltat_xp = zeros(n); 
deriv_deltat_yp = zeros(n); 

deriv_deltat_xp(U==1)= FDOM(U==1).^2./vs^2.* (XVEC_P(U==1)-XVEC_Q(U==1))./DELTA_t(U==1);
deriv_deltat_yp(U==1) = FDOM(U==1).^2./vs^2.* (YVEC_P(U==1)-YVEC_Q(U==1))./DELTA_t(U==1);



%%%%%%%%%%%%%

load('/export/storage/davidr/sandpit/davidr/thesis_version2/diags/revised_2dacoustics/foldederrorbar_fit3_normtype1.mat')
% loads: tw, xout_mean_twp75,  xout_sigma_twp75
MU1 = make_errorbar_mean(xout_mean_twp75, DELTA_t);
deriv_MU1_deltat = make_errorbar_mean_derivative(xout_mean_twp75, DELTA_t);
deriv_MU1_xp = deriv_MU1_deltat.*deriv_deltat_xp;
deriv_MU1_yp = deriv_MU1_deltat.*deriv_deltat_yp;
SIGMA1 = make_foldederrorbar_sigma(xout_sigma_twp75, DELTA_t);
SIGMA1 = U.*SIGMA1;
deriv_SIGMA1_deltat = make_errorbar_mean_derivative(xout_sigma_twp75, DELTA_t);
deriv_SIGMA1_xp = deriv_SIGMA1_deltat.*deriv_deltat_xp;
deriv_SIGMA1_yp = deriv_SIGMA1_deltat.*deriv_deltat_yp;

Phi_MU1_SIGMA_1_ZERO = zeros(n);
tmp_matrix =  normcdf(zeros(n),MU1,SIGMA1);
Phi_MU1_SIGMA_1_ZERO(U==1) = tmp_matrix(U==1);

% Phi_MUN_SIGMA_N_ZERO =  normcdf(zeros(n),MUN,SIGMAN);
% Phi_MUN_SIGMA_N_ZERO = U.*Phi_MUN_SIGMA_N_ZERO;


s = [-10^7 -10^6, -10^5, -10^4,-10^3,-10^2,-10:0.05:-1.05,-1:0.01:0];
%s = [-10:0.05:-1.05,-1:0.01:0];
deriv_Phi_xp = zeros(n,n);
deriv_Phi_yp = zeros(n,n);

delta_CWI = [0:0.01:1.2];

deriv_lnP_xp = zeros(n,n);
deriv_lnP_yp = zeros(n,n);

P_CWIN_deltat = zeros(n,n);
A = zeros(n); 
deriv_A_xp = zeros(n);
deriv_A_yp = zeros(n);

%g = zeros(n);
for i = 1:n-1
    for j = i+1:n
        g = -(s-MU1(i,j)).^2/(2*SIGMA1(i,j).^2);
%         deriv_g_xp_tmp = ( 2*SIGMA1(i,j).^2 *(-2*(s-MU1(i,j)).*(-deriv_MU1_xp(i,j))) + ...
%                             4*SIGMA1(i,j).*deriv_SIGMA1_xp(i,j)*(s-MU1(i,j)).^2)  ./ ...
%                             ( 4*SIGMA1(i,j).^4 );
        Phi_MU1_SIGMA_1_ZERO(i,j) = 1/(SIGMA1(i,j)*sqrt(2*pi))*trapz(s,exp(g));                
        deriv_g_xp_tmp = (4.*SIGMA1(i,j).^2.*(s-MU1(i,j)).*deriv_MU1_xp(i,j) + ...
                            +4.*SIGMA1(i,j).*(s-MU1(i,j)).^2.*deriv_SIGMA1_xp(i,j) )./ ...
                             ( 4*SIGMA1(i,j).^4 );
                        
        deriv_g_yp_tmp = (4.*SIGMA1(i,j).^2.*(s-MU1(i,j)).*deriv_MU1_yp(i,j) + ...
                            +4.*SIGMA1(i,j).*(s-MU1(i,j)).^2.*deriv_SIGMA1_yp(i,j) )./ ...
                             ( 4*SIGMA1(i,j).^4 );                 
                         
       
        deriv_Phi_xp(i,j) = 1/sqrt(2*pi).* (...
                             (SIGMA1(i,j).*trapz(s,exp(g).*deriv_g_xp_tmp) - ...
                             deriv_SIGMA1_xp(i,j).*trapz(s,exp(g))) ./ ...
                             (SIGMA1(i,j).^2));
                         
        deriv_Phi_yp(i,j) = 1/sqrt(2*pi).* (...
                             (SIGMA1(i,j).*trapz(s,exp(g).*deriv_g_yp_tmp) - ...
                             deriv_SIGMA1_yp(i,j).*trapz(s,exp(g))) ./ ...
                             (SIGMA1(i,j).^2));                
                         
        
        
        % The following can be computed in a vectorized fashion
        % However they needed below so keeping these outside the loop will 
        % mean that another loop must be created. Doing these here is
        % viewed as more efficient. 
        deriv_A_xp(i,j) = -(deriv_Phi_xp(i,j).*SIGMA1(i,j)*sqrt(2*pi) + ...
                (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*deriv_SIGMA1_xp(i,j)*sqrt(2*pi) ) ./   ...
                ( (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*SIGMA1(i,j)*sqrt(2*pi)).^2;
        deriv_A_yp(i,j) = -(deriv_Phi_yp(i,j).*SIGMA1(i,j)*sqrt(2*pi) + ...
                (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*deriv_SIGMA1_yp(i,j)*sqrt(2*pi) ) ./   ...
                ( (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*SIGMA1(i,j)*sqrt(2*pi)).^2;
       
        
        
        
        h = -(delta_CWI - MU1(i,j)).^2./(2*SIGMA1(i,j).^2);
        
        deriv_h_xp = (  -4.*SIGMA1(i,j).^2.*(delta_CWI - MU1(i,j)).*(deriv_MU1_xp(i,j)) - ...
                        -4.*(delta_CWI - MU1(i,j)).^2.*SIGMA1(i,j).*deriv_SIGMA1_xp(i,j))./ ...
                        (4*SIGMA1(i,j).^4);
        deriv_h_yp = (  -4.*SIGMA1(i,j).^2.*(delta_CWI - MU1(i,j)).*(deriv_MU1_yp(i,j)) - ...
                        -4.*(delta_CWI - MU1(i,j)).^2.*SIGMA1(i,j).*deriv_SIGMA1_yp(i,j))./ ...
                        (4*SIGMA1(i,j).^4);
                    
 
                    
        deriv_B_xp = -exp(h).*deriv_h_xp;
        deriv_B_yp = -exp(h).*deriv_h_yp;
        
        % Compute A
        A(i,j) = 1./((1-Phi_MU1_SIGMA_1_ZERO(i,j)).*SIGMA1(i,j)*sqrt(2*pi));
        
        % Compute B
        B = exp( -(delta_CWI- MU1(i,j)).^2 ./ (2*SIGMA1(i,j).^2));
        
        % Compute C
        [rubish,Phi_mu2_sigma2_0] = ...
            folded_normal_pdf(0,delta_CWI,SIGMAN(i,j));
        C = 1./((1-Phi_mu2_sigma2_0).*SIGMAN(i,j)*sqrt(2*pi));
        
        % Compute D
        D = exp( -(MUN(i,j)-delta_CWI).^2 ./ (2.*SIGMAN(i,j).^2));
        
        integrand = B.*C.*D;
        P_CWIN_deltat(i,j) = A(i,j).* trapz(delta_CWI,B.*C.*D);   
        
        deriv_P_xp =    deriv_A_xp(i,j).*trapz(delta_CWI,B.*C.*D) + ...
                        A(i,j).*trapz(delta_CWI,deriv_B_xp.*C.*D);
        deriv_P_yp =    deriv_A_yp(i,j).*trapz(delta_CWI,B.*C.*D) + ...
                        A(i,j).*trapz(delta_CWI,deriv_B_yp.*C.*D);
                    
        deriv_lnP_xp(i,j) = 1/P_CWIN_deltat(i,j).*deriv_P_xp;
        deriv_lnP_yp(i,j) = 1/P_CWIN_deltat(i,j).*deriv_P_yp;
          
        
       
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f = -sum(sum(A));
% grad_xp = -sum(deriv_A_xp - deriv_A_xp',2);  % note minus in the middle because we swap the order of p and q
% grad_yp = -sum(deriv_A_yp - deriv_A_yp',2);

f = -sum(sum(Phi_MU1_SIGMA_1_ZERO));
grad_xp = -sum(deriv_Phi_xp - deriv_Phi_xp',2);  % note minus in the middle because we swap the order of p and q
grad_yp = -sum(deriv_Phi_yp - deriv_Phi_yp',2);


grad = zeros(size(x));
grad(1)  =  grad_xp(2);  % i.e. becuase x1 = y1 = 0 (first event fixed at origin) and y2
count = 2;
for i = 3:n 
    grad(count) = grad_xp(i);
    grad(count+1) = grad_yp(i);
    count = count+2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%
% %Compute Lstar 
% f = -sum(sum(P_CWIN_deltat));
%             
% grad_xp = -sum(deriv_lnP_xp - deriv_lnP_xp',2);
% grad_yp = -sum(deriv_lnP_yp - deriv_lnP_yp',2);
% 
% grad = zeros(size(x));
% grad(1)  =  grad_xp(2);
% count = 2;
% for i = 3:n
%     grad(count) = grad_xp(i);
%     grad(count+1) = grad_yp(i);
%     count = count+2;
% end
            

            

            
            
            
            
            
            

      
      
      
            
            
            
                       
