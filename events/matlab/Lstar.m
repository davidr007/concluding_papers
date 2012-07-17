function [f,grad] = Lstar(x)
% This function takes a bunch of locations and returns the value of L* and 
% its gradient at x (the location of N events with some minor
% modification - see below). 
%
% INPUTS: 
% x         [double - vector] A column vector containing the coordinates 
%           for the events of interest. Note that if there are N events
%           then length(x) = N-5 and has form 
%             x = [x2, x3,y3,z3,x4,y4,z4,...,xN,yN,zN]'
%           Note that x1 = y1 = z1 = 0 sets the first event to be the
%           origin and y2 = z2 =0 sets the second event to be on the
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


%vs = 3300;

x = x(:);


% setup
xvec = [0; x(1);x(2); x(4:3:end)];
yvec = [0; 0; x(3); x(5:3:end)];
zvec = [0; 0; 0;x(6:3:end)];
n = length(xvec);


load CWI_stat.mat  % loads in CWI_stat and events_oi
MUN = zeros(n,n);
SIGMAN = zeros(n,n);
FDOM = zeros(n,n);
U = zeros(n,n);   %An upper triangular matrix
VS = zeros(n,n);


for i = 1:length(events_oi)
   ind1 = find(CWI_stat(:,1) == events_oi(i));
   for j = 1:length(ind1)
       row = i;
       col = find(events_oi==CWI_stat(ind1(j),2));
       MUN(row,col) =  CWI_stat(ind1(j),6);
       SIGMAN(row,col) = CWI_stat(ind1(j),7);
       FDOM(row,col) = CWI_stat(ind1(j),8);
       U(row,col) = 1;
       VS(row,col) = CWI_stat(ind1(j),9);
   end
end


% count = 1;
% for i = 1:n-1
%     disp(i)
%     MUN(i,i+1:n) = CWI_stat(count:count+n-i-1,6);
%     SIGMAN(i,i+1:n) = CWI_stat(count:count+n-i-1,7);
%     FDOM(i,i+1:n) = CWI_stat(count:count+n-i-1,8);
%     U(i,i+1:n) = ones(size([i+1:n]));
%     count = count +n-i;
% end

whos
XVEC_P = repmat(xvec,1,length(xvec));
XVEC_P = XVEC_P.*U;
XVEC_Q = repmat(xvec',length(xvec),1);
XVEC_Q = XVEC_Q.*U;
YVEC_P = repmat(yvec,1,length(yvec));
YVEC_P = YVEC_P.*U;
YVEC_Q = repmat(yvec',length(yvec),1);
YVEC_Q = YVEC_Q.*U;
ZVEC_P = repmat(zvec,1,length(zvec));
ZVEC_P = ZVEC_P.*U;
ZVEC_Q = repmat(zvec',length(zvec),1);
ZVEC_Q = ZVEC_Q.*U;
DELTA_t = FDOM./VS.*sqrt((XVEC_P - XVEC_Q).^2 + (YVEC_P - YVEC_Q).^2 + (ZVEC_P - ZVEC_Q).^2);
%DELTA_t = DELTA_t.*U;

deriv_deltat_xp = FDOM.^2./VS.^2.* (XVEC_P-XVEC_Q)./DELTA_t;
deriv_deltat_yp = FDOM.^2./VS.^2.* (YVEC_P-YVEC_Q)./DELTA_t;
deriv_deltat_zp = FDOM.^2./VS.^2.* (ZVEC_P-ZVEC_Q)./DELTA_t;

load('/export/storage/davidr/sandpit/davidr/thesis_version2/diags/revised_2dacoustics/foldederrorbar_fit3_normtype1.mat')
% loads: tw, xout_mean_twp75,  xout_sigma_twp75
MU1 = make_errorbar_mean(xout_mean_twp75, DELTA_t);
deriv_MU1_deltat = make_errorbar_mean_derivative(xout_mean_twp75, DELTA_t);
deriv_MU1_xp = deriv_MU1_deltat.*deriv_deltat_xp;
deriv_MU1_yp = deriv_MU1_deltat.*deriv_deltat_yp;
deriv_MU1_zp = deriv_MU1_deltat.*deriv_deltat_zp;
SIGMA1 = make_foldederrorbar_sigma(xout_sigma_twp75, DELTA_t);
SIGMA1 = U.*SIGMA1;
deriv_SIGMA1_deltat = make_errorbar_mean_derivative(xout_sigma_twp75, DELTA_t);
deriv_SIGMA1_xp = deriv_SIGMA1_deltat.*deriv_deltat_xp;
deriv_SIGMA1_yp = deriv_SIGMA1_deltat.*deriv_deltat_yp;
deriv_SIGMA1_zp = deriv_SIGMA1_deltat.*deriv_deltat_zp;

deriv_SIGMA1_deltat
deriv_MU1_deltat


Phi_MU1_SIGMA_1_ZERO =  normcdf(zeros(n),MU1,SIGMA1);
Phi_MU1_SIGMA_1_ZERO = U.*Phi_MU1_SIGMA_1_ZERO;


% Phi_MUN_SIGMA_N_ZERO =  normcdf(zeros(n),MUN,SIGMAN);
% Phi_MUN_SIGMA_N_ZERO = U.*Phi_MUN_SIGMA_N_ZERO;


s = [-10^7, -10^6, -10^5,-10^4,-10^3,-10^2,-10:0.05:-1-0.05,-1:0.01:0];

% Here we redefine Phi_MU1_SIGMA_1_ZERO so it is done the same way as the
% fortran code
for i = 1:n
    for j= 1:n
        Tester(i,j) = (1/(SIGMA1(i,j)*sqrt(2*pi)))*trapz(s,exp(-(s-MU1(i,j)).^2/(2*SIGMA1(i,j)^2)));
    end
end
Phi_MU1_SIGMA_1_ZERO =Tester;

deriv_Phi_xp = zeros(n,n);
deriv_Phi_yp = zeros(n,n);
deriv_Phi_zp = zeros(n,n);

delta_CWI = [0:0.005:0.51-0.01,0.51:0.01:1.2];

deriv_lnP_xp = zeros(n,n);
deriv_lnP_yp = zeros(n,n);
deriv_lnP_zp = zeros(n,n);

P_CWIN_deltat = zeros(n,n);
deriv_A_xp = zeros(n,n);
deriv_A_yp = zeros(n,n);
deriv_A_zp = zeros(n,n);
deriv_P_xp = zeros(n,n);
deriv_P_yp = zeros(n,n);
deriv_P_zp = zeros(n,n);
A= zeros(n,n);

for i = 1:n-1
    for j = i+1:n
        if MUN(i,j) ~=-99999  % i.e. do not consider pairs which do not have CWI constraint
            g = - exp( log((s-MU1(i,j)).^2) - log(2*SIGMA1(i,j).^2) );
            %        g = - exp(log( (s-MU1(i,j))**2) - log(2*SIGMA1(i,j)**2))
            deriv_g_xp_tmp = ( 4*SIGMA1(i,j).^2 *(s-MU1(i,j)).*deriv_MU1_xp(i,j) + ...
                4*SIGMA1(i,j)*deriv_SIGMA1_xp(i,j)*(s-MU1(i,j)).^2) ./ ...
                ( 4*SIGMA1(i,j).^4 );
            deriv_g_yp_tmp = ( 4*SIGMA1(i,j).^2 *(s-MU1(i,j))*deriv_MU1_yp(i,j) + ...
                4*SIGMA1(i,j)*deriv_SIGMA1_yp(i,j)*(s-MU1(i,j)).^2) ./...
                (4*SIGMA1(i,j).^4);
            deriv_g_zp_tmp = ( 4*SIGMA1(i,j).^2 *(s-MU1(i,j))*deriv_MU1_zp(i,j) + ...
                4*SIGMA1(i,j)*deriv_SIGMA1_zp(i,j)*(s-MU1(i,j)).^2) ./...
                (4*SIGMA1(i,j).^4);

            %deriv_Phi_xp(i,j) = 1/(SIGMA1(i,j)^2*sqrt(2*pi)).*trapz(s,exp(g).*deriv_g_xp_tmp);
            %deriv_Phi_yp(i,j) = 1/(SIGMA1(i,j)^2*sqrt(2*pi)).*trapz(s,exp(g).*deriv_g_yp_tmp);
            %deriv_Phi_zp(i,j) = 1/(SIGMA1(i,j)^2*sqrt(2*pi)).*trapz(s,exp(g).*deriv_g_zp_tmp);

            deriv_Phi_xp(i,j) = 1/sqrt(2*pi).*(...
                SIGMA1(i,j)*trapz(s,exp(g).*deriv_g_xp_tmp) - ...
                deriv_SIGMA1_xp(i,j)*trapz(s,exp(g)))/SIGMA1(i,j)^2;
            deriv_Phi_yp(i,j) = 1/sqrt(2*pi).*(...
                SIGMA1(i,j)*trapz(s,exp(g).*deriv_g_yp_tmp) - ...
                deriv_SIGMA1_yp(i,j)*trapz(s,exp(g)))/SIGMA1(i,j)^2;
            deriv_Phi_zp(i,j) = 1/sqrt(2*pi).*(...
                SIGMA1(i,j)*trapz(s,exp(g).*deriv_g_zp_tmp) - ...
                deriv_SIGMA1_zp(i,j)*trapz(s,exp(g)))/SIGMA1(i,j)^2;



            %         deriv_Phi_xp(i,j) = 1/sqrt(2*pi)* ( &
            %                       (SIGMA1(i,j)*trap_int_unevenspacing(integrand_tmp1,s) - &
            %                       deriv_SIGMA1_xp(i,j)*trap_int_unevenspacing(exp_g,s)) / &
            %                       (SIGMA1(i,j)**2))

            % The following can be computed in a vectorized fashion
            % However they needed below so keeping these outside the loop will
            % mean that another loop must be created. Doing these here is
            % viewed as more efficient.
            deriv_A_xp(i,j) = -(-deriv_Phi_xp(i,j).*SIGMA1(i,j)*sqrt(2*pi) + ...
                (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*deriv_SIGMA1_xp(i,j)*sqrt(2*pi) ) ./   ...
                ( (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*SIGMA1(i,j)*sqrt(2*pi)).^2;
            deriv_A_yp(i,j) = -(-deriv_Phi_yp(i,j).*SIGMA1(i,j)*sqrt(2*pi) + ...
                (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*deriv_SIGMA1_yp(i,j)*sqrt(2*pi) ) ./   ...
                ( (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*SIGMA1(i,j)*sqrt(2*pi)).^2;
            deriv_A_zp(i,j) = -(-deriv_Phi_zp(i,j).*SIGMA1(i,j)*sqrt(2*pi) + ...
                (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*deriv_SIGMA1_zp(i,j)*sqrt(2*pi) ) ./   ...
                ( (1-Phi_MU1_SIGMA_1_ZERO(i,j)).*SIGMA1(i,j)*sqrt(2*pi)).^2;




            h = -(delta_CWI - MU1(i,j)).^2./(2*SIGMA1(i,j).^2);

            deriv_h_xp = (  4.*SIGMA1(i,j).^2.*(delta_CWI - MU1(i,j)).*(deriv_MU1_xp(i,j)) + ...
                4.*(delta_CWI - MU1(i,j)).^2.*SIGMA1(i,j).*deriv_SIGMA1_xp(i,j))./ ...
                (4*SIGMA1(i,j).^4);
            deriv_h_yp = (  4.*SIGMA1(i,j).^2.*(delta_CWI - MU1(i,j)).*(deriv_MU1_yp(i,j)) + ...
                4.*(delta_CWI - MU1(i,j)).^2.*SIGMA1(i,j).*deriv_SIGMA1_yp(i,j))./ ...
                (4*SIGMA1(i,j).^4);

            deriv_h_zp = (  4.*SIGMA1(i,j).^2.*(delta_CWI - MU1(i,j)).*(deriv_MU1_zp(i,j)) + ...
                4.*(delta_CWI - MU1(i,j)).^2.*SIGMA1(i,j).*deriv_SIGMA1_zp(i,j))./ ...
                (4*SIGMA1(i,j).^4);


            deriv_B_xp = exp(h).*deriv_h_xp;
            deriv_B_yp = exp(h).*deriv_h_yp;
            deriv_B_zp = exp(h).*deriv_h_zp;

            % Compute A
            A(i,j) = 1./( (1-Phi_MU1_SIGMA_1_ZERO(i,j)) .*SIGMA1(i,j)*sqrt(2*pi) );

            % Compute B
            B = exp( -(delta_CWI- MU1(i,j)).^2 ./ (2*SIGMA1(i,j).^2));

            % Compute C
            [rubish,Phi_mu2_sigma2_0] = ...
                folded_normal_pdf(0,MUN(i,j),SIGMAN(i,j));
            %Phi_mu2_sigma2_0 = Phi_mu2_sigma2_0_tmp(1);
            C = 1./((1-Phi_mu2_sigma2_0).*SIGMAN(i,j)*sqrt(2*pi));

            % Compute D
            D = exp( -(MUN(i,j)-delta_CWI).^2 ./ (2.*SIGMAN(i,j).^2));

            integrand = B.*C.*D;
            P_CWIN_deltat(i,j) = A(i,j).* trapz(delta_CWI,B.*C.*D);

            deriv_P_xp(i,j) =    deriv_A_xp(i,j).*C*trapz(delta_CWI,B.*D) + ...
                A(i,j).*C*trapz(delta_CWI,deriv_B_xp.*D);
            deriv_P_yp(i,j) =    deriv_A_yp(i,j).*C*trapz(delta_CWI,B.*D) + ...
                A(i,j).*C*trapz(delta_CWI,deriv_B_yp.*D);
            deriv_P_zp(i,j) =    deriv_A_zp(i,j).*C*trapz(delta_CWI,B.*D) + ...
                A(i,j).*C*trapz(delta_CWI,deriv_B_zp.*D);

            deriv_lnP_xp(i,j) = 1/P_CWIN_deltat(i,j).*deriv_P_xp(i,j);
            deriv_lnP_yp(i,j) = 1/P_CWIN_deltat(i,j).*deriv_P_yp(i,j);
            deriv_lnP_zp(i,j) = 1/P_CWIN_deltat(i,j).*deriv_P_zp(i,j);
            
        else
            P_CWIN_deltat(i,j) = 1
            deriv_P_xp = 0
            deriv_P_yp = 0
            deriv_P_zp = 0

        end

    end

end


% Compute Lstar 

%f = -sum(sum(P_CWIN_deltat));

% for i = 1:n
%     for j = i+1:n
%         disp(['i = ', num2str(i), ' j = ', num2str(j)])
%         disp(['P_CWIN_deltat(i,j) = ', num2str(P_CWIN_deltat(i,j))])
%         disp(['log(P_CWIN_deltat(i,j)) = ', num2str(log(P_CWIN_deltat(i,j)))])
%     end
% end

f = 0;
for i = 1:n
    for j = i+1:n
      %disp(['i =', num2str(i), ' j = ', num2str(j)])  
      %if P_CWIN_deltat(i,j) ~=0
        f = f-log(P_CWIN_deltat(i,j));
      %end
    end
end

grad_xp = -sum(deriv_lnP_xp - deriv_lnP_xp',2);
grad_yp = -sum(deriv_lnP_yp - deriv_lnP_yp',2);
grad_zp = -sum(deriv_lnP_zp - deriv_lnP_zp',2);


grad = zeros(size(x));
grad(1)  =  grad_xp(2);
grad(2) = grad_xp(3);
grad(3) = grad_yp(3);
count = 4;
for i = 4:n
    grad(count) = grad_xp(i);
    grad(count+1) = grad_yp(i);
    grad(count+2) = grad_zp(i);   
    count = count+3;
end
            

            

            
            
            
            
            
            

      
      
      
            
            
            
                       
