function [med1]= gen_random_medium(dim,deltas,std_gwn, mean_medium,a, diag_switch,mtype,rseed)

% INPUTS:
% dim           [vector 3x1]: [dimx, dimy dimz]. Note that if dimy = 1 a 2d random
%               medium will be created.
% deltas        [vector 3x1]: grid spacing [deltax, deltay deltaz]
% st_gwn        [scalar]: standard deviation of the gaussian taken of all points
%               in the grid after the filter is applied
% mean_medium   [scalar]: mean value for medium
% a             [scalar]: correlation length (m)
% diag_switch   [scalar - switch] 1=> draw diagnostics
% mtype         [str] type of media
%                   'grm'  Gaussian Random Medium
%                   'exp'   Exponential Medium
% rseed         [scalar]
%                   1 => reset random number seeds
%                   0 => do not reset seeds
% USEAGE:
% 3D GRM: [med9]= gen_random_medium([192,192,192],[40,40,40],200, 5000,150,0,'grm');
% 2D vk:  [med] = gen_random_medium([192,1,192],[40,40,40],200, 5000,150,1,'vk');


% cswitch   case switch defining the type of rm to create
%           'original'
%           'spectral'
%
% saveswitch  case switch to define whether or not to save the random media
%           1 => save data
%           0 => do not save the data

if length(dim)~=3
    error('dim must have length 3')
end
nx = dim(1);
ny = dim(2);
nz = dim(3);
dx = deltas(1);
dy = deltas(2);
dz = deltas(3);

if ny ==1 % force dy =1
    dy = 1
end

if ~exist('rseed') | rseed ==1 | isempty('rseed') % reset seed (otherwise don't)
    randn('state',0)
    rand('state',0)
end

% compute the positive wave numbers.
kx = repmat([0:nx/2]./(nx*dx),[nz/2+1,1,floor(ny/2)+1]);
kz = repmat([0:nz/2]'./(nz*dz),[1,nx/2+1, floor(ny/2)+1]);
ky = zeros(nz/2+1,nx/2+1,floor(ny/2)+1);
for s = 1:floor(ny/2)+1
    tmp = [0:ny/2]./(ny*dy);
    ky(1,1,s) = tmp(s);
end
ky = repmat(ky(1,1,:),[nz/2+1,nx/2+1, 1]);
%ky = repmat([0:ny/2]'./(ny*dy),[nz/2+1,nx/2+1, 1]);
k = sqrt(kx.^2+kz.^2+ky.^2);

eps = calc_eps(std_gwn, mean_medium);
disp('   ')
disp(['Total Fractional Variance is: ', num2str(eps)])
disp([    'for mean: ',num2str(mean_medium), ' and std: ', num2str(std_gwn)])
disp('   ')

switch mtype
    case 'grm'

        if ny ==1 % the 2d case
            fG = (eps^2*a^2)/2*exp(-a^2.*k.^2/4);  % From Frankel eith eps term from Baig
        elseif ny >1
            %fG = sqrt(nx*dx*nz*dz*ny*dy)*a^3/(2*sqrt(2))*exp(-a^2.*k.^2/4);
            fG = (eps^2*a^3)/(2*sqrt(2))*exp(-a^2.*k.^2/4);  % From Baig
        else
            error('problem with value of ny=dim(2)')
        end
        %
    case 'exp'
        if ny ==1 %2D
            fG =(eps^2*a^2*2*pi)./(1+k.^2*a^2).^(3/2); % From Frankel or Hong +eps term from Baig
        elseif ny >1
            fG = (eps^2*a^3*2*sqrt(2))./(sqrt(pi)*(1+k.^2*a^2).^2);  % From Baig
        end
end


clear k kx ky kz
pack
fG = extend2negwave(fG,dim);

%fft_med1 = make_fft2_ofrm(dim);
med1 = zeros([nz,nx,ny],'single');
med1(:,:,:) = rand([nz,nx,ny])-0.5;
fft_med1 = fftn(med1);

fft_med1 = fG.*fft_med1;

if diag_switch ==1
    % Figure 1: plot med1 and fG
    figure

    if ny>1, subplot(4,2,1), else, subplot(2,2,1), end
    imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],med1(:,:,ceil(ny/2)))
    title('med1: x-z slice at y=ceil(ny/2)')
    xlabel('x')
    ylabel('z')
    if ny>1
        subplot(4,2,2)
        b = permute(med1,[3,2,1]);
        imagesc([dx/2:dx:nx*dx], [dy/2:dy:ny*dy],b(:,:,nz/2))
        title('med1: x-y slice at z=nz/2')
        xlabel('x')
        ylabel('y')
        subplot(4,2,3)
        b = permute(med1,[1,3,2]);
        imagesc([dy/2:dy:ny*dy] ,[dz/2:dz:nz*dz],b(:,:,nx/2))
        title('med1: y-z slice at x=nx/2')
        xlabel('y')
        ylabel('z')
    end
    if ny>1, subplot(4,2,4), else, subplot(2,2,2), end
    hist(med1(:))
    title('hist(med1(:))')


    if ny>1, subplot(4,2,5), else, subplot(2,2,3), end
    imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],fG(:,:,ceil(ny/2)))
    title('fG: x-z slice at y=ceil(ny/2)')
    xlabel('x')
    ylabel('z')
    if ny>1
        subplot(4,2,6)
        b = permute(fG,[3,2,1]);
        imagesc([dx/2:dx:nx*dx], [dy/2:dy:ny*dy],fG(:,:,nz/2))
        title('fG: x-y slice at z=nz/2')
        xlabel('x')
        ylabel('y')
        subplot(4,2,7)
        b = permute(fG,[1,3,2]);
        imagesc([dy/2:dy:ny*dy] ,[dz/2:dz:nz*dz],fG(:,:,nx/2))
        title('fG: y-z slice at x=nx/2')
        xlabel('y')
        ylabel('z')
    end
end
clear fG b
pack

med1 = ifftn(fft_med1);
clear fft_med1
pack

% now we re-scale med2 to ensure eq(12) of Baig03 is obeyed
eps_med1 = tot_fracvar(med1,mean_medium);
med1 = eps*med1/eps_med1;

%check1 = sum(med2(:).^2)/(nx*ny*nz)
check1 = mean(med1(:).^2);
check2 = eps^2*mean_medium^2;
if (check1-check2)/check1 > 0.05 % error greater than 5%
    error(['Problem with eps in gen_random_medium: rel error = ',num2str((check1-check2)/check1)] )
end

if diag_switch ==1
    % Figure 2: plots med2 and med3 (done later)
    figure

    if ny>1, subplot(4,2,1), else, subplot(2,2,1), end
    imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],med1(:,:,ceil(ny/2)))
    title('med2: x-z slice at y=ceil(ny/2)')
    xlabel('x')
    ylabel('z')
    if ny>1
        subplot(4,2,2)
        b = permute(med1,[3,2,1]);
        imagesc([dx/2:dx:nx*dx], [dy/2:dy:ny*dy],b(:,:,nz/2))
        title('med2: x-y slice at z=nz/2')
        xlabel('x')
        ylabel('y')
        subplot(4,2,3)
        b = permute(med1,[1,3,2]);
        imagesc([dy/2:dy:ny*dy] ,[dz/2:dz:nz*dz],b(:,:,nx/2))
        title('med2: y-z slice at x=nx/2')
        xlabel('y')
        ylabel('z')
    end
    if ny>1, subplot(4,2,4), else, subplot(2,2,2), end
    hist(med1(:))
    title('hist(med2(:))')
end


med1 = med1 + mean_medium;


if diag_switch ==1

    % Figure 2 (continued): plots med2 (done earlier) and med3

    if ny>1, subplot(4,2,5), else, subplot(2,2,3), end
    imagesc([dx/2:dx:nx*dx] ,[dz/2:dz:nz*dz],med1(:,:,ceil(ny/2)))
    title('med3: x-z slice at y=ceil(ny/2)')
    xlabel('x')
    ylabel('z')
    if ny>1
        subplot(4,2,6)
        b = permute(med1,[3,2,1]);
        imagesc([dx/2:dx:nx*dx], [dy/2:dy:ny*dy],b(:,:,nz/2))
        title('med3: x-y slice at z=nz/2')
        xlabel('x')
        ylabel('y')
        subplot(4,2,7)
        b = permute(med1,[1,3,2]);
        imagesc([dy/2:dy:ny*dy] ,[dz/2:dz:nz*dz],b(:,:,nx/2))
        title('med3: y-z slice at x=nx/2')
        xlabel('y')
        ylabel('z')
    end
    if ny>1, subplot(4,2,8), else, subplot(2,2,4), end
    hist(med1(:))
    title('hist(med3(:))')

end

%