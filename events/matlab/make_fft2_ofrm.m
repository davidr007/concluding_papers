function a = make_fft2_ofrm(dim)
% the second value in dim must be 1 for the 2d case
% the top row is all real
% the nz/2+1 is all positive
% reverse conjugate symmetry (eg flipud) exists across the nz/2+1 row

if length(dim)~=3
    error('dim must have length 3')
end

nx = dim(1);
ny = dim(2);  % must be 1 for the 2d case
nz = dim(3);


a = exp(i*pi*(rand([nz,nx,ny])-pi));

% for k=1:ny
%     a(2:nz/2,1:nx,k) = exp(i*2*pi*rand(nz/2-1,nx));
%     a(nz/2+2:end,1:nx,k) = conj(flipud(a(2:nz/2,1:nx,k)));
% end


% if dim(2) == 0 % then it is the 2d case
%     a = zeros([nz,nx])
%     a(2:nz/2,1:nx) = exp(i*2*pi*rand(nz/2-1,nx));
%     a(nz/2+2:end,1:nx) = conj(flipud(a(2:nz/2,1:nx)));
% else dim(2)  % the it is the 3d case
%     a = zeros([nz,nx,ny])
%     for k=1:ny
%         a(2:nz/2,1:nx,k) = exp(i*2*pi*rand(nz/2-1,nx));
%         a(nz/2+2:end,1:nx,k) = conj(flipud(a(2:nz/2,1:nx,k)));
%     end
%             
% end
    

