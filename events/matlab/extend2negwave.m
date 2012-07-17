function fG = extend2negwave(fG,dim)
% Now we can extend the filter into the negative wavenumbers using:
% 1) the spatial representation of the filter is real and therefore H(-f) = conj(H(f))
% 2) the fft stores information for wavenumbers as follows
%       f(0), f(1), f(2), ..., f(c-1), f(c) or -f(c), -f(c-1), ..., -f(2), -f(1) 
% Note here that f refers to 1/wavenumber (it would be frequency in the time domain)
nx = dim(1);
ny = dim(2);
nz = dim(3);


fG(nz/2+2:nz,:,:) = conj(flipdim(fG(2:nz/2,:,:),1));
fG(:,nx/2+2:nx,:) = conj(flipdim(fG(:,2:nx/2,:),2));
fG(:,:,floor(ny/2)+2:ny) = conj(flipdim(fG(:,:,2:floor(ny/2)),3));