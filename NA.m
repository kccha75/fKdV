% Function for linear operator Au=f for equation -u_xx+au_x+bu using
% Fourier spectral methods

function A=NA(u,k,a,b,c)

A=-ifft(-k'.^2.*fft(u))+ifft(a.*1i.*k'.*fft(u))+b.*u+c.*u.^2;

end