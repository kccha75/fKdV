% Function for linear operator Au=f for equation -u_xx+au_x+bu using
% Fourier spectral methods

function A=LA(u,k,a,b)

A=-ifft(-k'.^2.*fft(u))+ifft(a.*1i.*k'.*fft(u))+b.*u;

end