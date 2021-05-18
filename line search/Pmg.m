% Prolongation operator for 1D Fourier Multigrid

% Inputs:
% vc - solution on coarse grid
% Nf - N in fine grid

% Outputs:
% vf - correction to previous guess

function vf=Pmg(vc,Nf)

% vf=ifft(fft(vc),Nf);
vc_hat=fft(vc);
% Pad high frequency with 0's
vf=ifft([vc_hat(1:Nf/4);zeros(Nf/2+1,1);vc_hat(Nf/4+2:end)]);

end
