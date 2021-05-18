% Restriction operator for 1D Fourier Multigrid

% Inputs:
% rf - residual for fine grid
% Nc - N in coarse grid

% Outputs:
% fc - RHS of coarse grid

function fc=Rmg(rf,Nc)

% fc=ifft(fft(rf),Nc);
r_hat=fft(rf);
% Remove high frequency components
fc=ifft([r_hat(1:Nc/2);0;r_hat(3*Nc/2+2:end)]);

end
