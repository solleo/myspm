function [w,gi,lag,ri] = myspm_est_res(SPM)
% code robbed from spm_est_V.m:

dt = SPM.xY.RT;
X  = SPM.xX.X;
X  = [X SPM.xX.K.X0];                         % add low-freq cosine bases
m  = size(X,1);
R  = speye(m,m) - X*spm_pinv(X);              % residual forming matrix

C  = SPM.xVi.Cy;                              % sample covariance
C  = C*trace(R*R')/trace(R*C*R');             % scaling

% % covariance components (a mixture of exponentials)
% %==========================================================================
% dt    = SPM.xY.RT;                            % TR (seconds)
% T     = (0:(m - 1))*dt;                       % time
% d     = 2.^(floor(log2(dt/4)):log2(64));      % time constants (seconds)
% QQ    = {};                                   % dictionary of components
% for i = 1:length(d)
%     for j = 0:1
%         QQ{end + 1} = toeplitz((T.^j).*exp(-T/d(i)));
%     end
% end
% Q     = QQ(1:end);                            % full (exactly)
% [V,h] = spm_reml(C,X,Q,1,1,0,4);              % ReML to est V
V  = SPM.xVi.V;                               % intrinsic cov after HPF
W  = spm_inv(spm_sqrtm(V));                   % W*W' = V^-1
W  = W*sqrt(trace(R*R')/trace(R*W*C*W'*R'));  % scaling

% Power spectra:
S  = spm_sqrtm(C);
g  = [  sum(abs(fft(full(R)).^2),2)];         % residual forming matrix
g  = [g sum(abs(fft(full(R*S)).^2),2)];       % residuals unwhitened
g  = [g sum(abs(fft(full(R*W*S)).^2),2)] ;    % residuals whitened
i  = fix(2:m/2);
w  = (1:length(i))/(2*length(i)*dt);
gi = g(i,:);

% Autocovariance:
r     = [];
for i = 1:size(g,2)
  f      = ifft(g(:,i));
  r(:,i) = real(fftshift(f));
end
lag   = (-32:32);
i     = lag + fix(m/2);
ri    = r(i,:);
lag   = lag*dt;

if ~nargout
  figure
  subplot(121)
  plot(w,gi)
  legend({'ideal','unwhitened','whitened'})
  xlabel('Freq [Hz]'); ylabel('Power')
  title('Spectral density')
  subplot(122)
  plot(lag,ri)
  xlabel('Lag [s]'); ylabel('Cov');
  title('Autocovariance func')
end

end