function JOB = myspm_conv(JOB)
% convolute a given timeseries JOB.u and returns JOB.uc
%
% JOB requires:
%  .u
%  .TR_sec
% (.hrfname') can be 'hrf' (default), or
%                 'hrf (with time derivative)'
%                 'hrf (with time and dispersion derivatives)'
%                 'Fourier set'
%                 'Fourier set (Hanning)'
%                 'Gamma functions'
%                 'Finite Impulse Response'
% 
% [TODO:] non-default options for HRFs
% xBF.length  - window length (seconds)
% xBF.order   - order
% xBF.T       - microtime resolution (for 'hrf*')
% xBF.bf      - array of basis functions
%
% ref: http://spm.martinpyka.de/?p=41
% (cc) 2017, sgKIM.
if ~isfield(JOB,'hrfname'), JOB.hrfname='hrf'; end

xBF = struct('dt',JOB.TR_sec, 'name',JOB.hrfname);
xBF = spm_get_bf(xBF);
U.u = reshape(JOB.u, [numel(JOB.u) 1]);
U.name = {''};
JOB.uc = spm_Volterra(U, xBF.bf);
JOB.xBF = xBF;
end

%% test:
% % ref: http://spm.martinpyka.de/?p=41
% % % with a tr of 2.2 and HRF as basis function you receive
% % dt: 2.2000
% % name: 'hrf'
% % length: 33
% % order: 1
% TR_sec=1;
% bf=struct('dt',TR_sec, 'name','hrf');
% bf = spm_get_bf(bf);
% clf;
% subplot(311)
% plot(bf.bf); xlim([0 80])
% % bf: [33x1 double]
% % these are necessary definitions for the spm_Volterra-function
% U.u = [zeros(1,10), ones(1,10),zeros(1,10), ones(1,10) zeros(1,10), ones(1,10),zeros(1,10), ones(1,10)]';
% subplot(312)
% plot(U.u);
% U.name = {'reg'}; % U needs a name, but the string has no meaning
% % % convolve reg with the hrf
% convreg = spm_Volterra(U, bf.bf);
% subplot(313)
% plot(convreg)
