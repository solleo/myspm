function [stdout,stderr] = myunix(cmd,echo)
% [stdout,stderr] = myunix(cmd, echo)
%
% (cc) 2017, sgKIM

% if strcmp(cmd(1:7),'tkmedit') || strcmp(cmd(1:8),'tksurfer')
%  setenv('LD_LIBRARY_PATH','glnx86:/afs/cbs.mpg.de/software/matlab/9.1/sys/os/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/bin/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/extern/lib/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/runtime/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/matlabfix/ubuntu-xenial-amd64')
% else
%  setenv('LD_LIBRARY_PATH','/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/matlab/9.1/extern/lib/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/runtime/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/matlabfix/ubuntu-xenial-amd64');
% end

% setenv('LD_LIBRARY_PATH','/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/matlab/9.1/extern/lib/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/runtime/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/matlabfix/ubuntu-xenial-amd64');
if ~exist('echo','var'), echo=1; end
if echo
 disp(['$ ',cmd]);
end
% if strcmp(cmd(1:2),'fl')
%  cmd=['/afs/cbs.mpg.de/software/scripts/FSL --version 5.0 ',cmd(3:end)];
%  setenv('FSLOUTPUTTYPE','NIFTI')
% elseif strcmp(cmd(1:2),'fs')
%  cmd=['/afs/cbs.mpg.de/software/scripts/FREESURFER --version 6.0.0 ',cmd(3:end)];
% end
[stdout,stderr]=unix(cmd);
if ~isempty(stderr), disp(stderr); end
if nargout == 0
 clear stdout stderr
end
end
