function myspm_fmriprep12_par(JOBs, num_workers)
% myspm_fmriprep12_par(JOBs, num_workers)
% (cc) 2019, sgKIM


delete(gcp('nocreate'))
if ~exist('num_workers','var')
  num_workers = 4;
end
parpool(num_workers)
parfor i = 1:numel(JOBs)
  JOB = JOBs{i};
  try
    JOB = myspm_fmriprep12(JOB)
  catch ME
    warning(ME.identifier, '##### ERROR DURING PROCESSING: %s\n',JOB.fname_epi);
    warning(ME.identifier, '%s', ME.message);
    for iStack=1:numel(ME.stack)
      fprintf('#####  In %s (line %i):\n', ...
        ME.stack(iStack).file, ME.stack(iStack).line);
    end
    JOB.ME = ME;
    JOB.error_with = JOB.fname_epi;
  end
  JOBs{i} = JOB;
end
delete(gcp('nocreate'))
end