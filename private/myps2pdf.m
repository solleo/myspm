function myps2pdf(fn_ps)
% myps2pdf(fn_ps)

gsbin='/usr/bin/ghostscript';
if ~exist(gsbin,'file')
  gsbin='/opt/local/bin/gs';
end
[p1,f1,e1] = myfileparts(fn_ps);
addpath(
ps2pdf('psfile',fn_ps, 'pdffile',[p1,'/',f1,'.pdf'],'gscommand',gsbin);
delete(fn_ps)
end