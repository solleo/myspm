function EXP=myspm_compare_denosing(EXP)
% EXP=myspm_compare_denosing(EXP)
%
% (cc) 2015, sgKIM.

desc={'orig+trnd+cc+gs+scrb','orig+trnd+rp+cc+gs'};
cov_names1={'orig+trnd','orig+trnd','orig+trnd+cc','orig+trnd+cc+gs',};
cov_names2={'orig+trnd','orig+trnd+rp','orig+trnd+rp+cc','orig+trnd+rp+cc+gs'};
J1=[1 1 2 3];
J2=[1 2 3 4];
figure('position',[1963 55 894 1089]);
subjID = fsss_subjID(EXP.subjID);
path0=pwd;
for i=1:numel(EXP.subjID)
  clf
  subjid = subjID{i};
  cd (['/scr/vatikan3/APConn/rest12.410/',subjid]);
  load('resy_wmcsf99_n16d1v1b0.00-Inf_z3m0.6_b0.01-0.1orig+trnd+cc+gs+scrb.mat','sumstat');
  sumstat1=sumstat;
  load('resy_wmcsf99_n16d1v1b0.00-Inf_z3m0.6_b0.01-0.1orig+trnd+rp+cc+gs.mat','sumstat');
  sumstat2=sumstat;
  load art_regression_outliers_and_movement_uarest410.mat;
  subplot(5,1,1)
  plot(R(:,end),'color',[.2 .2 1]);
  title([subjid,':bpf[0.01,0.10]Hz'])
  xlabel('TR'); ylabel('||dm/dt||_2')
  xlim([1 410]); grid on; box on;
  
  for j=1:4
    subplot(5,1,j+1); hold on;
    if isfield(EXP,'absZ')&&EXP.absZ
      plot(sumstat2.meanabsz(J2(j),:)','color',[1 .2 .2]);
      plot(sumstat1.meanabsz(J1(j),:)','color',[0 .8 0]);
      fname_png=['/scr/vatikan3/APConn/rest12.410/fig_denoising/meanabsZ_',subjid,'.png'];
      ylabel('mean |Z|')
    else
      plot(sumstat2.meanz(J2(j),:)','color',[1 .2 .2]);
      plot(sumstat1.meanz(J1(j),:)','color',[0 .8 0]);
      fname_png=['/scr/vatikan3/APConn/rest12.410/fig_denoising/meanZ_',subjid,'.png'];
      ylabel('mean Z')
    end
    xlabel('TR'); 
    xlim([1 410]); grid on; box on;
    legend({cov_names1{j}, cov_names2{j}},'location','NorthWest');
  end
  
  screen2png(fname_png);
end
cd(path0)
end