%-------------------------------------------------------------------------------------%
%  NoHiMDO                                                                            %
%                                                                                     %
%  A solver for Multi-Disciplinary Optimization, based on Non-Hierarchical Analytical %
%  Target Cascading                                                                   %
%  Version 2.0.1                                                                      %
%                                                                                     %
%  Copyright (C) 2012-2016  Bastien Talgorn - McGill University, Montreal             %
%                                                                                     %
%  Author: Bastien Talgorn                                                            %
%  email: bastientalgorn@fastmail.com                                                 %
%                                                                                     %
%  This program is free software: you can redistribute it and/or modify it under the  %
%  terms of the GNU Lesser General Public License as published by the Free Software   %
%  Foundation, either version 3 of the License, or (at your option) any later         %
%  version.                                                                           %
%                                                                                     %
%  This program is distributed in the hope that it will be useful, but WITHOUT ANY    %
%  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    %
%  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   %
%                                                                                     %
%  You should have received a copy of the GNU Lesser General Public License along     %
%  with this program. If not, see <http://www.gnu.org/licenses/>.                     %
%                                                                                     %
%  You can find information on NoHiMDO at https://github.com/bastientalgorn/NoHiMDO   %
%-------------------------------------------------------------------------------------%
close all
clear all

% ================== %
%    PARAM LISTS     %
% ================== %
w0 = 1;
tab_gamma = [0.0 0.2 0.4 0.6 0.8 1.0];
tab_beta = [1.4 1.8 2.2 2.6];

tab_NI = [   8  16  32   64  128  256  512  ];
tab_NO = [ 512 256 128   64   32   16    8  ];

tab_pb = {'geo','sac','biquad','sbj'};

% ================== %
%   BUILD RUN LIST   %
% ================== %
list_run = cell(0);
  
  pb = 'sbj';
  gamma = 0.4;
  beta = 2.2;
  NI = 640;
  NO = 640;
  list_run{end+1} = {pb,NI,NO,gamma,beta,w0,'mads'};


if false
  %Gamma & beta comparison
  NI = 64;
  NO = 64;
  for i_gamma = 1:length(tab_gamma)
      gamma = tab_gamma(i_gamma);
      for i_beta = 1:length(tab_beta)
          beta  = tab_beta (i_beta );
          for i_pb = 1:length(tab_pb)
              pb = tab_pb{i_pb};
              list_run{end+1} = {pb,NI,NO,gamma,beta,w0,'mads'};
          end
      end
  end

  %NI & NO comparison
  gamma = 0.4;
  beta = 2.2;
  for i = 1:length(tab_NI)
      NI = tab_NI(i);
      NO = tab_NO(i);
      for i_pb = 1:length(tab_pb)
          pb = tab_pb{i_pb};
          list_run{end+1} = {pb,NI,NO,gamma,beta,w0,'mads'};
          list_run{end+1} = {pb,NI,NO,gamma,beta,w0,'madsx'};
      end
  end
end


disp('===========================================================');
disp('===========================================================');
disp('===========================================================');


% ============= %
%     TEST      %
% ============= %
for i=1:length(list_run)
    pb   = list_run{i}{1};
    NI   = list_run{i}{2};
    NO   = list_run{i}{3};
    gamma= list_run{i}{4};
    beta = list_run{i}{5};
    w0   = list_run{i}{6};
    optim_method = list_run{i}{7};

    out_file = [pb '_' num2str(NO) '_' num2str(NI) '_' num2str(gamma) '_' num2str(beta) '_' optim_method '.mat'];
    disp(out_file)
    if exist([out_file '.mat'],'file')
        disp('EXISTS ! (skip)');
    else
        disp(['MISSING : ' out_file]);
    end
end

list_run = list_run(randperm(length(list_run)));


disp('===========================================================');
disp('===========================================================');
disp('===========================================================');

disp(['Nb of runs : ' num2str(length(list_run))]);

disp('===========================================================');
disp('===========================================================');
disp('===========================================================');


% ============= %
%     RUN       %
% ============= %
for i=1:length(list_run)
    pb   = list_run{i}{1};
    NI   = list_run{i}{2};
    NO   = list_run{i}{3};
    gamma= list_run{i}{4};
    beta = list_run{i}{5};
    w0   = list_run{i}{6};
    optim_method = list_run{i}{7};
    out_file = [pb '_' num2str(NO) '_' num2str(NI) '_' num2str(gamma) '_' num2str(beta) '_' optim_method '.mat'];
    disp('===========================================================');
    disp(out_file);
    for nexp=1:40
      try
        run_nexp(pb,NI,NO,gamma,beta,w0,optim_method,nexp);
      catch exception
          display_exception(exception);
      end
    end
end

disp('===========================================================');
disp('===========================================================');
disp('===========================================================');



disp('=====================EXIT==================================');
exit;


