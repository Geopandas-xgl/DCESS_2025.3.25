function [LL,HL,AT,LB,LLCO3,HLCO3,LLCO3s,HLCO3s,pCalLL,pCalHL,pOrgLL,pOrgHL,dwpCalLL,dwpCalHL,...
    fimLL,fimHL,dwpCalssLL,dwpCalssHL,dwpOrgLL,dwpOrgssLL,dwpOrgHL,dwpOrgssHL,wsedLL,wsedHL]=InitCnd_M


ParVal_M; % Activate and define Global parameters once

global GAw GAc

GAw = load('GAw_Eo.txt');
GAc = load('GAc_Eo.txt');

GAw=GAw/GAw(1);
GAc=GAc/GAc(1);

load ResThilda_M.mat

% LL(8,:) = LL(8,:) -0.04;


return
