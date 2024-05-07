function Tab = EventExtract(tabinfos,event)
% function Tab = TimeStepCompute(tabinfos)
Tabinfos_New = struct();
for n_trial = 1:length(tabinfos)
    Tabinfos_New(n_trial).Start       = tabinfos(n_trial,2);
    Tabinfos_New(n_trial).TargetON    = tabinfos(n_trial,3);
    Tabinfos_New(n_trial).TargetOFF   = tabinfos(n_trial,3) + tabinfos(n_trial,4);
    Tabinfos_New(n_trial).GO          = tabinfos(n_trial,3) + tabinfos(n_trial,5);
    Tabinfos_New(n_trial).Touch       = tabinfos(n_trial,3) + tabinfos(n_trial,6);
    Tabinfos_New(n_trial).RT          = tabinfos(n_trial,3) + tabinfos(n_trial,5) + tabinfos(n_trial,9);
end
Tab = struct();
for n_trial = 1:length(tabinfos)
    Tab(n_trial).Start       = round((1.0000e-03*(Tabinfos_New(n_trial).Start - Tabinfos_New(n_trial).(event))),3);
    Tab(n_trial).TargetON    = round((1.0000e-03*(Tabinfos_New(n_trial).TargetON - Tabinfos_New(n_trial).(event))),3);
    Tab(n_trial).TargetOFF   = round((1.0000e-03*(Tabinfos_New(n_trial).TargetOFF - Tabinfos_New(n_trial).(event))),3);
    Tab(n_trial).GO          = round((1.0000e-03*(Tabinfos_New(n_trial).GO - Tabinfos_New(n_trial).(event))),3);
    Tab(n_trial).Touch       = round((1.0000e-03*(Tabinfos_New(n_trial).Touch - Tabinfos_New(n_trial).(event))),3);
    Tab(n_trial).RT          = round((1.0000e-03*(Tabinfos_New(n_trial).RT - Tabinfos_New(n_trial).(event))),3);
end
