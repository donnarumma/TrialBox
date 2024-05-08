function [valmin,t_point] = dpca_pCompute(icomp,data_trials,par)
%function [valmin,t_point] = dpca_pCompute(icomp,data_trials,par)

try
    it_start = par.it_start;
    it_stop = par.it_stop;
catch
    par.it_start        = data_trials(1).timedpca(1) + 0.5;
    par.it_stop         = data_trials(1).timedpca(end) - 0.5;
end
classes = unique([data_trials.trialType]);
time = data_trials(1).timedpca;
ind_time = time>it_start & time<it_stop;
new_time = time(ind_time); 
dpca_ic = nan(length(data_trials([data_trials.trialType]==1)),length(classes));
p = nan(length(new_time),1);
for it=1:length(new_time)
    for ic=1:length(classes)
        data_class=data_trials([data_trials.trialType]==classes(ic));
        for inc=1:length(data_class)
            data_time = data_class(inc).dpca(:,ind_time);
            dpca_ic(inc,ic)=data_time(icomp,it);
        end
    end
    p(it)=anova1(dpca_ic,[],'off');
end

% %% Plot log(p)
% figure
% plot(time,log(p))

[valmin,indmin]= min(log(p));
t_point = new_time(indmin);