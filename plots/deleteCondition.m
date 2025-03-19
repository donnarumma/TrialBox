function   deleteCondition(out,ind)
% function deleteCondition(out,ind)
xl=get(gca,'xlim');
yl=get(gca,'ylim');
zl=get(gca,'zlim');
for i=1:length(out.h{ind})
    delete(out.h  {ind}(i));
    delete(out.hpe{ind}(i));
    delete(out.hp0{ind}(i));
    delete(out.hpc{ind}(i));
end

xlim(xl);
ylim(yl);
zlim(zl);