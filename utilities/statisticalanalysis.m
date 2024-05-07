function res = statisticalanalysis(labels,predictions)

% Confusion matrix
Cmatrxix = confusionmat(labels, predictions);
% kappa value
res.kappaValue = kappaModel(Cmatrxix);

res.Accuracy                = sum(y_pred == labs)/length(labs)*100;
res.Accuracy_class(1,:)     = accuracy4classes(labs,y_pred);

