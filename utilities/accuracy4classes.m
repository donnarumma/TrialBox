function   class_accuracy = accuracy4classes(y_true,y_pred)
% function class_accuracy = accuracy4classes(y_true,y_pred)
% Accuracy for class evaluation

% Number of classes
classes                     = unique(y_true);
num_classes                 = length(classes);
% Initialize a vector to store class-wise accuracy for each class
class_accuracy              = nan(1, num_classes);
% Calculate class-wise accuracy for each class
for ic = 1:num_classes
    class                   = classes(ic);
    samples_in_class        = y_true == class;
    correct_predictions     = sum(y_pred(samples_in_class) == y_true(samples_in_class));
    total_samples_in_class  = sum(samples_in_class);
    class_accuracy(ic)      = 100*correct_predictions / total_samples_in_class;
end
