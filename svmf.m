
data = load('C:\Users\Umair Khan\Desktop\watershed (2)\watershed\tfeat.txt');
%loading the feature vectors from the file in which they were saved 
%in this case they were saved in the file tfeat.txt in wtrshed.m

feat= data(:,1:4);%separating the features from the data
tgt = data(:,5);%separating the target file - labels 

svmodel = fitcsvm(feat,tgt); %training the SVM classifier using cross validation of 10 folds by default....

cvsvm = crossval(svmodel);

predlabels = kfoldPredict(cvsvm); % testing our SVM model. 

cm = confusionmat(tgt, res) % confusion matrix created for known labels (tgt) with predicted labels in predlabels

%since the feature vectors are separated very distinctively so accuracy
%comes out to be 100%

