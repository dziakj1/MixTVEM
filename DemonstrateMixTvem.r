source("MixTvem.r");
obsData <- read.table("MixTvemSampleObservationLevel.txt",sep="\t",header=TRUE);
subData <- read.table("MixTvemSampleSubjectLevel.txt",sep="\t",header=TRUE);
obsData$Intercept <- 1;
obsData$CenteredNegAffect <- obsData$NegAffect - mean(obsData$NegAffect);
windows();
model1 <- TVEMMixNormal(dep=obsData$Urge,
                        id=obsData$ID,
                        numInteriorKnots=6,
                        numClasses=2,
                        numStarts=10,
                        tcov=obsData$Intercept,
                        time=obsData$Time);
windows();

set.seed(42); # random seed;
class1 <- which(model1$bestFit$postProbsBySub[,1]>
                              model1$bestFit$postProbsBySub[,2]);
   # Gets the indexes (as integers from 1 through #subjects, in data order)
   # of all the subjects estimated to be more likely to belong to class 1 than
   # to class 2.
class1exemplars <- sample(class1,4);
class2 <- which(model1$bestFit$postProbsBySub[,2]>
                              model1$bestFit$postProbsBySub[,1]);
   # Gets the indexes (as integers from 1 through #subjects, in data order)
   # of all the subjects estimated to be more likely to belong to class 1 than
   # to class 2.
class2exemplars <- sample(class2,4);
par(mfrow=c(2,4));
for (i in class1exemplars) {
plot(model1$time[which(model1$intId==i)],
     model1$dep[which(model1$intId==i)],type="l",
     main=paste("Class",1,": Subject",i),
     ylim=c(0,10),xlab="Time",ylab="Urge");
points(model1$time[which(model1$intId==i)],
     model1$dep[which(model1$intId==i)],pch=16);
}
for (i in class2exemplars) {
plot(model1$time[which(model1$intId==i)],
     model1$dep[which(model1$intId==i)],type="l",
     main=paste("Class",2,": Subject",i),
     ylim=c(0,10),xlab="Time",ylab="Urge");
points(model1$time[which(model1$intId==i)],
     model1$dep[which(model1$intId==i)],pch=16);
}



windows();
model1b <- TVEMMixNormal(dep=obsData$Urge,
                        id=obsData$ID,
                        numInteriorKnots=6,
                        numClasses=2,
                        numStarts=10,
                        tcov=obsData$Intercept,
                        xcov=obsData$CenteredNegAffect,
                        time=obsData$Time);
windows();
model2 <- TVEMMixNormal(dep=obsData$Urge,
                        id=obsData$ID,
                        numInteriorKnots=6,
                        numClasses=2,
                        numStarts=10,
                        tcov=cbind(obsData$Intercept,obsData$CenteredNegAffect),
                        time=obsData$Time);
windows();
obsData <- merge(obsData,subData,by="ID");
model3 <- TVEMMixNormal(dep=obsData$Urge,
                        id=obsData$ID,
                        numInteriorKnots=6,
                        numClasses=2,
                        numStarts=10,
                        tcov=cbind(obsData$Intercept,obsData$CenteredNegAffect),
                        scov=cbind(obsData$BaselineCigarettesPerDay,obsData$MinutesToFirstCigarette),
                        time=obsData$Time);
windows();

model4 <- TVEMMixNormal(dep=obsData$Urge,
                        id=obsData$ID,
                        numInteriorKnots=6,
                        numClasses=2,
                        numStarts=10,
                        tcov=cbind(obsData$Intercept,obsData$CenteredNegAffect),
                        scov=obsData$RelapseAtOneMonth,
                        time=obsData$Time);