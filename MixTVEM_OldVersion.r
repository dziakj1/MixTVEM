TVEMMixNormal <- function( dep,     # The dependent variable as a vector, one
                                    # entry per assessment per subject
                     id,            # The subject ID as a vector, one
                                    # entry per assessment per subject
                                    # (could be integers or strings)
                     numInteriorKnots,   # The number of interior knots
                                         # for the splines representing
                                         # the time-varying coefficients.
                                         # Assumed to be the same
                                         # for each coefficient.
                     numClasses,    # Number of classes in the mixture model
                     tcov,          # Time-varying covariates as a matrix, one
                                    # row per assessment per person.  Usually,
                                    # the first column should be a column of
                                    # ones representing an intercept.
                     time,          # Assessment time as a vector
   # These arguments have default values, given after the equal sign for each:
                     convergenceCriterion=1e-8,  # Convergence criterion for the
                                    # maximum absolute deviation of parameter
                                    # estimates between successive
                                    # steps of the EM algorithm
                     deg=3,         # Degree of the polynomial between
                                    # successive knots for the time-
                                    # varying coefficient functions.
                     doPlot=TRUE,   # Whether to draw a plot of the time-
                                    # varying coefficient functions.
                     equalVariance=FALSE, # Whether to assume equal error
                                    # variance between classes
                     gammaPriorConstant=1,  # Optional weight of a Bayesian
                                    # prior used to keep class prevalences
                                    # away from zero or one; similar to those
                                    # used in LatentGOLD (Vermunt and Magidson,
                                    # 2005a,b) and in PROC LCA (Lanza et al.,
                                    # 2011)
                     getSEs=TRUE,   # Whether to calculate standard errors.
                                    # Setting this to FALSE would save
                                    # computational time.
                     gridSize=1000, # The number of values of time for which
                                    # to obtain estimates of the coefficient
                                    # functions
                     maxIterations=5000, # Maximum number of EM iterations
                                    # to attempt;
                     maxVarianceRatio=20, # Maximum ratio between estimated
                                    # variances of different classes;
                     numStarts=50,  # Number of random starting values to use.;
                     referenceClass=NA,  # Reference class.  If it is not
                                    # specified here, then the last class
                                    # (i.e., NumClasses) will be selected
                                    # as the reference class.  The logistic
                                    # regression model for the class membership
                                    # will predict the probability of being
                                    # in each other class relative to this
                                    # reference class as a baseline.
                     seed=NA,       # The initial random seed.  It will be used
                                    # to fit the model (if numStarts=1) or to
                                    # generate more seeds (if numStarts>1).
                     scov=NULL,     # Subject-level class membership prediction
                                    # covariates, as a matrix.  A column of 1's
                                    # should not be provided for scov because,
                                    # unlike in the case of tcov, an intercept
                                    # column is included automatically
                                    # by the code.  This should have the same
                                    # number of rows as tcov, but should be
                                    # identical on all observations within the
                                    # same person.
                     useRoughnessPenalty=TRUE, # Whether or not to use
                                    # a penalty function to make the estimated
                                    # time-varying coefficients more smooth
                     xcov=NULL      # Optional matrix of covariates assumed
                                    # to have time-invariant effects.
                     ) {

    ##################################################################
    # MixTVEM macro Version 1.1
    # By John DZIAK, Xianming TAN, and Runze LI
    # Fits a mixture of nonparametric trajectories to longitudinal data.
    #
    # Copyright:
    # (c) 2013 The Pennsylvania State University
    #
    # License:
    # This program is free software; you can redistribute it and/or
    # modify it under the terms of the GNU General Public License as
    # published by the Free Software Foundation; either version 2 of
    # the License, or (at your option) any later version.
    #
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    # General Public License for more details.
    #
    # Acknowledgments and references:
    # We fit a mixture of nonparametric varying-coefficient models, using
    # a penalized B-spline approach.  See
    #    Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing using
    #        B-splines and penalized likelihood. Statistical Science
    #        11(2): 89-121.
    #    Hastie, T., & Tibshirani, R. (1993). Varying-coefficient models. Journal
    #        of the Royal Statistical Society, Series B, 55, 757-796.
    #    Shiyko, M. P., Lanza, S. T., Tan, X., Li, R., Shiffman, S. (2012). Using
    #        the Time-Varying Effect Model (TVEM) to Examine Dynamic Associations
    #        between Negative Affect and Self Confidence on Smoking Urges:
    #        Differences between Successful Quitters and Relapsers. Prevention
    #        Science, 13, 288-299.
    #    Ramsay, J., Hooker, G., & Graves, S. (2009). Functional Data Analysis with
    #        R and MATLAB. New York: Springer.
    #    Tan, X., Shiyko, M. P., Li, R., Li, Y., & Dierker, L. (2011, November 21).
    #        A time-varying effect model for intensive longitudinal data.
    #        Psychological Methods. Advance online publication. doi: 10.1037/a0025814.
    # Estimation is done using the EM algorithm for finite mixtures.
    #    McLachlan, G. J., and Peel, D. (2000). Finite mixture models. New York: Wiley.
    #    Dempster, A. P., Laird, N. M., and Rubin, D. B. (1977). Maximum likelihood
    #        from incomplete data via the EM algorithm. Journal of the Royal
    #        Statistical Society, B, 39, 1-38.
    # The standard error calculations for the mixture approach are based on
    # those used by Turner (2000) and Turner's mixreg R package, which are based on the
    # ideas of Louis (1982).
    #    Louis, T. A. (1982). Finding the Observed Information Matrix when Using
    #        the EM Algorithm. Journal of the Royal Statistical Society, B, 44, 226-233.
    #    Turner, T. R. (2000) Estimating the rate of spread of a viral infection of
    #        potato plants via mixtures of regressions. Applied Statistics, 49,
    #        pp. 371-384.
    #    Turner, R. (2009). mixreg: Functions to fit mixtures of regressions.
    #        R package version 0.0-3. http://CRAN.R-project.org/package=mixreg
    # The use of a sandwich formula to adjust for within-subject correlation
    # when calculating the standard errors is inspired by
    #    Liang, K.-Y., and Zeger, S. L. (1986). Longitudinal data analysis
    #        using generalized linear models. Biometrika, 73, 13-22.
    # Clustering functional data modeled with splines is described in
    #    James, G., and Sugar, C. (2003) Clustering for sparsely sampled
    #        functional data.  Journal of the American Statistical
    #        Association 98, 397-408.
    # The optional gamma stabilization prior is similar to that used
    # in Latent GOLD and PROC LCA.  See
    #    Lanza, S. T., Dziak, J. J., Huang, L., Xu, S., & Collins, L. M. (2011).
    #        PROC LCA & PROC LTA users' guide (Version 1.2.7). University Park: The
    #        Methodology Center, Penn State. Retrieved from http://methodology.psu.edu.
    #    Chung, H., Flaherty, B. P., and Schafer, J. L. (2006). Latent class
    #        logistic regression: application to marijuana use and attitudes
    #        among high school seniors.  Journal of the Royal Statitical
    #        Society, A, 169, 723-743.
    #    Clogg, C. C. and Goodman, L. A. (1984) Latent structure analysis of a
    #        set of multidimensional contingency tables.  JASA, 79, 762–771.
    # The model fit criteria used are adapted versions of the standard AIC, BIC
    # and GCV of:
    #    Akaike, H. (1973). Information theory and an extension of the maximum
    #         likelihood principle. In B. N. Petrov & F. Csaki (Eds.), Second
    #         international symposium on information theory (p. 267-281).
    #         Budapest, Hungary: Akademai Kiado.
    #     Schwarz, G. (1978). Estimating the dimension of a model. Annals of
    #         Statistics, 6, 461-464.
    #     Craven, P., and Wahba, G. (1978). Smoothing noisy data with spline
    #         functions: Estimating the correct degree of smoothing by the
    #         method of generalized cross-validation. Numerische
    #         Mathematik, 31, 377–403.
    ##################################################################
  library(splines);
  library(nlme);
  ## Define required helper functions:
          MixTVEMFitInner <- function(convergenceCriterion,
                                    equalVariance=TRUE,
                                    gammaPriorConstant=0,
                                    intId,
                                    maxIterations=5000,
                                    maxVarianceRatio,
                                    modelOrder,
                                    numClasses,
                                    penaltyMatrix=NULL,
                                    referenceClass,
                                    roughnessPenalty=0,
                                    S,
                                    useRoughnessPenalty=FALSE,
                                    X,
                                    Y  ) {
           ## Prepare to begin loop;
               theta <- matrix(0,ncol(X),numClasses);
               oldGamma <- Inf;
               oldtheta <- Inf;
               maxAbsDev <- Inf;
               numSubjects <- max(intId);
               numTotal <- length(intId);
               nonreferenceClasses <- (1:numClasses)[-referenceClass];
               iteration <- 1;
               numObsBySub <- table(intId); # assumes that intId
                       # consists of consecutive integers starting at 1;
               stopifnot(length(Y)==numTotal);
               stopifnot(length(unique(intId))==numSubjects);
               stopifnot(all.equal(unique(intId),1:numSubjects));
               stopifnot(as.integer(rownames(numObsBySub))==1:numSubjects);
               stopifnot(identical(intId,sort(intId)));
               stopifnot(length(numObsBySub)==numSubjects);
               stopifnot(length(intId)==numTotal);
               stopifnot(identical(intId,sort(intId)));
           ## Initial E step (generate random posterior probabilities;
               temp <- matrix(rexp(numSubjects*numClasses),
                                   numSubjects,
                                   numClasses);
               postProbsBySub <- (temp+gammaPriorConstant/numSubjects)/
                              (apply(temp,1,sum)+
                              numClasses*gammaPriorConstant/numSubjects);
           ## Begin loop;
               while ((iteration<maxIterations)&
                       (maxAbsDev>convergenceCriterion)) {
                  ## Initial work;
                      iteration <- iteration+1;
                      postProbsByAssessment <- matrix(0,numTotal,numClasses);
                      for (k in 1:numClasses) {
                         postProbsByAssessment[,k] <-
                                 rep(postProbsBySub[,k],numObsBySub);
                      }
                  ## M Step;
                      ## Fit the model;
                          modelList <- list();
                          for (k in 1:numClasses) {
                             w <- postProbsByAssessment[,k];
                             if (!useRoughnessPenalty) {
                                 theta[,k] <- solve(t(X)%*%(w*X),t(X)%*%(w*Y))
                             } else {
                                 theta[,k] <- solve(t(X)%*%(w*X)+
                                   roughnessPenalty*penaltyMatrix,t(X)%*%(w*Y));
                             }
                          }
                          fittedY <- X%*%theta; # one column for each class;
                          residY <- as.matrix(Y-fittedY);
                      ## Get new sigma estimates;
                          f <- function(i) {
                             return(apply(residY[which(intId==i),,drop=FALSE]^2,
                                    2,sum))};
                          rssBySubjectAndClass <- t(sapply(1:numSubjects,f));
                          if (equalVariance==TRUE) {
                             sigsq <- rep(sum(postProbsBySub*
                                   rssBySubjectAndClass/numTotal),numClasses);
                          } else {
                             sigsq <- rep(0,numClasses);
                             for (k in 1:numClasses) {
                                sigsq[k] <- sum(rssBySubjectAndClass[,k]*
                                                       postProbsBySub[,k]) /
                                            sum(numObsBySub*postProbsBySub[,k]);
                            }
                          }
                          if (!is.na(maxVarianceRatio)) {
                              if (max(sigsq)/min(sigsq)>maxVarianceRatio) {
                                 sigsq[which(sigsq<(max(sigsq)
                                            /maxVarianceRatio))] <-
                                             max(sigsq)/maxVarianceRatio;
                              }
                          }
                      ## Get new gamma estimates;
                          SExtended <- kronecker(S,rep(1,numClasses));
                          colnames(SExtended) <- colnames(S);
                          classForLR <- rep(1:numClasses,times=numSubjects);
                          weightForLR <- as.vector(t(postProbsBySub));
                          dataForLR <- cbind(rep(1:numSubjects,each=numClasses),
                                                 classForLR,
                                                 weightForLR,
                                                 SExtended);
                          numGamma <- ncol(SExtended);
                          gamma <- matrix(0,numGamma,numClasses);
                          eta <- matrix(0,numSubjects,numClasses);
                          for (k in nonreferenceClasses) {
                             outcomeForLogisticRegression <- 1*(classForLR==k);
                             outcomeForLogisticRegression[(classForLR!=k)&
                                         (classForLR!=referenceClass)] <- NA;
                             warnOption <- getOption("warn");
                             options(warn=-1);
                             thisModel <- glm(outcomeForLogisticRegression~
                                                           SExtended+0,
                                              family=binomial,
                                              weights=weightForLR);
                             options(warn=warnOption);
                             gamma[,k] <- thisModel$coef;
                             eta[,k] <- S%*%gamma[,k];
                          }
                          fittedProb <- (exp(eta)+gammaPriorConstant/numSubjects)/
                                         (apply(exp(eta),1,sum)+
                                          numClasses*gammaPriorConstant/
                                                 numSubjects);
                  ## E Step (get new likelihood contributions
                  ## and posterior probabilities);
                      logProbability <- matrix(0,numSubjects,numClasses);
                      for (k in 1:numClasses) {
                         logProbability[,k] <- -(numObsBySub/2)*
                                     log(2*3.1415926535*(sigsq[k]+1e-30)) -
                                     rssBySubjectAndClass[,k]/(2*sigsq[k]+1e-30);
                      }
                      tempMatrix <- logProbability;
                      for (k in 1:numClasses) {
                         tempMatrix[,k] = tempMatrix[,k] + log(fittedProb[,k]+1e-30);
                      }
                      logLikelihoodBySubject <- log(apply(exp(tempMatrix),1,sum));
                      logLik <- sum(logLikelihoodBySubject);
                      tempMax <- apply(tempMatrix,1,max);
                      for (k in 1:numClasses) {
                         tempMatrix[,k] <- tempMatrix[,k] - tempMax;
                      }
                      expTempMatrix <- exp(tempMatrix);
                      for (k in 1:numClasses) {
                         postProbsBySub[,k] <- (expTempMatrix[,k]+
                                             gammaPriorConstant/numSubjects)/
                                               (apply(expTempMatrix,1,sum)+
                                             numClasses*gammaPriorConstant/
                                                numSubjects);
                      }
                      maxAbsDev <- max(abs(c(as.vector(theta)-oldtheta,
                                             as.vector(gamma)-oldGamma)));
                      oldGamma <- as.vector(gamma);
                      oldtheta <- as.vector(theta);
               }
               if (useRoughnessPenalty) {
                   enrpByClass <- rep(NA,numClasses);
                   for (k in 1:numClasses) {
                       w <- postProbsByAssessment[,k];
                       enrpByClass[k] <- sum(diag(as.matrix((X*sqrt(w))%*%
                                              solve(t(X)%*%(w*X)+
                                              roughnessPenalty*penaltyMatrix)%*%
                                              t(X*sqrt(w)))));
                   }
                   enrp <- sum(enrpByClass);
                   enp <- enrp +
                          ncol(S)*(numClasses-1)+
                          ifelse(equalVariance,1,numClasses);
                   np <- ncol(X)*numClasses + ncol(S)*(numClasses-1)+
                          ifelse(equalVariance,1,numClasses);
               } else {
                   enrp <- ncol(X)*numClasses;
                   enp <- enrp + ncol(S)*(numClasses-1)+
                          ifelse(equalVariance,1,numClasses);
                   np <- enp;
               }
               converged <- maxAbsDev<=convergenceCriterion;
               hardRSS <- 0;
               weightedRSS <- sum(postProbsBySub*rssBySubjectAndClass);
               stopifnot(nrow(postProbsByAssessment)==nrow(fittedY));
               stopifnot(length(Y)==nrow(fittedY));
               for (i in 1:nrow(fittedY)) {
                   hardRSS <- hardRSS + ((Y[i]-
                              fittedY[i,which.max(postProbsByAssessment[i,])])^2); 
                   weightedRSS <- weightedRSS + sum(postProbsByAssessment[i,]*((Y[i]-fittedY[i,])^2));
               }
               return(list(aic=-2*logLik+2*enrp,
                                  # A statistic similar to AIC (Akaike, 1973)
                                  # but using the effective (penalized) number of
                                  # parameters instead of a count of parameters
                                  # (see Eilers and Marx, 1996)
                           hardRSS=hardRSS,
                                  # The hard-classified residual sum of squares
                                  # error measure.  It is like the classic
                                  # regression sum of squared residuals
                                  # (sum(y-yhat)^2), but it uses yhat for the
                                  # class with the highest posterior probability
                                  # for each subject.
                           weightedRSS=weightedRSS,  
                                  # The weighted residual sum of squares
                                  # error measure.  It is like the classic
                                  # regression sum of squared residuals
                                  # (sum(y-yhat)^2), but it is weighted
                                  # by the posterior probabilities estimated
                                  # for each individual and class. 
                           bic=-2*logLik+log(numTotal)*enrp,
                                  # A statistic somewhat analogous to BIC
                                  # (Schwarz, 1978).  The effective (penalized)
                                  # number of parameters is used instead of a
                                  # count of parameters (see Eilers and Marx,
                                  # 1996).  Also, in order to select a simpler
                                  # model, the total number of observations is
                                  # used as the sample size rather than the
                                  # number of subjects.
                           converged=converged,
                                  # Indicates whether the EM algorithm converged
                           enp=enp,
                                  # The effective number of regression
                                  # parameters (enrp) plus the number of free
                                  # gamma parameters plus the number of free
                                  # sigma parameters, summed across classes.
                           enrp=enrp,
                                  # The effective number of regression
                                  # parameters (see Eilers and Marx, 1996)
                           fittedProb=fittedProb,
                                  # The fitted values for class probabilities
                                  # given the S covariates for each subject,
                                  # as given by the logistic regression model
                                  # for class membership.
                           fittedY=fittedY,
                                  # The fitted values for the dependent
                                  # variable at each assessment.
                           gamma=gamma,
                                  # The estimated logistic regression
                                  # parameters in the model for class
                                  # membership.                                     
                           hardGCV=hardRSS/(numTotal*((1-(enrp/numTotal))**2)),     
                                  # A statistic similar to GCV (see Wahba, 1990;
                                  # Eilers and Marx, 1996) but calculated
                                  # with the hardRSS instead of an ordinary
                                  # residual sum of squares.
                           iteration=iteration,
                                  # The number of iterations run by the EM
                                  # algorithm.
                           lambda=roughnessPenalty,
                                  # The current candidate value of the penalty 
                                  # weight for the roughness penalty
                           logLik=logLik,
                                  # The fitted log-likelihood
                           np=np,
                                  # The counted number of parameters,
                                  # ignoring shrinkage imposed by the penalty
                           postProbsBySub=postProbsBySub,
                                  # Each subject's estimated posterior
                                  # probability of belonging to each class.
                                  # As a caveat, these posterior probabilities
                                  # assume local independence within
                                  # subject conditional upon the class mean.
                                  # They are likely to exaggerate differences
                                  # between classes and overstate certainty
                                  # of class assignment.
                           postProbsByAssessment=postProbsByAssessment,
                                  # The same as postProbsBySub, except that
                                  # they are repeated within subject so that
                                  # there is one row per assessment rather than
                                  # only one row per subject.
                           sigsq=sigsq,
                                  # Estimated sigma squared (variance)
                                  # parameters for each class.
                           theta=theta,
                                  # Estimated regression coefficients for each
                                  # of the regression coefficients, including
                                  # spline basis coefficients     
                           weightedGCV=weightedRSS/(numTotal*((1-(enrp/numTotal))**2))
                                  # A statistic similar to GCV (see Wahba, 1990;
                                  # Eilers and Marx, 1996) but calculated
                                  # with the weightedRSS instead of an ordinary
                                  # residual sum of squares.
                        ));
        }
        MixTVEMCovMats <- function(  equalVariance,
                                     fittedProb,
                                     gamma,
                                     intId,
                                     mu,
                                     numSubjects,
                                     numTotal,
                                     postProb,
                                     referenceClass,
                                     S,
                                     sigsq,
                                     theta,
                                     X,
                                     Y) {
            numthetas <- nrow(theta);
            numClasses <- ncol(theta);
            stopifnot(numthetas==ncol(X));
            numGammas <- nrow(gamma);
            stopifnot(numClasses==ncol(mu));
            nonreferenceClasses <- (1:numClasses)[-referenceClass];
            if (!equalVariance) {
                numParams <- (numClasses*numthetas) +
                                (numClasses-1)*numGammas + numClasses;
            } else {
                numParams <- (numClasses*numthetas) +
                                 (numClasses-1)*numGammas + 1;
            }
            thetaIndex <- matrix(0,numthetas,numClasses);
            if (numClasses>1) {
                gammaIndex <- matrix(0,numGammas,numClasses-1);
            }
            if (equalVariance) {
               for (k in 1:numClasses) {
                  thetaIndex[,k] <- (((k-1)*numthetas+(k-1)*numGammas) +
                                        (1:numthetas));
                  if (k < numClasses) {
                      gammaIndex[,k] <- (k*numthetas+(k-1)*numGammas) +
                                              (1:numGammas);
                  }
               }
               sigmaIndex <- matrix(numClasses*numthetas+
                               (numClasses-1)*numGammas+1,1,numClasses);
                               # same index for each class;
            } else {
               sigmaIndex <- matrix(0,1,numClasses);
               for (k in 1:numClasses) {
                  thetaIndex[,k] <- ( ((k-1)*numthetas+(k-1)*numGammas+(k-1)) +
                                              (1:numthetas) );
                  if (k < numClasses) {
                     gammaIndex[,k] <- ( (k*numthetas+(k-1)*numGammas+(k-1)) +
                                           (1:numGammas) );
                     sigmaIndex[k] <- (k*numthetas+k*numGammas+(k-1)) + 1;
                  } else {
                     sigmaIndex[k] <- (k*numthetas+(k-1)*numGammas+(k-1)) + 1;
                  }
               }
            }
            # These calculations are based on the mixreg package in Turner
            # (2009) with some modifications as described in our paper;
            I1 <- matrix(0,numParams,numParams);
            I2 <- matrix(0,numParams,numParams);
            I3 <- matrix(0,numParams,numParams);
            # Calculate matrix I1
            for (k in 1:numClasses) {
               for (i in 1:numSubjects) {
                  xi <- X[which(intId==i),,drop=FALSE];
                  ni <- sum(intId==i);
                  contrib <- (1/sigsq[k])*postProb[i,k]*
                                     t(as.matrix(xi))%*%(as.matrix(xi));
                  I1[thetaIndex[,k],thetaIndex[,k]] <-
                                I1[thetaIndex[,k],thetaIndex[,k]] + contrib;
                  I1[sigmaIndex[,k],sigmaIndex[,k]] <-
                                I1[sigmaIndex[k],sigmaIndex[k]] +
                                postProb[i,k]*ni/(2*sigsq[k]*sigsq[k]);
               }
            }
            if (numClasses>1) {
               for (k in 1:(numClasses-1)) {
                  for (m in 1:(numClasses-1)) {
                     for (i in 1:numSubjects) {
                        nonrefk <- nonreferenceClasses[k];
                        nonrefm <- nonreferenceClasses[m];
                        Si <- S[i,,drop=FALSE];
                        contrib <- fittedProb[i,nonrefk]*
                                  (1*(nonrefk==nonrefm)-fittedProb[i,nonrefm])*
                                   crossprod(Si);
                        I1[gammaIndex[,k],gammaIndex[,m]] <-
                                  I1[gammaIndex[,k],gammaIndex[,m]] + contrib;
                     }
                  }
               }
            }
            # Calculate matrix I2
            for (k in 1:numClasses) {
               #print(paste("k=",k));print(Sys.time());
               for (m in 1:numClasses) {
                  #print(paste("m=",m));print(Sys.time());
                  for (i in 1:numSubjects) {
                     for (j in 1:numSubjects) {
                        ni <- sum(intId==i);
                        nj <- sum(intId==j);
                        xi <- X[which(intId==i),,drop=FALSE];
                        yi <- Y[which(intId==i),drop=FALSE];
                        xj <- X[which(intId==j),,drop=FALSE];
                        yj <- Y[which(intId==j),drop=FALSE];
                        Si <- S[i,,drop=FALSE];
                        Sj <- S[j,,drop=FALSE];
                        muik <- mu[which(intId==i),k,drop=FALSE];
                        mujm <- mu[which(intId==j),m,drop=FALSE];
                        hik <- (1/sigsq[k])*t(xi)%*%(yi-muik);
                        hjm <- (1/sigsq[m])*t(xj)%*%(yj-mujm);
                        tik <- -ni/(2*sigsq[k])+
                                     (crossprod(yi-muik))/(2*sigsq[k]*sigsq[k]);
                        tjm <- -nj/
                                     (2*sigsq[m])+(crossprod(yj-mujm))/
                                      (2*sigsq[m]*sigsq[m]);
                        Gkmij <- postProb[i,k]*postProb[j,m] +
                                  postProb[i,k]*(i==j)*((k==m)-postProb[j,m]);
                        I2[thetaIndex[,k],thetaIndex[,m]] <-
                                    I2[thetaIndex[,k],thetaIndex[,m]] +
                                      Gkmij * hik%*%t(hjm);
                        I2[thetaIndex[,k],sigmaIndex[m]] <-
                                     I2[thetaIndex[,k],sigmaIndex[m]] +
                                         Gkmij * hik%*%tjm;
                        I2[sigmaIndex[m],thetaIndex[,k]] <-
                                          t(I2[thetaIndex[,k],sigmaIndex[m]]);
                        I2[sigmaIndex[k],sigmaIndex[m]] <-
                                         I2[sigmaIndex[k],sigmaIndex[m]] +
                                          Gkmij*tik*tjm;
                        if ((k < numClasses)&(m < numClasses)) {
                           nonrefk <- nonreferenceClasses[k];
                           nonrefm <- nonreferenceClasses[m];
                           Gnonrefkmij <- postProb[i,nonrefk]*postProb[j,m] +
                                           postProb[i,nonrefk]*(i==j)*
                                           ((nonrefk==m)-postProb[j,m]);
                           Gknonrefmij <- postProb[i,k]*postProb[j, nonrefm] +
                                             postProb[i,k]*(i==j)*
                                             ((k==nonrefm)-postProb[j,nonrefm]);
                           Gnonrefknonrefmij <- postProb[i,nonrefk]*
                                                 postProb[j,nonrefm]+
                                                 postProb[i,nonrefk]*(i==j)*
                                                 ((nonrefk==nonrefm)-
                                                 postProb[j, nonrefm]);
                           contrib <- (Gnonrefknonrefmij-postProb[i,nonrefk]*
                                                 fittedProb[j,nonrefm]-
                                                  postProb[j,nonrefm]*
                                                  fittedProb[i,nonrefk] +
                                                   fittedProb[i,nonrefk]*
                                                   fittedProb[j,nonrefm])*
                                                   (t(Si)%*%Sj);
                           I2[gammaIndex[,k],gammaIndex[,m]] <-
                                        I2[gammaIndex[,k],gammaIndex[,m]] +
                                          contrib;
                        }
                        if (k < numClasses) {
                           nonrefk <- nonreferenceClasses[k];
                           Gnonrefkmij <- postProb[i,nonrefk]*postProb[j,m]+
                                       postProb[i,nonrefk]*(i==j)*
                                       ((nonrefk==m)-postProb[j,m]);
                           contrib <- t(Gnonrefkmij%*%tjm%*%Si);
                           I2[gammaIndex[,k],sigmaIndex[m]] <-
                                           I2[gammaIndex[,k],sigmaIndex[m]];
                           I2[sigmaIndex[m],gammaIndex[,k]] <-
                                     t(I2[gammaIndex[,k],sigmaIndex[m]]);
                        }
                        if (m < numClasses) {
                           nonrefm <- nonreferenceClasses[m];
                           Gknonrefmij <- postProb[i,k]*postProb[j,nonrefm]+
                                       postProb[i,k]*(i==j)*
                                       ((k==nonrefm)-postProb[j,nonrefm]);
                           contrib = Gknonrefmij*hik%*%Sj-
                                        postProb[i,k]*fittedProb[j,m]*hik%*%Sj;
                           I2[thetaIndex[,k],gammaIndex[,m]] <-
                                       I2[thetaIndex[,k],gammaIndex[,m]] +
                                        contrib;
                           I2[gammaIndex[,m],thetaIndex[,k]] <-
                                       t(I2[thetaIndex[,k],gammaIndex[,m]]);
                        }
                     }
                  }
               }
            }
            # Calculate matrix I3
            #print("Starting I3");print(Sys.time());
            for (i in 1:numSubjects) {
               scorei <- matrix(0,numParams,1);
               for (k in 1:numClasses) {
                  ni <- sum(intId==i);
                  xi <- X[which(intId==i),,drop=FALSE];
                  yi <- Y[which(intId==i),drop=FALSE];
                  Si <- S[i,,drop=FALSE];
                  muik <- mu[which(intId==i),k];
                  hik <- (1/sigsq[k])*t(xi)%*%(yi-muik);
                  scorei[thetaIndex[,k]] <- postProb[i,k]*hik;
                  tik <- -ni/(2*sigsq[k]) + crossprod(yi-muik)/
                                     (2*sigsq[k]*sigsq[k]);
                  scorei[sigmaIndex[k]] <- scorei[sigmaIndex[k]] +
                                          postProb[i,k]%*%tik;
                  if (k < numClasses) {
                      nonrefk <- nonreferenceClasses[k];
                      scorei[gammaIndex[,k]] <- (postProb[i,nonrefk]-
                                           fittedProb[i,nonrefk])*Si;
                  }
               }
               I3 <- I3 + crossprod(t(scorei));
            }
            # Calculate naiveCovarianceMatrix (in theta, gamma and sigma)
            if (kappa(I1-I2)<Inf) {
               naiveCovarianceMatrix <- solve(I1-I2);
            } else {
               naiveCovarianceMatrix <- matrix(NA,nrow(I1),ncol(I1));
            }
            # Calculate covarianceMatrix (in theta, gamma and sigma)
            if (det(I3)>0) {
               covarianceMatrix <- naiveCovarianceMatrix%*%
                                        I3%*%naiveCovarianceMatrix;
            }
            stopifnot(min(eigen(covarianceMatrix)$values)>0);
            # Calculate covarianceForPrevalences
            if (numClasses > 1) {
               allGammas <- as.vector(gammaIndex);
               naiveCovGamma <- naiveCovarianceMatrix[allGammas,allGammas];
               covGamma <- covarianceMatrix[allGammas,allGammas];
               jacobian <- matrix(0,numClasses,(numClasses-1)*numGammas);
               for (m in 1:(numClasses-1)) {
                  for (k in 1:numClasses) {
                      for (q in 1:numGammas) {
                        w <- (m-1)*numGammas + q;
                        nonrefm <- nonreferenceClasses[m];
                        # Jacobian[m,i] = derivative of lambda[k] in gamma[m,q], ;
                        #   i.e., gamma[w] when the free gammas are all listed in a ;
                        #   vector, one class on top of another, excluding the
                        #   reference class. ;
                        temp <- sum(((1*(k==nonrefm))-fittedProb[,k])*
                                    fittedProb[,nonrefm]*S[,q]);
                        jacobian[k,w] <- temp/numSubjects;
                       }
                   }
               }
                covForPrevalence <- jacobian %*% covGamma %*% t(jacobian);
                naiveCovForPrevalence <- jacobian %*% naiveCovGamma  %*%
                                                    t(jacobian);
            }
            return(list( naiveCovarianceMatrix=naiveCovarianceMatrix,
                            # The estimated covariance matrix for all of
                            # the model parameters, ignoring within-subject
                            # correlation entirely
                         covarianceMatrix=covarianceMatrix,
                            # The sandwich covariance matrix for the model
                            # parameters, which may be more robust to
                            # misspecification of the within-subject
                            # correlation (see Liang and Zeger 1986)
                         naiveCovForPrevalence=naiveCovForPrevalence,
                            # A naive covariance matrix for the estimated
                            # overall average class prevalences
                         covForPrevalence=covForPrevalence,
                            # A sandwich-based covariance matrix for the
                            # estimated overall average class prevalences
                         thetaIndex=thetaIndex,
                            # Tells which rows and columns of the covariance
                            # matrix refer to thetas (linear regression
                            # parameters)
                         gammaIndex=gammaIndex,
                            # Tells which rows and columns of the covariance
                            # matrix refer to gammas (logistic regression
                            # parameters)
                         sigmaIndex=sigmaIndex,
                            # Tells which rows and columns of the covariance
                            # matrix refer to sigmas (variance
                            # parameters)
                         I1=I1,
                         I2=I2,
                         I3=I3)); # I1, I2, and I3 are described in our paper
                                  # (Dziak et al) and are based on combining
                                  # the approach of Turner (2000) with a
                                  # sandwich estimation approach
        }
  #####################################################################
  ## Main body of MixTVEM function
  ## Process the subject ID's;
       if (!is.integer(id)) {
           temporary.id <- as.integer(as.factor(id));
       } else {
           temporary.id <- id;
       }
       numTotal <- length(id);
       intId <- rep(0,numTotal);
       stopifnot(is.integer(unique(temporary.id)));
       numSubjects <- length(unique(temporary.id));
       for (i in 1:numSubjects) {
          intId[which(temporary.id==unique(temporary.id)[i])] <- i;
       }
       stopifnot(length(unique(intId))==numSubjects);
       stopifnot(all.equal(unique(intId),1:numSubjects));
   ## Decide on a reference class;
       if (is.na(referenceClass)) {referenceClass <- numClasses;}
       stopifnot(referenceClass>0);
       stopifnot(referenceClass<=numClasses);
   ## Process the class membership predictor variables S;
       if (!is.null(scov)) {
          scov <- as.matrix(scov);
          stopifnot(nrow(scov)==numTotal);
          S <- matrix(0,numSubjects,1+ncol(scov));
          for (i in 1:numSubjects) {
             theseObs <- which(intId==i);
             S[i,] <- c(1, scov[min(theseObs),]);
          }
       } else {
          S <- matrix(1,numSubjects,1);
       }
       if (is.null(colnames(scov))) {
           colnames(S) <- paste("S",0:(ncol(S)-1),sep="");
       };
       ## Process the time variable;
       stopifnot(!is.null(time));
       stopifnot((deg==1)|(deg==2)|(deg==3));
       time <- as.matrix(time);
       stopifnot(nrow(time)==numTotal);
       stopifnot(ncol(time)==1);
       stopifnot(is.numeric(time));
       if (is.null(colnames(time))) {
          colnames(time) <- "Time";
       }
       timeBasis <- list();
       designMatrix <- NULL;
       whichBeta <-  NULL;
   ## Construct the first part of the design matrix including ... ;
       ## ... the intercept column if there are no time-varying covariates;
           if (is.null(tcov)) {        interceptColumn <- matrix(1,numTotal,1);
                colnames(interceptColumn) <- "Intercept";
                designMatrix <- cbind(designMatrix,interceptColumn);
                whichBeta <- c(whichBeta,0);
           }
       ## ... and the covariates without time-varying effects;
           if (!is.null(xcov)) {
              if (min(apply(xcov,2,var))<1e-10) {
                  stop("Please do not include a constant column in xcov.");
              }
              xcov <- as.matrix(xcov);
              stopifnot(nrow(xcov)==numTotal);
              numXCov <- ncol(xcov);
              if (is.null(colnames(xcov))) {
                 colnames(xcov) <- paste("X",1:ncol(xcov),sep="");
              };
              designMatrix <- cbind(designMatrix,xcov);
           } else {
              numXCov <- 0;
           }
           whichBeta <- c(whichBeta,rep(0,numXCov));
       ## Get the tcov matrix ready;
           if (is.null(tcov)) {
              numTCov <- 0;
              stop(paste("tcov is null.",
                          "No time-varying betas or effects are in the model."));
           } else {
              tcov <- as.matrix(tcov);
              stopifnot(nrow(tcov)==numTotal);
              numTCov <- ncol(tcov);
              if (is.null(colnames(tcov))) {
                 colnames(tcov) <- paste("TV",1:ncol(tcov),sep="");
              };
   ## Create the scale vector to be used with the penalty for the time-varying
   ## covariates;
      if (!identical(as.vector(tcov[,1]),rep(1,numTotal))) {
          warning("tcov does not seem to contain an intercept (trajectory) column.");
      }
      tcov.scale <- apply(tcov,2,sd);
      if (length(tcov.scale)>1) {
          if (min(tcov.scale[-1])<1e-10) {
              stop(paste("Please include the intercept column as the first",
                         "column in tcov, and do not include any other columns",
                         "which are constant across all",
                         "subjects and assessments."));
          }
      }
      if (identical(as.vector(tcov[,1]),rep(1,numTotal))) {
          tcov.scale[1] <- 1; # don't scale intercept column;
      }
   ## Create the basis functions;
      if (length(numInteriorKnots)!=1) {stop();}
      stopifnot(numInteriorKnots>0);
      num.intervals <- numInteriorKnots+1;
      dx <- (max(time)-min(time))/num.intervals;
      all.knot.locations <- seq(min(time)-(deg)*dx,
                                max(time)+(deg)*dx,
                                by=dx);
      all.knot.locations[1] <- all.knot.locations[2]-1e-8;
      all.knot.locations[length(all.knot.locations)] <-
                      all.knot.locations[length(all.knot.locations)-1]+1e-8;
      interior.knot.locations <- seq(min(time)+dx,
                                 max(time)-dx,
                                 by=dx);
      timeGrid <- seq(min(time),max(time),length=gridSize);
      timeBasis <- spline.des(all.knot.locations,
                              time,
                              deg+1,
                              rep(0,length(time)))$design;
                         # Adapted from Eilers and Marx, 1996;
      colnames(timeBasis) <- paste(colnames(time),".Spline.",
                                             1:ncol(timeBasis),
                                             sep="");
      timeBasisByGrid <- spline.des(all.knot.locations,
                              timeGrid,
                              deg+1,
                              rep(0,length(timeGrid)))$design;
                         # Adapted from Eilers and Marx, 1996;
      colnames(timeBasisByGrid) <- paste(colnames(time),".Spline.",
                                             1:ncol(timeBasisByGrid),
                                             sep="");
      # Now generate regression matrix;
      for (j in 1:numTCov) {
          covariateTimesTimeBasis <- tcov[,j,drop=TRUE]*timeBasis;
                                        # Elementwise product;
          colnames(covariateTimesTimeBasis) <- paste(colnames(tcov)[j],
                                                     "times",
                                                     colnames(timeBasis),
                                                     sep=".");
          designMatrix <- cbind(designMatrix, covariateTimesTimeBasis);
          whichBeta <- c(whichBeta,rep(j,ncol(covariateTimesTimeBasis)));
      }
      stopifnot(sum(is.nan(designMatrix))==0);
      stopifnot(sum(is.null(designMatrix))==0);
      if(sum(is.na(designMatrix))>0) {
          stop("Missing data is not yet supported in this software.");
      }
   }
   ## Create the penalty weight matrix (see Eilers and Marx, 1996)
   diffMatrix <- crossprod(diff(diff(diag(rep(1,ncol(timeBasis))))));
   penaltyMatrix <- matrix(0,ncol(designMatrix),ncol(designMatrix));
   if (useRoughnessPenalty) {
       for (j in 1:numTCov) {
           indices <- numXCov+(j-1)*ncol(timeBasis)+(1:ncol(timeBasis));
           stopifnot(length(indices)==nrow(diffMatrix));
           stopifnot(length(indices)==ncol(diffMatrix));
           penaltyMatrix[indices,indices] <- diffMatrix/(tcov.scale[j]^2);
       }
   }
   #######################################################################
   ## Run a loop to find the best seed;
   logLikBySeed <- rep(NA,numStarts);
   hardRSSBySeed <- rep(NA,numStarts);
   weightedRSSBySeed <- rep(NA,numStarts);
   if (numStarts > 1) {
       if (!is.na(seed)) {
       };
       seeds <- round(runif(numStarts)*1e8);
   } else {
       if (!is.na(seed)) {
           seeds <- seed;
       } else {
           seeds <- round(runif(1)*1e8);
       }
   }
   postProbsList <- list();
   for (i in 1:numStarts) {
       set.seed(seeds[i]);
       thisFit <- MixTVEMFitInner(convergenceCriterion=convergenceCriterion,
                                  equalVariance=equalVariance,
                                  intId=intId,
                                  gammaPriorConstant=gammaPriorConstant,
                                  maxIterations=maxIterations,
                                  maxVarianceRatio=maxVarianceRatio,
                                  modelOrder=modelOrder,
                                  numClasses=numClasses,
                                  penaltyMatrix=penaltyMatrix,
                                  referenceClass=referenceClass,
                                  S=S,
                                  useRoughnessPenalty=FALSE,
                                  X=designMatrix,
                                  Y=dep);
      logLikBySeed[i] <- thisFit$logLik;
      hardRSSBySeed[i] <- thisFit$hardRSS;
      weightedRSSBySeed[i] <- thisFit$weightedRSS;
      postProbsList[[i]] <- thisFit$postProbsBySub;
   }
   bestSeed <- seeds[which.max(logLikBySeed)];
   #######################################################################
   ## Find the best tuning parameter;
   if (useRoughnessPenalty) {
       f <- function(p) {
           set.seed(bestSeed);
           lambda <- exp(p);
           thisFit <- MixTVEMFitInner(convergenceCriterion=convergenceCriterion,
                                      equalVariance=equalVariance,
                                      intId=intId,
                                      gammaPriorConstant=gammaPriorConstant,
                                      maxIterations=maxIterations,
                                      maxVarianceRatio=maxVarianceRatio,
                                      modelOrder=modelOrder,
                                      numClasses=numClasses,
                                      penaltyMatrix=penaltyMatrix,
                                      referenceClass=referenceClass,
                                      roughnessPenalty=lambda,
                                      S=S,
                                      useRoughnessPenalty=TRUE,
                                      X=designMatrix,
                                      Y=dep);
           return(thisFit$weightedGCV);
       }
       opt1 <- optimize(f,interval=c(log(1e-10),log(1e10)),tol=.01);
       lambda <- exp(opt1$minimum);
   } else {
       lambda <- 0;
   }
   #######################################################################
   ## Finally do the analysis;
       set.seed(bestSeed);
       bestFit <- MixTVEMFitInner(convergenceCriterion=convergenceCriterion,
                                  equalVariance=equalVariance,
                                  intId=intId,
                                  gammaPriorConstant=gammaPriorConstant,
                                  maxIterations=maxIterations,
                                  maxVarianceRatio=maxVarianceRatio,
                                  modelOrder=modelOrder,
                                  numClasses=numClasses,
                                  penaltyMatrix=penaltyMatrix,
                                  referenceClass=referenceClass,
                                  roughnessPenalty=lambda,
                                  S=S,
                                  useRoughnessPenalty=useRoughnessPenalty, 
                                  X=designMatrix,
                                  Y=dep);
   ## Now record the fitted values on the original assessment times;
       fitted <- designMatrix%*%bestFit$theta;
       fittedCoefficients <- list();
       for (j in 1:numTCov) {
           fittedCoefficients[[j]] <- timeBasis%*%
                                      bestFit$theta[which(whichBeta==j),];
       }
       estimate <- list();
       for (which.coef in 1:numTCov) {
           estimate[[which.coef]] <- list();
           for (which.class in 1:numClasses) {
               estimate[[which.coef]][[which.class]] <-
                   fittedCoefficients[[which.coef]][order(time),which.class];
           }
        }
   ## Record the fitted values on the fit grid;
       fittedCoefficientsByGrid <- list();
       for (j in 1:numTCov) {
           fittedCoefficientsByGrid[[j]] <- timeBasisByGrid%*%
                                      bestFit$theta[which(whichBeta==j),];
       }
       estimate <- list();
       for (which.coef in 1:numTCov) {
           estimate[[which.coef]] <- list();
           for (which.class in 1:numClasses) {
               estimate[[which.coef]][[which.class]] <-
                   fittedCoefficientsByGrid[[which.coef]][order(timeGrid),which.class];
           }
        }
   ## Get the standard errors;
       if (getSEs) {
           covarMats <- MixTVEMCovMats(equalVariance=equalVariance,
                         fittedProb=bestFit$fittedProb,
                         gamma=bestFit$gamma,
                         intId=intId,
                         mu=bestFit$fittedY,
                         numSubjects=numSubjects,
                         numTotal=numTotal,
                         postProb=bestFit$postProbsBySub,
                         referenceClass=referenceClass,
                         S=S,
                         sigsq=bestFit$sigsq,
                         theta=bestFit$theta,
                         X=designMatrix,
                         Y=dep);
           fittedCoefficientsSE <- list();
           for (j in 1:numTCov) {
               fittedCoefficientsSE[[j]] <- matrix(NA,numTotal,numClasses);
               for (class in 1:numClasses) {
                   this.covmat <-
                        covarMats$covarianceMatrix[
                             covarMats$thetaIndex[which(whichBeta==j),class],
                             covarMats$thetaIndex[which(whichBeta==j),class]];
                   fittedCoefficientsSE[[j]][,class] <- sqrt(diag(
                                                         as.matrix(timeBasis %*%
                                                         this.covmat %*%
                                                         t(timeBasis))));
               }
           }
           fittedCoefficientsSEByGrid <- list();
           for (j in 1:numTCov) {
               fittedCoefficientsSEByGrid[[j]] <- matrix(NA,gridSize,numClasses);
               for (class in 1:numClasses) {
                   this.covmat <-
                        covarMats$covarianceMatrix[
                             covarMats$thetaIndex[which(whichBeta==j),class],
                             covarMats$thetaIndex[which(whichBeta==j),class]];
                   fittedCoefficientsSEByGrid[[j]][,class] <- sqrt(diag(
                                                         as.matrix(timeBasisByGrid %*%
                                                         this.covmat %*%
                                                         t(timeBasisByGrid))));
               }
           }
       } else {
           covarMats <- NULL;
           fittedCoefficientsSE <- NULL;
           fittedCoefficientsSEByGrid <- NULL;
       }
       if (doPlot==TRUE&length(time)>10) {
            stopifnot(numClasses<10);
            stopifnot(numTCov<4);
            theColors <- c("red","blue","darkgreen","purple","pink",
                           "cyan","magenta","green","darkgray");
            if (numTCov==2) {par(mfrow=c(2,1));}
            if (numTCov>2) {par(mfrow=c(2,2));}
            for (which.coef in 1:numTCov) {
              for (which.class in 1:numClasses) {
                if (which.class==1) {
                    plot(x=timeGrid[order(timeGrid)],
                         y=fittedCoefficientsByGrid[[which.coef]][
                                         order(timeGrid),which.class],
                         ylim=c(min(unlist(fittedCoefficientsByGrid[[which.coef]])),
                                max(unlist(fittedCoefficientsByGrid[[which.coef]]))),
                         col=theColors[which.class],type="l",
                         xlab="Time",
                         ylab=paste("Coefficient",which.coef));
                    markers <- round(quantile(1:length(timeGrid),(1:9)/10));
                    text(x=timeGrid[order(timeGrid)][markers],
                         y=fittedCoefficientsByGrid[[which.coef]]
                                             [order(timeGrid),which.class][markers],
                         col=theColors[which.class],
                         labels=which.class);
                } else {
                    lines(x=timeGrid[order(timeGrid)],
                         y=fittedCoefficientsByGrid[[which.coef]]
                                              [order(timeGrid),which.class],
                         col=theColors[which.class]);
                    markers <- round(quantile(1:length(time),(1:9)/10));
                    text(x=timeGrid[order(timeGrid)][markers],
                         y=fittedCoefficientsByGrid[[which.coef]][
                                    order(timeGrid),which.class][markers],
                         col=theColors[which.class],
                         labels=which.class);
                }
               }
            }
       }
       prop.best.logLik <- mean(logLikBySeed>max(logLikBySeed)-.01);
       prop.best.hardRSS <- mean(hardRSSBySeed>max(hardRSSBySeed)-.01);
       prop.best.weightedRSS <- mean(weightedRSSBySeed>max(weightedRSSBySeed)-.01);
       best.agree <- (hardRSSBySeed[which.max(logLikBySeed)] - min(hardRSSBySeed) < .01 ) & 
                     (hardRSSBySeed[which.min(weightedRSSBySeed)] - min(hardRSSBySeed) < .01);
       cat("MixTVEM R Function \n");
       cat(sprintf("Number of subjects: %29.0f \n",numSubjects));
       cat(sprintf("Total number of observations: %19.0f \n",numTotal));
       if(deg==1) {cat("Effect of time between knots treated as linear \n")};
       if(deg==2) {cat("Effect of time between knots treated as quadratic \n")};
       if(deg==3) {cat("Effect of time between knots treated as cubic \n")};
       if(equalVariance) {print("Equal variances assumed for all classes");}
       if(min(apply(bestFit$fittedProb,2,mean))<.025) {
           warning("The smallest fitted class is very small");
       }
       if (!bestFit$converged) {
           warning("EM algorithm did not converge");
       }
       cat(sprintf("Hard-classified squared error (RSS): %12.2f \n",
                   bestFit$hardRSS));
       cat(sprintf("Weighted RSS statistic: %25.2f \n", bestFit$weightedRSS));
       cat(sprintf("Weighted GCV statistic: %26.2f \n", bestFit$weightedGCV));
       cat(sprintf("Log-likelihood: %33.2f \n",bestFit$logLik));
       cat(sprintf("AIC: %44.2f \n",bestFit$aic));
       cat(sprintf("BIC: %44.2f \n",bestFit$bic));
       if (numStarts>1) {
           cat("Proportion of starting values giving approximately the best obtained ...\n");
           cat(sprintf(" ... log-likelihood: %28.2f \n",prop.best.logLik));
           cat(sprintf(" ... hard-classified squared error (RSS): %7.2f \n",
                       prop.best.hardRSS));
           cat(sprintf(" ... weighted RSS: %29.2f \n",prop.best.weightedRSS));
           if (!best.agree) {
               print("These criteria for the best starting value may not agree.");
               print("The starting value that gave the best log-likelihood was used.");
           }
       }
       cat(sprintf("Count number of parameters: %21.2f \n", bestFit$np));
       cat(sprintf("Smoothed number of parameters: %18.2f \n", bestFit$enp));
       cat("Class proportions: \n");
       cat(round(apply(bestFit$fittedProb,2,mean),4)); cat("\n");
       cat("Standard deviations: \n");
       cat(round(sqrt(bestFit$sigsq),4)); cat("\n");
   # Provide the logistic regression output;
    logistic.reg.output <- NULL;
    for (which.class in 1:numClasses) {
        if (which.class==referenceClass) {next;}
        this.class.est <- as.vector(bestFit$gamma[,which.class]);
        this.class.se <- as.vector(sqrt(diag(as.matrix(
                         covarMats$covarianceMatrix
                         [covarMats$gammaIndex[,1,drop=FALSE],
                          covarMats$gammaIndex[,1,drop=FALSE]]))));
        this.class.z <- this.class.est/this.class.se;
        this.class.p <- 2*(1-pnorm(abs(this.class.z)))
        this.class.output <- data.frame(Column=paste("S",1:length(this.class.est),
                                                      sep=""),
                                        Class=which.class,
                                        Estimate=this.class.est,
                                        SE=this.class.se,
                                        z=round(this.class.z,8),
                                        p=round(this.class.p,8));
        logistic.reg.output <- rbind(logistic.reg.output,this.class.output);
    }
    cat("Logistic regression for class membership: \n");
    print(logistic.reg.output);
    nontvem.reg.output <- NULL;
    if (!is.null(xcov)) {
        for (which.class in 1:numClasses) {
            if (which.class==referenceClass) {next;}
            this.class.est <- bestFit$theta[which(whichBeta==0),which.class];
            these.indices <- covarMats$thetaIndex[which(whichBeta==0),which.class];
            this.class.se <- as.vector(sqrt(diag(covarMats$covarianceMatrix[
                                      these.indices, these.indices,drop=FALSE])));
            this.class.z <- this.class.est/this.class.se;
            this.class.p <- 2*(1-pnorm(abs(this.class.z)))
            this.class.output <- data.frame(Column=paste("X",
                                                   1:length(this.class.est),sep=""),
                                            Class=which.class,
                                            Estimate=this.class.est,
                                            SE=this.class.se,
                                            z=round(this.class.z,8),
                                            p=round(this.class.p,8));
            nontvem.reg.output <- rbind(nontvem.reg.output,this.class.output);
        }
        cat("Non-time-varying coefficients: \n");
        print(nontvem.reg.output);
    }
   # Return the answers;
       return(list( bestFit=bestFit,
                           # Answer from the MixTVEMFitInner function,
                           # as obtained for the best starting value.
                    allBSplineKnots=all.knot.locations,
                           # List of all knot locations for the spline
                    bestSeed=bestSeed,
                           # Selected random number seed for the best
                           # solution obtained.
                    covarMats=covarMats,
                           # Answer from MixTVEMCovMats for the
                           # estimates obtained by the best starting values.
                    dep=dep,
                           # The dependent variable as a vector; if there
                           # are no NA's (missing data) in the data provided
                           # in the input, then this is the same as the dep
                           # vector that was provided.  It is provided again
                           # here in the output list for convenience in
                           # comparing observed to fitted values.
                    beta=fittedCoefficients,
                           # The fitted values for the coefficient functions
                           # at each assessment time value
                    betaByGrid=fittedCoefficientsByGrid,
                           # The fitted values for the coefficient functions
                           # at time value on a regular grid of GridSize points
                    betaSE=fittedCoefficientsSE,
                           # The standard errors for the coefficient functions
                           # at each assessment time value
                    betaSEByGrid=fittedCoefficientsSEByGrid,
                           # The standard errors for the coefficient functions
                           # at time value on a regular grid of GridSize points
                    fittedValues=fitted,
                           # The fitted values for the dependent variable
                           # for each class, at each assessment time value
                    hardRSSBySeed=hardRSSBySeed,
                           # Fitted hard-classified RSS for each random seed
                           # tried.  Hard-classified RSS means that each subject
                           # is assumed to be in his or her most likely class for
                           # purposes of obtaining a predicted value, and then
                           # the sum of squared residuals is calculated.
                    knotLocations=interior.knot.locations,
                           # Interior knots of the spline.  This may differ
                           # from allBSplineKnots because a B-spline may have
                           # invisible "exterior" knots which help to make the
                           # math work.
                    logisticRegOutput=logistic.reg.output,
                           # Information on the logistic regression for
                           # predicting class memberships.
                    logLikBySeed=logLikBySeed,
                           # Fitted log-likelihood for each random seed which was
                           # tried.
                    referenceClass=referenceClass,
                           # The reference class for the logistic regression for
                           # predicting class memberships.
                    S=S,
                           # The matrix of subject-level covariates.  For
                           # convenience, they are presented in a one-row-per
                           # subject format without duplication; an intercept
                           # column is also included.  Multiplying these
                           # by the fitted gammas from bestFit obtains the
                           # generalized linear model linear predictors
                           # for class membership.
                    time=as.vector(time),
                           # The vector of assessment times.  It is provided again
                           # here in the output list for convenience.
                    timeBasis=timeBasis,
                           # The basis matrix used to express the effect
                           # of time.
                    weightedRSSBySeed=weightedRSSBySeed,
                    whichBeta=whichBeta,
                           # Used for interpreting the results of bestFit and
                           # covarMats.  It tells which time-varying coefficient
                           # is represented by each column of the design matrix.
                    X=designMatrix
                           # The design matrix used in the regression, including
                           # all of the columns for all of the splines.
                    ));
}

##################################################################
