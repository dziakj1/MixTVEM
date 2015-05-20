%MACRO TVEM_Mix_Normal( MyData,
                        Time,
                        Dep,
                        ID,
                        Latent_Classes,
                        Cov=,
                        SCov=,
                        TCov=,
                        Knots=,
                        Num_Starts=50,
                        Convergence=1e-8 ,
                        Coverage=0.95 ,
                        Deg=3, 
                        Initial_Seed=100000,
                        Max_Iterations=1000,
                        Max_Variance_Ratio=10, 
                        Penalty_Max=.,
                        Ref=1,
                        Scale=1000,
                        Std_Err_Option=yes);
    /***************************************************************
    | MixTVEM macro Version 1.1
    | By John DZIAK, Xianming TAN, and Runze LI
    | Fits a mixture of nonparametric trajectories to longitudinal data.
    |
    | Copyright:
    | (c) 2015 The Pennsylvania State University
    |
    | License:
    | This program is free software; you can redistribute it and/or
    | modify it under the terms of the GNU General Public License as
    | published by the Free Software Foundation; either version 2 of
    | the License, or (at your option) any later version.
    |
    | This program is distributed in the hope that it will be useful,
    | but WITHOUT ANY WARRANTY; without even the implied warranty of
    | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    | General Public License for more details.
    |
    | Acknowledgments and references:
    | We fit a mixture of nonparametric varying-coefficient models, using
    | a penalized B-spline approach.  See
    |    de Boor, C. (1977). Package for calculating with B-splines. 
    |        SIAM Journal on Numerical Analysis, 14, 441-72.
    |    Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing using
    |        B-splines and penalized likelihood. Statistical Science
    |        11(2): 89-121.
    |    Hastie, T., & Tibshirani, R. (1993). Varying-coefficient models. Journal
    |        of the Royal Statistical Society, Series B, 55, 757-796.
    |    Shiyko, M. P., Lanza, S. T., Tan, X., Li, R., Shiffman, S. (2012). Using
    |        the Time-Varying Effect Model (TVEM) to Examine Dynamic Associations
    |        between Negative Affect and Self Confidence on Smoking Urges:
    |        Differences between Successful Quitters and Relapsers. Prevention
    |        Science, 13, 288-299.
    |    Ramsay, J., Hooker, G., & Graves, S. (2009). Functional Data Analysis with
    |        R and MATLAB. New York: Springer.
    |    Tan, X., Shiyko, M. P., Li, R., Li, Y., & Dierker, L. (2011, November 21).
    |        A time-varying effect model for intensive longitudinal data.
    |        Psychological Methods. Advance online publication. doi: 10.1037/a0025814.
    | Estimation is done using the EM algorithm for finite mixtures.
    |    McLachlan, G. J., and Peel, D. (2000). Finite mixture models. New York: Wiley.
    |    Dempster, A. P., Laird, N. M., and Rubin, D. B. (1977). Maximum likelihood
    |        from incomplete data via the EM algorithm. Journal of the Royal
    |        Statistical Society, B, 39, 1-38.
    | The standard error calculations for the mixture approach are based on
    | those used by Turner (2000) and Turner's mixreg R package, which are based on the
    | ideas of Louis (1982).
    |    Louis, T. A. (1982). Finding the Observed Information Matrix when Using
    |        the EM Algorithm. Journal of the Royal Statistical Society, B, 44, 226-233.
    |    Turner, T. R. (2000) Estimating the rate of spread of a viral infection of
    |        potato plants via mixtures of regressions. Applied Statistics, 49,
    |        pp. 371-384.
    |    Turner, R. (2009). mixreg: Functions to fit mixtures of regressions.
    |        R package version 0.0-3. http://CRAN.R-project.org/package=mixreg
    | The use of a sandwich formula to adjust for within-subject correlation
    | when calculating the standard errors is inspired by
    |    Liang, K.-Y., and Zeger, S. L. (1986). Longitudinal data analysis
    |        using generalized linear models. Biometrika, 73, 13-22.
    | Clustering functional data modeled with splines is described in
    |    James, G., and Sugar, C. (2003) Clustering for sparsely sampled
    |        functional data.  Journal of the American Statistical
    |        Association 98, 397-408.
    | Some web resources which were helpful in writing this macro:
    |    SAS FAQ: How can I recode my ID variable to be short and numeric?.
    |       UCLA: Academic Technology Services, Statistical Consulting Group.
    |       Accessed at http://www.ats.ucla.edu/stat/sas/faq/enumerate_id.htm.
    |       Accessed Feb. 22, 2012 (and before).
    |    Statistical Computing Seminar: Proc Logistic and Logistic Regression Models.
    |       UCLA: Academic Technology Services, Statistical Consulting Group.
    |       Accessed at http://www.ats.ucla.edu/stat/sas/seminars/sas_logistic/logistic1.htm.
    |       Accessed Feb. 22, 2012 (and before).
    |    Deng, C. Q. (June 27, 2009). Spaghetti plot. Blog entry on "On Biostatistics
    |       and Clinical Trials: CQ's web blog on the issues in biostatistics and
    |       clinical trials."  Accessed at
    |       http://onbiostatistics.blogspot.com/2009/06/spaghetti-plot.html.
    | Some of the code is adapted from the TVEM macro suite:
    |    TVEM SAS Macro Suite (Version 2.1.0) [Software]. (2012). University
    |       Park: The Methodology Center, Penn State.
    |       Retrieved from http://methodology.psu.edu.
    |    Yang, J., Tan, X., Li, R., & Wagner, A. (2012). TVEM (time-varying
    |       effect model) SAS macro suite users' guide (Version 2.1.0). University
    |       Park: The Methodology Center, Penn State.
    |       Retrieved from http://methodology.psu.edu.
    | The model fit criteria used are adapted versions of the standard AIC, BIC 
    | and GCV of:
    |    Akaike, H. (1973). Information theory and an extension of the maximum 
    |         likelihood principle. In B. N. Petrov & F. Csaki (Eds.), Second 
    |         international symposium on information theory (p. 267-281). 
    |         Budapest, Hungary: Akademai Kiado.
    |     Schwarz, G. (1978). Estimating the dimension of a model. Annals of 
    |         Statistics, 6, 461-464.
    |     Craven, P., and Wahba, G. (1978). Smoothing noisy data with spline 
    |         functions: Estimating the correct degree of smoothing by the 
    |         method of generalized cross-validation. Numerische 
    |         Mathematik, 31, 377–403.
    /***************************************************************/
    OPTIONS NONOTES;
    %LET MixTVEMDebug = 1;
    %LET SEIgnoreInfoLoss=0; /* do not try to make calculations faster by leaving I2 matrix at zero */
    /*      First process the user input values.  If these trigger an error
         right away (i.e., they appear to be un-implementable or nonsensical)
         then we do not do anything else.  This part of the code goes through
         the user-provided input arguments
         and checks whether the values are unreasonable or impossible.  If
         it finds impossible values, it gives an error message.  Otherwise,
         the checked input values are processed and saved in new macro
         variables and in a dataset named MixTVEMSettings. */
    PROC IML;
        MacroError = 0;
        /* MyData */
        MyData= "&MyData";
        IF (LENGTH("_%TRIM(&MyData)")<2) THEN DO;
            PRINT("Error --  The input dataset name must be provided.");
            MacroError=1;
        END;
        /* Time */
        Time = "&Time";
        IF (LENGTH("_%TRIM(&Time)")<2) THEN DO;
            PRINT("Error --  The time variable name must be provided.");
            MacroError=1;
        END;
        IF (LENGTH("_%TRIM(&Time)")>25) THEN DO;
            PRINT("Error --  The time variable name is too long.");
            MacroError=1;
        END;
        /* Dep */
        Dep = "&Dep";
        IF (LENGTH("_%TRIM(&Dep)")<2) THEN DO;
            PRINT("Error --  The response variable name must be provided.");
            MacroError=1;
        END;
        /* ID */
        ID = "&ID";
        IF (LENGTH("_%TRIM(&ID)")<2) THEN DO;
            PRINT("Error --  The subject ID variable name must be provided.");
            MacroError=1;
        END;
        /* Latent_Classes */
        IF (LENGTH("_%TRIM(&Latent_Classes)")<2) THEN DO;
            PRINT("Error --  The number of classes must be provided.");
            MacroError=1;
        END;
        NumClasses = &Latent_Classes;
        IF (NumClasses < 1) THEN DO;
            PRINT("Error --  There must be at least one class.");
            MacroError = 1;
        END;
        IF (NumClasses > 19) THEN DO;
            PRINT("Error --  Cannot fit that many classes.");
            MacroError = 1;
        END;
        %GLOBAL MixTVEMNumClasses;
        CALL SYMPUT("MixTVEMNumClasses",CHAR(NumClasses));
        /* Ref */
        IF (LENGTH("_%TRIM(&Ref)")<2) THEN RefClass = .;
            ELSE RefClass = &Ref;
        IF (RefClass < 1) THEN DO;
            PRINT("Error --  Ref must be at least 1.");
            MacroError = 1;
        END;
        IF (RefClass > NumClasses) THEN DO;
            PRINT("Error --  Ref cannot be greater than the number of classes.");
            MacroError = 1;
        END;
        /* Cov */
        CovInput = " &Cov"; /* note the trailing space; it makes things work*/
        Cov = TRANSLATE(STRIP(COMPBL(CovInput))," ",",");
        IF Cov ="" THEN NumCov = 0; ELSE NumCov = 1 + COUNT(Cov," ");
        %GLOBAL MixTVEMNumCov;
        CALL SYMPUT("MixTVEMNumCov",CHAR(NumCov));
        %GLOBAL MixTVEMCov;
        IF NumCov>0 THEN DO;
            CALL SYMPUT("MixTVEMCov",Cov);
        END; ELSE DO;
            CALL SYMPUT("MixTVEMCov","DummyPlaceholder");
        END;
        IF Cov ="" THEN Cov = .;
        /* SCov */
        SCovInput = " &SCov"; /* note the trailing space; it makes things work*/
        SCov = TRANSLATE(STRIP(COMPBL(SCovInput))," ",",");
        IF SCov ="" THEN NumSCov = 0; ELSE NumSCov = 1 + COUNT(SCov," ");
        %GLOBAL MixTVEMNumSCov;
        CALL SYMPUT("MixTVEMNumSCov",CHAR(NumSCov));
        %GLOBAL MixTVEMSCov;
        IF SCov ="" THEN SCov = .;
        IF NumSCov>0 THEN DO;
            CALL SYMPUT("MixTVEMSCov",SCov);
        END; ELSE DO;
            CALL SYMPUT("MixTVEMSCov","DummyPlaceholder");
        END;
        /* TCov */
        TCovInput = " &TCov"; /* note the trailing space; it makes things work*/
        TCov = TRANSLATE(STRIP(COMPBL(TCovInput))," ",",");
        IF TCov ="" THEN NumTCov = 0; ELSE NumTCov = 1 + COUNT(TCov," ");
        %GLOBAL MixTVEMNumTCov;
        CALL SYMPUT("MixTVEMNumTCov",CHAR(NumTCov));
        %GLOBAL MixTVEMTCov;
        IF NumTCov>0 THEN DO;
            CALL SYMPUT("MixTVEMTCov",TCov);
        END; ELSE DO;
            CALL SYMPUT("MixTVEMTCov","DummyPlaceholder");
        END;
        IF TCov ="" THEN TCov = .;
        /* Knots */
        Knots= "&Knots";
        IF (NumTCov>0) THEN DO;
                 IF (LENGTH("_%TRIM(&Knots)")<2) THEN DO;
                      PRINT("Error --  Knot information must be provided.");
                      MacroError = 1;
                 END;
             /* Have to create dataset of knots */
             USE &MyData;
                  READ ALL VAR {&Time} INTO Time;
             CLOSE &MyData;
             KnotsString = TRANSLATE("'&Knots'"," ",",");
             CALL SYMPUT("MixTVEMTemp",KnotsString);
             NumKnotsVec = {. &Knots}`;  /* The extra entry is here to prevent a macro error */
             NumKnotsVec=NumKnotsVec[2:NROW(NumKnotsVec)];
             IF NumKnotsVec[><]<1 THEN DO;
                      PRINT ("Error -- Please use at least one knot per time-varying coefficient.");
                       MacroError = 1;
             END;
             StartTime = Time[><];
             StopTime = Time[<>];
             IF (NumTCov^=NROW(NumKnotsVec)) THEN DO;
                  PRINT ("Error -- Please check the knots string.");
                   MacroError = 1;
             END;
                 InteriorKnotsMatrix = J(NumKnotsVec[<>],NumTCov,.);
                 DO k = 1 TO NumTCov;
                     IF NumKnotsVec[k]>0 THEN DO;
                         InteriorKnotsMatrix[1:NumKnotsVec[k],k] = (StartTime + ((1:NumKnotsVec[k]))*
                                                         (StopTime-StartTime)/(NumKnotsVec[k]+1))`;
                     END;
                 END;
                 CREATE MixTVEMInteriorKnotsMatrix FROM InteriorKnotsMatrix;
                 APPEND FROM InteriorKnotsMatrix;
                 CLOSE MixTVEMInteriorKnotsMatrix;
                 KnotsDatasetName = "MixTVEMInteriorKnotsMatrix";
        END;
        /* Deg */
        IF (LENGTH("_%TRIM(&Deg)")<2) THEN Deg = .;
        Deg = &Deg;
        IF (Deg = .) THEN DO;
            Deg = 1;
        END;
        IF (^((Deg=1)|(Deg=2)|(Deg=3))) THEN DO;
            PRINT("Error --  Cannot fit that order of penalized spline.");
            MacroError = 1;
        END;
        /* Standard error option */
        StdErrOptionInput = "&Std_Err_Option";
        IF ( SUBSTR(UPCASE(STRIP(StdErrOptionInput )),1,1) = "Y" |
                  SUBSTR(UPCASE(STRIP(StdErrOptionInput )),1,1) = "S" |
                  SUBSTR(UPCASE(STRIP(StdErrOptionInput )),1,1) = "1" ) THEN DO;
            StdErrOption= 1;
        END; ELSE IF ( SUBSTR(UPCASE(STRIP(StdErrOptionInput )),1,1) = "N" |
                  SUBSTR(UPCASE(STRIP(StdErrOptionInput )),1,1) = "0" ) THEN DO;
            StdErrOption= 0;
        END;
        ELSE DO;
            PRINT("Error --  Please provide yes or no for Std_Err_Option.");
            MacroError = 1;
        END;
        IF (^((StdErrOption=0)|(StdErrOption=1))) THEN DO;
            PRINT("Error --  This standard error option is not currently supported.");
            MacroError = 1;
        END;
         /* Coverage */
        /* "Coverage" here means nominal coverage.  It is not to be taken too
           literally, because the confidence intervals which are constructed are
           only pointwise intervals and not a familywise confidence band for the
           whole curve, and because the confidence intervals are conditional on
           the regression model being correct (e.g., knots are in the right place,
           response is normal with variance not a function of time).  In reality
           these assumptions will be at least slightly inaccurate. Also,
           concerns about under-coverage of the sandwich matrix (see Kauermann
           and Carroll (2001, JASA) should be considered.*/
        IF (LENGTH("_%TRIM(&Coverage)")<2) THEN DO;
            Coverage = .95;
        END; ELSE DO;
            Coverage = &Coverage;
            IF (Coverage < 0.50) THEN DO;
                PRINT("Error --  The specified coverage level is too low.");
                MacroError = 1;
            END;
            IF (Coverage > 0.99999) THEN DO;
                PRINT("Error --  The specified coverage level is too high.");
                MacroError = 1;
            END;
        END;
        /* Convergence */
        IF (LENGTH("_%TRIM(&Convergence)")<2) THEN DO;
            Convergence = 1e-8;
        END; ELSE DO;
            Convergence = &Convergence;
            IF (Convergence < 1e-15) THEN DO;
                PRINT("Error --  The specified convergence criterion is too small.");
                MacroError = 1;
            END;
            IF (Convergence > 1e-4) THEN DO;
                PRINT("Error --  The specified convergence criterion is too large.");
                MacroError = 1;
            END;
        END;
        /* MaxIterations */
        MaxIterations = &Max_Iterations;
        IF (MaxIterations < 50) THEN DO;
            PRINT("Error --  The specified maximum number of iterations is too small.");
                MacroError = 1;
        END;
        /* Max_Variance_Ratio */
        MaxVarianceRatio = &Max_Variance_Ratio;
        IF (MaxVarianceRatio = .) THEN MaxVarianceRatio = 1000;
        IF (MaxVarianceRatio < 0) THEN DO;
            PRINT("Error --  The value supplied for Max_Variance_Ratio should not be negative.");
            MacroError = 1;
        END;
        /* Penalty_Max */
        Penalty_Max = &Penalty_Max;
        IF ((Penalty_Max < 0)&(Penalty_Max ^= .)) THEN DO;
            PRINT("Error --  The value supplied for Penalty_Max should not be negative.");
            MacroError = 1;
        END;
        /* Scale */
        ScaleForGraph = &Scale;
        IF (ScaleForGraph < 0) THEN DO;
            PRINT("Error --  The value supplied for Scale should not be negative.");
            MacroError = 1;
        END;
        /* Initial_Seed */
        IF (LENGTH("_%TRIM(&Initial_Seed)")<2) THEN DO;
            InitialSeed = .;
        END;
        InitialSeed = &Initial_Seed;
        IF (InitialSeed = .) THEN DO;
            InitialSeed = 0;
        END;
        IF (InitialSeed < 0) THEN DO;
            PRINT("Error --  Need a nonnegative random seed.");
            MacroError = 1;
        END;
        SelectedSeed = .;
        /* Num_Starts */
        IF (LENGTH("_%TRIM(&Num_Starts)")<2) THEN DO;
            PRINT("Error --  The number of starting values must be provided.");
            MacroError=1;
        END;
        NumStarts = &Num_Starts;
        IF (NumStarts < 1) THEN DO;
            PRINT("Error --  There must be at least one starting value.");
            MacroError = 1;
        END;
        IF ((NumStarts=1)&(InitialSeed=0)) THEN DO;
            PRINT("Error --  Please either specify a higher Num_Starts or a nonzero Initial_Seed.");
            MacroError = 1;
        END;
        Iterations = .; /* Not filled in until later */
        Lowest_Num_Observations = .; /* Not filled in until later */
        Highest_Num_Observations = .; /* Not filled in until later */
        CREATE MixTVEMSettings VAR {
            Convergence
            Cov
            Coverage
            Deg
            Dep
            Highest_Num_Observations
            ID
            InitialSeed
            Iterations
            Knots
            Lowest_Num_Observations
            MacroError
            MaxIterations
            MaxVarianceRatio
            NumClasses
            NumCov
            NumSCov
            NumStarts
            NumTCov
            Penalty_Max
            RefClass
            RoughnessPenalty
            ScaleForGraph
            SCov
            StdErrOption
            Tcov};
        APPEND;
        CLOSE MixTVEMSettings;
    QUIT;
    DATA MixTVEMData;
        SET &MyData;
        KEEP &ID &Time
             %IF %EVAL(&MixTVEMNumCov>0) %THEN %DO; &MixTVEMCov %END;
             %IF %EVAL(&MixTVEMNumSCov>0) %THEN %DO; &MixTVEMSCov %END;
             %IF %EVAL(&MixTVEMNumTCov>0) %THEN %DO; &MixTVEMTCov %END;
             &Dep;
        WHERE ((&ID IS NOT MISSING) &
               (&Time IS NOT MISSING) &
               %DO MixTVEMThisVarIndex = 1 %TO &MixTVEMNumCov;
                   (%QSCAN(&MixTVEMCov,&MixTVEMThisVarIndex," ") IS NOT MISSING) &
               %END;
               %DO MixTVEMThisVarIndex = 1 %TO &MixTVEMNumSCov;
                   (%QSCAN(&MixTVEMSCov,&MixTVEMThisVarIndex," ") IS NOT MISSING) &
               %END;
               %DO MixTVEMThisVarIndex = 1 %TO &MixTVEMNumTCov;
                   (%QSCAN(&MixTVEMTCov,&MixTVEMThisVarIndex," ") IS NOT MISSING) &
               %END;
               (&Dep ^= .));
    RUN;
    PROC IML;
        USE MixTVEMData; READ ALL INTO MixTVEMData; CLOSE MixTVEMData;
        USE MixTVEMSettings; READ VAR {InitialSeed
                                       MacroError
                                       MaxIterations
                                       NumStarts
                                       StdErrOption };
        CLOSE MixTVEMSettings;
        %GLOBAL MixTVEMInitialSeed; CALL SYMPUT("MixTVEMInitialSeed",CHAR(InitialSeed));
        %GLOBAL MixTVEMMacroError; CALL SYMPUT("MixTVEMMacroError",CHAR(MacroError));
        %GLOBAL MixTVEMMaxIterations; CALL SYMPUT("MixTVEMMaxIterations",CHAR(MaxIterations));
        %GLOBAL MixTVEMNumStarts; CALL SYMPUT("MixTVEMNumStarts",CHAR(NumStarts));
        %GLOBAL MixTVEMNumTotal; CALL SYMPUT("MixTVEMNumTotal",CHAR(NROW(MixTVEMData)));
        %GLOBAL MixTVEMStdErrOption; CALL SYMPUT("MixTVEMStdErrOption",CHAR(StdErrOption));
    QUIT;
    PROC SORT DATA=MixTVEMData TAGSORT NOTHREADS;
        BY &ID;
    RUN;
    DATA MixTVEMData; SET MixTVEMData;
        BY &ID;
        RETAIN MixTVEMIntID 0;
        IF first.&ID THEN MixTVEMIntID = MixTVEMIntID + 1;
            /* This creates a recoded version of the subject ID variable,
               called "MixTVEMIntID," that is composed of consecutive
               integers starting with 1.  This makes it possible for some
               of the code in the macro to be simpler than if subjects could
               have anything for ID's.
             See http://www.ats.ucla.edu/stat/sas/faq/enumerate_id.htm; */
        OUTPUT;
        %GLOBAL MixTVEMNumSubjects;
        CALL SYMPUT("MixTVEMNumSubjects",MixTVEMIntID);
    RUN;
    DATA MixTVEMSettings;
        SET MixTVEMSettings;
        NumTotal = &MixTVEMNumTotal;
        NumSubjects = &MixTVEMNumSubjects;
        IF (Penalty_Max = .) THEN Penalty_Max = NumTotal;
    RUN;
    /* Now create a datafile, MixTVEMSBySub, which contains
    the S predictors if there are any.  Unlike MixTVEMData and MixTVEMDataX,
    which have one record per observation, this dataset will be much smaller,
    having only one record per subject.  This is because class membership is
    considered to be a stable characteristic that does not change.  Thus, its
    predictors in the model are not expected to be time-varying either. */
    %IF %EVAL(&MixTVEMDebug>0) %THEN %PUT Now prepare subject variables matrix;
    PROC IML;
        USE MixTVEMSettings;
            READ VAR{ NumSCov } INTO  NumSCov ;
            READ VAR{ SCov } INTO  SCov ;
            READ VAR{ NumSubjects } INTO  NumSubjects ;
            READ VAR{ NumTotal } INTO  NumTotal ;
            READ VAR{ MacroError} INTO  MacroError;
        CLOSE MixTVEMSettings;
        USE MixTVEMData;
            READ ALL VAR{MixTVEMIntID} INTO IntID;
        CLOSE MixTVEMData;
        AvailableCov = CONTENTS(MixTVEMData) ;
        SCovArray = "Intercept";
        IF NumSCov>0 THEN DO;
            DO i = 1 TO NumSCov;
                ThisVariable = SCANQ(SCov,i);
                VarFound = (COMPARE(ThisVariable,AvailableCov,"li")=0)[+];
                IF VarFound = 0 THEN DO;
                    PRINT "Error -- Could not find this variable in the dataset:" ThisVariable;
                    MacroError = 1;
                END;
                IF (COMPARE(ThisVariable,"Intercept","li")=0) THEN DO;
                    PRINT("Error --  An intercept variable should not be included in the S variables.");
                    MacroError = 1;
                    /* The user should not have provided an intercept column for S,
                       since one will be provided automatically and we do not want
                       redundant columns. */
                END;
                SCovArray = SCovArray // ThisVariable;
            END;
            USE MixTVEMData; READ ALL VAR {&MixTVEMSCov} INTO SCov ; CLOSE MixTVEMSCov ;
            IF (((SCov[<>]=SCov[><])[+])>0) THEN DO;
                PRINT("Error --  An intercept or constant column should not be included in the S variables.");
                MacroError = 1;
            END;
            IF (NROW(SCov)^=NumTotal) THEN DO;
                PRINT("Error --  Unexpected amount of data in S variables.");
                MacroError = 1;
            END;
            SByObs = J(NumTotal,1,1) || SCov;
        END; ELSE DO;
            SByObs = J(NumTotal,1,1);
        END;
        SNum = (1:(NumSCov+1))`;
        SCovName = SCovArray;
        CREATE MixTVEMSCovNames VAR{SNum SCovName};
            APPEND;
        CLOSE MixTVEMSCovNames;
        IF (NROW(SByObs) ^= NumTotal)|(NCOL(SByObs) ^= NumSCov+1) THEN DO;
            PRINT("Error --  Unexpected amount of data in S variables.");
            MacroError = 1;
        END;
        IF (NROW(IntID) ^= NumTotal)|(NCOL(IntID) ^= 1) THEN DO;
            PRINT("Error --  Unexpected amount of data in ID variable.");
            MacroError = 1;
        END;
        IF (IntID[><] ^= 1)|(IntID[<>]^=NumSubjects) THEN DO;
            PRINT("Error --  Unexpected values of data in ID variable.");
            MacroError = 1;
        END;
        SBySub = J(NumSubjects,NCOL(SByObs),0);
        ThisRowOkayForThisSubject = J(NumTotal,1,0);
        sub = 1;
        SBySub[sub,] = SByObs[1,];
        GaveErrorMessage = 0;
        /* The following loop is somewhat tricky but involves compressing
        the dataset which the user provides (which has NumTotal rows,
        one for each observation on each subject), into a new dataset
        with NumSubject rows (having only NumSubject rows, i.e., one row
        for each subject, where NumSubject is less than or equal to NumTotal.)
        The loop also checks to try to make sure that this is meaningful and
        does not cause confusion due to misspecified ID or S variables.  In other
        words, the subjects in MixTVEMData should have a one-to-one mapping
        onto the subjects in MixTVEMSBySub, and the observations in MixTVEMData
        should have a many-to-one mapping to the values in MixTVEMData.*/
        DO obs = 2 TO NumTotal;
            IF IntID[obs]<IntID[obs-1] THEN DO;
                PRINT("Error --  Possible bug -- Internal ID variable not sorted.");
                MacroError = 1;
            END;
            IF IntID[obs]>IntID[obs-1] THEN DO;
                /* Start new subject */
                sub = sub + 1;
                SBySub[sub,] = SByObs[obs,];
            END; ELSE DO;
                /* Continuing old subject.   Check to make sure that the S variables
                   are constant within-subject. */
                IF SBySub[sub,] ^= SByObs[obs,] THEN DO;
                    IF (GaveErrorMessage=0) THEN DO;
                        PRINT("Error -- The S variables do not appear to be constant within subject.");
                        PRINT "Check observation # " obs " and subject # " sub "in the dataset.";
                        GaveErrorMessage = 1;
                    END;
                    MacroError = 1;
                END;
            END;
        END;
        IF((SBySub[<>,]=SBySub[><,])[+]>1) THEN DO;
            PRINT("Error -- One of the S variables (not the intercept) seems to be a constant.");
            MacroError = 1;
        END;
        IF(NCOL(SBySub)^=NROW(SCovArray)) THEN DO;
            PRINT("Error -- Difficulty deciding number of S variable columns.");
            MacroError = 1;
        END;
        IF MacroError>0 THEN DO;
            EDIT MixTVEMSettings VAR{MacroError}; MacroError = MacroError; REPLACE; CLOSE MixTVEMSettings;
        END;
        CREATE MixTVEMSBySub FROM SBySub [ COLNAME = SCovArray ]; APPEND FROM SBySub; CLOSE MixTVEMSBySub;
    QUIT;
    %IF %EVAL(&MixTVEMDebug>0) %THEN %PUT Finished preparing S matrix;
    %IF %EVAL(&MixTVEMDebug>0) %THEN %PUT Now prepare random seeds;
    DATA MixTVEMRandomSeeds;
        DO SeedID = 1 TO &MixTVEMNumStarts;
            Seeds = ROUND(RANUNI(&MixTVEMInitialSeed)*1000000);
            LogLik = .;
            RSS = .;
            GCV = .;
            AIC = .;
            BIC = .;
            SmallestPrevalence = .;
            LargestPrevalence = .;
            Estim_ProportionNugget = .;
            Estim_Rho = .;
            OUTPUT;
        END;
    RUN;
    %IF %EVAL(&MixTVEMDebug>0) %THEN %PUT Finished preparing random seeds;
    PROC IML;
        %IF %EVAL(&MixTVEMDebug>0) %THEN %PUT Start preparing regressor matrix;
        START SafeExp(x);
            y = .*x;
            IF ((x< 30)[+]>0) THEN y[LOC(x< 30)] = EXP(x[LOC(x<30)]);
            IF ((x>=30)[+]>0) THEN y[LOC(x>=30)] = EXP(30);
            RETURN(y);
        FINISH;
        START nls_function(p) GLOBAL(crossprods,lags,weights);
            propnonnugget = SafeExp(p[1])/(1+SafeExp(p[1]));
            logrho = p[2];
            f = (weights#((crossprods-propnonnugget*SafeExp(logrho*lags))##2))[+]; 
            RETURN(f);
        FINISH nls_function; 
        START CreateWeightingMatrices;                        
            Inv_Y_CorrMats_Stacked = J(NumTotal,Highest_Num_Observations,0);
            Det_Y_Covmats = J(NumSubjects,1,0);
            DO i = 1 TO IntID[<>];
                these=LOC(IntID=i);
                ni = SUM(IntID=i);
                ti = Time[these];
                ti_matrix = REPEAT(ti,1,NROW(ti));
                IF (Working_Rho>1E-5) THEN DO;
                    y_covmat_these = (1-Working_ProportionNugget)*(Working_Rho**(ABS(ti_matrix - ti_matrix`))) + 
                                     Working_ProportionNugget*I(NROW(ti));                 
                    Inv_Y_CorrMats_Stacked[these,1:ni] = INV(y_covmat_these);
                    Det_Y_Covmats[i] = DET(y_covmat_these);
                END; ELSE DO;
                    Inv_Y_CorrMats_Stacked[these,1:ni] = I(NROW(ti));
                    Det_Y_Covmats[i] = 1;
                END;    
            END; 
            CREATE Inv_Y_CorrMats_Stacked FROM Inv_Y_CorrMats_Stacked;
                APPEND FROM Inv_Y_CorrMats_Stacked;
            CLOSE Inv_Y_CorrMats_Stacked;
            CREATE Det_Y_Covmats FROM Det_Y_Covmats;
                APPEND FROM Det_Y_Covmats;
            CLOSE Det_Y_Covmats;
        FINISH CreateWeightingMatrices;
        START my_Bspline(x, knots, d);
            /*    evaluate the values of (num of knots) + d +1 B-spline basis functions at x, x is a real number*/
            /*    output a vector of (num of knots) + 1 - d real numbers */
            /*  knots include both exterior and inner knots */
            n_ext = ncol(knots);
            all_knots = (1:(n_ext+2));
            all_knots[1] = knots[1] - (1e-12);
            all_knots[2:(n_ext+1)]= knots;
            all_knots[n_ext+2] = knots[n_ext] + (1e-12);
            n_ext =n_ext+2;
            tmp=0.0*(1:(n_ext-1));
            DO i = 1 TO (n_ext-1);
                IF (x>=all_knots[i] & x<all_knots[i+1]) THEN tmp[i]=1.0;
            END;
            j=1;
            DO WHILE (j<=d); /* the De Boor formula, refer to Eilers & Marx (1996) */
                DO i = 1 TO (n_ext-j-1);
                    w1 = (x-all_knots[i])/(all_knots[i+j] - all_knots[i]);
                    w2 = (all_knots[i+j+1]-x)/(all_knots[i+j+1] - all_knots[i+1]);
                    tmp[i] = w1*tmp[i] +w2*tmp[i+1];
                END;
                j=j+1;
            END;
            RETURN(tmp[1:(n_ext-d-1)]);
        FINISH my_Bspline;
        START my_vec_Bspline(xx, knots, d);
            /*    evaluate the values of the (num of knots) + d +1 B-spline basis functions at xx */
            /*    xx is a vector of real numbers*/
            /*    output a real matrix with dim(xx) rows,  and (num of knots) + d +1 columns */
            n_row = ncol(xx);
            n_col = ncol(knots)+2-d-1;
            out = J(n_row, n_col, .);
            DO i=1 TO n_row;
                out[i, ] = t(my_Bspline(xx[i], knots, d));
            END;
            RETURN(out);
        FINISH my_vec_Bspline;
        START GenerateTimeBasis(TimeBasis,Time, InteriorKnots,deg);
            * B-Spline basis;
            numIntervals = NCOL(InteriorKnots)+1;
            dx = (time[<>]-time[><])/numIntervals;
            EarlyKnots = InteriorKnots[><]-dx*((deg):1);
            LateKnots = InteriorKnots[<>]+dx*(1:(deg));
            AllKnots = EarlyKnots || InteriorKnots || LateKnots;
            CREATE MixTVEMAllBSplineKnots FROM AllKnots;
               APPEND FROM AllKnots;
            CLOSE MixTVEMAllBSplineKnots;
            TimeBasis = my_vec_Bspline(Time`, AllKnots, deg);
        FINISH GenerateTimeBasis;
        START PrepareRegressionMatrix;
            USE MixTVEMData;
                READ ALL VAR {&time} INTO Time;
            CLOSE MixTVEMData;
            USE MixTVEMSettings;
                READ ALL VAR {deg MacroError NumCov
                              NumTCov RoughnessPenalty};
            CLOSE MixTVEMSettings;
            USE MixTVEMInteriorKnotsMatrix;
                READ ALL INTO InteriorKnotsMatrix;
            CLOSE MixTVEMInteriorKnotsMatrix;
            IF (NumCov > 0) THEN DO;
                USE MixTVEMData;
                    READ ALL VAR{&MixTVEMCov} INTO Cov;
                CLOSE MixTVEMData;
                Regressors = Cov;
                whichBeta = J(NCOL(Cov),1,.);
                whichThetaWithinBeta = (1:NCOL(Cov))`;
            END;
            IF (NumTCov > 0) THEN DO;
                USE MixTVEMData;
                    READ ALL VAR{&MixTVEMTCov} INTO TCov;
                CLOSE MixTVEMData;
                tcovStdDev = J(NCOL(TCov),1,1);
                DO thisTCov = 1 TO NCOL(TCov);
                    IF NROW(TCov)>1 THEN DO;
                        TcovStdDev[thisTCov] = SQRT(((TCov[,thisTCov]-TCov[:,thisTCov])[##])/(NROW(TCov[,thisTCov])-1));
                    END;
                    IF (ABS(TcovStdDev[thisTCov])<1e-20) THEN DO;
                        IF tcov[<>,thisTCov]=1 THEN DO;
                            TcovStdDev[thisTCov] = 1; /* do not rescale intercept term */
                        END;
                    END;
                    InteriorKnots = InteriorKnotsMatrix[LOC(InteriorKnotsMatrix[,thisTCov]^=.), thisTCov]`;
                    CALL GenerateTimeBasis(TimeBasis,
                                           Time,
                                           InteriorKnots,
                                           deg);
                    datasetName = COMPRESS(CONCAT("MixTVEMTimeBasis",CHAR(thisTCov)));
                    commandText = CONCAT("CREATE ",
                                         datasetName,
                                         " FROM TimeBasis;",
                                         " APPEND FROM TimeBasis;",
                                         " CLOSE ",
                                         datasetName,
                                         ";");
                    CALL EXECUTE(commandText);
                    covariateTimesTimeBasis = TCov[,thisTCov]#TimeBasis;
                    IF ((NumCov=0)&(thisTCov =1)) THEN DO;
                        Regressors = covariateTimesTimeBasis;
                        whichBeta = J(NCOL(TimeBasis),1,thisTCov);
                        whichThetaWithinBeta = (1:NCOL(TimeBasis))`;
                    END; ELSE DO;
                        Regressors = Regressors || covariateTimesTimeBasis;
                        whichBeta = whichBeta // J(NCOL(TimeBasis),1,thisTCov);
                        whichThetaWithinBeta = whichThetaWithinBeta // (1:NCOL(TimeBasis))`;
                    END;
                END;
            END;
            PenaltyMatrix = J(NROW(whichBeta),NROW(whichBeta),0);
            IF whichBeta[<>]>0 THEN DO;
                DO i = 1 TO whichBeta[<>];
                    IF (whichBeta[i] >0) THEN DO;
                       n_theta = NCOL(LOC(whichBeta=i));
                       PenaltyMatrix[LOC(whichBeta=i),LOC(whichBeta=i)] =
                          CreatePenaltyMatrix(n_theta);
                    END;
                END;
            END;
            TVScales = J(nrow(whichBeta),1,1);
            DO i = 1 TO whichBeta[<>];
                IF tcovStdDev[i]>0 THEN TVScales[LOC(whichBeta=i)]=1/(tcovStdDev[i]+1e-20);
            END;
            PenaltyMatrix = TVScales#TVScales#PenaltyMatrix;
            CREATE MixTVEMRegressors FROM Regressors; APPEND FROM Regressors; CLOSE MixTVEMRegressors;
            CREATE MixTVEMWhichBeta FROM whichBeta; APPEND FROM whichBeta; CLOSE MixTVEMWhichBeta;
            CREATE MixTVEMWhichThetaWithinBeta FROM whichThetaWithinBeta; APPEND FROM whichThetaWithinBeta; CLOSE MixTVEMWhichThetaWithinBeta;
            CREATE MixTVEMPenaltyMatrix FROM PenaltyMatrix; APPEND FROM PenaltyMatrix; CLOSE MixTVEMPenaltyMatrix;
            EDIT MixTVEMSettings VAR{MacroError};
                MacroError = MacroError;
                REPLACE;
            CLOSE MixTVEMSettings;
        FINISH PrepareRegressionMatrix;
        START CreatePenaltyMatrix(n);
           /* The R version is
                 crossprod(diff(diff(diag(as.vector(rep(1,n))))));
              see Eilers & Marx (1996) p. 100 */
           M = J(n,n,0);
           IF n>2 THEN DO;
               DO r = 1 TO n;
                    DO c = 1 TO n;
                        IF r = c THEN DO;
                           M[r,c] = 6;
                           IF ((r=1)|(r=n)) THEN M[r,c] = 1;
                           IF ((r=2)|(r=n-1)) THEN DO;
                                IF n>3 THEN M[r,c] = 5; ELSE M[r,c] = 4;
                           END;
                        END;
                        IF (ABS(r-c)=1) THEN DO;
                           M[r,c] = -4;
                           IF ((r=1)|(c=1)|(r=n)|(c=n)) THEN M[r,c] = -2;
                        END;
                        IF (ABS(r-c)=2) THEN DO;
                            M[r,c] = 1;
                        END;
                    END;
               END;
            END;
            RETURN(M);
        FINISH CreatePenaltyMatrix;
        START WeightedLogistic(b,
                               w, /* w,x,y, are all for observations */
                               x,
                               y);
            localw = w[LOC(y^=.)];
            localx = x[LOC(y^=.),];
            localy = y[LOC(y^=.)];
            b = J(NCOL(localx),1,0);
            MaxAbsDev = 10000000;
            CriterionMaxAbsDev = 1e-8;
            MaxIter=1000;
            iter=0;
            DO UNTIL((iter=maxIter) | (maxAbsDev<=criterionMaxAbsDev));
              iter = iter + 1;
                 bOld = b;
                 expEta = EXP(((localx*b)><20)<>-20);
                 mu = expEta/(1+expEta);
                 b = b + SOLVE(localx`*(mu#(1-mu)#localw#localx),
                               localx`*(localw#(localy-mu)));
                 maxAbsDev = MAX(ABS(b-bOld));
            END;
        FINISH WeightedLogistic;
        START WeightedRidgeRegression(b,
                                      w,  /* weights for subjects, not observations! */
                                      x,  /* x,y are for observations */ 
                                      y,
                                      sigsqtotal_thisclass,
                                      RoughnessPenalty,
                                      PenaltyMatrix,
                                      intID,
                                      Time) GLOBAL(Inv_Y_CorrMats_Stacked); 
            denominator = J(NCOL(x),NCOL(x),0);
            numerator = J(NCOL(x),1,0);
            DO i = 1 TO intID[<>];
                these = LOC(IntID=i);
                ti = Time[these,];
                ti_matrix = REPEAT(ti,1,NROW(ti));
                ni = SUM(IntID=i);
                inv_y_covmat_these = Inv_Y_CorrMats_Stacked[these,1:ni]/sigsqtotal_thisclass;
                xi = x[these,];
                yi = y[these,];
                denominator = denominator + w[i]*xi`*inv_y_covmat_these*(xi);
                numerator = numerator + w[i]*xi`*inv_y_covmat_these*(yi);
            END;
            b = SOLVE( denominator + RoughnessPenalty*PenaltyMatrix,
                       numerator ); 
        FINISH WeightedRidgeRegression;
        START WeightedHatMatrixTrace(enrp_this_class,
                                     wbysub,
                                     wbyobs,
                                     x,
                                     y,
                                      sigsqtotal_thisclass,
                                     RoughnessPenalty,
                                     PenaltyMatrix,
                                     intID,
                                     Time) GLOBAL(Inv_Y_CorrMats_Stacked);
                /* Called once per class */ 
            XTWX = J(NCOL(x),NCOL(x),0);
            DO i = 1 TO intID[<>];
                these = LOC(intID=i);
                ni = SUM(IntID=i);
                inv_y_covmat_these = Inv_Y_CorrMats_Stacked[these,1:ni]/sigsqtotal_thisclass;
                xi = x[these,];
                yi = y[these,];
                XTWX = XTWX + wbysub[i]*xi`*inv_y_covmat_these*(xi);
            END;
            denominator = XTWX + RoughnessPenalty*PenaltyMatrix ; 
            inv_denominator = INV(denominator); 
            enrp_this_class = 0;
                /* Effective number of regression parameters for this class */
            DO i = 1 TO intID[<>];
                these = LOC(intID=i);
                ni = SUM(IntID=i);
                inv_y_covmat_these = Inv_Y_CorrMats_Stacked[these,1:ni]/sigsqtotal_thisclass;
                xi = x[these,]; 
                temp = xi*inv_denominator*(xi`)*inv_y_covmat_these;
                enrp_this_class = enrp_this_class + wbysub[i]*TRACE(temp);
            END;                    
        FINISH;
        START CopyToLongForm(longTarget,shortSource,intID);
            longTarget = J(NROW(intID),NCOL(shortSource),0);
            DO I = 1 TO NROW(shortSource);
                longTarget[LOC(intID=i),] = shortSource[i,] @ J(SUM(intID=i),1,1);
            END;
        FINISH CopyToLongForm;
        START InitialRandomEStep;
            postProbsBySubject = J(NumSubjects,NumClasses,0);
            CALL RANDGEN(postProbsBySubject,"exponential");
            postProbsBySubject = (1/postProbsBySubject[,+])#postProbsBySubject;
            SigSqTotal = J(NumClasses,1,1); /* arbitrary initial values of 1 */
        FINISH InitialRandomEStep;
        START MStep;
            /* ... for Theta ... */
            Theta = J(NCOL(X),numClasses,.); 
            DO class = 1 TO numClasses; 
                CALL WeightedRidgeRegression(b,
                                             postProbsBySubject[,class],
                                             X,
                                             Y,
                                             SigSqTotal[class],
                                             RoughnessPenalty,
                                             PenaltyMatrix,
                                             intID,
                                             Time);
                Theta[,class] = b;
            END;
            fittedY = X*Theta;
            /* ... for sigma squareds ... */ 
            sqdResidY = (fittedY - Y@J(1,numClasses,1))##2;
            rssBySubjectAndClass = J(numSubjects,numClasses,0);
            numObsBySubject = J(numSubjects,1,0);
            DO i = 1 TO numSubjects;
                rssBySubjectAndClass[i,] = (sqdResidY[LOC(intId=i),])[+,];
                numObsBySubject[i] = NCOL(LOC(intId=i));
            END;
            SigSqTotal = ((rssBySubjectAndClass#postProbsBySubject)[+,])/
                           ((numObsBySubject#postProbsBySubject)[+,]);
            IF ((SigSqTotal[<>]/SigSqTotal[><])>maxVarianceRatio) THEN DO;
                 SigSqTotal[LOC(SigSqTotal<(SigSqTotal[<>]/maxVarianceRatio))]=SigSqTotal[<>]/maxVarianceRatio;
            END;
            /* ... and for gammas ... */ 
            USE MixTVEMSBySub; READ ALL INTO S; CLOSE MixTVEMSCov;
            SForLogisticRegression = S @ J(numClasses,1,1);
            ClassForLogisticRegression = SHAPE((1:numClasses) @ J(numSubjects,1,1),
                                               numClasses*numSubjects);
            WeightsForLogisticRegression = SHAPE(PostProbsBySubject,numClasses*numSubjects);
            IF (numClasses>1) THEN DO;
                NonRefClasses = (1:numClasses)[LOC((1:numClasses)^=refClass)];
            END;
            gamma = J(NCOL(S),numClasses,0);
            IF (NumClasses>1) THEN DO;
                expEta = J(numSubjects,numClasses,1);
                DO NonRefClassIndex = 1 TO NROW(NonRefClasses);
                    ThisNonRefClass = NonRefClasses[NonRefClassIndex];
                    OutcomeForLogisticRegression = 1*(ClassForLogisticRegression=ThisNonRefClass);
                    IF (NROW(LOC((ClassForLogisticRegression^=ThisNonRefClass)
                                &(ClassForLogisticRegression^=RefClass)))>0) THEN DO;
                          OutcomeForLogisticRegression[LOC((ClassForLogisticRegression^=ThisNonRefClass)
                                                          &(ClassForLogisticRegression^=RefClass))] = .;
                    END;
                    CALL WeightedLogistic(b,
                                          WeightsForLogisticRegression,
                                          SForLogisticRegression,
                                          OutcomeForLogisticRegression);
                    gamma[,ThisNonRefClass] = b;
                    expEta[,ThisNonRefClass] = SafeExp(S*b);
                END;
                fittedProb = expEta / expEta[,+];
           END; ELSE DO;
                fittedProb = J(numSubjects,1,1);
           END;
        FINISH MStep;
        START EStep;
            /* Get new likelihood contributions and posterior probabilities; */
            logProbabilityBySubjectAndClass = J(numSubjects,numClasses,0);
            IF (Working_Rho<1e-20) THEN Working_Rho =1e-20;  
                                /* This is mainly to avoid problems with raising 0 to the 0th power */
            DO class = 1 TO numClasses; 
                DO i = 1 TO intID[<>];
                    these = LOC(intID=i);
                    numObs = NCOL(these);
                    ti = Time[these,];
                    mui = fittedY[these,class];
                    yi = y[these,];
                    ni = SUM(IntID=i);
                    inv_y_covmat_these = Inv_Y_CorrMats_Stacked[these,1:ni]/SigSqTotal[class]; 
                    thisWeightedSumSquares = (yi-mui)`*inv_y_covmat_these*(yi-mui);
                    determinant = (SigSqTotal[class]**ni)*Det_Y_Covmats[i];
                    logProbabilityBySubjectAndClass[i,class] = -(numObs/2)*LOG(2*3.1415926535) -
                                                                0.5*(LOG(determinant+1e-30)) -
                                                                0.5*thisWeightedSumSquares;    
                END;
            END;
            tempMatrix = logProbabilityBySubjectAndClass + LOG(fittedProb+1e-30);
            logLik = (LOG((EXP(tempMatrix))[,+]))[+];
            tempMax = tempMatrix[,<>];
            tempMatrix = tempMatrix-tempMax;
            expTempMatrix = EXP(tempMatrix);
            postProbsBySubject = expTempMatrix / expTempMatrix[,+]; 
        FINISH EStep;
        START EMLoop; 
            OldTheta = 1e20;
            OldGamma = 1e20;
            maxAbsDev = 1e20;
            iteration = 0;
            Converged = 0;
            DO UNTIL ((maxAbsDev<Convergence)|(iteration > MaxIterationsToUse));
                iteration = iteration + 1;
                postProbsByObservation = J(NumTotal,NumClasses,0);
                CALL CopyToLongForm(postProbsByObservation,postProbsBySubject,intID); 
                CALL MStep; 
                CALL EStep; 
                maxAbsDev = MAX( (ABS(gamma-oldGamma))[<>,<>], (ABS(theta-oldTheta))[<>,<>]);
                oldGamma = Gamma;
                oldTheta = Theta;
            END;
            IF (maxAbsDev<Convergence) THEN DO;
                Converged = 1;
            END;
            ClassPrevalences = fittedProb[+,]/numSubjects;
        FINISH EMLoop;
        START FitTheModel;
            CALL InitialRandomEStep;
            CALL EMLoop;
            RSS = (postProbsBySubject#rssBySubjectAndClass)[+,+];
            BestClass = postProbsBySubject[,<:>];
            sumTrace = 0;
            DO class = 1 TO NumClasses; 
                CALL WeightedHatMatrixTrace(enrp_this_class,
                                   postProbsBySubject[,class],
                                   postProbsByObservation[,class],
                                   X,
                                   y,
                                   SigSqTotal[class],
                                   RoughnessPenalty,
                                   PenaltyMatrix,
                                   intID,
                                   Time);
                sumTrace = sumTrace + enrp_this_class;
            END;
            EffectiveNumParams = sumTrace  + ncol(S)*(numClasses-1)+ numClasses+2;
            CountNumParams = NCOL(X)*numClasses + ncol(S)*(numClasses-1)+ numClasses+2;
            AIC = -2*logLik + 2*EffectiveNumParams;
            BIC = -2*logLik + LOG(NumSubjects)*EffectiveNumParams;
            GCV = RSS/(numTotal*((1-(EffectiveNumParams/numTotal))**2));
            smallestPrevalence = ClassPrevalences[,><];
            largestPrevalence = ClassPrevalences[,<>];
        FINISH FitTheModel;
        START SaveFittedValues;
            fittedValues = X * Theta;
            temp = Time || IntID || fittedValues;
            CREATE MixTVEMFittedValues FROM temp
               [COLNAME=("&Time" || "MixTVEMIntID" || 
                         COMPRESS(CONCAT("Class",CHAR(1:NumClasses))))];
               APPEND FROM temp;
            CLOSE MixTVEMFittedValues;
        FINISH SaveFittedValues;
        START SaveFittedValuesSEs;
            IF StdErrOption = 1 THEN DO;
                fittedValuesSEs = J(NumTotal,NumClasses,0);
                DO class = 1 TO NumClasses;
                    DO i = 1 TO NumSubjects;
                        CovMatOfFittedValues = X[LOC(intId=i),] *
                                               CovarianceMatrix[thetaIndex[,class],
                                                                thetaIndex[,class]] *
                                               X[LOC(intId=i),]`;
                        fittedValuesSEs[LOC(intId=i),class] = SQRT(VECDIAG(CovMatOfFittedValues));
                    END;
                END;
                temp = Time || IntID || fittedValuesSEs;
                CREATE MixTVEMFittedValueSEs FROM temp
                   [COLNAME=("&Time" || "MixTVEMIntID" || 
                         COMPRESS(CONCAT("SE_Class",CHAR(1:NumClasses))))];
                   APPEND FROM temp;
                CLOSE MixTVEMFittedValueSEs;
            END;
        FINISH SaveFittedValuesSEs;
        START SaveFittedBetasAndSEs;
            DatasetIn = COMPRESS(CONCAT("MixTVEMTimeBasis",CHAR(thisTCov)));
            CALL EXECUTE(CONCAT("USE ",DatasetIn,";"));
            CALL EXECUTE("READ ALL INTO TimeBasis;");
            CALL EXECUTE(CONCAT("CLOSE ",DatasetIn,";"));
            /* Calculate estimates of time-specific regression coefficients */
            fittedBetas = TimeBasis*Theta[LOC(WhichBeta=thisTCov),];
            DatasetOut = COMPRESS(CONCAT("MixTVEMFittedBeta",CHAR(thisTCov)));
            temp = Time || fittedBetas;
            CALL EXECUTE(CONCAT("CREATE ",
                                DatasetOut,
                                " FROM temp",
                                "[COLNAME=('&Time' || COMPRESS(CONCAT('Class',CHAR(1:NumClasses))))];"));
            CALL EXECUTE("APPEND FROM temp;");
            CALL EXECUTE(CONCAT("CLOSE ",DatasetOut,";"));
            /* Calculate estimates of their SEs */
            IF StdErrOption = 1 THEN DO;
                fittedBetaSEs = J(NumTotal,NumClasses,0);
                DO class = 1 TO NumClasses;
                    DO i = 1 TO NumSubjects;
                        CovMatOfFittedBetas = TimeBasis[LOC(intId=i),] *
                                               CovarianceMatrix[thetaIndex[LOC(WhichBeta=thisTCov),class],
                                                                thetaIndex[LOC(WhichBeta=thisTCov),class]] *
                                               TimeBasis[LOC(intId=i),]`;
                        fittedBetaSEs[LOC(intId=i),class] = SQRT(VECDIAG(CovMatOfFittedBetas));
                    END;
                END;
                DatasetOut = COMPRESS(CONCAT("MixTVEMFittedBetaSE",CHAR(thisTCov)));
                temp = Time || fittedBetaSEs;
                CALL EXECUTE(CONCAT("CREATE ",
                                    DatasetOut,
                                    " FROM temp",
                                    "[COLNAME=('&Time' || COMPRESS(CONCAT('SE_Class',CHAR(1:NumClasses))))];"));
                CALL EXECUTE("APPEND FROM temp;");
                CALL EXECUTE(CONCAT("CLOSE ",DatasetOut,";"));
            END;
        FINISH SaveFittedBetasAndSEs;
        START SaveGridFittedBetasAndSEs;
            HasInteriorKnots = 0;
            InteriorKnots = .;
            IF (SUM(InteriorKnotsMatrix[,thisTCov]^=.)>0) THEN DO;
                HasInteriorKnots = 1;
                InteriorKnots = InteriorKnotsMatrix[LOC(InteriorKnotsMatrix[,thisTCov]^=.), thisTCov]`;
            END;
            CALL GenerateTimeBasis(GridTimeBasis ,
                                   GridTime,
                                   InteriorKnots,
                                   deg);
            GridFittedBetas = GridTimeBasis*Theta[LOC(WhichBeta=thisTCov),];
            DatasetOut = COMPRESS(CONCAT("MixTVEMGridFittedBeta",CHAR(thisTCov)));
            temp = GridTime || GridFittedBetas;
            CREATE MixTVEMGridTimeBasis FROM GridTimeBasis;
                APPEND FROM GridTimeBasis;
            CLOSE MixTVEMGridTimeBasis;
            CALL EXECUTE(CONCAT("CREATE ",
                                DatasetOut,
                                " FROM temp",
                                "[COLNAME=('&Time' || COMPRESS(CONCAT('Class',CHAR(1:NumClasses))))];"));
            CALL EXECUTE("APPEND FROM temp;");
            CALL EXECUTE(CONCAT("CLOSE ",DatasetOut,";"));
            /* Calculate estimates of their SEs */
            IF StdErrOption = 1 THEN DO;
                GridFittedBetaSEs = J(NROW(GridTimeBasis),NumClasses,0);
                DO class = 1 TO NumClasses;
                    DO i = 1 TO NROW(GridTimeBasis);
                        CovMatOfGridFittedBetas = GridTimeBasis[i,] *
                                               CovarianceMatrix[thetaIndex[LOC(WhichBeta=thisTCov),class],
                                                                thetaIndex[LOC(WhichBeta=thisTCov),class]] *
                                               GridTimeBasis[i,]`;
                        GridFittedBetaSEs[i,class] = SQRT(VECDIAG(CovMatOfGridFittedBetas));
                    END;
                END;
                DatasetOut = COMPRESS(CONCAT("MixTVEMGridFittedBetaSE",CHAR(thisTCov)));
                temp = GridTime || GridFittedBetaSEs;
                CALL EXECUTE(CONCAT("CREATE ",
                                    DatasetOut,
                                    " FROM temp",
                                    "[COLNAME=('&Time' || COMPRESS(CONCAT('SE_Class',CHAR(1:NumClasses))))];"));
                CALL EXECUTE("APPEND FROM temp;");
                CALL EXECUTE(CONCAT("CLOSE ",DatasetOut,";"));
            END;
        FINISH SaveGridFittedBetasAndSEs;
        START DoParameterIndexing;
            numThetas = NROW(theta);
            numClasses = NCOL(theta);
            numGammas = NROW(gamma);
            IF (numClasses>1) THEN DO;
                NonRefClasses = (1:numClasses)[LOC((1:numClasses)^=refClass)];
            END;
            numParams = (numClasses*numThetas) + (numClasses-1)*numGammas;
            thetaIndex = J(numThetas,numClasses,0);
            IF (numClasses>1) THEN DO;
                gammaIndex = J(numGammas,numClasses-1,0);
            END;
            DO class = 1 TO numClasses;
                thetaIndex[,class] = ( ((class-1)*numThetas+(class-1)*numGammas)+(1:numThetas)`);
                IF (class < numClasses) THEN DO;
                    gammaIndex[,class] = ((class*numThetas+(class-1)*numGammas)+(1:numGammas)`);
                END; 
            END;
        FINISH DoParameterIndexing;
        START GetCovarianceMatrices;
            penalty = RoughnessPenalty*PenaltyMatrix;
            I1 = J(numParams,numParams,0);
            I2 = J(numParams,numParams,0);
            I3 = J(numParams,numParams,0);
            /* Calculate matrix I1 */
            DO c = 1 TO numClasses;
                DO i = 1 TO numSubjects;
                   xi = X[LOC(IntId=i),];
                   ni = SUM(IntId=i);
                   inv_y_covmat_i = Inv_Y_CorrMats_Stacked[LOC(intId=i),1:ni]/SigSqTotal[c];  
                   contrib = postProbsBySubject[i,c]*xi`*(inv_y_covmat_i)*xi;
                   I1[thetaIndex[,c],thetaIndex[,c]] = I1[thetaIndex[,c],thetaIndex[,c]] + contrib;
                END; 
                I1[thetaIndex[,c],thetaIndex[,c]] = I1[thetaIndex[,c],thetaIndex[,c]] + penalty; 
            END;
            IF (numClasses>1) THEN DO;
                DO c = 1 TO (numClasses-1);
                    DO k = 1 TO (numClasses-1);
                        DO i = 1 TO numSubjects;
                           nonrefc = nonrefClasses[c];
                           nonrefk = nonrefClasses[k];
                           Si = S[i,];
                           contrib = fittedProb[i,nonrefc]*(1*(nonrefc=nonrefk)-fittedProb[i,nonrefk])* (Si`*Si);
                           I1[gammaIndex[,c],gammaIndex[,k]] = I1[gammaIndex[,c],gammaIndex[,k]] + contrib;
                        END;
                    END;
                END;
            END;
            /* Calculate matrix I2 */
            DO c = 1 TO numClasses;
                DO k = 1 TO numClasses;
                    DO i = 1 TO numSubjects;
                        ni = SUM(intId=i);
                        xi = X[LOC(intId=i),];
                        yi = Y[LOC(intId=i)];
                        Si = S[i,];
                        muic = fittedValues[LOC(intId=i),c];
                        inv_y_covmat_i = Inv_Y_CorrMats_Stacked[LOC(intId=i),1:ni]/SigSqTotal[c]; 
                        hic = (xi`)*inv_y_covmat_i*(yi-muic);    
                        DO j = 1 TO numSubjects;
                            * Information loss in linear regression coefficients:;
                            nj = SUM(intId=j);
                            xj = X[LOC(intId=j),];
                            yj = Y[LOC(intId=j)];
                            Sj = S[j,];
                            mujk = fittedValues[LOC(intId=j),k];
                            inv_y_covmat_j = Inv_Y_CorrMats_Stacked[LOC(intId=j),1:nj]/SigSqTotal[k];
                            hjk = (xj`)*inv_y_covmat_j*(yj-mujk);    
                            Gckij = postProbsBySubject[i,c]*postProbsBySubject[j,k]*(i^=j) +
                                     postProbsBySubject[i,c]*(i=j)*(c=k);
                            I2[thetaIndex[,c],thetaIndex[,k]] = I2[thetaIndex[,c],thetaIndex[,k]] + 
                                     Gckij * hic*hjk` - 
                                     (1/(numSubjects**2))*postProbsBySubject[i,c]*hic*(theta[,k]`)*penalty -
                                     (1/(numSubjects**2))*postProbsBySubject[j,k]*hjk*(theta[,c]`)*penalty +
                                     (1/(numSubjects**2))*penalty*theta[,c]*(theta[,k])`*penalty;
                            * Information loss in logistic regression coefficients:;
                            IF ((c < numClasses)&(k < numClasses)&(numClasses>1)) THEN DO;
                               nonrefc = nonrefClasses[c];
                               nonrefk = nonrefClasses[k]; 
                               Gcnonrefkij = postProbsBySubject[i,c]*postProbsBySubject[j,nonrefk]*(i^=j) +
                                   postProbsBySubject[i,c]*(i=j)*(c=nonrefk);
                               Gnonrefcnonrefkij = postProbsBySubject[i,nonrefc]*postProbsBySubject[j,nonrefk]*(i^=j) +
                                   postProbsBySubject[i,nonrefc]*(i=j)*(nonrefc=nonrefk);
                               contrib = (Gnonrefcnonrefkij-
                                   postProbsBySubject[i,nonrefc]*fittedProb[j,nonrefk]-
                                   postProbsBySubject[j,nonrefk]*fittedProb[i,nonrefc] +
                                   fittedProb[i,nonrefc]*fittedProb[j,nonrefk])* (Si`*Sj);
                               I2[gammaIndex[,c],gammaIndex[,k]] = I2[gammaIndex[,c],gammaIndex[,k]] + contrib;
                            END;
                            * Information loss in covariance between linear and;
                            * logistic regression coefficients:;
                            IF (k < numClasses) THEN DO; 
                                DO kprime = 1 TO numClasses;
                                    nonrefk = nonrefClasses[k];
                                    Gcnonrefkij = postProbsBySubject[i,c]*postProbsBySubject[j,nonrefk]+
                                          postProbsBySubject[i,c]*(i=j)*
                                          ((c=nonrefk)-postProbsBySubject[j,nonrefk]);
                                    contrib = Gcnonrefkij*((k=kprime)-fittedProb[j,k])*hic*(Sj) +
                                              ( (1/numSubjects)*postProbsBySubject[j,kprime] *
                                                ((k=kprime)-fittedProb[j,k])*penalty*theta[,c]*(Sj) );
                                    I2[thetaIndex[,c],gammaIndex[,k]] =
                                          I2[thetaIndex[,c],gammaIndex[,k]] +
                                         contrib;
                                    I2[gammaIndex[,k],thetaIndex[,c]] =
                                           (I2[thetaIndex[,c],gammaIndex[,k]])`;
                                END;
                            END;
                        END;
                    END;
                 END;
             END;
            /* Calculate matrix I3 */
            DO i = 1 TO numSubjects;
                scorei = J(numParams,1,0);
                DO c = 1 TO numClasses;
                    ni = SUM(intId=i);
                    xi = X[LOC(intId=i),];
                    yi = Y[LOC(intId=i)];
                    Si = S[i,];
                    muic = fittedValues[LOC(intId=i),c];
                    inv_y_covmat_i = Inv_Y_CorrMats_Stacked[LOC(intId=i),1:ni]/SigSqTotal[c]; 
                    hic = (xi`)*inv_y_covmat_i*(yi-muic);    
                    scorei[thetaIndex[,c]] = postProbsBySubject[i,c]*hic-
                       (1/numSubjects)*penalty*theta[,c];
                    IF c < numClasses THEN DO;
                        nonrefc = nonrefClasses[c];
                        scorei[gammaIndex[,c]] = (postProbsBySubject[i,nonrefc]-fittedProb[i,nonrefc])*Si;
                    END;
                END;
                I3 = I3 + scorei*scorei`;
            END;
            /* Calculate model-based covariance matrix (in beta, gamma and sigma) */
            naiveCovarianceMatrix = INV(I1-I2);
            /* Calculate sandwich covariance matrix (in beta, gamma and sigma) */
            covarianceMatrix = naiveCovarianceMatrix*I3*naiveCovarianceMatrix;
            IF (EIGVAL(covarianceMatrix)[><] < (-.000001)) THEN DO;
                PRINT("Warning - The covariance matrix is not positive definite");
                PRINT("Minimum eigenvalue:");
                PRINT(EIGVAL(covarianceMatrix)[><]);
                MacroError = 1;
            END;
        FINISH GetCovarianceMatrices;
        START PrintOutput;
            FILE PRINT;
            PRINT("MixTVEM Macro");
            PUT "Time variable name:            &time";
            PUT "Response variable name:        &dep";
            PUT "Number of subjects:            " NumSubjects;
            PUT "Total number of observations:  " NumTotal;
            IF NumCov>0 THEN PUT  "Constant-Coefficient Covariates:     &Cov";
            IF NumTCov>0 THEN PUT "Varying-Coefficient Covariates:      &TCov";
            IF NumSCov>0 THEN PUT "Subject-Level Trajectory Covariates: &SCov";
            IF Deg=1 THEN PUT "Effect of time between knots treated as: linear";
            IF Deg=2 THEN PUT "Effect of time between knots treated as: quadratic";
            IF Deg=3 THEN PUT "Effect of time between knots treated as: cubic";
            IF MacroError>0 THEN PUT "<<<< Macro error reported! >>>>";
            IF SmallestPrevalence<.025 THEN PUT "Warning: The smallest fitted class is very small";
            IF Converged=0 THEN PUT "Warning: EM algorithm did not converge";
            IF NROW(LogLik)=1 THEN  PUT "Log-likelihood:                       " LogLik;
            IF NROW(RSS)=1 THEN PUT "Residual squared error (RSS) weighted by posterior probability:  " RSS;
            IF NROW(GCV)=1 THEN     PUT "GCV:                                  " GCV;
            IF NROW(AIC)=1 THEN     PUT "AIC:                                  " AIC;
            IF NROW(BIC)=1 THEN     PUT "BIC:                                  " BIC;
            IF NROW(BestProportionNugget)=1 THEN PUT "Estimated proportion nugget: " BestProportionNugget;
            IF NROW(BestRho)=1 THEN PUT "Estimated autocorrelation parameter rho: " BestRho;
            IF NumStarts>1 THEN DO;
                PUT "Proportion of starting values giving approximately the best obtained log-likelihood:         "
                       PropBestLL;
            END;
            PUT "Count number of parameters:     " countNumParams;
            PUT "Smoothed number of parameters:  " effectiveNumParams;
        FINISH PrintOutput;
        /*****************************************************************/
        /* The following is the main body of the IML part of the macro   */
        /*****************************************************************/
        /******** Make beginning preparations. *****************/ 
        USE MixTVEMSettings;
            READ ALL VAR { Convergence 
                           MacroError MaxIterations
                           MaxVarianceRatio NumClasses NumStarts
                           NumSCov NumTCov NumSubjects
                           NumTotal Penalty_Max RefClass ScaleForGraph StdErrOption };
        CLOSE MixTVEMSettings;
        CALL PrepareRegressionMatrix;
        X = Regressors;
        USE MixTVEMData;
            READ ALL VAR{&dep} INTO Y;
            READ ALL VAR{&time} INTO Time;
            READ ALL VAR{MixTVEMIntID} INTO IntID;
        CLOSE MixTVEMData;
        USE MixTVEMRandomSeeds;
            READ ALL VAR{Seeds} INTO Seeds;
            READ ALL VAR{RSS} INTO RSSBySeed; /* will be all missing values */
            READ ALL VAR{GCV} INTO GCVBySeed; /* will be all missing values */
            READ ALL VAR{AIC} INTO AICBySeed; /* will be all missing values */
            READ ALL VAR{BIC} INTO BICBySeed; /* will be all missing values */
            READ ALL VAR{LogLik} INTO LogLikBySeed; /* will be all missing values */
            READ ALL VAR{SmallestPrevalence} INTO SmallestPrevalenceBySeed; /* will be all missing values */
            READ ALL VAR{LargestPrevalence} INTO LargestPrevalenceBySeed; /* will be all missing values */
            READ ALL VAR{Estim_Rho} INTO Estim_Rho_BySeed; /* will be all missing values */
            READ ALL VAR{Estim_ProportionNugget} INTO Estim_ProportionNugget_BySeed; /* will be all missing values */
        CLOSE MixTVEMRandomSeeds;
        USE MixTVEMWhichBeta;
            READ ALL INTO WhichBetaVector;
        CLOSE MixTVEMWhichBeta;
        USE MixTVEMWhichThetaWithinBeta;
            READ ALL INTO WhichThetaWithinBetaVector;
        CLOSE MixTVEMWhichThetaWithinBeta;
        Lowest_Num_Observations = .;
        Highest_Num_Observations = .;
        DO i = 1 TO IntID[<>];
            ni = SUM(IntID=i);
            IF (Lowest_Num_Observations=.) THEN Lowest_Num_Observations = ni;
            IF (Highest_Num_Observations=.) THEN Highest_Num_Observations = ni;
            IF (ni < Lowest_Num_Observations) THEN Lowest_Num_Observations = ni;
            IF (ni > Highest_Num_Observations) THEN Highest_Num_Observations = ni;        
        END; 
        EDIT MixTVEMSettings VAR{Lowest_Num_Observations Highest_Num_Observations}; 
            Lowest_Num_Observations = Lowest_Num_Observations; 
            Highest_Num_Observations = Highest_Num_Observations;
        REPLACE; CLOSE MixTVEMSettings; 
        /******** Loop over possible seeds to try to find the best seed. *****************/  
        MaxIterationsToUse = MaxIterations/2;  /* saves time for initial searching among seeds */ 
        DO SeedIndex = 1 TO NROW(Seeds); 
            /* We initially assume a very strong penalty in order to get an 
               initial parametric fit. We also initially assume uncorrelated 
               data because we have not yet estimated the correlation parameters. */
            Working_ProportionNugget = 0.5;
            Working_Rho = 0; 
            CALL CreateWeightingMatrices;            
            sdY = SQRT(((Y - Y[:])##2)[+]/(NROW(Y)-1));
            RoughnessPenalty=1e9*sdY;
            /* Do the initial fit of the data. */
            CALL RANDSEED(Seeds[SeedIndex],1);
            CALL FitTheModel;
            /* Now estimate rho and proportion nugget.*/
            crossprods = .;
            lags = .;
            weights = .; 
            DO k = 1 TO numClasses;
                DO i = 1 TO numSubjects;
                   rows_i = (LOC(intId=i)`);
                   ni = SUM(intId=i); 
                   IF (ni>1) THEN DO;
                       these_lags = J(ni*(ni-1)/2,1,.);
                       these_crossprods = J(ni*(ni-1)/2,1,.);
                       these_weights = J(ni*(ni-1)/2,1,.);
                       temp_start = 0;
                       time_i = Time[rows_i,]; 
                       resids_ik = Y[rows_i,] - fittedY[rows_i,k]; 
                       DO j = 1 TO (ni-1);
                            these_lags[(temp_start+1):(temp_start+ni-j)] = ABS(time_i[(j+1):ni]-time_i[j]); 
                            these_crossprods[(temp_start+1):(temp_start+ni-j)] = (resids_ik[j]*resids_ik[(j+1):ni])/SigSqTotal[k];
                            these_weights[(temp_start+1):(temp_start+ni-j)] = postProbsBySubject[i,k];
                            temp_start = temp_start + (ni-j);
                       END;
                       lags = lags // these_lags;
                       crossprods = crossprods // these_crossprods;
                       weights = weights // these_weights;
                   END;
                END;
            END;
            lags = lags[2:NROW(lags)]; * drop initial empty placeholder;
            crossprods = crossprods[2:NROW(crossprods)]; * drop initial empty placeholder;
            weights = weights[2:NROW(weights)]; * drop initial empty placeholder;
            CREATE MixTvemNLSData VAR {lags crossprods weights};
                APPEND;
            CLOSE MixTvemNLSData;
            /* Now use the residuals from the model that was fitted under independence 
               to estimate the dependence parameters.*/
            average_crossprod = ((crossprods#weights)[+])/(weights[+]); 
            average_lag = ((lags#weights)[+])/(weights[+]); 
            initial_logitnonnugget = 0;
            initial_logrho = (LOG(average_crossprod)-LOG(.5))/average_lag; 
            initialvals = J(2,1,.);
            initialvals[1] = initial_logitnonnugget;
            initialvals[2] = initial_logrho; 
            CALL NLPNRR(code,nls_result,"nls_function",initialvals);
            IF nls_result[1] < 30 THEN DO;
                Estim_ProportionNugget = 1-EXP(nls_result[1])/(1+EXP(nls_result[1]));
            END; ELSE DO;
                Estim_ProportionNugget = 1;
            END;
            Estim_Rho = EXP(nls_result[2]);
            Working_ProportionNugget = Estim_ProportionNugget;
            Working_Rho = Estim_Rho;  
            CALL CreateWeightingMatrices;            
            /* Now fit the data with the estimated correlation structure. */
            CALL InitialRandomEStep;
            CALL RANDSEED(Seeds[SeedIndex],1);
            CALL FitTheModel;
            RSSBySeed[SeedIndex] = RSS;
            GCVBySeed[SeedIndex] = GCV;
            AICBySeed[SeedIndex] = AIC;
            BICBySeed[SeedIndex] = BIC;
            logLikBySeed[SeedIndex] = logLik;
            smallestPrevalenceBySeed[SeedIndex] = ClassPrevalences[,><];
            largestPrevalenceBySeed[SeedIndex] = ClassPrevalences[,<>];
            Estim_ProportionNugget_BySeed[SeedIndex] = Estim_ProportionNugget;
            Estim_Rho_BySeed[SeedIndex] = Estim_Rho;
        END;
        BestSeed = Seeds[LogLikBySeed[<:>]]; 
        BestRho = Estim_Rho_BySeed[LogLikBySeed[<:>]]; 
        BestProportionNugget = Estim_ProportionNugget_BySeed[LogLikBySeed[<:>]]; 
        PropBestLL = (LogLikBySeed>=LogLikBySeed[<>]-.1)[+]/NROW(LogLikBySeed); 
        EDIT MixTVEMRandomSeeds VAR{Seeds
                                    RSS GCV AIC BIC logLik
                                    SmallestPrevalence
                                    LargestPrevalence 
                                    Estim_ProportionNugget Estim_Rho};
            Seeds=Seeds; 
            RSS=RSSBySeed;
            GCV=GCVBySeed;
            AIC=AICBySeed;
            BIC=BICBySeed;
            LogLik = LogLikBySeed;
            SmallestPrevalence=SmallestPrevalenceBySeed;
            LargestPrevalence=LargestPrevalenceBySeed;
            Estim_ProportionNugget = Estim_ProportionNugget_BySeed;
            Estim_Rho = Estim_Rho_BySeed;
            REPLACE ALL;
        CLOSE MixTVEMRandomSeeds;
        USE MixTVEMPenaltyMatrix;
            READ ALL INTO PenaltyMatrix;
        CLOSE MixTVEMPenaltyMatrix;
        Working_ProportionNugget = BestProportionNugget;
        Working_Rho = BestRho;  
        MaxIterationsToUse = MaxIterations;
        CALL CreateWeightingMatrices;
        /****** Found the best seed.  Now find the best penalty weight (lambda) .*********/
        /* First pass of search */
        sdY = SQRT(((Y-Y[:])[##])/(numTotal-1));
        penVector = (10##((-3):(FLOOR(LOG10(numTotal)))))`*sdY;
        enrpVector = J(NROW(penVector),1,.); 
        bicVector = J(NROW(penVector),1,.);
        DO vectorIndex = 1 TO NROW(bicVector);
            RoughnessPenalty = penVector[vectorIndex];
            CALL RANDSEED(BestSeed,1);
            CALL FitTheModel; 
            penVector[vectorIndex] = RoughnessPenalty; 
            enrpVector[vectorIndex] = EffectiveNumParams; 
            bicVector[vectorIndex] = BIC;
        END;
        IF (bicVector[>:<] = 1) THEN DO;
            PRINT "Warning: Penalty_Max may need to be made smaller.";
        END;
        IF (bicVector[>:<] = NROW(bicVector)) THEN DO;
            PRINT "Warning: Penalty_Max may need to be made larger.";
        END;
        bestLambdaSoFar = penVector[bicVector[>:<]]; 
         /* Second pass of search */
        penVector = ((10**((0:4)*.3-.6)`))*bestLambdaSoFar;
        enrpVector = J(NROW(penVector),1,.);
        bicVector = J(NROW(penVector),1,.);
        DO vectorIndex = 1 TO NROW(bicVector);
            RoughnessPenalty = penVector[vectorIndex];
            CALL RANDSEED(BestSeed,1);
            CALL FitTheModel;
            penVector[vectorIndex] = RoughnessPenalty;
            enrpVector[vectorIndex] = EffectiveNumParams;
            bicVector[vectorIndex] = BIC;
        END;
        *PRINT penVector bicVector;
        bestLambdaSoFar = penVector[bicVector[>:<]];
        /* Third pass of search */
        penVector = ((10**((0:4)*.1-.2)`))*bestLambdaSoFar;
        enrpVector = J(NROW(penVector),1,.);
        bicVector = J(NROW(penVector),1,.);
        DO vectorIndex = 1 TO NROW(bicVector);
            RoughnessPenalty = penVector[vectorIndex];
            CALL RANDSEED(BestSeed,1);
            CALL FitTheModel;
            penVector[vectorIndex] = RoughnessPenalty;
            enrpVector[vectorIndex] = EffectiveNumParams;
            bicVector[vectorIndex] = BIC;
        END;
        *PRINT penVector bicVector;
        bestLambdaSoFar = penVector[bicVector[>:<]];
        /****** Found the best lambda.  Now fit the final model. *********/
        RoughnessPenalty = bestLambdaSoFar;
		PRINT "Using lambda = " RoughnessPenalty;
        CALL RANDSEED(BestSeed,1);
        MaxIterationsToUse = MaxIterations;
        CALL FitTheModel; 
        /****** Now do follow-up calculations, such as standard errors. ******/        
        CALL DoParameterIndexing;
        CALL SaveFittedValues;
        GridTime = Time[><]+(Time[<>]-Time[><])*(0:(ScaleForGraph-1))`/(ScaleForGraph-1);
        IF StdErrOption = 1 THEN DO;
            CALL GetCovarianceMatrices;
            CALL SaveFittedValuesSEs;
        END;
        DO thisTCov = 1 TO NumTCov;
            CALL SaveFittedBetasAndSEs;
            CALL SaveGridFittedBetasAndSEs;
        END;
        colnames = COMPRESS(CONCAT("Class",CHAR(1:NumClasses)));
        CREATE MixTVEMPi FROM ClassPrevalences[COLNAME=colnames];
            APPEND FROM ClassPrevalences;
        CLOSE MixTVEMPi;
        USE MixTVEMWhichBeta; READ ALL INTO WhichBetaVector; CLOSE MixTVEMWhichBeta;
        USE MixTVEMWhichThetaWithinBeta;
            READ ALL INTO WhichThetaWithinBetaVector;
        CLOSE MixTVEMWhichThetaWithinBeta;
        OutClass = COLVEC((1:NumClasses)`@J(NROW(Theta),1,1));
        WhichBeta = COLVEC(WhichBetaVector`@J(NumClasses,1,1));
        WhichThetaWithinBeta = COLVEC(WhichThetaWithinBetaVector`@J(NumClasses,1,1));
        OutTheta = OutClass || WhichBeta || WhichThetaWithinBeta || COLVEC(Theta`);
        ThetaColNames = "Class" || "WhichBeta" || "WhichThetaWithinBeta" || "Theta";
        CREATE MixTVEMTheta FROM OutTheta[COLNAME=ThetaColNames];
            APPEND FROM OutTheta;
        CLOSE MixTVEMTheta;        
        OutClass = COLVEC((1:NumClasses)`@J((NumSCov+1),1,1));
        OutSNum = COLVEC(J(NumClasses,1,1)@(1:(NumSCov+1))`);
        OutGamma = OutClass || OutSNum || COLVEC(Gamma`);
        GammaColNames = "Class" || "SNum" || "Gamma";
        CREATE MixTVEMGamma FROM OutGamma[COLNAME=GammaColNames];
            APPEND FROM OutGamma;
        CLOSE MixTVEMGamma;
        colnames = COMPRESS(CONCAT("Class",CHAR(1:NumClasses)));
        CREATE MixTVEMSigSqTotal FROM SigsqTotal[COLNAME=colnames];
            APPEND FROM SigsqTotal;
        CLOSE MixTVEMSigSqTotal;
        IF StdErrOption>0 THEN DO;
            CREATE MixTVEMMatrixI1 FROM I1;
                APPEND FROM I1;
            CLOSE MixTVEMMatrixI1;
            CREATE MixTVEMMatrixI2 FROM I2;
                APPEND FROM I2;
            CLOSE MixTVEMMatrixI2;
            CREATE MixTVEMMatrixI3 FROM I3;
                APPEND FROM I3;
            CLOSE MixTVEMMatrixI3;
            CREATE MixTVEMNaiveCovMatrix FROM naiveCovarianceMatrix;
                APPEND FROM naiveCovarianceMatrix;
            CLOSE MixTVEMNaiveCovMatrix;
            CREATE MixTVEMCovMatrix FROM CovarianceMatrix;
                APPEND FROM CovarianceMatrix;
            CLOSE MixTVEMCovMatrix;
            CREATE MixTVEMThetaIndex FROM ThetaIndex;
                APPEND FROM ThetaIndex;
            CLOSE MixTVEMThetaIndex;
            SETheta = J(NROW(ThetaIndex),numClasses,0);
            IF NumClasses > 1 THEN DO;
                CREATE MixTVEMGammaIndex FROM GammaIndex;
                    APPEND FROM GammaIndex;
                CLOSE MixTVEMGammaIndex;
                SEGamma = J(NROW(GammaIndex),numClasses,0);
            END;
            colnames = COMPRESS(CONCAT("Class",CHAR(1:NumClasses)));
            ParameterStdErrs = SQRT(VECDIAG(CovarianceMatrix));
            IF (numClasses>1) THEN DO;
                NonRefClasses = (1:numClasses)[LOC((1:numClasses)^=refClass)];
            END;
            DO k = 1 TO numClasses;
                SETheta[,k] = ParameterStdErrs[ThetaIndex[,k]];
                IF (k<numClasses) THEN DO;
                    SEGamma[,NonRefClasses[k]] = ParameterStdErrs[GammaIndex[,k]];
                END;
            END;
            OutClass = COLVEC((1:NumClasses)`@J((NumSCov+1),1,1));
            OutSNum = COLVEC(J(NumClasses,1,1)@(1:(NumSCov+1))`);
            IF (NumClasses > 1) THEN DO;
                OutSEGamma = OutClass || OutSNum || COLVEC(SEGamma`);
                GammaColNames = "Class" || "SNum" || "SE_Gamma";
                CREATE MixTVEMStdErrGamma FROM OutSEGamma[COLNAME=GammaColNames];
                    APPEND FROM OutSEGamma;
                CLOSE MixTVEMStdErrGamma;
            END;
            OutClass = COLVEC((1:NumClasses)`@J(NROW(Theta),1,1));
            OutThetaNum = COLVEC(J(NumClasses,1,1)@(1: NROW(Theta))`);
            WhichBeta = COLVEC(WhichBetaVector`@J(NumClasses,1,1));
            WhichThetaWithinBeta = COLVEC(WhichThetaWithinBetaVector`@J(NumClasses,1,1));
            OutSETheta = OutClass || WhichBeta || WhichThetaWithinBeta || COLVEC(SETheta`);
            ThetaColNames = "Class" || "WhichBeta" || "WhichThetaWithinBeta" || "SE_Theta";
            CREATE MixTVEMStdErrTheta FROM OutSETheta[COLNAME=ThetaColNames];
                APPEND FROM OutSETheta;
            CLOSE MixTVEMStdErrTheta;
        END;
        temp = ((1:numSubjects)`) || postProbsBySubject;
        CREATE MixTVEMPostProbs FROM temp
           [COLNAME=("MixTVEMIntID" || COMPRESS(CONCAT("PostProbClass",CHAR(1:NumClasses))))];
            APPEND FROM temp;
        CLOSE MixTVEMPostProbs;
        EDIT MixTVEMSettings VAR{Iterations MacroError RoughnessPenalty};
            MacroError = MacroError;
            Iterations = Iteration;
            RoughnessPenalty = RoughnessPenalty;
            REPLACE;
        CLOSE MixTVEMSettings;
        deg = &deg;  
        AdjustedN = (NumTotal+2)/24;
        Cov = "' &Cov '";
        SCov = "' &SCov '";
        TCov = "' &TCov '";
        Knots = "' &Knots '";
        CREATE MixTVEMFitStatistics VAR{RSS GCV 
                                        AIC BIC
                                        logLik
                                        iteration
                                        Converged
                                        PropBestLL
                                        numSubjects
                                        numTotal
                                        deg 
                                        SmallestPrevalence
                                        LargestPrevalence
                                        RoughnessPenalty
                                        Knots Cov SCov TCov
                                        numClasses BestProportionNugget BestRho
                                        effectiveNumParams
                                        countNumParams};
            APPEND;
        CLOSE MixTVEMFitStatistics;
        RUN;
        CALL PrintOutput;
        RUN;
   QUIT;
   PROC SORT DATA=MixTVEMGamma; BY SNum Class; RUN;
        DATA MixTVEMGamma;
           MERGE MixTVEMGamma MixTVEMSCovNames;
           BY SNum;
        RUN;
   PROC SORT DATA=MixTVEMGamma; BY Class SNum; RUN;
   PROC IML;
        USE MixTVEMPostProbs;
            READ ALL INTO PostProbs;
        CLOSE MixTVEMPostProbs;
        MixTVEMIntID = PostProbs[,1];
        PostProbs = PostProbs[,2:NCOL(PostProbs)];
        MaxIndex = PostProbs[,<:>];
        Output = MixTVEMIntID || MaxIndex;
        ColNames = "MixTVEMIntID" || "BestAssign";
        CREATE MixTVEMHardAssign FROM Output [COLNAME=ColNames];
            APPEND FROM Output;
        CLOSE MixTVEMHardAssign;
        RUN;
   QUIT;
   DATA MixTVEMIDKey;
       SET MixTVEMData;
       KEEP MixTVEMIntID &id;
   RUN;
   PROC SORT DATA=MixTVEMIDKey NODUPKEY;
       BY MixTVEMIntID;
   RUN;
   PROC SORT DATA=MixTVEMPostProbs;
       BY MixTVEMIntID;
   RUN;
   PROC SORT DATA=MixTVEMHardAssign;
       BY MixTVEMIntID;
   RUN;
   DATA MixTVEMPostProbs; MERGE MixTVEMIDKey MixTVEMPostProbs MixTVEMHardAssign;
       BY MixTVEMIntID;
   RUN;
   PROC DATASETS NOLIST NOWARN; DELETE MixTVEMIDKey MixTVEMHardAssign; RUN;
   PROC SORT DATA=MixTVEMRandomSeeds;
       BY LogLik;
   RUN; 
   PROC SORT DATA=MixTVEMFittedValues;
       BY MixTVEMIntID &Time;
   RUN;
   PROC SORT DATA=MixTVEMFittedValueSEs;
       BY MixTVEMIntID &Time;
   RUN;
   %IF %EVAL(&MixTVEMStdErrOption>0) %THEN %DO;
       DATA MixTVEMFittedValues;
           MERGE MixTVEMFittedValues MixTVEMFittedValueSEs;
           BY MixTVEMIntID &Time;
       RUN;
       DATA MixTVEMFittedValues;
           SET MixTVEMFittedValues;
           %DO MixTVEMThisClass = 1 %TO &Latent_Classes;
                    Lower_Class&MixTVEMThisClass =
                            Class&MixTVEMThisClass +
                            QUANTILE('NORMAL',(1-((1-&Coverage)*.5)))*SE_Class&MixTVEMThisClass;
                    Upper_Class&MixTVEMThisClass =
                            Class&MixTVEMThisClass +
                            QUANTILE('NORMAL',(1-((1-&Coverage)*.5)))*SE_Class&MixTVEMThisClass;
           %END;
       RUN;
       %IF %EVAL(&MixTVEMNumTCov>0) %THEN %DO;
           %DO MixTVEMTemp = 1 %TO &MixTVEMNumTCov;
               DATA MixTVEMFittedBeta&MixTVEMTemp;
                   MERGE MixTVEMFittedBeta&MixTVEMTemp MixTVEMFittedBetaSE&MixTVEMTemp;
               RUN;
               DATA MixTVEMFittedBeta&MixTVEMTemp;
                   SET MixTVEMFittedBeta&MixTVEMTemp;
                   %DO MixTVEMThisClass = 1 %TO &Latent_Classes;
                            Lower_Class&MixTVEMThisClass =
                                    Class&MixTVEMThisClass -
                                    QUANTILE('NORMAL',(1-((1-&Coverage)*.5)))*SE_Class&MixTVEMThisClass;
                            Upper_Class&MixTVEMThisClass =
                                    Class&MixTVEMThisClass +
                                    QUANTILE('NORMAL',(1-((1-&Coverage)*.5)))*SE_Class&MixTVEMThisClass;
                   %END;
               RUN;
               PROC SORT DATA=MixTVEMGridFittedBetaSE&MixTVEMTemp; BY &Time; RUN;
               DATA MixTVEMGridFittedBeta&MixTVEMTemp;
                   MERGE MixTVEMGridFittedBeta&MixTVEMTemp MixTVEMGridFittedBetaSE&MixTVEMTemp;
                   BY &Time;
               RUN;
               DATA MixTVEMGridFittedBeta&MixTVEMTemp;
                   SET MixTVEMGridFittedBeta&MixTVEMTemp;
                   %DO MixTVEMThisClass = 1 %TO &Latent_Classes;
                            Lower_Class&MixTVEMThisClass =
                                    Class&MixTVEMThisClass -
                                    QUANTILE('NORMAL',(1-((1-&Coverage)*.5)))*SE_Class&MixTVEMThisClass;
                            Upper_Class&MixTVEMThisClass =
                                    Class&MixTVEMThisClass +
                                    QUANTILE('NORMAL',(1-((1-&Coverage)*.5)))*SE_Class&MixTVEMThisClass;
                   %END;
               RUN;
           %END;
       %END;
   %END;
   %IF %EVAL(&MixTVEMStdErrOption=1) %THEN %DO;
        %IF %EVAL(&MixTVEMNumClasses>1) %THEN %DO;
            DATA MixTVEMGamma;
               MERGE MixTVEMGamma MixTVEMStdErrGamma;
               BY Class SNum;
               z = .;
               IF (ABS(Gamma)>0 & ABS(SE_Gamma) > 0) THEN z = Gamma/SE_Gamma;
               IF (z ^= .) THEN p = 2*(1-PROBNORM(ABS(z)));
            RUN;
            PROC DATASETS NOLIST NOWARN;
               DELETE MixTVEMStdErrGamma
                      MixTVEMSCovNames;
            RUN;
       %END;
       DATA MixTVEMTheta;
          MERGE MixTVEMTheta MixTVEMStdErrTheta;
          BY Class WhichBeta WhichThetaWithinBeta;
             z = Theta/SE_Theta;
          p = 2*(1-PROBNORM(ABS(z)));
       RUN;
    %END; %ELSE %DO;
        PROC SORT DATA=MixTVEMGamma; BY SNum Class; RUN;
        DATA MixTVEMGamma;
           MERGE MixTVEMGamma MixTVEMSCovNames;
           BY SNum;
        RUN;
        PROC SORT DATA=MixTVEMGamma; BY Class SNum; RUN;
       * PROC DATASETS NOLIST NOWARN;
       *    DELETE MixTVEMSCovNames;
       * RUN;
    %END;
    DATA MixTVEMGamma;
        SET MixTVEMGamma;
        DROP SNum;
    RUN;
    %IF %EVAL(&MixTVEMNumTCov>0) %THEN %DO;
        PROC SORT DATA=MixTVEMTheta;
            BY WhichBeta;
        RUN;
        PROC IML;
            BetaNamesSpaced = "&MixTVEMTCov";
            USE MixTVEMSettings; READ VAR {NumTCov}; RUN;
            Name = J((NumTCov),1,"                                ");
            DO CovIndex = 1 TO (NumTCov);
                Name[CovIndex] = SCANQ(BetaNamesSpaced,CovIndex);
            END;
            WhichBeta = (1:(NumTCov))`;
            CREATE MixTVEMNames VAR {WhichBeta Name}; APPEND;
            CLOSE MixTVEMNames;
        RUN;
        DATA MixTVEMTheta;
           MERGE MixTVEMNames MixTVEMTheta;
           BY WhichBeta;
        RUN;
    %END;
    PROC SORT DATA=MixTVEMTheta;
        BY Class WhichBeta WhichThetaWithinBeta;
    RUN;
    %IF %EVAL(&MixTVEMNumCov>0) %THEN %DO;
        PROC IML;
            CovNamesSpaced = "&MixTVEMCov";
            USE MixTVEMSettings; READ VAR {NumCov}; CLOSE MixTVEMSettings;
            VarName = J((NumCov),1,"                                ");
            DO i = 1 TO NumCov;
                VarName[i] = SCANQ(CovNamesSpaced,i);
            END;
            WhichThetaWithinBeta = (1:NumCov)`;
            CREATE MixTVEMNonTVCovNames VAR {WhichThetaWithinBeta VarName};
                APPEND;
            CLOSE MixTVEMNonTVCovNames;
        QUIT;
        DATA MixTVEMNonTimeVarying;
            SET MixTVEMTheta;
            WHERE WhichBeta IS MISSING;
        RUN;
        PROC SORT DATA=MixTVEMNonTimeVarying;
            BY WhichThetaWithinBeta;
        RUN;
        DATA MixTVEMNonTimeVarying;
            MERGE MixTVEMNonTVCovNames MixTVEMNonTimeVarying ;
            BY WhichThetaWithinBeta;
            DROP Name WhichBeta;
        RUN;
        DATA MixTVEMNonTimeVarying;
            SET MixTVEMNonTimeVarying;
            Beta = Theta;
            SE_Beta = SE_Theta;
            DROP Theta SE_Theta;
        RUN;
        PROC SORT DATA=MixTVEMNonTimeVarying;
            BY Class WhichThetaWithinBeta;
        RUN;
        %IF %EVAL(&MixTVEMStdErrOption=1) %THEN %DO;
            PROC PRINT DATA=MixTVEMNonTimeVarying NOOBS;
               TITLE2 "Constant Regression Effects";
               VAR Class VarName Beta SE_Beta z p;
            RUN;
            TITLE2;
        %END; %ELSE %DO;
            PROC PRINT DATA=MixTVEMNonTimeVarying NOOBS;
               TITLE2 "Constant Regression Effects";
               VAR Class VarName Beta ;
            RUN;
            TITLE2;
        %END;
    %END;
    PROC PRINT DATA=MixTVEMPi NOOBS;
       TITLE2 "Class Proportions";
    RUN;
    TITLE2;
    DATA MixTVEMSigmas;
        SET MixTVEMSigSqTotal;
        ARRAY sigmasquared(&Latent_Classes) Class1 - Class&Latent_Classes;
        ARRAY sigma(&Latent_Classes) Sigma_Class1 - Sigma_Class&Latent_Classes;
        DO i = 1 TO &Latent_Classes;
            sigma(i)=SQRT(sigmasquared(i));
          END;
        KEEP Sigma_Class1 - Sigma_Class&Latent_Classes;
    RUN;
    PROC PRINT DATA=MixTVEMSigmas NOOBS;
       TITLE2 "Standard Deviations";
    RUN;
    TITLE2;
    %IF %EVAL((&MixTVEMNumClasses>1)&(&MixTVEMStdErrOption=1)) %THEN %DO;
        PROC PRINT DATA=MixTVEMGamma NOOBS;
            TITLE2 "Logistic Regression for Class Membership";
            VAR Class SCovName Gamma SE_Gamma z p;
        RUN;
        TITLE2;
    %END; %ELSE %DO;
        PROC PRINT DATA=MixTVEMGamma NOOBS;
            TITLE2 "Logistic Regression for Class Membership";
            VAR Class SCovName Gamma ;
        RUN;
        TITLE2;
    %END;
/*    PROC DATASETS NOLIST NOWARN;
        DELETE  FittedValueses
                NaiveCovMatrix
                Names
                Names
                NonTVCovNames
                StdErrTheta
                WhichBeta
                WhichThetaWithinBeta
                %DO MixTVEMThisVarIndex = 1 %TO &MixTVEMNumTCov;
                    FittedBetaSE&MixTVEMThisVarIndex
                    GridFittedBetaSE&MixTVEMThisVarIndex
               %END;
    RUN;QUIT;*/
%MEND;
