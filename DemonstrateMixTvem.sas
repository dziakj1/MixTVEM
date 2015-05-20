%INCLUDE "C:\Users\jjd264\Documents\Sims-MixTvem\Tvem_Mix_Normal.sas";

PROC IMPORT OUT=SubjectLevelData 
            DATAFILE= "C:\Users\jjd264\Documents\Sims-MixTvem\MixTvemSampleSubjectLevel.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

PROC IMPORT OUT=ObservationLevelData 
            DATAFILE= "C:\Users\jjd264\Documents\Sims-MixTvem\MixTvemSampleObservationLevel.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
 
DATA Combined; 
	MERGE SubjectLevelData ObservationLevelData;
	BY ID;
RUN;

PROC MEANS DATA=Combined;
RUN;

DATA Combined;
	SET Combined;
	CenteredNA = NegAffect - 1.5688092;
	Intercept = 1;
RUN;
 
%TVEM_Mix_Normal(mydata = Combined,
	                time = Time,
	                dep = urge,
	                id = ID ,             
	                tcov = Intercept  , 
	                latent_classes = 2, 
	                Knots = 6 , 
	                Num_Starts = 5,
	                Std_Err_Option  = yes );

SYMBOL1 COLOR="BLUE" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL2 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL3 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL4 COLOR="RED" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL5 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL6 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
PROC GPLOT DATA=MixTVEMGridFittedBeta1;
    PLOT Class1*Time Upper_Class1*Time Lower_Class1*Time 
         Class2*Time Upper_Class2*Time Lower_Class2*Time  
	   / OVERLAY;
RUN;QUIT; 
SYMBOL1; SYMBOL2; SYMBOL3; SYMBOL4; SYMBOL5; SYMBOL6;
   
 
   
%TVEM_Mix_Normal(mydata = Combined,
	                time = Time,
	                dep = urge,
	                id = ID ,             
	                tcov = Intercept  ,
                    cov = CenteredNA, 
	                latent_classes = 2, 
	                Knots = 6 , 
	                Num_Starts = 5,
	                Std_Err_Option  = yes );
SYMBOL1 COLOR="BLUE" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL2 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL3 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL4 COLOR="RED" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL5 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL6 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
PROC GPLOT DATA=MixTVEMGridFittedBeta1;
    PLOT Class1*Time Upper_Class1*Time Lower_Class1*Time 
         Class2*Time Upper_Class2*Time Lower_Class2*Time  
	   / OVERLAY;
RUN;QUIT; 
   
   

  
%TVEM_Mix_Normal(mydata = Combined,
	                time = Time,
	                dep = urge,
	                id = ID ,
	                deg = 3,                  
	                tcov = Intercept CenteredNA, 
	                latent_classes = 2, 
	                Knots = 6 6 , 
	                Num_Starts = 5,
	                Std_Err_Option  = yes );
SYMBOL1 COLOR="BLUE" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL2 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL3 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL4 COLOR="RED" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL5 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL6 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
PROC GPLOT DATA=MixTVEMGridFittedBeta1;
    PLOT Class1*Time Upper_Class1*Time Lower_Class1*Time 
         Class2*Time Upper_Class2*Time Lower_Class2*Time  
	   / OVERLAY;
RUN;QUIT; 
PROC GPLOT DATA=MixTVEMGridFittedBeta2;
    PLOT Class1*Time Upper_Class1*Time Lower_Class1*Time 
         Class2*Time Upper_Class2*Time Lower_Class2*Time  
	   / OVERLAY;
RUN;QUIT; 

 
%TVEM_Mix_Normal(mydata = Combined,
	                time = Time,
	                dep = urge,
	                id = ID ,   
	                tcov = Intercept CenteredNA, 
	                latent_classes = 2, 
					scov = BaselineCigarettesPerDay  MinutesToFirstCigarette, 
	                Knots = 6 6 ,  
	                Num_Starts = 5,
	                Std_Err_Option  = yes );
SYMBOL1 COLOR="BLUE" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL2 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL3 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL4 COLOR="RED" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL5 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL6 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
PROC GPLOT DATA=MixTVEMGridFittedBeta1;
    PLOT Class1*Time Upper_Class1*Time Lower_Class1*Time 
         Class2*Time Upper_Class2*Time Lower_Class2*Time  
	   / OVERLAY;
RUN;QUIT; 
PROC GPLOT DATA=MixTVEMGridFittedBeta2;
    PLOT Class1*Time Upper_Class1*Time Lower_Class1*Time 
         Class2*Time Upper_Class2*Time Lower_Class2*Time  
	   / OVERLAY;
RUN;QUIT; 
   
   

%TVEM_Mix_Normal(mydata = Combined,
	                time = Time,
	                dep = urge,
	                id = ID ,   
	                tcov = Intercept CenteredNA, 
	                latent_classes = 2, 
					scov = RelapseAtOneMonth, 
	                Knots = 6 6 ,  
	                Num_Starts = 5,
	                Std_Err_Option  = yes );
SYMBOL1 COLOR="BLUE" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL2 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL3 COLOR="BLUE" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL4 COLOR="RED" INTERPOL=JOIN LINE=1 VALUE=NONE WIDTH=2; 
SYMBOL5 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
SYMBOL6 COLOR="RED" INTERPOL=JOIN LINE=2 VALUE=NONE WIDTH=2; 
PROC GPLOT DATA=MixTVEMGridFittedBeta1;
    PLOT Class1*Time Upper_Class1*Time Lower_Class1*Time 
         Class2*Time Upper_Class2*Time Lower_Class2*Time  
	   / OVERLAY;
RUN;QUIT; 
PROC GPLOT DATA=MixTVEMGridFittedBeta2;
    PLOT Class1*Time Upper_Class1*Time Lower_Class1*Time 
         Class2*Time Upper_Class2*Time Lower_Class2*Time  
	   / OVERLAY;
RUN;QUIT; 
   
   
