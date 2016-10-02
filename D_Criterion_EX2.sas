/*----------------------------------------------------------------------------------------

This module performs an itterative analysis of a set of observations:

1) Ingest Orbits
2) Calculate a mean orbit
3) Compare orbits against mean orbit
4) If there are orbits that exceed D Criteria match, remove worst match 
5) Repeat (2-4) untill all orbits are within threshold
6) Compare rejected orbits with new mean Orbit, add back to the main set the best fit within threshold
7) Recalculate mean orbits
8) Repeat (6-7) until no more rejected orbits are within threshold.

NOTE: 	Input data set must cntain orbits that are thought to comprise a single stream.
	ONLY DSH and DH are calculated in this example.
 
 INPUT DATASET: CSV file (UFO Orbit UNIFIED observations in csv format)
 
 OUTPUT DATASET: WORK.D_CRITERION
 
------------------------------------------------------------------------------
*/


/*****************************************************************************************
CONFIGURATION
/*****************************************************************************************
;
/* Setup Threshold value and input filepath macro variables*/
%let THRESHOLD = 0.1;   /* D Criterion threshold value     */
%let FullFilePath = /folders/myshortcuts/SAS_Uni/myfolders/UNIFIED.csv;
%let MaxIter = 5;       /* Max iterations per pass         */

/*****************************************************************************************
 MACRO TILITY DEFINITION
/*****************************************************************************************
;

/* Macro to convert degrees to radians for trig functions */
%macro radians(x);
%global PI;
%let PI = 3.141592653589793;
&x. * &PI/ 180
%mend;

/* Macro to calculate mean orbits of a dataset and return orbital elements */
%macro get_mean_orbit (dset=,outtab=N);
	proc means data = &dset. 
	%IF (%UPCASE(&outtab.) NE Y) %THEN noprint;
	;
	output 	out=work.means  
		mean(_e) = mean_e
		mean(_peri) = mean_peri
		mean(_incl) = mean_incl
		mean(_node) = mean_node
		mean(_q) = mean_q;
	var _incl _node _peri _e _q;
	run;
	/* Put results into Macro variables */
	proc sql noprint;
		select mean_e, mean_peri, mean_node, mean_incl, mean_q
		into :x_e2, :x_peri2, :x_node2, :x_incl2, :x_q2
		from work.means;
	Quit;

/* Make as GLOBAL variables */
%global _e2 _peri2 _node2 _incl2 _q2;	
%let _e2     	= &x_e2.;
%let _peri2	= &x_peri2.;
%let _node2	= &x_node2.;
%let _incl2   	= &x_incl2.;
%let _q2 	= &x_q2.;
	
%mend;

/* Macro to count pass fail against criterion threshold */ 
%macro Filter_D0 (dset=, D0=0.2);
DATA _NULL_;
set &dset. end = EOF;
If D_SH >= &D0. THEN FailCount + 1;
ELSE PassCount + 1;
IF EOF THEN DO;
	If PassCount = . THEN PassCount = 0;
	If FailCount = . THEN FailCount = 0;
	CALL SYMPUT("x_D0_Pass",PassCount);
	CALL SYMPUT("x_D0_Fail",FailCount);
	END;
RUN;

%global D0_fail D0_Pass;
%let D0_pass = &x_D0_pass;
%let D0_Fail = &x_D0_Fail;
%mend;

/* Macro to calculate D_Criterion for each observation against mean orbit */
%macro calculate_D(Dset=);

	DATA &dset. ;
		
	/* The prevents the seed variables being reinitialised to null */
	RETAIN _incl2 _q2 _e2 _peri2 _node2;
	
	FORMAT D_SH 12.2;
	FORMAT D_H  12.2;
	
	/* Specify input data set */
	SET &dset. ;
	
	/* Prevent temorary variables being written to output data set */
	DROP   _incl1 _peri1 _node1  _incl2 _q2 _e2 _peri2 _node2 I II ;
	
	/* If first row of input dataset, set the seed orbit variables
	   from values held outside the data step in MACRO variables */
	IF _N_ = 1 THEN DO
		_incl2    = %radians(&_incl2); 
		_q2    = &_q2;
		_e2    = &_e2;
		_peri2 = %radians(&_peri2);
		_node2 = %radians(&_node2);
		END;
		
	/* Convert all angles to RADIANS (input variables) */
		_node1 = %radians(_node);
		_incl1 = %radians(_incl);
		_peri1 = %radians(_peri);
		
	/* Calculate I and II (Southworth & Hawkins) formula */		
		I = arcos(cos(_incl1) * cos(_incl2) + sin(_incl1) * sin(_incl2) * cos(_node1 - _node2) );
		ADIFF = ABS(_node1-_node2);
		IF ADIFF > &PI. THEN ADIFF = ABS(ADIFF - 2 * &PI.);
		IF ADIFF <= &PI.
		THEN II = _peri2 - _peri1 + 2 * arsin(cos( (_incl2 + _incl1)/2 ) * sin( (_node2 + _node1)/2 ) * sec(I/2));
		ELSE II = _peri2 - _peri1 - 2 * arsin(cos( (_incl2 + _incl1)/2 ) * sin( (_node2 + _node1)/2 ) * sec(I/2));
		
	/* Calculate the Southworth & Hawkins Dsh parameter */	
		D_SH = sqrt ( (_q - _q2)**2 + (_e - _e2)**2 + (2 * sin(I/2))**2 + ( (_e - _e2)/2 * 2 * sin(II /2) )**2 );
		D_H  = sqrt ( ((_q - _q2)/(_q + _q2) )**2 + (_e - _e2)**2 + (2 * sin(I/2))**2 + ( (_e - _e2)/2 * 2 * sin(II /2) )**2 );		
	run;

%mend;


/*****************************************************************************************
 IMPORT UfO Orbit Data 
/*****************************************************************************************
;
/* Import UFO Orbit data and create temporary dataset UA_IMPORT in the SAS work library */

FILENAME REFFILE "&FullFilePath";
PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.UA_IMPORT REPLACE;
	GETNAMES=YES;
RUN;

/* Select the unified observations; split Unified observations     */
/* into two datasets - passed QA criterion and failed QA criretion */                
DATA WORK.UNIFIED_QA_PASS WORK.UNIFIED_QA_FAILED;
	SET WORK.UA_IMPORT; 
	IF substr(_ID1,1,9) = "_(UNIFIED" THEN DO;
	   IF _e <= 1 and abs(_dv12_ <=10.0) and abs(_Gm_ > 80.0) THEN OUTPUT UNIFIED_QA_PASS;
	   ELSE OUTPUT UNIFIED_QA_FAILED;
	   END;
RUN;


/*****************************************************************************************
 MAIN PROGRAM
/*****************************************************************************************
;
/* Copy data to working table */
DATA WORK.D_CRITERION;
SET work.unified_qa_pass;
IF _stream = "_J5_Per"; /* Filter for test only */
RUN;

/* Main routine is run as a macro to allow conditional processing and loops */

%macro Iteration;
Options Nonotes;

/* Create empty work tables */

DATA WORK.ALL_REJECTS;
SET WORK.D_CRITERION;
STOP;
RUN;

DATA WORK.ACCEPTS;
SET WORK.D_CRITERION;
STOP;
RUN;

%let Itercount = 0;
	
/* Keep processing untill no more observations fail against criteria threshold */
%DO %UNTIL (&D0_Fail. = 0 or &IterCount. = &MaxIter.);

	%let Itercount = %eval(&Itercount + 1);

	%get_mean_orbit(dset=WORK.D_CRITERION );

	/* Calculate D criterion */
	%calculate_D(dset=WORK.D_CRITERION );

	/* Count pass fail against D Criterion threshold */
	%Filter_D0(dset=WORK.D_CRITERION, D0=&THRESHOLD.);
	%PUT Pass 1, Iter: &Itercount.. &D0_PASS Meteors below D_criterion threshold. ;

	/* If any rejects, remove the one with greatest variance and add to rejects*/
	%IF (&D0_Fail. NE 0) %THEN %DO;
	
		PROC SORT DATA = WORK.D_CRITERION;
		BY DESCENDING D_SH;
		RUN;
	
		DATA WORK.D_CRITERION WORK.REJECTS;
		SET WORK.D_CRITERION;
		IF _N_ NE 1 THEN OUTPUT WORK.D_CRITERION;
		ELSE OUTPUT WORK.REJECTS;
		RUN;
		
		PROC APPEND BASE=WORK.ALL_REJECTS DATA = WORK.REJECTS FORCE NOWARN;
		RUN;
		
	%END;
%END;

%if (&IterCount. = &MaxIter.) %then %put Pass 1 terminated early, max iterations reached (&MaxIter);
%let Itercount = 0;

/* Keep processing untill no more observations fail against criteria threshold */
%DO %UNTIL (&D0_PASS. = 0 or &IterCount. = &MaxIter.);

	%let Itercount = %eval(&Itercount + 1);

	/* Calculate mean orbital elements */
	%get_mean_orbit(dset=WORK.D_CRITERION );
	
	/* Calculate D criterion for rejects */
	%calculate_D(dset=WORK.ALL_REJECTS );

	/* Count passes and fails against threshold */
	%Filter_D0(dset=ALL_REJECTS, D0=&THRESHOLD.);
	%PUT Pass 2, Iter: &Itercount.. &D0_PASS Rejects below D_criterion threshold. ;
	
	/* If any passes, extract first and add back to dataset */
	%IF (&D0_PASS. NE 0) %THEN %DO;
	
		PROC SORT DATA = WORK.ALL_REJECTS;
		BY D_SH;
		RUN;
	
		DATA WORK.D_CRITERION WORK.REJECTS;
		SET WORK.D_CRITERION;
		IF _N_ EQ 1 THEN OUTPUT WORK.ACCEPTS;
		ELSE OUTPUT WORK.REJECTS;
		RUN;
		
		PROC APPEND BASE=WORK.D_CRITERION DATA = WORK.ACCEPTS FORCE NOWARN;
		RUN;
		
	%END;
		
%END; 

%if (&IterCount. = &MaxIter.) %then %put Pass 1 terminated early, max iterations reached (&MaxIter);

/* Display summary of results */
%get_mean_orbit(dset=WORK.D_CRITERION, Outtab=Y);

Options Notes;		
%mend;	
 
%Iteration;