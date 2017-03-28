
/*
================================================================================================
 
 OVERVIEW:

 This programme performs Orbital Simimilarity analysis of meteor showers using one of three
 selected D-criterion algorithms and compares the resulting stream identification with the 
 identity assigned in initial processing by UFO ORBIT.  The three D -criterion algorithms
 imple,emted in this code are:
 
 - Drummond (DD, see Reference 1)
 - Southworth and Hawkins (Dsh, see Reference 2)
 - Jopek (DH, see Reference 3)
 
 orbital elements of observed meteors ARE COMPARED against those of known meteor showers 
 (reference orbits).  Observations are then assigned to Stream based on a best fit to the 
 reference orbits.  
 
 Best fit is defined as the lowest D criterion value when comparing streams to the 
 reference orbits.  However, D criterion values greater than a threashold value are
 not assigned to a shower; these are candidate sporadic (non shower) meteors.

 The process steps are as follows:
 
 1) Read reference orbits into a table
 2) Read observations 
 3) Remove observed orbits which do not meet EDMOND quality criteria
 4) For each reference orbit:
 	a) 	Calculate D criterion value for each observation against the reference orbit.
 	b) 	If D value is less than threshold AND D value is lower than a previous D value
 		then assign observation to shower (i.e. find best fit)
 
 NOTE: This is an iterative algorithm - it can take a while to process large datasets
 
 INPUTS (files): 	

 	- UFO Orbit UNIFIED observations in csv format 
 	  (see http://sonotaco.com/soft/UO2/UO21Manual_EN.pdf)
 	- Reference stream orbital elements in csv format
 
 OUTPUTS (SAS Dataset): 

 	- Correlation summary (WORK.CORRELATION)
 
 REFERENCES:
 
 1) Drummond, J. D., 1981, Icarus 45, 545
 2) Southworth, R. B., Hawkins, G. S., 1963, Smithson. Contrib. Astrophys. 7, 261
 3) Jopek, T. J., 1993, Icarus 106, 603
 
 AUTHOR:

 Peter Campbell-Burns, UK Meteor Observation Network (UKMON) 
 
 VERSION CONTROL:
 
 Date	     Version		Notes
 -------------------------------------------------------
 26/03/2017  1.0			First release
 
 LICENSE:
 
 This work is licensed under Creative Commons Attribution-ShareAlike 4.0 International License.
 
================================================================================================
*/

/* -------------------------------------------------------------------------------------
   CONFIGURATION
   ------------------------------------------------------------------------------------- */

/* Select algorithm and set threshold value got orbital similarity */
%let THRESHOLD    		= 0.8;   	/* D Criterion threshold value  */   
%let D_CRITERION  		= DD;		/* D Criterion algorithm */
	
/* Raw Data */
%let Reference_Orbits 	= /folders/myfolders/shower_list.csv;
%let Observation_Data 	= /folders/myfolders/unified.csv;

/* Quality Criterion (see Sonotaco UFO Orbit user manual for full definitions) */
%let QA_dv12     		= 10.0;	/*Percentage velocity difference */
%let QA_GM      		= -100;	/* Percentage overlap of two observed trajectories */
%let QA_Dur      		= 0.1;	/* Event duration (seconds) */
%let QA_QA       		= 0.15;	/* Total quality assessment */
%let QA_Qo       		= 1.0;	/* Minimum observed trajectory angle (degrees) */
%let QA_Qc       		= 10.0;	/* Cross angle of two observed planes */
%let QA_Delta_GP 		= 0.5;	/* Difference of pole of two round tragectories */
%let QA_H1      		= 200;	/* Altitude of start of ablation (in km) */
%let QA_H2       		= 20;	/* Altitude of end of ablation   (in km) */

/* Define Pi */
%let PI = 3.141592653589793;

/* Define macro to convert degrees to radians */
%macro radians(x);
	&x. * &PI/ 180
%mend;

/* -------------------------------------------------------------------------------------
   Import reference orbit data  
   Build lookup macro variables from orbit reference data in table (STREAM.LIS
   ------------------------------------------------------------------------------------- */
FILENAME REFFILE "&Reference_Orbits";
PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.STREAM_LIST REPLACE;
	GETNAMES=YES;
RUN;

DATA _NULL_;
	SET WORK.STREAM_LIST;
	CALL SYMPUT("shower_name_" || strip(put(_N_,3.)),shower_name);
	CALL SYMPUT("_peri_"       || strip(put(_N_,3.)), put(%radians(_peri),best16.));
	CALL SYMPUT("_node_"       || strip(put(_N_,3.)), put(%radians(_node),best16.));
	CALL SYMPUT("_incl_"       || strip(put(_N_,3.)), put(%radians(_incl),best16.));
	CALL SYMPUT("_e_"          || strip(put(_N_,3.)), put(_e,best16.));
	CALL SYMPUT("_q_"          || strip(put(_N_,3.)), put(_q,best16.));
	CALL SYMPUT("Candidates",put(_N_,8.)); /* Count of the number of rows */
	IF d_sh_threshold NE . 
	 	THEN CALL SYMPUT("_D_THRESHOLD_" || strip(put(_N_,3.)),d_sh_threshold);
	 	ELSE CALL SYMPUT("_D_THRESHOLD_" || strip(put(_N_,3.)),&Threshold.);
RUN;


/* -------------------------------------------------------------------------------------
	Import raw observations:    
	- Select only unified observations 
	- Apply QA filter and split observations into passed QA / failed QA  
	- Transform _STREAM to 3 UC characters, set initial stream assignment to SPO (Sporadic)
	  and convert angle to radians
   ------------------------------------------------------------------------------------- */

FILENAME OBSFILE"&Observation_Data";
PROC IMPORT DATAFILE=OBSFILE
	DBMS=CSV
	OUT=WORK.RAW_OBS REPLACE;
	GETNAMES=YES;
RUN;	

/* Select only the the unified observations and apply QA filter */

DATA WORK.FILTERED_OBS WORK.FILTERED_OBS_QA_FAIL ;
	
	DROP QA_PASS QA_FAIL;
	SET WORK.RAW_OBS END = EOF; 
	
	/* Select only UNIFIED observations, exclude single statin observations */
	WHERE substr(_ID1,1,9) = "_(UNIFIED";
	
	IF  
	 _QA >= &QA_QA.
	 and abs(_dv12_ <= &QA_dv12.) 
	 and abs(_Gm_ >= &QA_GM.) 
	 and _dur >= &QA_dur.
	 and _H1 <= &QA_H1.
	 and _H2 >= &QA_H2.
	 and _Qo >= &QA_Qo.
	 and _QC >= &QA_QC.
	THEN DO /* Passes QA */;
		QA_PASS +1;
		OUTPUT FILTERED_OBS;
		END;
	ELSE DO /* Fails QA */;
		QA_FAIL +1;
		OUTPUT FILTERED_OBS_QA_FAIL;
		END;
	
	/* If last record then output summary stats */
	IF EOF THEN DO;
		PUTLOG "QA Filter Criteria:"
		     / "-------------------"
		     / "_QA    >= &QA_QA.   "
		     / "_dv12_ <= &QA_dv12. "
		     / "_Gm_   >= &QA_GM.   "
		     / "_dur   >= &QA_dur.  "
		     / "_H1    <= &QA_H1.   "
		     / "_H2    >= &QA_H2.   "
		     / "_Qo    >= &QA_Qo.   "
		     / "_QC    >= &QA_QC.   "
		     / ". "	
		     / "QA Filter Results:"
		     / "------------------"				
		     / "Accepted" QA_PASS  COMMA10.
		     / "Rejected" QA_FAIL  COMMA10.
		     / ". ";	
		END;
RUN;

/* Transformation steps */

DATA WORK.TRANSFORMED_OBS (RENAME = (__ = ID));
	SET WORK.FILTERED_OBS (KEEP   = __ _LOCALTIME _NODE _PERI _INCL _E _Q _STREAM ) ;
	FORMAT 	MATCH_D 12.4;
	LENGTH 	Stream_D $50.;	
	
	/* add coumns and set initial match result as NO MATCH */
	RETAIN 	MATCH_D STREAM_D;
	IF _N_ = 1 THEN DO;
		MATCH_D  = .;
		STREAM_D = "SPO";
		END;

	/* Convert angles from degrees to radians */
	_node = %radians(_node);
	_peri = %radians(_peri);
	_incl = %radians(_incl);

	/* Refformat _Stream column */
	IF UPCASE(_STREAM) = "_SPO" THEN _STREAM = "SPO";
	ELSE IF UPCASE(SUBSTR(_STREAM,1,2)) = "_J" THEN _STREAM = UPCASE(SUBSTR(_STREAM,5,3));
RUN;



/* -------------------------------------------------------------------------------------
   Define D analysis macros
   ------------------------------------------------------------------------------------- */

/* D analysis calculation */

%macro D_Calc(D_TYPE=DD);

%*Algorithms for DD, DSH and DH below follow references (1), (2) and(3) respectively;
	
	I21 = arcos(cos(_i1) * cos(_i2) + sin(_i1) * sin(_i2) * cos(_n1 - _n2) );
	ADIFF = ABS(_n1-_n2);
	IF ADIFF > &PI. THEN ADIFF = ABS(ADIFF - 2 * &PI.);
	IF ADIFF <=&PI.
		THEN II21 = _p2 - _p1 + 2 * arsin(cos((_i2 + _i1)/2) * sin((_n2 + _n1)/2) * sec(I21/2));
	    ELSE II21 = _p2 - _p1 - 2 * arsin(cos((_i2 + _i1)/2) * sin((_n2 + _n1)/2) * sec(I21/2));

	/* Drummond */
	%if &D_TYPE EQ DD %THEN
		%DO;
			DROP THETA B1 B2 G1 G2;
			B1 = arsin(sin(_i1) * sin(_p1));
			B2 = arsin(sin(_i2) * sin(_p2));
			G1 = _n1 + atan(cos(_i1) * tan(_p1));
			G2 = _n2 + atan(cos(_i2) * tan(_p2));
			IF cos(_P1) LT 0 THEN G1 = G1 + &PI;
			IF cos(_P2) LT 0 THEN G2 = G2 + &PI;
			Theta = arcos(sin(B1) * sin(B2) + cos(B1) * cos(B2) * cos(G2-G1) );
			D_Result = sqrt(((_q - _q2)/(_q +_q2))**2 + ((_e - _e2)/(_e + _e2))**2 + (I21 
				/ &PI)**2 +((_e2+_e)/2)**2 * (Theta/&PI)**2);
		%end;

	/* Southworth and Hawkns */
	%if &D_TYPE EQ DSH %THEN
		%DO;
			D_Result = sqrt ((_q - _q2)**2 + (_e - _e2)**2 + (2 * 
				sin(I21/2))**2 + ((_e - _e2)/2 * 2 * sin(II21/2) )**2);
		%end;
		
	/* Jopek */
	%if &D_TYPE EQ DH %THEN
		%DO;
			D_Result = sqrt (((_q - _q2)/(_q + _q2))**2 + (_e - _e2)**2 + (2 * 
				sin(I21/2))**2 + ((_e - _e2)/2 * 2 * sin(II21/2) )**2);
		%end;
%mend D_calc;

/* tterative application of D analysis to input data set  */

%macro D_Similarity (dset = );

OPTIONS NONOTES; /* Keeps log to reasonable size, no rows will be added or dropped */

%DO Reference_row = 1 %to &Candidates.;
	
	DATA &dset.;
	
		/* Ensure that ref_orbits are not reinitialised (to null) */
		RETAIN _i2 _q2 _e2 _p2 _n2;
		
		/* Specify input data set */
		SET &dset. ;
		
		/* Prevent intermediate values being written to output data set */
		DROP	_i1 _p1 _n1 _i2 _q2 _e2 _p2 _n2 I21 ADIFF II21 D_RESULT;
		
		/* If first row of input dataset, set the reference orbit variables */
		IF _N_ = 1 THEN DO
			_i2 = &&_incl_&Reference_row.; 
			_q2 = &&_q_&Reference_row.;
			_e2 = &&_e_&Reference_row.;
			_p2 = &&_peri_&Reference_row.;
			_n2 = &&_node_&Reference_row.;
			END;
				
			/* Convert all angles to RADIANS (input variables) */
			_n1 = _node;
			_i1 = _incl;
			_p1 = _peri;
			
			/* Calculate D value */
			%D_Calc(D_TYPE=&D_CRITERION);
	
			/* Assign new stream IF we are below threshold AND D_NEW is below the 
			   D value for at the last stream assignment (i.e. can we better the match) */
			IF D_RESULT LE &&_D_THRESHOLD_&Reference_row. AND (D_RESULT LT MATCH_D OR MATCH_D = .) THEN DO;
				Stream_D = "&&shower_name_&Reference_row.";
				MATCH_D = D_RESULT ;
				END;  
	RUN;
	
%END;	

OPTIONS NOTES;

%mend D_Similarity;

/* -------------------------------------------------------------------------------------
   Perform Orbital Similarity analysis
   ------------------------------------------------------------------------------------- */

%D_Similarity(dset = WORK.TRANSFORMED_OBS);



/* -------------------------------------------------------------------------------------
   Summarise correlation between UFO ORBIT and D CRITERION Correlation
   Note nested case: "NO MATCH" is equivalent to a sporadic meteor ("SPO")
   ------------------------------------------------------------------------------------- */
PROC SQL NOPRINT;
	CREATE TABLE WORK.CORRELATION AS
	SELECT 	_STREAM, 
	        COUNT(*) AS TOTAL_OBS, 
			SUM(CORRELATION) AS CORRELATED_STREAMS, 
			CALCULATED CORRELATED_STREAMS / CALCULATED TOTAL_OBS AS CORRELATION_RATE  FORMAT PERCENT6.0
		FROM (	SELECT _STREAM,
				CASE WHEN _STREAM = Stream_D 
				     THEN 1 
				     ELSE 0 
				     END AS CORRELATION
			FROM WORK.TRANSFORMED_OBS)
	GROUP BY _STREAM
	ORDER BY CORRELATION_RATE DESCENDING;
QUIT;

TITLE "Correlation of stream identification between UFO Orbit and Orbital Similarity analysis";
PROC PRINT DATA=WORK.CORRELATION label;
	LABEL 	 _STREAM = "UFO Identification"
			 TOTAL_OBS = "Total Meteor Observations"
			 CORRELATED_STREAMS= "&D_CRITERION. criterion matches"
			 CORRELATION_RATE = "Correlation Rate";
RUN;

/* -------------------------------------------------------------------------------------
   Diagnostic plots
   ------------------------------------------------------------------------------------- */

%let diag = PER;
ods noproctitle;

TITLE "Detailed analysis of correlation between UFO Orbit and Orbital Similarity analysis for Perseids (PER)";
PROC FREQ DATA = WORK.TRANSFORMED_OBS;
	WHERE _STREAM = "&diag" or STREAM_D = "&diag";
	TABLE _STREAM*STREAM_D / nopct nocol norow;
	RUN;

TITLE "Frequency distribution of &D_CRITERION. values for Perseid meteors";
PROC UNIVARIATE DATA = WORK.TRANSFORMED_OBS;
	WHERE STREAM_D = "&diag";
	ODS SELECT Histogram cdfplot;
	VAR MATCH_D;
	HISTOGRAM MATCH_D / VSCALE = COUNT ENDPOINTS = 0.0 to 0.3 by 0.005;
	CDFPLOT  MATCH_D ;
RUN;

TITLE;


