
/*
================================================================================================
 
 ORBITAL SIMILARITY ANALYSIS (SINGLE SHOWER)
 
 OVERVIEW:

 This programme performs an Orbital Simimilarity analysis of a defined meteor shower using one 
 of three D-criterion algorithms and plots a frequency distribution and cumulative distribution 
 function for the calculated D values.  
 
 The three D-criterion algorithms implemented by thi code are:
 
 - Drummond (DD, see Reference 1)
 - Southworth and Hawkins (Dsh, see Reference 2)
 - Jopek (DH, see Reference 3)
  
 D criterion values greater than the threashold value are not assigned to a the shower; these are 
 candidate sporadic (non shower) meteors.
 
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
   Reference orbit data  
   ------------------------------------------------------------------------------------- */

%let SHOWER_NAME = Peseids;
%let _PERI 		 = 149.513875;		/* Argument of perihelion (in Degrees) */
%let _NODE		 = 139.194266; 		/* Longitude of the ascending node (in Degrees) */
%let _INCL		 = 113.038548;		/* Inclination (in Degrees) */
%let _E    		 = 0.912969;		/* Eccetricity */
%let _Q          = 0.946507;		/* Perihelion distance (in AU) */


/* -------------------------------------------------------------------------------------
   CONFIGURATION
   ------------------------------------------------------------------------------------- */

/* Select algorithm and set threshold value got orbital similarity */
%let THRESHOLD    		= 0.8;   	/* D Criterion threshold value  */   
%let D_CRITERION  		= DD;		/* D Criterion algorithm */
	
/* Raw Data file */
%let Observation_Data 	= /folders/myfolders/unified.csv; /* UFO Orbit output file (Unified) */

/* Quality Criterion (see Sonotaco UFO Orbit user manual for full definitions) */
%let QA_dv12     		= 10.0;	/*Percentage velocity difference */
%let QA_GM      		= -100;	/* Percentage overlap of two observed trajectories */
%let QA_Dur      		= 0.1;		/* Event duration (seconds) */
%let QA_QA       		= 0.15;	/* Total quality assessment */
%let QA_Qo       		= 1.0;		/* Minimum observed trajectory angle (degrees) */
%let QA_Qc       		= 10.0;	/* Cross angle of two observed planes */
%let QA_Delta_GP 		= 0.5;		/* Difference of pole of two round tragectories */
%let QA_H1      		= 200;		/* Altitude of start of ablation (in km) */
%let QA_H2       		= 20;		/* Altitude of end of ablation   (in km) */

/* Define Pi */
%let PI = 3.141592653589793;

/* Define macro to convert degrees to radians */
%macro radians(x);
	&x. * &PI/ 180
%mend;


/* -------------------------------------------------------------------------------------
	Import raw observations:    
	- Select only unified observations 
	- Apply QA filter and split observations into passed QA / failed QA  
	- Transform _STREAM to 3 UC characters, set initial stream assignment to SPO (Sporadic)
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

	/* Refformat _Stream column */
	IF UPCASE(_STREAM) = "_SPO" THEN _STREAM = "SPO";
	ELSE IF UPCASE(SUBSTR(_STREAM,1,2)) = "_J" THEN _STREAM = UPCASE(SUBSTR(_STREAM,5,3));
RUN;



/* -------------------------------------------------------------------------------------
   Define D analysis macros
   ------------------------------------------------------------------------------------- */
  
/* Define Pi */
%let PI = 3.141592653589793;

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

	
DATA WORK.ANALYSED_OBS;

	/* Ensure that ref_orbits are not reinitialised (to null) */
	RETAIN _i2 _q2 _e2 _p2 _n2;
	
	/* Specify input data set */
	SET WORK.TRANSFORMED_OBS ;
	
	/* Prevent intermediate values being written to output data set */
	DROP	_i1 _p1 _n1 _i2 _q2 _e2 _p2 _n2 I21 ADIFF II21 D_RESULT;
	
	/* If first row of input dataset, set the reference orbit variables */
	IF _N_ = 1 THEN DO
		_q2 = &_q.;
		_e2 = &_e.;
		_i2 = %radians(&_incl.); 
		_p2 = %radians(&_peri.);
		_n2 = %radians(&_node.);
		END;
			
		/* Convert all angles to RADIANS (input variables) */
		_n1 = %radians(_node);
		_i1 = %radians(_incl);
		_p1 = %radians(_peri);
		
		/* Calculate D value */
		%D_Calc(D_TYPE=&D_CRITERION);

		/* Assign new stream IF we are below threshold AND D_NEW is below the 
		   D value for at the last stream assignment (i.e. can we better the match) */
		IF D_RESULT LE &THRESHOLD. THEN DO;
			Stream_D = "&SHOWER_NAME.";
			MATCH_D = D_RESULT ;
			END;  
RUN;
		

/* -------------------------------------------------------------------------------------
   Plot Frequeny distribution and CDF
   ------------------------------------------------------------------------------------- */

ods noproctitle;

TITLE "Frequency distribution of &D_CRITERION. values for Perseid meteors";
PROC UNIVARIATE DATA = WORK.ANALYSED_OBS;
	WHERE STREAM_D = "&SHOWER_NAME";
	ODS SELECT Histogram cdfplot;
	VAR MATCH_D;
	HISTOGRAM MATCH_D / VSCALE = COUNT ENDPOINTS = 0.0 to &THRESHOLD. by 0.005;
	CDFPLOT  MATCH_D ;
RUN;

TITLE;


