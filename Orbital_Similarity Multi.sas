
/*
******************************************************************************************
 This programme performs Orbital Simimilarity analysis of meteor showers using one of three
 D-criterion algorithms:

 - Drummond (DD)
 - Southworth and Hawkins (Dsh)
 - Jopek (DH)
 
 Orbital elements of observed meteors ARE COMPARED against those of known meteor showers 
 (reference orbits) and the D-criterion calculated.  
 
 D-values up to THRESHOLD are then plotted as a histogram and CDF.  Only those D-Values
 are plotted where the solar longitude of the observation falls within the start and end
 solar longitudes of the reference orbit.
 
 INPUT DATASET: CSV file containing UFO Orbit UNIFIED observations
                CSV file containing list of reference orbits
                
 OUTPUT DATASET: WORK.D_OBSERVATIONS
 
******************************************************************************************
*/

/* -------------------------------------------------------------------------------------
   CONFIGURATION
   ------------------------------------------------------------------------------------- */

/* Select algorithm and set threshold value got orbital similarity */
%let THRESHOLD    = 0.6;   	/* D Criterion threshold value  */   
%let D_CRITERION  = DD;		/* D Criterion algorithm */
%let BIN_SIZE     = 0.01;	/* Plot bin size */

/* Input Data */
%let Reference_Orbits = /folders/myfolders/j8.csv;
%let Observation_Data = /folders/myfolders/Edmond05_v03.csv;

/* Quality Criterion (see Sonotaco UFO Orbit user manual for full definitions) */
%let QA_dv12     = 10.0;
%let QA_GM       = 80.0;
%let QA_Dur      = 0.1;		/* Event duration (seconds) */
%let QA_QA       = 0.15;
%let QA_Qo       = 1.0;
%let QA_Qc       = 10.0;
%let QA_Delta_GP = 0.5;
%let QA_H1       = 200;		/* Altitude of start of ablation (in km) */
%let QA_H2       = 20;		/* Altitude of end of ablation   (in km) */


/* -------------------------------------------------------------------------------------
   Import reference orbit data  
   ------------------------------------------------------------------------------------- */
FILENAME REFFILE "&Reference_Orbits";
PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.STREAM_LIST REPLACE;
	GETNAMES=YES;
RUN;


/* -------------------------------------------------------------------------------------
   Import raw observations    
   ------------------------------------------------------------------------------------- */
FILENAME OBSFILE"&Observation_Data";
PROC IMPORT DATAFILE=OBSFILE
	DBMS=CSV
	OUT=WORK.RAW_OBS REPLACE;
	GETNAMES=YES;
RUN;	

/* -------------------------------------------------------------------------------------
   FILTER:
   Select only the the unified observations
   Apply QA filter and split observations into two datasets - passed QA and failed QA  
   ------------------------------------------------------------------------------------- */
DATA WORK.FILTERED_OBS WORK.FILTERED_OBS_QA_FAIL ;
	
	DROP QA_PASS QA_FAIL;

	SET WORK.RAW_OBS END = EOF; 
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
	THEN DO;
		QA_PASS +1;
		OUTPUT FILTERED_OBS ;
		END;
	ELSE DO;
		QA_FAIL +1;
		OUTPUT FILTERED_OBS_QA_FAIL;
		END;
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


/* -------------------------------------------------------------------------------------
   Transform: 
   Reformat _STREAM to three upper case characters
   Add initial stream assignments (NO MATCH)
   ------------------------------------------------------------------------------------- */
DATA WORK.TRANSFORMED_OBS (RENAME = (__ = ID));
	SET WORK.FILTERED_OBS (KEEP   = __ _LOCALTIME _NODE _PERI _INCL _E _Q _STREAM _SOL ) ;
	
	/* Refformat _Stream column */
	IF UPCASE(_STREAM) = "_SPO" THEN _STREAM = "SPO";
	ELSE IF UPCASE(SUBSTR(_STREAM,1,2)) = "_J" THEN _STREAM = UPCASE(SUBSTR(_STREAM,5,3));
RUN;

/* -------------------------------------------------------------------------------------
   Build lookup macro variables from orbit reference data in table (STREAM.LIST)
   ------------------------------------------------------------------------------------- */
DATA _NULL_;
	SET WORK.STREAM_LIST;
	CALL SYMPUT("shower_name_" || strip(put(_N_,3.)),_name);
	CALL SYMPUT("_sol1_"       || strip(put(_N_,3.)),_sol1);
	CALL SYMPUT("_sol2_"       || strip(put(_N_,3.)),_sol2);
	CALL SYMPUT("_peri_"       || strip(put(_N_,3.)),_peri);
	CALL SYMPUT("_node_"       || strip(put(_N_,3.)),_node);
	CALL SYMPUT("_incl_"       || strip(put(_N_,3.)),_incl);
	CALL SYMPUT("_e_"          || strip(put(_N_,3.)),_e);
	CALL SYMPUT("_q_"          || strip(put(_N_,3.)),_q);
	CALL SYMPUT("Candidates",_N_); /* Count of the number of rows */
	IF d_sh_threshold NE . 
	 	THEN CALL SYMPUT("_D_THRESHOLD_" || strip(put(_N_,3.)),d_sh_threshold);
	 	ELSE CALL SYMPUT("_D_THRESHOLD_" || strip(put(_N_,3.)),&Threshold.);
RUN;

/* Define Pi */
%let PI = 3.141592653589793;

/* Define macro to convert degrees to radians */
%macro radians(x);
	&x. * &PI/ 180
%mend;

/* -------------------------------------------------------------------------------------
   Macro to perform D_Criterion calculation 
   ------------------------------------------------------------------------------------- */

%macro D_Calc(D_TYPE=DD);
	I21 = arcos(cos(_i1) * cos(_i2) + sin(_i1) * sin(_i2) * cos(_n1 - _n2) );;
	IF ABS(_n1-_n2) <=&PI.
		THEN II21 = _p2 - _p1 + 2 * arsin(cos((_i2 + _i1)/2) * sin((_n2 + _n1)/2) * sec(I21/2));
	    ELSE II21 = _p2 - _p1 - 2 * arsin(cos((_i2 + _i1)/2) * sin((_n2 + _n1)/2) * sec(I21/2));

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
			D_Result = sqrt(((_q2 - _q)/(_q2 +_q))**2 + ((_e2 - _e)/(_e2 + _e))**2 + (I21 
				/ &PI)**2 +((_e2+_e)/2)**2 * (Theta/&PI)**2);
		%end;

	%if &D_TYPE EQ DSH %THEN
		%DO;
			D_Result = sqrt ((_q - _q2)**2 + (_e - _e2)**2 + (2 * 
				sin(I21/2))**2 + ((_e - _e2)/2 * 2 * sin(II21/2) )**2);
		%end;

	%if &D_TYPE EQ DH %THEN
		%DO;
			D_Result = sqrt (((_q - _q2)/(_q + _q2))**2 + (_e - _e2)**2 + (2 * 
				sin(I21/2))**2 + ((_e - _e2)/2 * 2 * sin(II21/2) )**2);
		%end;
%mend D_calc;


%macro D_Similarity (dset = );

%DO Reference_row = 1 %to &Candidates.;

	DATA &dset._out;
	
		/* Ensure that ref_orbits are not reinitialised (to null) */
		RETAIN _i2 _q2 _e2 _p2 _n2;
		
		/* Specify input data set */
		SET &dset. ;
		WHERE _SOL >= &&_sol1_&Reference_row. AND _SOL <= &&_sol2_&Reference_row.;
		/* Prevent intermediate values being written to output data set */
		DROP	_i1 _p1 _n1 _i2 _q2 _e2 _p2 _n2 I21 ADIFF II21 ;
		
		/* If first row of input dataset, set the reference orbit variables */
		IF _N_ = 1 THEN DO
			_i2 = %radians(&&_incl_&Reference_row.); 
			_q2    = &&_q_&Reference_row.;
			_e2    = &&_e_&Reference_row.;
			_p2 = %radians(&&_peri_&Reference_row.);
			_n2 = %radians(&&_node_&Reference_row.);
			END;
							
		/* Convert all angles to RADIANS (input variables) */
			_n1 = %radians(_node);
			_i1 = %radians(_incl);
			_p1 = %radians(_peri);
			
			%D_Calc(D_TYPE=&D_CRITERION);
			
			&D_CRITERION = D_RESULT;
	
			run;
			

	ods noproctitle;
	
	TITLE "Frequency distribution of &D_CRITERION. values for &&shower_name_&Reference_row. meteors";
	PROC UNIVARIATE DATA = &dset._out;
		WHERE &D_CRITERION LE &THRESHOLD.;
		ODS SELECT Histogram cdfplot;
		VAR &D_CRITERION ;
		HISTOGRAM &D_CRITERION / GRID VSCALE = COUNT ENDPOINTS = 0.0 to &THRESHOLD. by &BIN_SIZE;
		CDFPLOT   &D_CRITERION ;
	RUN;

TITLE;
			
	
%END;	

%mend D_Similarity;

/* -------------------------------------------------------------------------------------
   Perform Orbital Similarity analysis
   ------------------------------------------------------------------------------------- */

%D_Similarity(dset = WORK.TRANSFORMED_OBS);

