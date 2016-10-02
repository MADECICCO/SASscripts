/*
-------------------------------------------------------------------------------------------
 This module will compare a set of observations against a set of reference orbits and will 
 assign observations to a Stream based on a best fit to the reference orbits.  Best fit is
 defined as the lowest D criterion value when comparing streamd to the reference orbits.
 
 1) 	Ingest and filter observations based on detailed Quality Criteria
 2)	Ingest reference orbits
 3) 	Compare all observations with each reference orbit.
 4)	Set Stream_D to reference stream if D value is less than threshold AND D value is 
 	lower than a previous D value (i.e. find best fit)
 
 NOTE:  Match can be pertformed usinf either DH, DD or DSH.
 	This is an itterative algorithm - it can take a while to process large datasets
 
 INPUT DATASET: CSV file (UFO Orbit UNIFIED observations in csv format)
 
 OUTPUT DATASET: WORK.D_OBSERVATIONS
 
------------------------------------------------------------------------------------------
*/

/*****************************************************************************************
CONFIGURATION
/*****************************************************************************************
;
/* Setup Threshold value and input filepath macro variables*/
%let THRESHOLD   = 0.1;   /* D Criterion threshold value     */
%let D_CRITERION = DD;

/* Quality Criterion */

%let QA_e = 1;
%let QA_dv12 = 10.0;
%let QA_GM = 80.0;
%let QA_Dur = 0.1;
%let QA_QA = 0.15;
%let QA_Qo = 1.0;
%let QA_Qc = 10.0;
%let QA_Delta_GP = 0.5;
%let QA_H1 = 200;
%let QA_H2 = 20;


%let FullFilePath_ref = /folders/myshortcuts/SAS_Uni/myfolders/shower_list.csv;
%let FullFilePath_obs = /folders/myshortcuts/SAS_Uni/myfolders/unified.csv;

/* Import observations */
FILENAME OBSFILE"&FullFilePath_obs";
PROC IMPORT DATAFILE=OBSFILE
	DBMS=CSV
	OUT=WORK.IMPORT REPLACE;
	GETNAMES=YES;
RUN;	

/* Select the unified observations; split Unified observations     */
/* into two datasets - passed QA criterion and failed QA criretion */  


DATA WORK.OBSERVATIONS WORK.FAIL_QA ;
	SET WORK.IMPORT END = EOF; 
	IF substr(_ID1,1,9) = "_(UNIFIED";
	IF  
/* 	_e <= &QA_e.  */
/* 	 and  */
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
		OUTPUT WORK.OBSERVATIONS ;
		END;
	ELSE DO;
		QA_FAIL +1;
		OUTPUT WORK.FAIL_QA;
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

/* Import reference orbits */
FILENAME REFFILE "&FullFilePath_ref";
PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.STREAM_LIST REPLACE;
	GETNAMES=YES;
RUN;

	
/* Turn reference orbit data into macro variables */
DATA _NULL_;
SET WORK.STREAM_LIST;
CALL SYMPUT("shower_name_" || strip(put(_N_,3.)),shower_name);
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

%let PI = 3.141592653589793;
%macro radians(x);
	&x. * &PI/ 180
%mend;


/* Macro to perform D_Criterion calculation */
%macro D_Calc(D_TYPE=DD);

	I21 = arcos(cos(_i1) * cos(_i2) + sin(_i1) * sin(_i2) * cos(_n1 - _n2) );
	ADIFF = ABS(_n1-_n2);
	IF ADIFF > &PI. THEN ADIFF = ABS(ADIFF - 2 * &PI.);
	IF ADIFF <= &PI.
	THEN II21 = _p2 - _p1 + 2 * arsin(cos( (_i2 + _i1)/2 ) * sin( (_n2 + _n1)/2 ) * sec(I21/2));
	ELSE II21 = _p2 - _p1 - 2 * arsin(cos( (_i2 + _i1)/2 ) * sin( (_n2 + _n1)/2 ) * sec(I21/2));
	
  	%if &D_TYPE EQ DD OR &D_TYPE = ALL %THEN %DO;
  		DROP THETA B1 B2 G1 G2;		
		B1 = arsin(sin(_i1) * sin(_p1));
		B2 = arsin(sin(_i2) * sin(_p2));
		G1 = _n1 + atan(cos(_i1) * tan(_p1));
		G2 = _n2 + atan(cos(_i2) * tan(_p2));
		IF cos(_P1) LT 0 THEN G1 = G1 + &PI;
		IF cos(_P2) LT 0 THEN G2 = G2 + &PI;
		Theta = arcos( sin(B1) * sin(B2) + cos(B1) * cos(B2) * cos(G2-G1) );
		D_Result_DD  = sqrt( ((_q - _q2)/(_q +_q2))**2 + ((_e - _e2)/(_e + _e2))**2 + (I21 / &PI)**2 +((_e2+_e)/2)**2 * (Theta/&PI)**2);
	%end;
		
  	%if &D_TYPE EQ DSH OR &D_TYPE = ALL %THEN %DO;	
		D_Result_DSH = sqrt ( (_q - _q2)**2 + (_e - _e2)**2 + (2 * sin(I21/2))**2 + ( (_e - _e2)/2 * 2 * sin(II21/2) )**2 );
	%end;
		
  	%if &D_TYPE EQ DH OR &D_TYPE = ALL %THEN %DO;	
		D_Result_DH  = sqrt ( ((_q - _q2)/(_q + _q2))**2 + (_e - _e2)**2 + (2 * sin(I21/2))**2 + ( (_e - _e2)/2 * 2 * sin(II21/2) )**2 );		
	%end;

%mend;


%macro D_Scan (dset = );

%DO I = 1 %to &Candidates.;

	DATA &dset.;

	FORMAT MATCH_D  12.2;
	LENGTH Stream_D $50.;
	
	/* The prevents the seed variables being reinitialised to null */
	RETAIN _i2 _q2 _e2 _p2 _n2;
	
	/* Specify input data set */
	SET &dset. ;
	
	/* Prevent temorary variables being written to output data set */
	DROP	_i1 _p1 _n1  
		_i2 _q2 _e2 _p2 _n2 
		I21 ADIFF II21 ;
	
	/* If first row of input dataset, set the seed orbit variables
	   from values held outside the data step in MACRO variables */
	IF _N_ = 1 THEN DO
		_i2 = %radians(&&_incl_&I.); 
		_q2    = &&_q_&I.;
		_e2    = &&_e_&I.;
		_p2 = %radians(&&_peri_&I.);
		_n2 = %radians(&&_node_&I.);
		END;
		
	If MATCH_D EQ . THEN DO;
		MATCH_D = 999999;
		Stream_D = "zz_No_match";
		END;	
			
	/* Convert all angles to RADIANS (input variables) */
		_n1 = %radians(_node);
		_i1 = %radians(_incl);
		_p1 = %radians(_peri);
		
		%D_Calc(D_TYPE=&D_CRITERION);

		/* Assign new stream IF we are below threshold AND D_NEW is below the 
		   D value for at the last stream assignment (i.e. can we better the match) */
		IF D_RESULT_&D_CRITERION. LE &&_D_THRESHOLD_&I AND D_RESULT_&D_CRITERION. LT MATCH_D THEN DO;
			Stream_D = "&&shower_name_&I.";
			MATCH_D = D_RESULT_&D_CRITERION. ;
			END;  
		run;
	
%END;	

%mend;

/* Run categorisation */
%D_SCAN(dset = WORK.OBSERVATIONS );


/* Sumarise */

options ORIENTATION=Landscape;
ods pdf file="/folders/myshortcuts/SAS_Uni/myfolders/Sample1.pdf";
TITLE1 "Summary";
TITLE2 "Orbit reclassification (using &D_CRITERION)";
Proc Tabulate data =work.observations;
Class _Stream Stream_D ;
Table (Stream_D ALL),(_Stream ALL) ;
Run;
ods pdf close;


/* Summary stream statistics */

TITLE1 "Summary";
TITLE2 "Mean orbits based on D_Criterion grouping";
proc sort data = work.OBSERVATIONS; by Stream_D; run;
proc means data = WORK.OBSERVATIONS noprint;
output 	out=work.means  
	mean(_peri) = _peri
	mean(_node) = _node
	mean(_incl) = _incl
	mean(_e) = _e
	mean(_q) = _q;
var _incl _node _peri _e _q;
by Stream_D;
run;

options ORIENTATION=Portrait;
ods pdf file="/folders/myshortcuts/SAS_Uni/myfolders/Sample2.pdf";
proc print data = work.means; run;
ods pdf close;

TITLE1;
TITLE2;



