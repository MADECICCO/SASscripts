
/*
------------------------------------------------------------------------------
 This is a simple D_Criterion example in which unified orbit data is filtered
 on two basic quality criteria; all orbits are then compared with a single 
 FIXED reference orbit.  Columns D_SH, and D_H are added to the output SAS 
 dataset.
 
 NOTE: ONLY DSH and DH are calculated in this example.
 
 INPUT DATASET: CSV file (UFO Orbit UNIFIED observations in csv format)
 
 OUTPUT DATASET: WORK.D_CRITERION
 
------------------------------------------------------------------------------
*/

FILENAME REFFILE '/folders/myshortcuts/SAS_Uni/myfolders/UNIFIED.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.UA_IMPORT;
	GETNAMES=YES;
RUN;

/* Split the data into unified and by station */
DATA WORK.UNIFIED_QA_PASS WORK.UNIFIED_QA_FAILED;
	SET WORK.UA_IMPORT; 
	IF substr(_ID1,1,9) = "_(UNIFIED" THEN DO;
		IF _e <= 1 and abs(_dv12_ <=10.0) and abs(_Gm_ > 80.0) THEN OUTPUT UNIFIED_QA_PASS;
		ELSE OUTPUT UNIFIED_QA_FAILED;
	END;
RUN;

/* Summarise the orbital data of meteors meeting Quality Criteria */
proc sort data = WORK.UNIFIED_QA_PASS;
	by _stream;
run; 


/* Macro to convert degrees to Radians */
%macro radians(x);
&x. * 3.141592653589793 / 180
%mend;

/* Process unified observations and compare orbit with reference orbit */
Data work.d_Criterion;

* Note: 
	Perihelion distance (_q) 
	Eccentricity (_e), 
	Inclination (_incl)
	Argument of pericenter / perihelion (_p)
	Longitude of ascending node (_node)
;

/* Test only - set a reference orbit against which to comare dataset */
RETAIN _incl2 _q2 _e2 _peri2 _node2;
DROP   _incl1 _q1 _e1 _peri1 _node1  _incl2 _q2 _e2 _peri2 _node2 I II;

FORMAT D_SH 12.2;
FORMAT D_H  12.2;

SET work.unified_qa_pass (rename = (__ = Meteor_ID) );

IF _N_ = 1 THEN DO
	_incl2    = %radians(110.2); 
	_q2    = 0.960;
	_e2    = 0.881;
	_peri2 = %radians(152.5);
	_node2 = %radians(139.7);
	END;

/* Convert all angles to RADIANS */

	_node1 = %radians(_node);
	_incl1 = %radians(_incl);
	_peri1 = %radians(_peri);
	
/* Calculate I and II (Southworth & Hawkins) formula */	
		
	I = arcos(cos(_incl1) * cos(_incl2) + sin(_incl1) * sin(_incl2) * cos(_node1 - _node2) );
	IF abs(_node1 -_node2)) <= 180 
	THEN II = _peri2 - _peri1 + 2 * arsin(cos( (_incl2 + _incl1)/2 ) * sin( (_node2 + _node1)/2 ) * sec(I/2));
	ELSE II = _peri2 - _peri1 - 2 * arsin(cos( (_incl2 + _incl1)/2 ) * sin( (_node2 + _node1)/2 ) * sec(I/2));
	
/* Calculate the Southworth & Hawkins Dsh parameter */

	D_SH = sqrt ( (_q - _q2)**2 + (_e - _e2)**2 + (2 * sin(I/2))**2 + ( (_e - _e2)/2 * 2 * sin(II /2) )**2 );
	D_H  = sqrt ( ((_q - _q2)/(_q + _q2) )**2 + (_e - _e2)**2 + (2 * sin(I/2))**2 + ( (_e - _e2)/2 * 2 * sin(II /2) )**2 );
	
run;

