/* See proc mix for difference*/
%let CurFile = &_SASPROGRAMFILE;
%let CurPath = %substr( &CurFile, 1, %eval(%length(&CurFile)-
%sysfunc(indexc(%sysfunc(reverse(&CurFile)),\/))) );
libname qclib "&CurPath/Data"; run;

DATA easi;
	set qclib.adea;
RUN;


data look;
              set easi;
              where PARAMCD='EASI' and avisitn>0;
              keep SUBJID AVISIT AVISITN TRTA TRTAN AVAL BASE ANL01FL CRIT2FN CRIT3FN dtype;
run;


data look1; set look;
	if (SUBJID='10011001' and AVISITN = 3 and ANL01FL='') or 
	   (SUBJID='10011001' and AVISITN = 12 and ANL01FL='')
	then delete;
run;

data datres; set look1;
	pchg = (aval - base)/base * 100;
	easi50 = (pchg <= -50) and pchg ne '.';
	easi75 = (pchg <= -75) and pchg ne '.';
	easi90 = (pchg <= -90) and pchg ne '.';
	easi100 = (pchg <= -100) and pchg ne '.';
run;

%macro dfv (ana_var, ana_visit);

/* select week */
data look2; set datres;
	where (trtan=19 or trtan = 20) and avisit = &ana_visit;
	if trtan = 19 then dose="trt";
    if trtan = 20 then dose="con";
run;

/* Chan and Zhang*/
 proc sort data=look2 out=look2;
 by dose;
run; 
PROC BINOMIAL DATA=look2 ALPHA=0.9 out=resultci;
 by dose;
 BI/BS;
 OU &ana_var;
 RUN;
data outputci2_trt;
	set resultci;
	where by_val ='trt' and item in ("TRIALS", "SUCCESS", "EST_PI", "B_LCI_PI", "B_UCI_PI");
	keep item value;
run;
data outputci2_pbo;
	set resultci;
	where by_val='con' and item in ("TRIALS", "SUCCESS", "EST_PI", "B_LCI_PI", "B_UCI_PI");
	keep item value;
run;

proc transpose  data=outputci2_trt out = wideoutputci2_trt;
    id item;
    var value;
run;
proc transpose  data=outputci2_pbo out = wideoutputci2_pbo;
    id item;
    var value;
run;

data wideci; set wideoutputci2_trt wideoutputci2_pbo;
run;

PROC BINOMIAL DATA=look2 GAMMA=0 ALPHA=0.9 out=result;
 PD/EX ONE STD;
 PO dose;
 OU &ana_var;
RUN;

data output2;
	set result;
	where item in ("STAT", "TLWR_CI", "TUPR_CI", "XCTPVAL2");
	keep item value;
run;

proc transpose  data= output2 out = wideoutput2;
    id item;
    var value;
run;

data tempoutput2;
	set wideoutput2;
	if stat <=0 then XCTPVAL1 = 1-xctpval2/2;
	if stat >0 then XCTPVAL1 = xctpval2/2;
run;

data finaldiff75; set tempoutput2;
	p_value_one_side=XCTPVAL1;
	keep stat TLWR_CI TUPR_CI p_value_one_side; 
run;

data outall; set wideci finaldiff75;
	length Visit Variable $ 10;
	Visit = &ana_visit;
run;
%mend dfv;

%dfv(easi50, 'Week 1')
data final_all;
 set outall;
%dfv(easi50, 'Week 2')
data final_all;
 set final_all outall;
%dfv(easi50, 'Week 4')
data final_all;
 set final_all outall;
%dfv(easi50, 'Week 6')
data final_all;
 set final_all outall;
%dfv(easi50, 'Week 8')
data final_all;
 set final_all outall;
%dfv(easi50, 'Week 12')
data final_all;
 set final_all outall;
%dfv(easi50, 'Week 16')
data final_all;
 set final_all outall;
run;

title "EASI50 response at Week 6, Chan and Zhang 1999, mITT (NRI)";
proc print data = final_all;
run;

%dfv(easi75, 'Week 1')
data final_all;
 set outall;
%dfv(easi75, 'Week 2')
data final_all;
 set final_all outall;
%dfv(easi75, 'Week 4')
data final_all;
 set final_all outall;
%dfv(easi75, 'Week 6')
data final_all;
 set final_all outall;
%dfv(easi75, 'Week 8')
data final_all;
 set final_all outall;
%dfv(easi75, 'Week 12')
data final_all;
 set final_all outall;
%dfv(easi75, 'Week 16')
data final_all;
 set final_all outall;
run;

title "EASI75 response at Week 6, Chan and Zhang 1999, mITT (NRI)";
proc print data = final_all;
run;

%dfv(easi90, 'Week 1')
data final_all;
 set outall;
%dfv(easi90, 'Week 2')
data final_all;
 set final_all outall;
%dfv(easi90, 'Week 4')
data final_all;
 set final_all outall;
%dfv(easi90, 'Week 6')
data final_all;
 set final_all outall;
%dfv(easi90, 'Week 8')
data final_all;
 set final_all outall;
%dfv(easi90, 'Week 12')
data final_all;
 set final_all outall;
%dfv(easi90, 'Week 16')
data final_all;
 set final_all outall;
run;

title "EASI90 response at Week 6, Chan and Zhang 1999, mITT (NRI)";
proc print data = final_all;
run;

%dfv(easi100, 'Week 1')
data final_all;
 set outall;
%dfv(easi100, 'Week 2')
data final_all;
 set final_all outall;
%dfv(easi100, 'Week 4')
data final_all;
 set final_all outall;
%dfv(easi100, 'Week 6')
data final_all;
 set final_all outall;
%dfv(easi100, 'Week 8')
data final_all;
 set final_all outall;
%dfv(easi100, 'Week 12')
data final_all;
 set final_all outall;
%dfv(easi100, 'Week 16')
data final_all;
 set final_all outall;
run;

title "EASI100 response at Week 6, Chan and Zhang 1999, mITT (NRI)";
proc print data = final_all;
run;

