/* See proc mix for difference*/
%let CurFile = &_SASPROGRAMFILE;
%let CurPath = %substr( &CurFile, 1, %eval(%length(&CurFile)-
%sysfunc(indexc(%sysfunc(reverse(&CurFile)),\/))) );
libname qclib "&CurPath/Data"; run;

DATA easi;
	set qclib.adea;
RUN;

PROC CONTENTS DATA=easi;
RUN;

data look;
              set easi;
              where PARAMCD='EASI' and avisitn>=0;
              keep  SUBJID AVISIT AVISITN VISIT VISITNUM TRTA TRTAN AVAL BASE dtype ANL01FL;
run;


data look1; set look;
	if (SUBJID='10011001' and AVISITN = 3 and ANL01FL='') or 
	   (SUBJID='10011001' and AVISITN = 12 and ANL01FL='')
	then delete;
run;

proc sort data=look1 out=look1;
              by subjid;
run;

data look1; set look1;
	if SUBJID='10011001' and AVISITN=12 then AVAL=.;
run;


proc means data= look1;
	class trta avisit;
run;

/* fit for the study: two treatment groups, week1 - week12 */
data look2; set look1;
	where (trtan=19 or trtan = 20) and avisitn <= 12;
	if trtan = 19 then dose=0;
    if trtan = 20 then dose=1000;
run;

data ex1;
              set look2;
              where avisitn = 12;
run;

proc sort data=ex1 out=ex2;
    *data must be sorted by regimen before usin proc mi;
              by dose;
run;

proc mi data=ex2 seed=5897 nimpute=1000 out=outimp  max=72 min=0;
              class dose;
              monotone regpmm(aval=base/details k=5);
              mnar model(aval/modelobs = (dose='0'));
              var base aval;
run;

data outimp1;
              set outimp;
              pchg = (aval - base)/base * 100;
run;

proc sort data = outimp1 out = outimp2;
    by _imputation_ dose subjid;
run;

proc mixed data=outimp2;
             by _imputation_;
             class dose(ref='0');
             model pchg = dose;
             lsmeans dose/ diff alpha=.1 ;
             ods output diffs=diffs lsmeans=lsmeans;
run;

data diffsout;
              set diffs;
              **only keep within regimen contrasts vs placebo;
              where _dose = 0;
run;

proc sort data=diffsout out=diffsout1;
              by dose _imputation_;
run;

proc mianalyze data=diffsout1  alpha=.1;  **specify alpha for 90% CIs here;
      by dose;
      modeleffects estimate;
      stderr stderr;
      ods output parameterestimates=parameterestimates_diff;
run;

proc sort data=parameterestimates_diff out = parameterestimates_diff1;
	by descending dose;
run;

data parameterestimates_diff2;
	set parameterestimates_diff1;
	if tvalue <=0 then one_side_Pvalue = Probt/2;
	if tvalue >0 then one_side_Pvalue = 1-Probt/2;
run;

proc transpose data= parameterestimates_diff2 out = transdiff;
	ID dose;
	var Estimate StdErr LCLMean UCLMean one_side_Pvalue;
run;

proc sort data=lsmeans out=lsmeans1;
              by dose _imputation_;
run;
proc mianalyze data=lsmeans1  alpha=.1;  **specify alpha for 90% CIs here;
      by dose;
      modeleffects estimate;
      stderr stderr;
      ods output parameterestimates=parameterestimates_lsmean;
run;

proc sort data=parameterestimates_lsmean out = parameterestimates_lsmean1;
	by descending dose;
run;

proc transpose data= parameterestimates_lsmean1 out = translsmean;
	ID dose;
	var Estimate StdErr LCLMean UCLMean;
run;

data dfinal;
	set translsmean transdiff;
run;

/* Output as csv file */
/* proc export data=dfinal outfile="&CurPath/test_ANCOVA_%CFBEASI.csv" dbms=csv  replace; run; */


/* Create the PDF file */
/* ods pdf file="&Curpath/test_ANCOVA_%CFBEASI.pdf";*/

title "Percent change from Baseline in EASI Score at Week 12 ANOVA, mITT (MI)";
proc print data=dfinal;
run;

/* ods pdf close;*/
