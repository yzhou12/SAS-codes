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

data look1;
              set easi;
              where PARAMCD='EASI' and avisitn>0 and anl01fl='Y';
              keep  SUBJID AVISIT AVISITN VISIT VISITNUM TRTA TRTAN AVAL BASE dtype ANL01FL;
run;

/* fit for the study: two treatment groups, week1 - week12 */
data look2; set look1;
	where (trtan=19 or trtan = 20);
	if trtan = 19 then dose=0;
    if trtan = 20 then dose=1000;
run;

data look3; set look2;
	keep trta avisit aval;
run;

proc means data= look3;
	class trta avisit;
run;

data look2;
              set look2;
              pchg = (aval - base)/base * 100;
run;

/* MMRM */
PROC MIXED DATA=look2 METHOD=REML;
 CLASS dose SUBJID avisitn;
 MODEL pchg = base dose avisitn dose*avisitn/DDFM=KR SOLUTION;
 REPEATED avisitn/ TYPE=UN SUBJECT=SUBJID;
 lsmeans avisitn*dose/CL alpha=0.1;
 ods output lsmeans=lsmeans;
 estimate 'xxxmg vs PLACEBO at week 1' dose -1 1 avisitn*dose -1 0 0 0 0 0 0 1 0 0 0 0 0 0 /cl alpha=0.1;
 estimate 'xxxmg vs PLACEBO at week 2' dose -1 1 avisitn*dose 0 -1  0 0 0 0 0 0 1 0 0 0 0 0 /cl alpha=0.1;
 estimate 'xxxmg vs PLACEBO at week 4' dose -1 1 avisitn*dose 0 0 -1 0 0 0 0 0 0 1 0 0 0 0 /cl alpha=0.1;
 estimate 'xxxmg vs PLACEBO at week 6' dose -1 1 avisitn*dose 0 0 0 -1 0 0 0 0 0 0 1 0 0 0 /cl alpha=0.1;
 estimate 'xxxmg vs PLACEBO at week 8' dose -1 1 avisitn*dose 0 0 0 0 -1 0 0 0 0 0 0 1 0 0 /cl alpha=0.1;
 estimate 'xxxmg vs PLACEBO at week 12' dose -1 1 avisitn*dose 0 0 0 0 0 -1 0 0 0 0 0 0 1 0 /cl alpha=0.1;
 estimate 'xxxmg vs PLACEBO at week 16' dose -1 1 avisitn*dose 0 0 0 0 0 0 -1 0 0 0 0 0 0 1 /cl alpha=0.1;
RUN;

/* longitudinal plot: need to modify (Nov 4) */
data w9; set lsmeans;
length window $ 7;
if AVISITN=7 then window="Week 1";
else if AVISITN=14 then window="Week 2";
else if AVISITN=28 then window="Week 4";
else if AVISITN=42 then window="Week 6";
else if AVISITN=56 then window="Week 8";
else if AVISITN=84 then window="Week 12";
else if AVISITN=112 then window="Week 16";
length Treatment $ 18;
if dose=0 then Treatment='Placebo';
else if dose=1000 then Treatment='PF-07242813 1000mg';
visit=1.0*avisitn;
mean1 = estimate;
week=avisitn/7;
if AVISITN=7 then labelx=0;
else labelx= week;
run;


proc sort data=w9; by Treatment visit;run;

proc sql;
	create table final_table as
	select * from w9 as x left join freqcnt as y
	on x.dose = y.dose and x.avisitn=y.avisitn;
quit;

proc sort data=final_table; by Treatment visit;run; 

data _orig1; set final_table;
by Treatment visit;
Mean=round(mean1, 0.001);
lcl=round(lower, 0.001);
ucl=round(upper, 0.001);
PCM1=compress(put(mean,8.1))||'['||compress(put(count,8.0))||']';
run;

proc template;
 define statgraph lipid_discrete_outer_stat_color;
  begingraph / designwidth=6in designheight=4.5in;
    entrytitle 'Mean percent change from baseline in EASI scores comparing xxx and Placebo - MMRM (mITT, OC)';
    entryfootnote halign=left '* xx.x[] is %CFB Mean [N]' ;
entryfootnote halign=left '** Referece studies are based on internal meta-analysis with simulated data (MBMA) in longitudinal models';

layout lattice / rows=2 rowweights=preferred  columndatarange=union;
      layout overlay / yaxisopts=(griddisplay=on label='mean %CFB in EASI with 90% CI') 
                       xaxisopts=(type=linear display=(TICKS TICKVALUES LINE LABEL) LABEL='Week'
 									linearopts=(tickvaluelist=(1 2 4 6 8 12 16)));
        scatterplot x=week y=Mean / group=Treatment groupdisplay=cluster clusterwidth=0.5 
          yerrorlower=lcl yerrorupper=ucl  
          markerattrs=(symbol=SquareFilled size=11) errorbarattrs=(thickness=2);
        seriesplot x=week y=Mean / group=Treatment groupdisplay=cluster clusterwidth=0.5 
          lineattrs=(pattern=solid thickness=3) name='s';
		referenceline y=-61.1 /
           lineattrs=(pattern=shortdash color=grey thickness=2) curvelabel="Dupilumab=-61.1 (Wk 6)" curvelabelposition=min;
		   		referenceline y=-66.2 /
           lineattrs=(pattern=shortdash color=grey thickness=2) curvelabel="Abrocitinib=-66.2 (Wk 6)" curvelabelposition=min;
		   referenceline y=-75.8 /
           lineattrs=(pattern=shortdash color=grey thickness=2) curvelabel="Upadacitinib=-75.8 (Wk 6)" curvelabelposition=min;
      endlayout;
      layout overlay  / walldisplay=none xaxisopts=(display=none) yaxisopts=(display=none);
        axistable x=labelx value=PCM1 / class=Treatment display=(label) 
          valueattrs=(size=10) 
          labelattrs=(size=10) colorgroup=Treatment;
      endlayout;
    endlayout;
  endgraph;
 end;
run;

%let dpi=300;
%let gpath="&CurPath/Output";
ods pdf file="&CurPath/Output/EASIpchg_all_LongiPlot&sysdate..pdf";
ods listing gpath=&gpath image_dpi=&dpi;
ods graphics /reset noscale attrpriority=color width=800px height=500px imagename='EASIpchg_all_LongiPlot';
proc sgrender data=_orig1 template=lipid_discrete_outer_stat_color;
run;
ods pdf close;

