data pre_diff;
length trt $41;
input trt $ trial success x col;
datalines;
act 26 0 0 26
act 26 0 10 0
pbo 19 0 10 0
pbo 19 0 0 19
;
run;
proc binomial data=pre_diff GAMMA=0 alpha=0.9 out=rel1; 
 RISKDIFF/EX ONE STD; 
 PO trt; 
 OU x;
 weight col;
run;

proc print data=rel1;
run;
