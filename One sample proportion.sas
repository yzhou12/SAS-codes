data testdata;
input trial x n alpha;
datalines;
1 81 263 0.05
2 0 20 0.05
;

data testdat1;
set testdata;
y = n - x;
p = x/n;
z = probit(1-alpha/2);
run;

proc transpose data = testdat1 out = x_data(rename=(_NAME_=outcome COL1=count));
by trial;
var x y;
run;

/* Clopper-Pearson (exact) method*/
ods output BinomialCLs = CLs1;
proc freq data = x_data ;
by trial;
tables outcome /binomial(cl=exact) alpha=0.1;
weight count / zero;
run;

/* koyama and chen (2008) for Simon's two-stage */
%macro Simon_p(num_success=, n=, n1=, r1=, p0=); 
*data pval; *uncomment to produce a data set along with run statement;
p_val = 0; 
do X = &num_success to &n; 
do X1=max((&r1+1),(X-(&n-&n1))) to min(X, &n1); 
X2 = X - X1; 
pr_x = PDF('BINOMIAL', X1, &p0, &n1) * PDF('BINOMIAL', X2, &p0, &n-&n1); 
p_val = p_val+pr_x; retain p_val; end; end; 
*run; *uncomment to produce a data set.; 
%mend;

%macro Simon_CI(num_success=, alpha=, n=, n1=, r1=, p0= ); 
data confint; length errors $100; *To catch generic errors; 
phat_estimate = &num_success/&n; *MLE for p; 
if &num_success > &r1 then do; *Only proceed if successes > r1; 
errors = "None!"; 
**Calculate the observed p-value; 
%Simon_p(num_success=&num_success, n=&n, n1=&n1, r1=&r1, p0=&p0); 
p_value = p_val;
**Calculate the lower confidence limit at alpha specified; 
p_val = 0; 
conf_L = 0; 
do while (p_val < &alpha); *start at 0. Increase until p>alpha; 
*Below sets the tolerance for the accuracy of the limit, defaulted to 0.0001. 
For more accuracy, choose a value closer to 0; 
conf_L = conf_L + 0.0001; 
%Simon_p(num_success=&num_success, n=&n, n1=&n1, r1=&r1, p0=conf_L); 
end;
***Calculate the upper confidence limit at alpha specified; 
p_val = 1; 
conf_U = 1; 
do while (p_val > (1-&alpha));*decrease from 1 until p<(1-alpha); 
conf_U = conf_U - 0.0001;*accuracy threshold: -0.0001; 
%Simon_p(num_success=&num_success, n=&n, n1=&n1, r1=&r1, p0=conf_U);
end; 
end;
***Generate report; 
if &num_success <= &r1 then do; 
p_value = .; 
conf_L = .; 
conf_U = .; 
errors = "Failed at stage 1, use exact binomial methods"; 
end; 
output; 
keep phat_estimate p_value conf_L conf_U errors; 
run;
%let conf_level=%SYSEVALF((1-&alpha)*100); 
%let conf_level2=%SYSEVALF((1-2*&alpha)*100);
proc print data=confint;
title "Simon's two-stage p-value and one-sided &conf_level% confidence limits"; 
title2 "N total: &n, Total successes: &num_success"; 
run; 
%mend;

%Simon_CI(num_success=3, alpha=0.1, n=24, n1=13, r1=1, p0=.1);
