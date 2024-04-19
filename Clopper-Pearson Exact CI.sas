data testdata;
input trtpn x n alpha;
datalines;
1 2 13 0.05 
2 0 1 0.05
;
data testdat1;
set testdata;
y = n - x;
p = x/n;
z = probit(1-alpha/2);
run;
proc transpose data = testdat1 out = x_data(rename=(_NAME_=outcome COL1=count));
by trtpn;
var x y;
run;
ods output BinomialCLs = CLs1;
proc freq data = x_data;
by trtpn;
tables outcome /binomial(cl=(exact)) alpha=0.1;
weight count / zero;
run;
