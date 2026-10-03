%let root =/home/u64552773/drugx_project ;
options dlcreatedir nodate nonumber;
libname sdtm "&root./sdtm";
libname adam "&root./adam";
libname outp "&root./output";

proc import datafile="&root./raw/DM.csv" out=work.raw_dm dbms=csv replace;
  getnames=yes; guessingrows=max;
run;

proc import datafile="&root./raw/EX.csv" out=work.raw_ex dbms=csv replace;
  getnames=yes; guessingrows=max;
run;

proc import datafile="&root./raw/LB.csv" out=work.raw_lb dbms=csv replace;
  getnames=yes; guessingrows=max;
run;

proc contents data=work.raw_dm varnum; run;
proc print data=work.raw_dm(obs=5); run;

data sdtm.dm;
  set work.raw_dm(rename=(RFSTDTC=_rfst));
  length DOMAIN $2 USUBJID $30 RFSTDTC $10 AGEU $5 ARMCD $8;
  DOMAIN  = "DM";
  USUBJID = catx("-", STUDYID, SUBJID);
  RFSTDTC = put(_rfst, yymmdd10.);
  AGEU    = "YEARS";
  if ARM = "Drug X" then ARMCD = "DRUGX";
  else if ARM = "Placebo" then ARMCD = "PBO";
  drop _rfst;
run;

data sdtm.ex;
  set work.raw_ex(rename=(EXSTDTC=_st EXENDTC=_en));
  length DOMAIN $2 USUBJID $30 EXSTDTC EXENDTC $10;
  DOMAIN  = "EX";
  USUBJID = catx("-", STUDYID, SUBJID);
  EXSTDTC = put(_st, yymmdd10.);
  EXENDTC = put(_en, yymmdd10.);
  drop _st _en;
run;

data sdtm.lb;
  set work.raw_lb(rename=(LBORRES=_res LBDTC=_dtc));
  length DOMAIN $2 USUBJID $30 LBORRES $12 LBDTC $10 LBSTRESU $10;
  DOMAIN   = "LB";
  USUBJID  = catx("-", STUDYID, SUBJID);
  if not missing(_res) then LBORRES = strip(put(_res, best12.));
  LBSTRESN = _res;
  LBSTRESU = LBORRESU;
  LBDTC    = put(_dtc, yymmdd10.);
  drop _res _dtc;
run;

proc sort data=sdtm.dm out=dm_s;                                  by USUBJID; run;
proc sort data=sdtm.ex out=ex_s(keep=USUBJID EXTRT EXSTDTC EXENDTC); by USUBJID; run;

data adam.adsl;
  merge dm_s(in=in_dm) ex_s;
  by USUBJID;
  if in_dm;

  length TRT01P TRT01A $10 SAFFL $1;
  TRT01P = ARM;
  TRT01A = EXTRT;
  if      TRT01P = "Placebo" then TRT01PN = 1;
  else if TRT01P = "Drug X"  then TRT01PN = 2;

  TRTSDT = input(EXSTDTC, yymmdd10.);
  TRTEDT = input(EXENDTC, yymmdd10.);
  format TRTSDT TRTEDT date9.;
  TRTDUR = TRTEDT - TRTSDT + 1;

  if not missing(TRTSDT) then SAFFL = "Y";
  else SAFFL = "N";

  keep STUDYID USUBJID SUBJID SITEID AGE AGEU SEX RACE
       TRT01P TRT01PN TRT01A TRTSDT TRTEDT TRTDUR SAFFL;
run;

data adlb0;
  set sdtm.lb;
  where LBTESTCD = "GLUC";
  length PARAMCD $8 PARAM $40 AVISIT $10 ABLFL $1;
  PARAMCD = LBTESTCD;
  PARAM   = "Glucose (mg/dL)";
  AVAL    = LBSTRESN;
  ADT     = input(LBDTC, yymmdd10.);
  format ADT date9.;
  AVISIT  = VISIT;
  select (AVISIT);
    when ("Screening") AVISITN = -1;
    when ("Baseline")  AVISITN = 0;
    when ("Week 4")    AVISITN = 4;
    when ("Week 8")    AVISITN = 8;
    when ("Week 12")   AVISITN = 12;
    otherwise          AVISITN = .;
  end;
  if AVISIT = "Baseline" then ABLFL = "Y";
  keep STUDYID USUBJID PARAMCD PARAM AVISIT AVISITN ADT AVAL ABLFL;
run;

data base;
  set adlb0;
  where ABLFL = "Y";
  BASE = AVAL;
  keep USUBJID BASE;
run;

proc sort data=adlb0;  by USUBJID AVISITN; run;
proc sort data=base;   by USUBJID;         run;
proc sort data=adam.adsl out=adsl_trt(keep=USUBJID TRT01P TRT01PN SAFFL);
  by USUBJID;
run;

data adam.adlb;
  merge adlb0(in=in_lb) base adsl_trt;
  by USUBJID;
  if in_lb;

  if AVISITN > 0 and nmiss(AVAL, BASE) = 0 then do;
    CHG  = AVAL - BASE;
    PCHG = 100 * CHG / BASE;
  end;

  label AVAL="Analysis Value" BASE="Baseline Value"
        CHG="Change from Baseline" PCHG="% Change from Baseline"
        AVISIT="Analysis Visit" TRT01P="Planned Treatment";
run;

proc print data=adam.adlb(obs=10) noobs; run;

data wk12;
  set adam.adlb;
  where AVISIT = "Week 12" and SAFFL = "Y";
run;

title1 "Table 14.1.1: Demographics";
proc means data=adam.adsl n mean std min max maxdec=1;
  class TRT01P;
  var AGE;
run;
proc freq data=adam.adsl;
  tables (SEX RACE)*TRT01P / nopercent norow nocum;
run;
title;

ods rtf file="&root./output/Table_Glucose_Wk12.rtf" style=journal;
title1 "Table 14.2.1: Change from Baseline in Fasting Glucose (mg/dL) at Week 12";
title2 "Safety Population";
footnote1 "Change = Week 12 value minus Baseline value. SD = standard deviation.";
proc means data=wk12 n mean std median min max maxdec=1;
  class TRT01P;
  var BASE AVAL CHG;
run;
title; footnote;
ods rtf close;

proc means data=wk12 noprint nway;
  class TRT01P;
  var CHG;
  output out=summary_table(drop=_type_ _freq_)
         n=N mean=MEAN std=SD median=MEDIAN min=MIN max=MAX;
run;

proc export data=summary_table outfile="&root./output/SUMMARY_TABLE.csv"
            dbms=csv replace;
run;

title1 "Two-sample t-test on change from baseline at Week 12";
proc ttest data=wk12;
  class TRT01P;
  var CHG;
run;

title1 "ANCOVA: CHG = Treatment + Baseline";
proc glm data=wk12;
  class TRT01P;
  model CHG = TRT01P BASE / solution;
  lsmeans TRT01P / pdiff cl stderr;
run;
quit;
title;

proc sort data=adam.adlb out=listing(keep=USUBJID TRT01P AVISIT ADT AVAL BASE CHG);
  by TRT01P USUBJID AVISITN;
run;

title1 "Listing 16.2.1: Fasting Glucose by Subject and Visit";
proc print data=listing noobs label;
run;
title;

proc export data=listing outfile="&root./output/LISTING.csv" dbms=csv replace;
run;

proc means data=adam.adlb noprint nway;
  where AVISITN > 0;
  class TRT01P AVISITN;
  var CHG;
  output out=vis_stats(drop=_type_ _freq_) n=N mean=MEAN stderr=SE;
run;

data vis_stats;
  set vis_stats;
  lower = MEAN - SE;
  upper = MEAN + SE;
run;

ods listing gpath="&root./output";
ods graphics / imagename="Fig_Glucose_Change" imagefmt=png width=7in height=4.5in;

title1 "Figure 14.2.1: Mean Change from Baseline in Fasting Glucose by Visit";
title2 "Error bars = +/- 1 standard error";
proc sgplot data=vis_stats;
  series  x=AVISITN y=MEAN / group=TRT01P name="s1" markers;
  scatter x=AVISITN y=MEAN / group=TRT01P yerrorlower=lower yerrorupper=upper;
  refline 0 / axis=y lineattrs=(pattern=dash);
  xaxis values=(4 8 12) label="Week";
  yaxis label="Mean change from baseline (mg/dL)";
  keylegend "s1" / title="Treatment";
run;
title;
