/*====================================================================
  PROJECT : Drug X vs Placebo in Type 2 Diabetes  (MOCK STUDY)
  FLOW    : raw CSV -> SDTM -> ADaM -> TLF (Tables, Listings, Figures)
  ENDPOINT: Change from baseline (CHG) in fasting glucose at Week 12
  TIP     : Run ONE STEP AT A TIME and read the Log after each step.
            Green/blue log = fine. Red "ERROR" = fix before moving on.
====================================================================*/


/*--------------------------------------------------------------------
  STEP 0 : SETUP
  WHY: SAS needs to know where your files live. A LIBNAME is a nickname
       for a folder; datasets saved in sdtm./adam. are permanent.
  BEFORE RUNNING:
   1. In SAS Studio's Files pane create a folder  drugx_project
   2. Inside it create a folder  raw  and upload DM.csv, EX.csv, LB.csv
   3. Right-click drugx_project > Properties > copy the path below
--------------------------------------------------------------------*/
%let root = /home/u12345678/drugx_project;      /* <-- CHANGE THIS */

options dlcreatedir nodate nonumber;            /* dlcreatedir = auto-create folders */
libname sdtm "&root./sdtm";
libname adam "&root./adam";
libname outp "&root./output";


/*--------------------------------------------------------------------
  STEP 1 : IMPORT RAW CSV FILES
  WHAT: PROC IMPORT reads a CSV and guesses each column's type.
  WHY : every project starts by bringing raw data into SAS.
--------------------------------------------------------------------*/
proc import datafile="&root./raw/DM.csv" out=work.raw_dm dbms=csv replace;
  getnames=yes; guessingrows=max;
run;

proc import datafile="&root./raw/EX.csv" out=work.raw_ex dbms=csv replace;
  getnames=yes; guessingrows=max;
run;

proc import datafile="&root./raw/LB.csv" out=work.raw_lb dbms=csv replace;
  getnames=yes; guessingrows=max;
run;

/* ALWAYS inspect what was imported. Look at the Type column:
   dates like 2025-01-25 are normally read as Num with format YYMMDD10.
   (If RFSTDTC shows as Char instead, delete the put() in Step 2.)      */
proc contents data=work.raw_dm varnum; run;
proc print data=work.raw_dm(obs=5); run;


/*--------------------------------------------------------------------
  STEP 2 : SDTM  (Study Data Tabulation Model)
  WHAT: a standard layout for the data AS COLLECTED, one dataset per
        "domain" (DM = demographics, EX = exposure/dosing, LB = labs).
  WHY : regulators (FDA etc.) expect this structure, so every company
        organises collected data the same way.
  KEY RULES used here:
   - USUBJID = unique subject ID across the whole study
   - DOMAIN  = two-letter domain code
   - dates are stored as ISO 8601 TEXT in variables ending in DTC
     (numeric SAS dates are created later, in ADaM)
--------------------------------------------------------------------*/

/* ---- SDTM.DM ---- */
data sdtm.dm;
  set work.raw_dm(rename=(RFSTDTC=_rfst));
  length DOMAIN $2 USUBJID $30 RFSTDTC $10 AGEU $5 ARMCD $8;
  DOMAIN  = "DM";
  USUBJID = catx("-", STUDYID, SUBJID);
  RFSTDTC = put(_rfst, yymmdd10.);          /* numeric date -> "2025-01-25" */
  AGEU    = "YEARS";
  if ARM = "Drug X" then ARMCD = "DRUGX";
  else if ARM = "Placebo" then ARMCD = "PBO";
  drop _rfst;
run;

/* ---- SDTM.EX ---- */
data sdtm.ex;
  set work.raw_ex(rename=(EXSTDTC=_st EXENDTC=_en));
  length DOMAIN $2 USUBJID $30 EXSTDTC EXENDTC $10;
  DOMAIN  = "EX";
  USUBJID = catx("-", STUDYID, SUBJID);
  EXSTDTC = put(_st, yymmdd10.);
  EXENDTC = put(_en, yymmdd10.);
  drop _st _en;
run;

/* ---- SDTM.LB ----
   LBORRES = result as originally collected (character)
   LBSTRESN = standardised numeric result, used for analysis           */
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


/*--------------------------------------------------------------------
  STEP 3A : ADaM - ADSL  (Subject-Level Analysis Dataset)
  WHAT: ONE row per subject; holds treatment group, population flags,
        treatment dates. Every other ADaM dataset borrows from it.
  WHY : analysis datasets are "analysis-ready" - the TLF programs
        should need almost no extra data cleaning.
--------------------------------------------------------------------*/
proc sort data=sdtm.dm out=dm_s;                                  by USUBJID; run;
proc sort data=sdtm.ex out=ex_s(keep=USUBJID EXTRT EXSTDTC EXENDTC); by USUBJID; run;

data adam.adsl;
  merge dm_s(in=in_dm) ex_s;
  by USUBJID;
  if in_dm;                                   /* keep every subject in DM */

  length TRT01P TRT01A $10 SAFFL $1;
  TRT01P = ARM;                               /* P = planned treatment  */
  TRT01A = EXTRT;                             /* A = actual treatment   */
  if      TRT01P = "Placebo" then TRT01PN = 1;
  else if TRT01P = "Drug X"  then TRT01PN = 2;

  TRTSDT = input(EXSTDTC, yymmdd10.);         /* text date -> numeric   */
  TRTEDT = input(EXENDTC, yymmdd10.);
  format TRTSDT TRTEDT date9.;
  TRTDUR = TRTEDT - TRTSDT + 1;               /* days on treatment      */

  if not missing(TRTSDT) then SAFFL = "Y";    /* Safety population flag */
  else SAFFL = "N";

  keep STUDYID USUBJID SUBJID SITEID AGE AGEU SEX RACE
       TRT01P TRT01PN TRT01A TRTSDT TRTEDT TRTDUR SAFFL;
run;


/*--------------------------------------------------------------------
  STEP 3B : ADaM - ADLB  (Lab Analysis Dataset, glucose only)
  WHAT: one row per subject per visit.
        AVAL = analysis value, BASE = baseline value,
        CHG  = AVAL - BASE (change from baseline)
  WHY : "change from baseline" is the standard way to measure whether
        a drug moved a lab value, since people start at different levels.
--------------------------------------------------------------------*/
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
  if AVISIT = "Baseline" then ABLFL = "Y";    /* baseline record flag */
  keep STUDYID USUBJID PARAMCD PARAM AVISIT AVISITN ADT AVAL ABLFL;
run;

/* pull out each subject's baseline value */
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
  merge adlb0(in=in_lb) base adsl_trt;        /* one-to-many merge */
  by USUBJID;
  if in_lb;

  /* only post-baseline visits get a change value; skip if anything missing */
  if AVISITN > 0 and nmiss(AVAL, BASE) = 0 then do;
    CHG  = AVAL - BASE;
    PCHG = 100 * CHG / BASE;                  /* % change from baseline */
  end;

  label AVAL="Analysis Value" BASE="Baseline Value"
        CHG="Change from Baseline" PCHG="% Change from Baseline"
        AVISIT="Analysis Visit" TRT01P="Planned Treatment";
run;

/* quick check: does the derivation look right? */
proc print data=adam.adlb(obs=10) noobs; run;


/*--------------------------------------------------------------------
  STEP 4 : TLF  (Tables, Listings, Figures)
--------------------------------------------------------------------*/
data wk12;                                    /* analysis set for the table */
  set adam.adlb;
  where AVISIT = "Week 12" and SAFFL = "Y";
run;

/* ---- TABLE 1 : who was in the study (demographics) ---- */
title1 "Table 14.1.1: Demographics";
proc means data=adam.adsl n mean std min max maxdec=1;
  class TRT01P;
  var AGE;
run;
proc freq data=adam.adsl;
  tables (SEX RACE)*TRT01P / nopercent norow nocum;
run;
title;

/* ---- TABLE 2 : main efficacy table, shown as a Word-readable RTF ---- */
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

/* same summary saved as a dataset, then exported to CSV (like the reference repo) */
proc means data=wk12 noprint nway;
  class TRT01P;
  var CHG;
  output out=summary_table(drop=_type_ _freq_)
         n=N mean=MEAN std=SD median=MEDIAN min=MIN max=MAX;
run;

proc export data=summary_table outfile="&root./output/SUMMARY_TABLE.csv"
            dbms=csv replace;
run;

/* ---- STATISTICAL TESTS ----
   (a) two-sample t-test : simple, good for learning
   (b) ANCOVA            : the standard in real trials - compares the
       treatment groups' change AFTER adjusting for baseline glucose.   */
title1 "Two-sample t-test on change from baseline at Week 12";
proc ttest data=wk12;
  class TRT01P;
  var CHG;
run;

title1 "ANCOVA: CHG = Treatment + Baseline";
proc glm data=wk12;
  class TRT01P;
  model CHG = TRT01P BASE / solution;
  lsmeans TRT01P / pdiff cl stderr;           /* adjusted means + difference + 95% CI */
run;
quit;
title;

/* ---- LISTING : subject-level glucose data ---- */
proc sort data=adam.adlb out=listing(keep=USUBJID TRT01P AVISIT ADT AVAL BASE CHG);
  by TRT01P USUBJID AVISITN;
run;

title1 "Listing 16.2.1: Fasting Glucose by Subject and Visit";
proc print data=listing noobs label;
run;
title;

proc export data=listing outfile="&root./output/LISTING.csv" dbms=csv replace;
run;

/* ---- FIGURE : mean change from baseline by visit, +/- 1 SE ---- */
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

ods listing gpath="&root./output";            /* saves the PNG into output folder */
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
