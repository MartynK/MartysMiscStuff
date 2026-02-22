use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
  
drop _merge
 
generate trt=.
replace trt=1 if Arm==1
replace trt=0 if Arm==0

generate Died2=.
replace Died2=1 if Died==1
replace Died2=0 if Died==0

generate HFH=.
replace HFH=1 if HF_anytime=="Yes"
replace HFH=0 if HF_anytime!="Yes"

***BNP change category

gen BNP_change = BNP_num_12 / BNP_num
gen BNP_change_5 = 0 if BNP_change > 0.95 & BNP_change !=.
replace BNP_change_5 = 1 if BNP_change <= 0.95 & BNP_change !=.

gen BNP_change_10 = 0 if BNP_change > 0.90 & BNP_change !=.
replace BNP_change_10 = 1 if BNP_change <= 0.90 & BNP_change !=.

***LVESV change category - response2

gen response2=.
replace response2= 0 if LVESV_change > 1 & LVESV_change !=. 
***non response
replace response2= 1 if LVESV_change > 0.85 & LVESV_change <= 1 
***0-15% közötti változás
replace response2= 2 if LVESV_change <= 0.85 & LVESV_change > 0.70
***15-30% közötti
replace response2= 3 if LVESV_change <= 0.70 
***30% feletti

save "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", replace


use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear

winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF)

winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_change c <) outcomes(BNP_num_12 c <) 

winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(responder c  >) outcomes(BNP_num_12 c <) 

winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(response c  >) outcomes(BNP_num_12 c <) 

winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(response2 c  >) outcomes(BNP_num_12 c <) 

winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(responder c  >) outcomes(BNP_diff_3 c <) 

winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(responder c  >) outcomes(BNP_change_5 c >) 

winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(responder c  >) outcomes(BNP_change_10 c >) 


***INCIDENCE death 10000 patient years

stset days_death, failure(Died==1) exit(time 365) id(Number) 
stptime, title(Person-years) per(10000)  dd(2) level(95)

*** at(1 (90) 365)                       

***INCIDENCE HFH 10000 patient years

stset days_HF, failure(HFH==1) exit(time 365) id(Number)
stptime, title(Person-years) per(10000) dd(2) level(95)

  
***SUBGROUPS

*** IMPORTANT: 1, RVpace_num <= 96.25 
keep if RVpace_num <= 96.25
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 


*** IMPORTANT: 1, RVpace_num > 96.25 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if RVpace_num > 96.25
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, Age <= 75 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if Age <= 75
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, Age > 75 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if Age > 75
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, LVEF_baseline <= 30 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if LVEF_baseline <= 30
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, LVEF_baseline > 30 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if LVEF_baseline > 30 & LVEF_baseline !=.
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, subgroup_ISCH_2 = 1
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if subgroup_ISCH_2==1
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, subgroup_ISCH_2 = 0
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if subgroup_ISCH_2==0
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, Diabetes == 0 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if Diabetes == 0
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, Diabetes == 1 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if Diabetes == 1
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, NYHA_baseline_num == 2 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if NYHA_baseline_num == 2
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1, NYHA_baseline_num == 3 | NYHA_baseline_num == 4 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if NYHA_baseline_num == 3 | NYHA_baseline_num == 4
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1,  BMI < 30 & BMI !=. 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if  BMI < 30 & BMI !=.
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1,  BMI >= 30 & BMI !=. 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if  BMI >= 30 & BMI !=.
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1,  BNP_num <= 2154  & BNP_num !=. 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if  BNP_num <= 2154  & BNP_num !=.
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1,  BNP_num > 2154 & BNP_num !=. 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if  BNP_num > 2154 & BNP_num !=.
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1,  egfr_ckdepi < 60 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if  egfr_ckdepi < 60
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 

*** IMPORTANT: 1,  egfr_ckdepi >= 60 
use "/Users/epi/Desktop/BUDAPEST CRT/Data\Stata Data\04_Outcomes_combined.dta", clear
keep if  egfr_ckdepi >= 60
winratio Number trt, outcomes(Died2 tf days_death) outcomes(HFH tf days_HF) outcomes(LVESV_12M c  <) outcomes(BNP_num_12 c <) 




