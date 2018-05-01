import numpy as np
import scr.MarkovClasses as MarkovCls

delta_t=1/12
discount_rate=0.03

#Part1: annual rate: well to next state
#Caculate rate well to pneumoniae
annual_incidence_pneumoniae=7952/100000
rate_well_to_pneumoniae=-np.log(1-annual_incidence_pneumoniae)

#Caculate rate well to bacteremia
annual_incidence_bactermia= 19.8/100000
rate_well_to_bactermia=-np.log(1-annual_incidence_bactermia)

#Rate well to meningitis
annual_incidence_meningitis= 14/100000
rate_well_to_meningitis=-np.log(1-annual_incidence_meningitis)

#rate to AOM (S and NS)
annual_incidence_aom= 11525/100000
rate_well_to_aom= -np.log(1-annual_incidence_aom)
probability_tympanostomy =0.059
rate_well_to_aomtympanostomy = rate_well_to_aom*probability_tympanostomy
rate_well_to_aomwithouttympanostomy= rate_well_to_aom*(1-probability_tympanostomy)


#Part2 disease to outcome
# pneumoniae
worklossday_pneumoniae=7/365
rate_pneumoniae_to_recovery=1/worklossday_pneumoniae
case_fatality_pneumoniae=0.00526
rate_pneumoniae_to_death=rate_pneumoniae_to_recovery*case_fatality_pneumoniae/(1-case_fatality_pneumoniae)

#meningitis
worklossday_meningitis=9/365
rate_meningitis_to_recovery=1/worklossday_meningitis
case_fatality_meningitis=0.083
deafness_prob_meningitis=0.13
disability_prob_meningitis=0.07
rate_meningitis_to_death = rate_meningitis_to_recovery*case_fatality_meningitis/(1-case_fatality_meningitis-deafness_prob_meningitis-disability_prob_meningitis)
rate_memingitis_to_deafness= rate_meningitis_to_recovery*deafness_prob_meningitis/(1-case_fatality_meningitis-deafness_prob_meningitis-disability_prob_meningitis)
rate_memingitis_to_disability= rate_meningitis_to_recovery*disability_prob_meningitis/(1-case_fatality_meningitis-deafness_prob_meningitis-disability_prob_meningitis)


#bacteremia
worklossday_bacteremia=8/365
rate_bacteremia_to_recovery=1/worklossday_bacteremia
case_fatality_bacteremia=0.046
rate_bacteremia_to_death = rate_bacteremia_to_recovery*case_fatality_bacteremia/(1-case_fatality_bacteremia)

#aom
worklossday_aom=5/365
rate_aom_to_recovery=1/worklossday_aom

rate_matrix=[
    [None, rate_well_to_pneumoniae, rate_well_to_meningitis, 0, 0, rate_well_to_bactermia, rate_well_to_aomtympanostomy,rate_well_to_aomwithouttympanostomy,0],
    [rate_pneumoniae_to_recovery, None, 0, 0, 0, 0, 0, 0, rate_pneumoniae_to_death],
    [rate_meningitis_to_recovery,0,None,rate_memingitis_to_disability,rate_memingitis_to_deafness,0,0,0,rate_meningitis_to_death],
    [0,0,0,None,0,0,0,0,0],
    [0,0,0,0,None,0,0,0,0],
    [rate_bacteremia_to_recovery,0,0,0,0,None,0,0,rate_bacteremia_to_death],
    [rate_aom_to_recovery,0,0,0,0,0,None,0,0],
    [rate_aom_to_recovery,0,0,0,0,0,0,None,0],
    [0,0,0,0,0,0,0,0,None]
]
prob_matrix = MarkovCls.continuous_to_discrete(rate_matrix,delta_t)
print(prob_matrix)

cost_matrix=[
    0,652,3851,0,0,2863,119,357,0
]

utility_matrix=[
    1,0.9941,0.9768,0.7393,0.8611,0.9941,0.995,0.82,0
]

work_loss_day=[0,worklossday_pneumoniae,worklossday_meningitis,0,0,worklossday_bacteremia,worklossday_aom,worklossday_aom,0]

salary=30

cohort_size= 200000

efficacy= [0,0.90, 0.974, 0,0,0.974,0.576,0.576,0]
coverage= [0,0.923,0.871,0,0,0.871,0.815,0.815,0]
proportion=[0,0.07,0.095,0,0,1,0.307,0.307,0]

shot_efficacy=[0,0.86,1]

rate_matrix_vaccine=[
    [None, rate_well_to_pneumoniae, rate_well_to_meningitis, 0, 0, rate_well_to_bactermia, rate_well_to_aomtympanostomy,rate_well_to_aomwithouttympanostomy,0],
    [rate_pneumoniae_to_recovery, None, 0, 0, 0, 0, 0, 0, rate_pneumoniae_to_death],
    [rate_meningitis_to_recovery,0,None,rate_memingitis_to_disability,rate_memingitis_to_deafness,0,0,0,rate_meningitis_to_death],
    [0,0,0,None,0,0,0,0,0],
    [0,0,0,0,None,0,0,0,0],
    [rate_bacteremia_to_recovery,0,0,0,0,None,0,0,rate_bacteremia_to_death],
    [rate_aom_to_recovery,0,0,0,0,0,None,0,0],
    [rate_aom_to_recovery,0,0,0,0,0,0,None,0],
    [0,0,0,0,0,0,0,0,None]
]

vaccine_cost=146
vaccine_administration=2

rr = [0, 1, 2, 3, 4, 5, 6, 7, 8]
new_rate_matrix = rate_matrix
for i in range(0, 9):
    rr[i] = (1 - proportion[i] * coverage[i] * efficacy[i] )
for i in range(1, 9):
    new_rate_matrix[0][i] = rate_matrix[0][i] * rr[i]
print(rr,new_rate_matrix)

prob_matrix_vaccine=MarkovCls.continuous_to_discrete(new_rate_matrix,delta_t)
print(prob_matrix,prob_matrix_vaccine)

print(range(1,8))
