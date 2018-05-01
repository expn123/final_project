import scr.MarkovClasses as MarkovCls
import Matrix as ma
import scr.RandomVariantGenerators as rndClasses
import scr.EconEvalClasses as EconCls
from enum import Enum



class HealthStats(Enum):
    """ health states of patients with risk of stroke """
    WELL = 0
    Pneumoniae = 1
    Meningitis = 2
    Disability = 3
    Deaf=4
    Bacteremia=5
    AOM_T=6
    AOM_NT=7
    DEATH=8

class PCV(Enum):
    without_vacc=0
    with_pcv=1

def gen_vaccine_prob_matrix(shot):
    ef=ma.shot_efficacy[shot-1]
    rr=[0,1,2,3,4,5,6,7,8]
    new_rate_matrix=ma.rate_matrix
    for i in range(0,8):
        rr[i]=(1-ma.proportion[i]*ma.coverage[i]*ma.efficacy[i]*ef)
    for i in range(1,8):
        new_rate_matrix[0][i]=ma.rate_matrix[0][i]*rr[i]
    return new_rate_matrix

class Patient:
    def __init__(self, id, vaccination):
        """ initiates a patient
        :param id: ID of the patient
        :param parameters: parameter object
        """

        self._id = id
        # random number generator for this patient
        self._rng = None
        self.healthstat=0
        self._ndisability=0
        self._ndeaf=0
        self._ndeath=0
        self.survival=0
        self.totalDiscountUtility=0
        self.totalDiscountCost=0
        self.vaccine=vaccination
        self.shot=0
        self.aom=0
        self.pneumonaie=0
        self.meningitis=0


    def simulate_fiveshort(self, sim_length_short):
        """ simulate the patient over the specified simulation length """

        # random number generator for this patient
        self._rng = rndClasses.RNG(self._id)


        if self.vaccine==0:
            k=0
            # while the patient is alive and simulation length is not yet reached
            while (self.healthstat != 8) and k * ma.delta_t < sim_length_short:
                # find the transition probabilities of the future states
                trans_probs = ma.prob_matrix[0][self.healthstat]
                # create an empirical distribution
                empirical_dist = rndClasses.Empirical(trans_probs)
                # sample from the empirical distribution to get a new state
                # (returns an integer from {0, 1, 2, ...})
                new_state_index = empirical_dist.sample(self._rng)
                # caculate cost and utality
                cost = ma.cost_matrix[self.healthstat] + ma.salary * ma.work_loss_day[self.healthstat]
                utility = ma.utility_matrix[self.healthstat]*ma.delta_t
                # update total discounted cost and utility (corrected for the half-cycle effect)
                self.totalDiscountCost += \
                    EconCls.pv(cost, ma.discount_rate * ma.delta_t, k + 1)
                self.totalDiscountUtility += \
                    EconCls.pv(utility, ma.discount_rate * ma.delta_t, k + 1)
                # update diseases:
                if self.healthstat==HealthStats.Pneumoniae.value:
                    self.pneumonaie+=1
                if self.healthstat==HealthStats.Meningitis.value:
                    self.meningitis+=1
                if self.healthstat==HealthStats.AOM_T.value or self.healthstat==HealthStats.AOM_NT.value:
                    self.aom+=1
                # update disability number
                if self.healthstat == HealthStats.Disability.value:
                    self._ndisability = 1
                # update deafness number
                if self.healthstat == HealthStats.Deaf.value :
                    self._ndeaf = 1
                # update health state
                self.healthstat = new_state_index[0]
                #update number of deahts
                if self.healthstat==HealthStats.DEATH.value:
                    self._ndeath=1
                # increment time step
                k += 1
        if self.vaccine==1:
            k=0
            while (self.healthstat!=8) and k * ma.delta_t<sim_length_short:
                # find the transition probabilities of the future states
                trans_probsv = ma.prob_matrix_vaccine[0][self.healthstat]
                # create an empirical distribution
                empirical_distv = rndClasses.Empirical(trans_probsv)
                # sample from the empirical distribution to get a new state
                # (returns an integer from {0, 1, 2, ...})
                new_state_indexv = empirical_distv.sample(self._rng)
                # caculate cost and utality
                cost = ma.cost_matrix[self.healthstat] + ma.salary * ma.work_loss_day[self.healthstat]
                utility = ma.utility_matrix[self.healthstat]*ma.delta_t
                # update total discounted cost and utility (corrected for the half-cycle effect)
                self.totalDiscountCost += \
                    EconCls.pv(cost, ma.discount_rate * ma.delta_t, k + 1)
                self.totalDiscountUtility += \
                    EconCls.pv(utility, ma.discount_rate * ma.delta_t, k + 1)
                # update diseases:
                if self.healthstat==HealthStats.Pneumoniae.value:
                    self.pneumonaie+=1
                if self.healthstat==HealthStats.Meningitis.value:
                    self.meningitis+=1
                if self.healthstat==HealthStats.AOM_T.value or self.healthstat==HealthStats.AOM_NT.value:
                    self.aom+=1
                # update disability number
                if self.healthstat == HealthStats.Disability.value:
                    self._ndisability = 1
                # update deafness number
                if self.healthstat == HealthStats.Deaf.value :
                    self._ndeaf = 1
                # update health state
                self.healthstat = new_state_indexv[0]
                #update number of deahts
                if self.healthstat==HealthStats.DEATH.value:
                    self._ndeath=1


                if k==3:
                    self.shot+=1
                    self.totalDiscountCost+= \
                        EconCls.pv(ma.vaccine_administration+ma.vaccine_cost, ma.discount_rate * ma.delta_t, k + 1)
                if k==5:
                    self.shot+=1
                    self.totalDiscountCost+= \
                        EconCls.pv(ma.vaccine_administration+ma.vaccine_cost, ma.discount_rate * ma.delta_t, k + 1)
                if k==11:
                    self.shot+=1
                    self.totalDiscountCost+= \
                        EconCls.pv(ma.vaccine_administration+ma.vaccine_cost, ma.discount_rate * ma.delta_t, k + 1)



                # increment time step
                k += 1
        if self.healthstat==3:
            for i in (6,16):
              self.totalDiscountCost+=EconCls.pv(2746,ma.discount_rate,i+1)




    def get_disability_number(self):
        """ returns the patient's survival time"""
        return self._ndisability

    def get_deaf_number(self):
        """ returns the patient's survival time"""
        return self._ndeaf

    def get_total_utility(self):
        return self.totalDiscountUtility

    def get_total_cost(self):
        return self.totalDiscountCost

class Cohort():
    def __init__(self,id,vaccine):
        self._initial_pop_size=ma.cohort_size
        self.id=id
        self._nudisability=0
        self._ndeaf=0
        self._ndeath=0
        self.pneumonaie=0
        self.meningitis=0
        self.aom=0
        self.totaldiscountedcost=[]
        self.totaldiscountedutility=[]
        self.vaccination=vaccine

    def simulate(self):
        for i in range(self._initial_pop_size):
            patient=Patient(self.id*self._initial_pop_size+i,self.vaccination)
            patient.simulate_fiveshort(5)
            self._nudisability+=patient.get_disability_number()
            self._ndeaf+=patient.get_deaf_number()
            self.pneumonaie+=patient.pneumonaie
            self.meningitis+=patient.meningitis
            self.aom+=patient.aom
            self._ndeath+=patient._ndeath
            self.totaldiscountedcost.append(patient.get_total_cost())
            self.totaldiscountedutility.append(patient.get_total_utility())
            print(i)
            if patient.healthstat==8:
                self._ndeath+=1



    def get_disability_number(self):
        return self._nudisability

    def get_deaf_number(self):
        """ returns the patient's survival time"""
        return self._ndeaf

    def get_death_number(self):
        return self._ndeath

    def get_pneumonaie_number(self):
        return self.pneumonaie

    def get_aom_number(self):
        return self.aom

    def get_meningitis_number(self):
        return self.meningitis

    def get_total_utility(self):
        return self.totaldiscountedutility

    def get_total_cost(self):
        return self.totaldiscountedcost



def report_CEA_CBA(simOutputs_none, simOutputs_anticoag):
    no_therapy_strategy=EconCls.Strategy(name="without vaccination", cost_obs=simOutputs_none.get_total_cost(),
                                      effect_obs=simOutputs_none.get_total_utility())
    anticoag_therapy_strategy=EconCls.Strategy(name="With vaccination", cost_obs=simOutputs_anticoag.get_total_cost(),
                                            effect_obs=simOutputs_anticoag.get_total_utility())

    listofStrategies = [no_therapy_strategy, anticoag_therapy_strategy]

    CEA = EconCls. CEA(listofStrategies, if_paired=False)

    CEA.show_CE_plane(
        title='Cost-Effectiveness Analysis',
        x_label='Additional discounted utility',
        y_label='Additional discounted cost',
        show_names=True,
        show_clouds=True,
        show_legend=True,
        figure_size=6,
        transparency=0.3
    )
    # report the CE table
    CEA.build_CE_table(
        interval=EconCls.Interval.CONFIDENCE,
        alpha=0.05,
        cost_digits=0,
        effect_digits=2,
        icer_digits=2,
    )

    CBA = EconCls.CBA(listofStrategies, if_paired=False)

    CBA.graph_deltaNMB_lines(
        min_wtp=0,
        max_wtp=50000,
        x_label="Willingness-to-pay for one additional QALY ($)",
        y_label="Incremental Net Monetary Benefit ($)",
        interval=EconCls.Interval.CONFIDENCE,
        transparency=0.4,
        show_legend=True,
        figure_size=6,
        title='cost benefit analysis')



a=Cohort(3,0)
b=Cohort(9,1)

a.simulate()
b.simulate()
print(a.get_pneumonaie_number(),a.get_aom_number(),a.get_meningitis_number(),a.get_death_number(),a.get_deaf_number(),a.get_disability_number())
print(b.get_pneumonaie_number(),b.get_aom_number(),b.get_meningitis_number(),b.get_death_number(),b.get_deaf_number(),b.get_disability_number())

report_CEA_CBA(a,b)
