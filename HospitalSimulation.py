import simpy
import numpy as np
import rngStream
import copy
import pandas as pd
import matplotlib.pyplot as plt
import sys
import time


SIM_TIME = 2500
PATIENTS = 100000
SCHEDULER = "first_min"
SCHEDULERS = ['first_min','last_min','uniform','custom']
SPECIALTIES = ["Home","Primary Care","Opthamology","Immunology","Cardiology","Other Dr Specialty","Urology","Gastroenterology","General Surgery","OB/GYN","Nephrology","Orthopedics","Psychiatry","Neurology","Oncology","Dermatology","Otorhinolaryngology"]
#SPECIALTY_COUNTS=[1,99,6,2,10,25,3,5,8,13,4,8,15,5,8,4,3]
#SPECIALTY_CAPS=[1,12,12,12,12,12,16,14,13,16,14,16,8,8,12,16,16]

# UNCAP settings
SPECIALTY_CAPS=[9999]*17
SPECIALTY_COUNTS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

TRANSITION_MATRIX = np.zeros([3,17,17])
TIME_LAG_MATRIX = np.zeros([3,17,17])
PATIENT_CLASSES = ['A','B','C']
PATIENT_CLASS_PROBS = [.7054,.1872,.1074]
trans_mat_raw = pd.read_csv("parameters/mrp_trans_probs.csv", sep=',',header=0,index_col='from_desc')
mtd_mat_raw = pd.read_csv("parameters/mrp_mtds.csv", sep=",",header=0,index_col='from_desc')

for i in range(3):
    TRANSITION_MATRIX[i] = trans_mat_raw.loc[trans_mat_raw["conditions"] == PATIENT_CLASSES[i]].iloc[:,1:]
    TIME_LAG_MATRIX[i] = mtd_mat_raw.loc[mtd_mat_raw["conditions"]==PATIENT_CLASSES[i]].iloc[:,1:]

DEBUG_LEVEL = 1
TRIAL = False

if TRIAL:
    # RNGSTREAM setup
    seed0 = [23456,23456,23456,23456,23456,23456]
    rngStream.SetPackageSeed(seed0)
unigens = {}
unigens["refer"] = rngStream.RngStream()
unigens["trans_exp"] = rngStream.RngStream()
unigens["init-class"] = rngStream.RngStream()
unigens["init-specialty"] = rngStream.RngStream()

def print_v(debug, *args):
    if debug <= DEBUG_LEVEL:
        print(*args)

class Specialty(object):
    scheduling: simpy.Resource
    specialty: str
    id: int
    appointments: [simpy.Resource]
    data: []
    daily_capacity: int
    def __init__(self, sched, specialty, id, appointments, net, daily_capacity):
        self.scheduling = sched
        self.specialty = specialty
        self.id = id
        self.appointments = appointments
        self.data=[]
        self.network=net
        self.daily_capacity= daily_capacity

    def log_schedule(self, env, patient, scheduling_failed, lag_time, delay, max_date, appt_date):
        call_date=env.now
        window_len = max_date - int(env.now + 1)
        day_full = self.appointments[appt_date].count
        day_remain = self.daily_capacity - day_full
        self.data.append((call_date, patient.id, patient.pt_class, scheduling_failed, lag_time, delay, window_len, max_date,appt_date,day_full,day_remain))

    def render_data(self, verbose=False):
        scheduler_data = pd.DataFrame(self.data, columns=["call_date","patient_id","patient_class","failed_schedule", "lag_time","delay","window_len","max_date","appt_date","date_filled","date_remaining"])
        scheduler_data["specialty"]=self.specialty
        scheduler_data["specialist_id"]=self.id
        return scheduler_data

    def graph_capacity(self):
        scheduler_data = self.render_data(False)
        class_colors = ['green' if patient_class == 0 else 'orange' if patient_class==1 else 'red' for patient_class in scheduler_data["patient_class"]]
        plt.scatter(scheduler_data["call_date"],scheduler_data["pct_10d_available"], color=class_colors,alpha=0.5)

    def graph_delays(self):
        scheduler_data = self.render_data(False)
        class_colors = ['green' if patient_class == 0 else 'blue' if patient_class==1 else 'red' for patient_class in scheduler_data["patient_class"]]
        plt.scatter(scheduler_data["call_date"],scheduler_data["lag_time"], color=class_colors,alpha=0.5)


class HealthNetwork(object):
    # Dict of list[Specialty] objects by type
    specialists = {}
    patient_vec = []

    def __init__(self, env):
        i=0
        for specialty in SPECIALTIES:
            self.specialists[specialty]=[]
            capacity= SPECIALTY_CAPS[SPECIALTIES.index(specialty)]
            for n in range(SPECIALTY_COUNTS[i]):
                self.specialists[specialty].append(
                Specialty(
                    sched = simpy.Resource(env, capacity=1),
                    specialty = specialty,
                    id = n,
                    appointments = [copy.copy(simpy.Resource(env, capacity=capacity)) for x in range(SIM_TIME + 14)],
                    net = self,
                    daily_capacity=capacity
                ))
            i+=1

    def populate_patients(self, env):
        for i in range(PATIENTS):

            cumul_class_probs = np.cumsum(PATIENT_CLASS_PROBS/np.sum(PATIENT_CLASS_PROBS))
            u1_class = unigens["init-class"].RandU01()
            patient_class = np.argmin(np.where(cumul_class_probs - u1_class < 0, 2, cumul_class_probs))


            cumul_specialties = np.cumsum(TRANSITION_MATRIX[patient_class,0])
            u1_specialty = unigens["init-specialty"].RandU01()
            specialty = SPECIALTIES[np.argmin(np.where(cumul_specialties - u1_specialty < 0, 2, cumul_specialties))]


            specialist = np.random.choice(self.specialists[specialty])


            pat = Patient(specialist, i, patient_class)
            self.patient_vec.append(pat)
            env.process(pat.patient(env))



class Patient():
    current_specialist: Specialty
    id: int
    data: []
    last_appt: float
    pt_class: int #[0,1,2] -> ["A","B","C"]
    transition_time: float

    def __init__(self, specialist, id, pt_class):
        self.current_specialist = specialist
        self.referring_specialist = specialist
        self.id = id
        self.pt_class = pt_class
        self.data = []

    def patient(self, env):
        # Wait for scheduling call
        self.last_appt = 0
        while True:
            # Patient requests an appointment
            print_v(3, f"Patient {self.id}'s Last Appointment: {self.last_appt}")
            with self.current_specialist.scheduling.request() as sched_appt:
                yield sched_appt
                # Logging
                req_date = env.now
                print_v(2, f"Patient {self.id} requests appointment at {self.current_specialist.specialty} on {env.now}")
            sched_lag = max(env.now - self.last_appt,0)
            print_v(3, f'Scheduling Lag: {sched_lag}')

            # Once the scheduler is available, the patient schedules the appointment
            patient_flexibility = min(sched_lag//5,7)
            appt,wait_time,delay = schedule(self.current_specialist, patient_flexibility, self, env)
            print_v(3, f"Wait time: {wait_time}\n")
            # Patient waits for the appointment
            yield env.timeout(wait_time)
            # Patient attends appointment
            yield appt
            self.last_appt = env.now
            print_v(2, f"Patient {self.id} attends appointment {self.current_specialist.specialty} at {env.now}")
            # Patient calculates their next specialty type and gets a referral
            new_specialist,self.transition_time = refer(self.current_specialist, self)
            print_v(4, f"Next Specialty: {new_specialist.specialty}")
            print_v(3, f"Transition Time: {self.transition_time}\n")
            yield env.timeout(self.transition_time)
            # Patient waits for the referral to process before returning to the beginning
            self.log_visit(req_date, self.last_appt, wait_time, delay, self.current_specialist.specialty, self.referring_specialist.specialty, self.transition_time, patient_flexibility)


            self.referring_specialist = self.current_specialist
            self.current_specialist = new_specialist


    def log_visit(self, req_date, appt_date, wait, delay, specialist, referrer, transition_time, flexibility):
        self.data.append((req_date, appt_date, wait, delay, specialist, referrer, transition_time, flexibility))

    def render_data(self, verbose=False):
        patient_data = pd.DataFrame(self.data, columns=["req",'appt','wait','delay','specialty','referring','transition','flex'])
        patient_data["patient_id"]=self.id
        patient_data["patient_class"]=self.pt_class
        if verbose:
            print_v(2, patient_data[["appt","wait","delay","flex","specialty","referring"]])
        return patient_data


def schedule(specialist, flex, patient, env, scheduler=SCHEDULER):
    if specialist.specialty=="Home":
        date = int(env.now + 1)
        # auto-releasing request to unbound home capacity
        specialist.log_schedule(env, patient, False, 0, 0, date,date)
        with specialist.appointments[date].request() as appointment:
            return (appointment, max(date - env.now,0), 0)
    min_date = int(env.now + 1)
    max_date = min(int(env.now + 2 + 2*flex), SIM_TIME+13)
    calendar = specialist.appointments[min_date:max_date]
    capacity = [appointment.count for appointment in calendar]
    date = min_date
    failed_schedule = False
    delaying = False
    if(env.now > 1):
        print_v(3,f"Flexibility: {flex}")
        print_v(4,f"Capacity: {capacity}")
    if scheduler == "first_min":
             date_raw = np.argmin(capacity)
             date = date_raw + int(env.now) + 1
             if capacity[date_raw] >= specialist.daily_capacity: # if appointment is full
                delaying = True
             else:
                 print_v(3,f"Date Selected: {date}")
    elif scheduler == "last_min":
             date_raw = len(capacity) - np.argmin(capacity[::-1])
             date = date_raw + int(env.now)

             if capacity[date_raw] >= specialist.daily_capacity: # if appointment is full
                delaying = True
             else:
                 print_v(3,f"Date Selected: {date}")
    elif scheduler == "uniform":
        for attempt in range(min(int(flex*4),4)): # if four times the flexibility window can't find an appointment there ain't one
            date = np.random.randint(int(env.now) + 1,max(int(env.now)+flex*2 + 1, SIM_TIME+13))
            if specialist.appointments[date].count < specialist.daily_capacity:
                print_v(3,f"Date Selected: {date}")
                break
        delaying = True
    else:
        print_v(2,f"Scheduler {scheduler} is unrecognized.")

    delay = 0

    if delaying:
        failed_schedule = True
        delayed_date = max_date

        while delayed_date < SIM_TIME + 13:
            if specialist.appointments[delayed_date].count < specialist.daily_capacity:
                date = delayed_date
                delay = delayed_date - max_date
                print_v(3, f"Date Selected: {date} with Delay {delay}.")
                break
            else:
                delayed_date += 1
        if delayed_date == SIM_TIME + 13:
            date = delayed_date
            delay = delayed_date - max_date
        print_v(3, f"Date Selected: {date} with Delay {delay}.")


    lag_time = max(date - env.now, 0)
    appointment = specialist.appointments[date].request() # NOT releasing the normal way
    if(env.now>1):
        print_v(4,f"Capacity - New: {[appointment.count for appointment in calendar]}")

    specialist.log_schedule(env, patient, failed_schedule, lag_time, delay, max_date, date)

    return (appointment, lag_time, delay)

def refer(specialist, patient) -> (Specialty, float):
    network=specialist.network
    #given an input specialty and a patient class, generate a new specialty
    specialty_id = SPECIALTIES.index(specialist.specialty)
    if(patient.pt_class > 2):
        print(patient.pt_class)
    trans_probs = TRANSITION_MATRIX[patient.pt_class,specialty_id] # [17, ]
    # Making sure to normalize, just in case
    cumul_trans_probs = np.cumsum(trans_probs/np.sum(trans_probs))
    u1_refer = unigens["refer"].RandU01()
    new_specialty_id = np.argmin(np.where(cumul_trans_probs - u1_refer < 0, 2, cumul_trans_probs))

    u1_trans = unigens["trans_exp"].RandU01()
    tau_ij = TIME_LAG_MATRIX[patient.pt_class,specialty_id,new_specialty_id]
    # Generating an exponential with mean tau_ij
    wait_time = -tau_ij * np.log(1-u1_trans)
    new_specialty = SPECIALTIES[new_specialty_id]
    #TODO choice -> rngstream
    new_specialist = np.random.choice(network.specialists[new_specialty])
    return(new_specialist, wait_time)


def main():
    # Try to use first CLI as the filename stem and second as debug, default otherwise
    try:
        filestub = sys.argv[1]
    except:
        filestub = "HospitalSim"

    try:
        debug = int(sys.argv[2])
        DEBUG_LEVEL=debug
    except:
        print_v(1,"Debug level unspecified, defaulting to %d" % DEBUG_LEVEL)
    # Initialize and run the simulation
    start=time.time()
    env = simpy.Environment()
    network = HealthNetwork(env)
    network.populate_patients(env)
    poptime = time.time() - start
    print_v(1,"Populated the network in %f seconds" % poptime)

    start=time.time()
    env.run(until=SIM_TIME)
    runtime=time.time() - start
    print_v(1,"Simulation run completed in %f seconds" % runtime)

    # Render the internal lists into data frames
    start=time.time()
    patient_data = pd.concat([patient.render_data() for patient in network.patient_vec])
    print(patient_data.groupby('patient_class')['delay'].mean())
    scheduler_data = pd.concat([specialist.render_data() for specialist_list in network.specialists.values() for specialist in specialist_list])
    rendertime = time.time() - start
    print_v(1,"Data rendering complete in %f seconds" % rendertime)


    # DATA WRITE OUT
    start=time.time()
    filename_stub = f"{filestub}T{SIM_TIME}N{PATIENTS}Trial{TRIAL}"
    patient_data.to_csv("data/"+filename_stub+"-Patient.csv")
    scheduler_data.to_csv("data/"+filename_stub+"-Scheduler.csv")
    writetime=time.time() - start
    print_v(1, "Data write complete in %f seconds" %writetime)

if __name__ == '__main__':
    main()
