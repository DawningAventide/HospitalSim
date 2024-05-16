import simpy
import numpy as np
# Make sure rngstream.so is in the folder if getting import issues
import rngStream
import copy
import pandas as pd
import matplotlib.pyplot as plt
import sys
import time

# Duration of the simulation in days (5-day week, 250d year)
SIM_TIME = 2500
# Number of patients
PATIENTS = 100000

# Scheduling heuristic (first_min or greedy recommended)
# Descriptions in original paper
SCHEDULER = "first_min"
SCHEDULERS = ['first_min','last_min','uniform','greedy']
# Specialty names - do not change, order-sensitive
SPECIALTIES = ["Home","Primary Care","Ophthalmology","Immunology", \
    "Cardiology","Other Dr Specialty","Urology","Gastroenterology", \
    "General Surgery","OB/GYN","Nephrology","Orthopedics", \
    "Psychiatry","Neurology","Oncology","Dermatology", "Otorhinolaryngology"]
# Count of each specialty to instantiate - based on rates per 100k patients

#manually optimized
#SPECIALTY_COUNTS=[1,120,14,2, 7,42,2,3, 2,10,2,7, 16,5,4,4,2]

# default
SPECIALTY_COUNTS=[1,99,6,2, 10,25,3,5, 8,13,4,8, 15,5,8,4,3]
# +2 Ophthalmologists
#SPECIALTY_COUNTS=[1,99,8,2,10,25,3,5,8,13,4,8,15,5,8,4,3]

# Daily capacity of each specialty, derived from NAMCS 2016
SPECIALTY_CAPS=[1,12,12,12,12,12,16,14,13,16,14,16,8,8,12,16,16]

# UNCAP settings - uncomment to set uncapacitated:
#SPECIALTY_CAPS=[99999]*17
#SPECIALTY_COUNTS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

# Loyalty Settings
# Set turnover_base >= 1 to disable preferential referral and randomly assign all.
# Yearly expected doctor turnover
GLOBAL_TURNOVER = 0.06
# Per-appointment turnover chance
TURNOVER_BASE = 0.01
#TURNOVER_BASE=1
# Rate modifiers on f(doctor-related turnover,t)
SPECIALIST_TURNOVER_MOD = [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]


# Verbosity, controls print_v calls
DEBUG_LEVEL = 1

# Set a specific replicable seed if True
TRIAL = False

# Number of progress updates to post
# (ideally divides SIM_TIME but you can be messy, it won't break,
#  I'll just judge you a little)
PROGRESS_BARS = 100

#===================================================================#
# END OF PARAMETER ZONE - PROCEED AT OWN RISK
#
#===================================================================#

# Initializing a few global counters
PROG_ITER = 0
LAST_CHKPT = 0
RUN_START = 0
# Initializing probabilistic matrices
TRANSITION_MATRIX = np.zeros([3,17,17])
TIME_LAG_MATRIX = np.zeros([3,17,17])
# Setting up patient classes (See original thesis for numerical breaks)
# TLDR: A- fewest chronic conditions, B- in between, C- most
PATIENT_CLASSES = ['A','B','C']
PATIENT_CLASS_PROBS = [.7054,.1872,.1074]

# Reading in transition parameters
trans_mat_raw = pd.read_csv("parameters/mrp_trans_probs.csv", sep=',',header=0,index_col='from_desc')
mtd_mat_raw = pd.read_csv("parameters/mrp_mtds.csv", sep=",",header=0,index_col='from_desc')
for i in range(3):
    TRANSITION_MATRIX[i] = trans_mat_raw.loc[trans_mat_raw["conditions"] == PATIENT_CLASSES[i]].iloc[:,1:]
    TIME_LAG_MATRIX[i] = mtd_mat_raw.loc[mtd_mat_raw["conditions"]==PATIENT_CLASSES[i]].iloc[:,1:]

# RNGSTREAM setup
if TRIAL:
    # adjust the seed if you want different replicable runs.
    # you know how seeds work. I hope.
    seed0 = [12345,23456,23456,23456,23456,23456]
    rngStream.SetPackageSeed(seed0)
unigens = {}
unigens["refer"] = rngStream.RngStream()
unigens["trans-exp"] = rngStream.RngStream()
unigens["init-class"] = rngStream.RngStream()
unigens["init-specialty"] = rngStream.RngStream()
unigens["init-specialist"] = rngStream.RngStream()
unigens["refer-specialist"] = rngStream.RngStream()
unigens["pref-referral"] = rngStream.RngStream()
unigens["uniform-choice"] = rngStream.RngStream()


# Globally tunable debug levels
def print_v(debug, *args):
    if debug <= DEBUG_LEVEL:
        print(*args)

# np.choice-like wrapper for RngStream for replicability
def rng_choice(generator,array):
    u = generator.RandU01()
    l = len(array)
    return array[int(u*l)]

# Specialist instance class
class Specialty(object):
    # Virtual point-time queue for scheduling calls
    scheduling: simpy.Resource
    # Type and numeric ID
    specialty: str
    id: int
    # List of resources representing calendar appointments
    appointments: [simpy.Resource]
    daily_capacity: int

    # Internal data buffer
    data: []

    def __init__(self, sched, specialty, id, appointments, net, daily_capacity):
        self.scheduling = sched
        self.specialty = specialty
        self.id = id
        self.appointments = appointments
        self.data=[]
        self.network=net
        self.daily_capacity= daily_capacity

    # Write scheduling details to internal buffer
    def log_schedule(self, env, patient, scheduling_failed, lag_time, delay, max_date, appt_date, first_avail, third_avail):
        call_date=env.now
        window_len = max_date - int(env.now + 1)
        day_full = self.appointments[appt_date].count
        day_remain = self.daily_capacity - day_full
        self.data.append((call_date, patient.id, patient.pt_class, scheduling_failed, lag_time, delay, window_len, max_date,appt_date,day_full,day_remain,first_avail,third_avail))

    # Render internal data buffer into a pandas dataframe
    def render_data(self, verbose=False):
        scheduler_data = pd.DataFrame(self.data, columns=["call_date","patient_id","patient_class","failed_schedule", "lag_time","delay","window_len","max_date","appt_date","date_filled","date_remaining","first_available","third_available"])
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

# Wrapper class that instantiates and manages specialties and patients
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

    # Create and initialize patients into the system
    def populate_patients(self, env):
        for i in range(PATIENTS):

            cumul_class_probs = np.cumsum(PATIENT_CLASS_PROBS/np.sum(PATIENT_CLASS_PROBS))
            u1_class = unigens["init-class"].RandU01()
            patient_class = np.argmin(np.where(cumul_class_probs - u1_class < 0, 2, cumul_class_probs))


            cumul_specialties = np.cumsum(TRANSITION_MATRIX[patient_class,0])
            u1_specialty = unigens["init-specialty"].RandU01()
            specialty = SPECIALTIES[np.argmin(np.where(cumul_specialties - u1_specialty < 0, 2, cumul_specialties))]


            specialist = rng_choice(unigens["init-specialist"],self.specialists[specialty])


            pat = Patient(specialist, i, patient_class)
            self.patient_vec.append(pat)
            env.process(pat.patient(env))


# Individual patient class
class Patient():
    current_specialist: Specialty
    id: int
    data: []
    last_appt: float
    pt_class: int #[0,1,2] -> ["A","B","C"]
    transition_time: float

    # Known practicioners
    care_team: {}

    def __init__(self, specialist, id, pt_class):
        self.current_specialist = specialist
        self.referring_specialist = specialist
        self.id = id
        self.pt_class = pt_class
        self.data = []
        # specialist id, date of last contact
        self.care_team = {sp:(-1,-1) for sp in SPECIALTIES}

    # Main patient loop - reference SimPy documentation for help
    def patient(self, env):

        # Wait for scheduling call
        self.last_appt = 0
        while True:
            global PROG_ITER
            global LAST_CHKPT
            #If we've hit a checkpoint, note it
            if(env.now > 0 and int(env.now)//(SIM_TIME//PROGRESS_BARS) >= PROG_ITER):
                PROG_ITER += 1
                now = time.time()
                dt = now-LAST_CHKPT
                elapsed=now-RUN_START
                LAST_CHKPT=now
                print_v(1,"%d/%d: %.3f sec, %.3f sec elapsed." %(PROG_ITER,PROGRESS_BARS,dt,elapsed))

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
            patient_flexibility = int(min(sched_lag//5,7))
            appt,wait_time,delay = schedule(self.current_specialist, patient_flexibility, self, env)
            print_v(3, f"Wait time: {wait_time}\n")
            # Patient waits for the appointment
            yield env.timeout(wait_time)
            # Patient attends appointment
            yield appt
            self.last_appt = env.now
            print_v(2, f"Patient {self.id} attends appointment {self.current_specialist.specialty} at {env.now}")
            # Patient calculates their next specialty type and gets a referral
            new_specialist,self.transition_time = refer(self.current_specialist, self, env)
            print_v(4, f"Next Specialty: {new_specialist.specialty}")
            print_v(3, f"Transition Time: {self.transition_time}\n")
            yield env.timeout(self.transition_time)
            # Patient waits for the referral to process before returning to the beginning

            # Log the visit details to internal data buffer
            self.log_visit(req_date, self.last_appt, wait_time, delay, self.current_specialist.specialty, self.referring_specialist.specialty, self.transition_time, patient_flexibility)


            self.referring_specialist = self.current_specialist
            self.current_specialist = new_specialist

    # Add visit data to internal storage
    def log_visit(self, req_date, appt_date, wait, delay, specialist, referrer, transition_time, flexibility):
        self.data.append((req_date, appt_date, wait, delay, specialist, referrer, transition_time, flexibility))

    # Write out the full patient trace for the individual
    def render_data(self, verbose=False):
        patient_data = pd.DataFrame(self.data, columns=["req",'appt','wait','delay','specialty','referring','transition','flex'])
        patient_data["patient_id"]=self.id
        patient_data["patient_class"]=self.pt_class
        if verbose:
            print_v(2, patient_data[["appt","wait","delay","flex","specialty","referring"]])
        return patient_data

# Schedules a patient with a specialist
def schedule(specialist, flex, patient, env, scheduler=SCHEDULER):
    if specialist.specialty=="Home":
        date = int(env.now + 1)
        # auto-releasing request to unbound home capacity
        specialist.log_schedule(env, patient, False, 0, 0, date,date,date,date)
        with specialist.appointments[date].request() as appointment:
            return (appointment, max(date - env.now,0), 0)
    min_date = int(env.now + 1)
    max_date = min(int(env.now + 2 + 2*flex), SIM_TIME+13)
    calendar = specialist.appointments[min_date:max_date]
    capacity = [appointment.count for appointment in calendar]
    capacity_unb = np.array([appointment.count for appointment in specialist.appointments[min_date:SIM_TIME+13]])
    first_avail = np.argmin(capacity_unb==specialist.daily_capacity) + min_date
    third_avail = np.argmax((np.cumsum(-1*(capacity_unb-specialist.daily_capacity)))>=3)+ min_date
    date = min_date
    failed_schedule = False
    delaying = False
    if(env.now > 1):
        print_v(3,f"Flexibility: {flex}")
        print_v(4,f"Capacity: {capacity}")
    if scheduler == "first_min":
             date_raw = np.argmin(capacity)
             date = date_raw + int(env.now) + 1
             if specialist.appointments[date].count >= specialist.daily_capacity: # if appointment is full
                delaying = True
             else:
                 print_v(3,f"Date Selected: {date}")
    elif scheduler == "last_min":
             date_raw = len(capacity) - np.argmin(capacity[::-1]) -1
             date = date_raw + int(env.now)

             if specialist.appointments[date].count >= specialist.daily_capacity: # if appointment is full
                delaying = True
             else:
                 print_v(3,f"Date Selected: {date}")
    elif scheduler == "greedy":
            date = first_avail
            if first_avail > max_date:
                delaying=True
            else:
                print_v(3,f"Date Selected:{date}")
    elif scheduler == "uniform":
        temp_calendar = np.array([a for a in range(min_date, min_date + len(capacity))])
        avail_win = temp_calendar[np.array(capacity)!=specialist.daily_capacity]
        if len(avail_win)!=0:
            date = rng_choice(unigens["uniform-choice"], avail_win)
        else:
            delaying = True

    # --- IF IMPLEMENTING NEW SCHEDULING HEURISTICS DO IT HERE ---

    else:
        print_v(0,f"Scheduler {scheduler} is unrecognized.")

    delay = 0

    # Get the first available
    if delaying:
        failed_schedule = True
        delayed_date = first_avail

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

    patient.care_team[specialist.specialty] = (specialist.id,date)
    specialist.log_schedule(env, patient, failed_schedule, lag_time, delay, max_date, date, first_avail, third_avail)

    return (appointment, lag_time, delay)

# Given an existing patient/specialist combo, refer to the next and calculate transition time
def refer(specialist, patient, env) -> (Specialty, float):
    network=specialist.network
    #given an input specialty and a patient class, generate a new specialty
    specialty_id = SPECIALTIES.index(specialist.specialty)
    trans_probs = TRANSITION_MATRIX[patient.pt_class,specialty_id] # [17, ]
    # Making sure to normalize, just in case
    cumul_trans_probs = np.cumsum(trans_probs/np.sum(trans_probs))
    u1_refer = unigens["refer"].RandU01()
    new_specialty_id = np.argmin(np.where(cumul_trans_probs - u1_refer < 0, 2, cumul_trans_probs))

    u1_trans = unigens["trans-exp"].RandU01()
    tau_ij = TIME_LAG_MATRIX[patient.pt_class,specialty_id,new_specialty_id]
    # Generating an exponential with mean tau_ij
    wait_time = -tau_ij * np.log(1-u1_trans)
    new_specialty = SPECIALTIES[new_specialty_id]

    # LOYAL/PREFERENTIAL REFERRAL:
    #if no cached specialist, refer as usual
    if(patient.care_team[new_specialty][0]>=0):
        #If cached specialist, autoselect with probability 1-f(time,spec)
        time_lag = int(env.now) + wait_time - patient.care_team[new_specialty][1]
        p_daily = 1-(1-GLOBAL_TURNOVER)**(wait_time/250)
        spec_mod = SPECIALIST_TURNOVER_MOD[new_specialty_id]
        transfer_prob = p_daily*spec_mod + TURNOVER_BASE
        u = unigens["pref-referral"].RandU01()
        if(u > transfer_prob):
            # Go back to the same specialist
            new_specialist = network.specialists[new_specialty][patient.care_team[new_specialty][0]]
        else:
            # clear the old specialist
            patient.care_team[new_specialty] = (-1,-1)
            new_specialist = rng_choice(unigens["refer-specialist"],network.specialists[new_specialty])
    else:
        # DEFAULT NON-PREFERENTIAL BEHAVIOR
        new_specialist = rng_choice(unigens["refer-specialist"],network.specialists[new_specialty])
    return(new_specialist, wait_time)


def main():
    # Try to use first CLI as the filename stem and second as debug, default otherwise
    global DEBUG_LEVEL
    global LAST_CHKPT
    global RUN_START
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
    LAST_CHKPT=start
    RUN_START=start
    env.run(until=SIM_TIME)
    runtime=time.time() - start
    print_v(1,"Simulation run completed in %f seconds" % runtime)

    # Render the internal lists into data frames
    start=time.time()
    patient_data = pd.concat([patient.render_data() for patient in network.patient_vec])
    #print(patient_data.groupby('patient_class')['delay'].mean())
    scheduler_data = pd.concat([specialist.render_data() for specialist_list in network.specialists.values() for specialist in specialist_list])
    util_tmp = scheduler_data.groupby(['specialty','appt_date'])['date_filled'].max()
    #print("-----------Utilization-----------")
    #print(util_tmp.groupby('specialty').mean()/util_tmp.groupby('specialty').max('appt_date'))
    print("----------- Delays -------------")
    print(scheduler_data.groupby('patient_class')['delay'].mean())
    rendertime = time.time() - start
    print_v(1,"Data rendering complete in %f seconds" % rendertime)


    # DATA WRITE OUT
    start=time.time()
    # feel free to modify the form to reflect the variables you intend to change often
    filename_stub = f"{filestub}T{SIM_TIME}N{PATIENTS}Trial{TRIAL}"
    patient_data.to_csv("data/"+filename_stub+"-Patient.csv")
    scheduler_data.to_csv("data/"+filename_stub+"-Scheduler.csv")
    writetime=time.time() - start
    print_v(1, "Data write complete in %f seconds" %writetime)

if __name__ == '__main__':
    main()
