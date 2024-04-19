# HospitalSim
GSMP simulation of an ambulatory specialty network serving a heterogeneous patient population for capacity planning.

Based on a Ph.D Thesis by Prashant Meckoni, "Capacity Planning in Heterogeneous Patient Populations." The original thesis described a simulation model based on longitudinal outpatient traces in the United States [recheck data years]. The model was largely focused on capacity planning for primary care providers, but the final chapter expanded the simulation to examine how patients move between various outpatient specialists. 

The model implemented here expands on that specialist model in several ways. Primarily, the original model does not cap the daily availability of any specialist, instead measuring the numbers of patients queued for each day as an approximation of the required capacity. Because of this, queues cannot build up as they often do in the real world and the amount of time spent delayed by insufficient specialist resources is difficult to assess. 
To more granularly model scheduling and delay, this simulation model introduces a set daily capacity for each specialty type and then creates instances of each specialty type to match the distributions of specialists per 100k population in the United States, per [grab the cite]. Both counts and daily capacities of specialties can be modified on a per-specialty basis, to more accurately track with the appointment lengths that different specialists tend towards. For example, Psychiatry and Orthopedics appointments tend to be longer than Primary Care appointments, and their daily capacities reflect that with help from appointment length data from NAMCS. For a more detailed discussion of this derivation, see the file 'namcs-analysis.Rmd' for the transformation code and an analysis of the decisions made during that transformation, as well as literature that analysis is based on.

The model is flexible, though, and uncapacitated behavior can be replicated simply by setting the vector of specialty caps to arbitrarily high values, as the commented code in lines 20 and 21 of the source file demonstrates. 

The second central innovation of the revised model is a preferential referral process, where patients will tend to default back to a specialist with whom they have a preexisting relationship. This was not relevant to the original, uncapacitated model, as there were only single instances of each specialist standing in for their fields. However, in the capacitated model, the assumption of independent random assignment of patients to specialists despite any existing relationship is not one that holds on the real-world system. Finding a new doctor is a pain, and it is rare that patients will go out of their way to switch providers unless given a reason. Importantly, though, this preferential behavior means that the actual specialist capacity required to satisfy an opinionated population may be higher than predicted, as this preferential selection can lead to waste of open appointments from non-preferred providers and increased delays at preferred ones. Patients' inital selection of a specialist from a field they have yet to contact remains random, as does reselection in the case that a patient decides to move providers. 
The behavior as implemented is a simple modification to the existing referral process, with the flexibility to model more complex trends if evidence is found for them. Unfortunately, hard numbers on patient-specialist discontinuity that can be split up by specialty type seem to be scarce. Due to this, though the model has the capability to adjust the rates of patient turnover by specialty type, it does not currently do so as no usable empirical data could be found that would provide a solid basis for those adjustments. Given this limitation, the model in its current form uses physician turnover rates between 2010 and 2020[^1] to model physician-related discontinuity as a probabilistic function of the time since the patient and physician were last in contact, derived from the average 6% yearly turnover rates during the relevant period of the cited study. Additionally, there is a fixed per-appointment probability (1% at current) that some patient-related externality triggers a search for a new physician. This rate is also adjustable, and is largely to acknowledge the amount of complexity and unpredictability there can be in navigating physician-patient relationships but also functions as a floor probability of changing providers that ensures patients move around at least somewhat, ensuring slightly more flow in the model that compensates for the worst of the preferential impacts. 

## Goals

The purpose of the model is to simulate patient flow of a heterogeneous population among a heterogeneous specialty network in a reproducible, flexible way and enable capacity planning and analysis on the network. This can be done in multiple ways: given a patient population and an existing capacity, the model can be used to simulate the effects of adjustments to that capacity on the delays and wait times experienced by patients in the system and expose any differential impacts on different classes of patients by those adjustments. On a larger scale, though, this kind of model would be usable for joint capacity optimization given an existing patient population to estimate the distribution of physicians across specialties that would be needed to meet the population's needs.

A secondary goal, though largely adjunct to the first, is to record traces of the full simulation such that the state of the system at any arbitrary point in the simulation run can be fully reconstructed from the data trace. This trace is written to disk as two .csv files, with a name prefix given optionally as a command line input and a summary of the most critical parameters to ensure traces are easily identifiable and distinct from each other. 

### Running the Model

The simulation is built to be run from the command line, but all non-verbosity/naming parameters must be adjusted within the file. 

Command line parameters are as follows:

`python HospitalSimulation.py [filenameprefix] [debug level 0-5]` 

i.e. `python HospitalSimulation.py uncapRun2 1`. 

Debug level 0 prints nothing, 1 shows percentile progress and timing information, 2 shows that AND patient-attend/patient-schedule events, and 3-5 show increasingly granular information as patients move. Levels 3-5 are not recommended for larger runs due to the performance hit and the console spam, but they are useful if you need specific debug information and are not reaching the trace file writes at the very end of the simulation.

All in-file parameter uses are documented at the top of the file in global declarations. 

A data dictionary for the trace format is coming soon. 


[^1]:Amelia M. Bond, Lawrence P. Casalino, Ming Tai-Seale, et al. Physician Turnover in the United States. Ann Intern Med.2023;176:896-903. [Epub 11 July 2023]. doi:10.7326/M22-2504


