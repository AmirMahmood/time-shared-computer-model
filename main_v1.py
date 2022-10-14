"""
TIME-SHARED COMPUTER MODEL
Simulation Modeling and Analysis 5th Edition

@author: AmirMahmood Ahrar
"""

import queue
import utils


# Define limits.
MAX_SVAR = 25  # Max number of sampst variables.
TIM_VAR  = 25  # Max number of timest variables.
MAX_TVAR = 50  # Max number of timest variables + lists.

# Define array sizes.
SVAR_SIZE = 26  # MAX_SVAR + 1.
TVAR_SIZE = 51  # MAX_TVAR + 1.

# Pre-define attribute numbers of transfer for event list.
EVENT_TIME = 0  # Attribute 0 in event list is event time.
EVENT_TYPE = 1  # Attribute 1 in event list is event type.

# Define some other values.
INFINITY = 1.E30  # Not really infinity, but a very large number.


################
STREAM_THINK = 1  # Random-number stream for think times.
STREAM_SERVICE = 2  # Random-number stream for service times.

EVENT_ARRIVAL = 1  # Event type for arrival of job to CPU.
EVENT_END_CPU_RUN = 2  # Event type for end of a CPU run.
EVENT_END_SIMULATION = 3  # Event type for end of the simulation.

SAMPST_RESPONSE_TIMES = 1  # sampst variable for response times.

_LIST_QUEUE = 1
_LIST_CPU = 2
_LIST_EVENTS = 25


def timest(value: float, variable: int) -> float:
    """
    Initialize, update, or report statistics on continuous-time processes:
    integral/average, max(default - 1E30), min(default 1E30)
    for timest variable "variable", where "variable":
        = 0 initializes counters
        > 0 updates area, min, and max accumulators with new level of variable
        < 0 reports stats on variable "variable" and returns them in transfer:
            [1] = time-average of variable updated to the time of this call
            [2] = maximum value variable has attained
            [3] = minimum value variable has attained
    Note that variables TIM_VAR + 1 through TVAR_SIZE are used for automatic
    record keeping on the length of lists 1 through MAX_LIST.
    """

    # static double area[TVAR_SIZE], max[TVAR_SIZE], min[TVAR_SIZE], preval[TVAR_SIZE], tlvc[TVAR_SIZE], treset;
    timest.area = getattr(timest, 'area', [0.0] * TVAR_SIZE)
    timest._max = getattr(timest, '_max', [0.0] * TVAR_SIZE)
    timest._min = getattr(timest, '_min', [0.0] * TVAR_SIZE)
    timest.preval = getattr(timest, 'preval', [0.0] * TVAR_SIZE)
    timest.tlvc = getattr(timest, 'tlvc', [0.0] * TVAR_SIZE)
    timest.treset = getattr(timest, 'treset', 0.0)

    # If the variable value is improper, stop the simulation.
    if (not(variable >= -MAX_TVAR) and (variable <= MAX_TVAR)):
        raise Exception(f"\n{variable} is an improper value for a timest variable at time {sim_time}\n")

    # Execute the desired option.

    if (variable > 0):
        # Update.
        timest.area[variable] += (sim_time - timest.tlvc[variable]) * timest.preval[variable]
        if (value > timest._max[variable]):
            timest._max[variable] = value
        if (value < timest._min[variable]):
            timest._min[variable] = value
        timest.preval[variable] = value
        timest.tlvc[variable] = sim_time
        return 0.0

    if (variable < 0):
        # Report summary statistics in transfer.
        transfer = [0] * 4
        ivar = -variable
        timest.area[ivar] += (sim_time - timest.tlvc[ivar]) * timest.preval[ivar]
        timest.tlvc[ivar] = sim_time
        transfer[1] = timest.area[ivar] / (sim_time - timest.treset)
        transfer[2] = timest._max[ivar]
        transfer[3] = timest._min[ivar]
        return transfer[1]

    # Initialize the accumulators.
    for ivar in range(1, MAX_TVAR+1):
        timest.area[ivar] = 0.0
        timest._max[ivar] = -INFINITY
        timest._min[ivar] = INFINITY
        timest.preval[ivar] = 0.0
        timest.tlvc[ivar] = sim_time

    timest.treset = sim_time


def filest(_list: int) -> float:
    """
    Report statistics on the length of list "list" in transfer:
        [1] = time-average of list length updated to the time of this call
        [2] = maximum length list has attained
        [3] = minimum length list has attained
    This uses timest variable TIM_VAR + list.
    """

    return timest(0.0, -(TIM_VAR + _list))


def sampst(value: float, variable: int) -> float:
    """ 
    Initialize, update, or report statistics on discrete-time processes:
    sum/average, max (default -1E30), min (default 1E30), number of observations
    for sampst variable "variable", where "variable":
        = 0 initializes accumulators
        > 0 updates sum, count, min, and max accumulators with new observation
        < 0 reports stats on variable "variable" and returns them in transfer:
            [1] = average of observations
            [2] = number of observations
            [3] = maximum of observations
            [4] = minimum of observations         
    """
    # static variable definition
    sampst.ivar = getattr(sampst, 'ivar', 0)
    sampst.num_observations = getattr(sampst, 'num_observations', [0] * SVAR_SIZE)

    sampst._max = getattr(sampst, '_max', [0.0] * SVAR_SIZE)
    sampst._min = getattr(sampst, '_min', [0.0] * SVAR_SIZE)
    sampst._sum = getattr(sampst, '_sum', [0.0] * SVAR_SIZE)

    # If the variable value is improper, stop the simulation.
    if (not (variable >= -MAX_SVAR) and (variable <= MAX_SVAR)):
        raise Exception(f"\n{variable} is an improper value for a sampst variable at time {sim_time}\n")

    # Execute the desired option.

    if variable > 0:
        # Update.
        sampst._sum[variable] += value
        if value > sampst._max[variable]:
            sampst._max[variable] = value
        if value < sampst._min[variable]:
            sampst._min[variable] = value
        sampst.num_observations[variable] += 1
        return 0.0

    # Report summary statistics in transfer.
    if (variable < 0):
        transfer = [0] * 5
        sampst.ivar = -variable
        transfer[2] = float(sampst.num_observations[sampst.ivar])
        transfer[3] = sampst._max[sampst.ivar]
        transfer[4] = sampst._min[sampst.ivar]
        if sampst.num_observations[sampst.ivar] == 0:
            transfer[1] = 0.0
        else:
            transfer[1] = sampst._sum[sampst.ivar] / transfer[2]
        return transfer[1]

    # Initialize the accumulators.
    for sampst.ivar in range(1, MAX_SVAR+1):
        sampst._sum[sampst.ivar] = 0.0
        sampst._max[sampst.ivar] = -INFINITY
        sampst._min[sampst.ivar] = INFINITY
        sampst.num_observations[sampst.ivar] = 0

    return 0.0


def add_to_list_event(transfer):
    LIST_EVENTS.append(transfer)
    LIST_EVENTS.sort()
    timest(float(len(LIST_EVENTS)), TIM_VAR + _LIST_EVENTS)


def event_schedule(time_of_event: float, type_of_event: int) -> None:
    """
    Schedule an event at time event_time of type event_type.  If attributes
    beyond the first two (reserved for the event time and the event type) are
    being used in the event list, it is the user's responsibility to place
    their values into the transfer array before invoking event_schedule.
    """
    add_to_list_event([time_of_event, type_of_event])


def timing():
    """
    Remove next event from event list, placing its attributes in transfer.
    Set sim_time (simulation time) to event time, transfer[1].
    Set next_event_type to this event type, transfer[2].
    """
    global sim_time

    # Remove the first event from the event list and put it in transfer[].
    transfer = LIST_EVENTS.pop(0)
    timest(float(len(LIST_EVENTS)), TIM_VAR + _LIST_EVENTS)

    # Check for a time reversal.
    if transfer[EVENT_TIME] < sim_time:
        raise Exception(
            f"\nAttempt to schedule event type {transfer[EVENT_TYPE]} for "
            f"time {transfer[EVENT_TIME]} at time {sim_time}\n"
        )

    # Advance the simulation clock and set the next event type.
    sim_time = transfer[EVENT_TIME]
    return transfer[EVENT_TYPE]


def start_CPU_run():
    """ Non-event function to start a CPU run of a job. """

    # Remove the first job from the queue.
    transfer = LIST_QUEUE.get()
    timest(float(LIST_QUEUE.qsize()), TIM_VAR + _LIST_QUEUE)

    # Determine the CPU time for this pass, including the swap time.
    if (quantum < transfer[1]):
        run_time = quantum + swap
    else:
        run_time = transfer[1] + swap

    # Decrement remaining CPU time by a full quantum.  (If less than a full
    # quantum is needed, this attribute becomes negative.  This indicates that
    # the job, after exiting the CPU for the current pass, will be done and is
    # to be sent back to its terminal.)
    transfer[1] -= quantum

    # Place the job into the CPU.
    LIST_CPU.put(transfer)
    timest(float(LIST_CPU.qsize()), TIM_VAR + _LIST_CPU)

    # Schedule the end of the CPU run.
    event_schedule(sim_time + run_time, EVENT_END_CPU_RUN)


def arrive():
    LIST_QUEUE.put([sim_time, rnd.expon(mean_service, STREAM_SERVICE)])
    timest(float(LIST_QUEUE.qsize()), TIM_VAR + _LIST_QUEUE)

    # If the CPU is idle, start a CPU run.
    if LIST_CPU.empty():
        start_CPU_run()


def end_CPU_run():
    """ Event function to end a CPU run of a job. """

    global num_responses

    # Remove the job from the CPU.
    transfer = LIST_CPU.get()
    timest(float(LIST_CPU.qsize()), TIM_VAR + _LIST_CPU)

    # Check to see whether this job requires more CPU time.
    if (transfer[1] > 0.0):

        # This job requires more CPU time, so place it at the end of the queue
        # and start the first job in the queue.
        LIST_QUEUE.put(transfer)
        timest(float(LIST_QUEUE.qsize()), TIM_VAR + _LIST_QUEUE)
        start_CPU_run()
    else:

        # This job is finished, so collect response-time statistics and send it
        # back to its terminal, i.e., schedule another arrival from the same
        # terminal.
        sampst(sim_time - transfer[0], SAMPST_RESPONSE_TIMES)

        event_schedule(sim_time + rnd.expon(mean_think, STREAM_THINK), EVENT_ARRIVAL)

        # Increment the number of completed jobs.
        num_responses += 1

        # Check to see whether enough jobs are done.
        if num_responses >= num_responses_required:

            # Enough jobs are done, so schedule the end of the simulation
            # immediately (forcing it to the head of the event list).
            event_schedule(sim_time, EVENT_END_SIMULATION)

        else:
            # Not enough jobs are done; if the queue is not empty, start
            # another job.
            if not LIST_QUEUE.empty():
                start_CPU_run()


def report():
    """Report generator function."""

    # Get and write out estimates of desired measures of performance.
    outfile.write(
        f"\n\n{num_terms:5d}{sampst(0.0, -SAMPST_RESPONSE_TIMES):16.3f}{filest(_LIST_QUEUE):16.3f}{filest(_LIST_CPU):16.3f}"
    )


# random number generation
rnd = utils.RandomNumberGenerator()

# Open output files.
outfile = open("tscomp.out", "w")

# Open input file & read input parameters.
with open("tscomp.in", "r") as infile:
    min_terms, max_terms, incr_terms, num_responses_required, mean_think, mean_service, quantum, swap \
        = infile.readline().split()
min_terms = int(min_terms)
max_terms = int(max_terms)
incr_terms = int(incr_terms)
num_responses_required = int(num_responses_required)
mean_think = float(mean_think)
mean_service = float(mean_service)
quantum = float(quantum)
swap = float(swap)

# Write report heading and input parameters.
outfile.write("Time-shared computer model\n\n")
outfile.write(f"Number of terminals{min_terms:9d} to{max_terms:4d} by {incr_terms:4d}\n\n")
outfile.write(f"Mean think time  {mean_think:11.3f} seconds\n\n")
outfile.write(f"Mean service time{mean_service:11.3f} seconds\n\n")
outfile.write(f"Quantum          {quantum:11.3f} seconds\n\n")
outfile.write(f"Swap time        {swap:11.3f} seconds\n\n")
outfile.write(f"Number of jobs processed{num_responses_required:12d}\n\n\n")
outfile.write("Number of      Average         Average")
outfile.write("       Utilization\n")
outfile.write("terminals   response time  number in queue     of CPU")

# Run the simulation varying the number of terminals.
for num_terms in range(min_terms, max_terms + 1, incr_terms):

    # Initialize simulation environment
    LIST_EVENTS = []
    LIST_QUEUE = queue.Queue()
    LIST_CPU = queue.Queue(maxsize=1)
    sim_time = 0.0
    sampst(0.0, 0)
    timest(0.0, 0)

    # Initialize the non-simulation statistical counter.
    num_responses = 0

    # Schedule the first arrival to the CPU from each terminal.
    for term in range(1, num_terms + 1):
        event_schedule(rnd.expon(mean_think, STREAM_THINK), EVENT_ARRIVAL)

    # Run the simulation until it terminates after an end-simulation event
    # (type EVENT_END_SIMULATION) occurs.
    while True:

        # Determine the next event.
        next_event_type = timing()

        # Invoke the appropriate event function.
        events_switch = {
            EVENT_ARRIVAL: arrive,
            EVENT_END_CPU_RUN: end_CPU_run,
            EVENT_END_SIMULATION: report
        }
        events_switch[next_event_type]()

        # If the event just executed was not the end-simulation event
        # (type EVENT_END_SIMULATION), continue simulating.  Otherwise, end
        # the simulation.
        if next_event_type == EVENT_END_SIMULATION:
            break

outfile.close()
