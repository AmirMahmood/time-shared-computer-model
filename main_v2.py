"""
TIME-SHARED COMPUTER MODEL
Simulation Modeling and Analysis 5th Edition

@author: AmirMahmood Ahrar
"""

from typing import List
import queue
import math
import argparse
from enum import Enum


class StatisticsMonitor:
    def __init__(self) -> None:
        self.sampst_num_observations = 0
        self.sampst_sum = 0.0

        self.timest_dict = {}

        self.timest_init()

    def timest_init(self):
        # Initialize the accumulators.
        for ivar in ['LIST_QUEUE', 'LIST_CPU']:
            self.timest_dict[ivar] = {'area': 0.0, 'preval': 0.0, 'tlvc': 0.0}

    def timest_update(self, value, variable, sim_time):
        """
        update statistics on continuous-time processes:
        """
        # Update.
        self.timest_dict[variable]['area'] += (sim_time - self.timest_dict[variable]['tlvc']) * \
            self.timest_dict[variable]['preval']
        self.timest_dict[variable]['preval'] = value
        self.timest_dict[variable]['tlvc'] = sim_time

    def timest_report(self, variable, sim_time) -> float:
        """
        report statistics on continuous-time processes:
        """
        # Report summary statistics in transfer.
        self.timest_dict[variable]['area'] += (sim_time - self.timest_dict[variable]
                                               ['tlvc']) * self.timest_dict[variable]['preval']
        self.timest_dict[variable]['tlvc'] = sim_time
        return self.timest_dict[variable]['area'] / sim_time

    def sampst_update(self, value: float) -> float:
        """ 
        update statistics on discrete-time processes:
        """
        self.sampst_sum += value
        self.sampst_num_observations += 1

    def sampst_report(self) -> float:
        """ 
        report statistics on discrete-time processes:
        """

        # Report summary statistics in transfer.
        if self.sampst_num_observations == 0:
            return 0.0
        else:
            return self.sampst_sum / float(self.sampst_num_observations)


class RandomNumberGenerator:
    """
    Prime modulus multiplicative linear congruential generator
    Z[i] = (630360016 * Z[i-1]) (mod(pow(2,31) - 1)), based on Marse and
    Roberts' portable FORTRAN random-number generator UNIRAN.  Multiple
    (100) streams are supported, with seeds spaced 100,000 apart.
    Throughout, input argument "stream" must be an int giving the
    desired stream number.
    Usage:
    To obtain the next U(0,1) random number from stream "stream,"
        execute
            u = lcgrand(stream);
        where lcgrand is a float function. The float variable u will
        contain the next random number.
    """

    def __init__(self) -> None:

        # Set the default seeds for all 100 streams.
        self.zrng = [
            1, 1973272912, 281629770, 20006270, 1280689831, 2096730329, 1933576050,
            913566091, 246780520, 1363774876, 604901985, 1511192140, 1259851944,
            824064364, 150493284, 242708531, 75253171, 1964472944, 1202299975,
            233217322, 1911216000, 726370533, 403498145, 993232223, 1103205531,
            762430696, 1922803170, 1385516923, 76271663, 413682397, 726466604,
            336157058, 1432650381, 1120463904, 595778810, 877722890, 1046574445,
            68911991, 2088367019, 748545416, 622401386, 2122378830, 640690903,
            1774806513, 2132545692, 2079249579, 78130110, 852776735, 1187867272,
            1351423507, 1645973084, 1997049139, 922510944, 2045512870, 898585771,
            243649545, 1004818771, 773686062, 403188473, 372279877, 1901633463,
            498067494, 2087759558, 493157915, 597104727, 1530940798, 1814496276,
            536444882, 1663153658, 855503735, 67784357, 1432404475, 619691088,
            119025595, 880802310, 176192644, 1116780070, 277854671, 1366580350,
            1142483975, 2026948561, 1053920743, 786262391, 1792203830, 1494667770,
            1923011392, 1433700034, 1244184613, 1147297105, 539712780, 1545929719,
            190641742, 1645390429, 264907697, 620389253, 1502074852, 927711160,
            364849192, 2049576050, 638580085, 547070247
        ]

    def lcgrand(self, stream: int) -> float:
        """ Generate the next random number. """

        # Define the constants.
        MULT1 = 24112
        MULT2 = 26143
        MODLUS = 2147483647

        zi = self.zrng[stream]
        lowprd = (zi & 65535) * MULT1
        hi31 = (zi >> 16) * MULT1 + (lowprd >> 16)
        zi = ((lowprd & 65535) - MODLUS) + ((hi31 & 32767) << 16) + (hi31 >> 15)
        if zi < 0:
            zi += MODLUS
        lowprd = (zi & 65535) * MULT2
        hi31 = (zi >> 16) * MULT2 + (lowprd >> 16)
        zi = ((lowprd & 65535) - MODLUS) + ((hi31 & 32767) << 16) + (hi31 >> 15)
        if zi < 0:
            zi += MODLUS
        self.zrng[stream] = zi
        return (zi >> 7 | 1) / 16777216.0

    def expon(self, mean: float, stream: int) -> float:
        """ Exponential variate generation function. """

        return -mean * math.log(self.lcgrand(stream))


class EventTypes(Enum):
    EVENT_ARRIVAL = 1
    EVENT_END_CPU_RUN = 3
    EVENT_END_SIMULATION = 2


class Event:
    def __init__(self, EVENT_TIME, EVENT_TYPE) -> None:
        self.EVENT_TIME = EVENT_TIME
        self.EVENT_TYPE = EVENT_TYPE


class Computer:
    def __init__(self, quantum, swap, environment, num_responses_required, STREAM_THINK, mean_service, mean_think) -> None:
        # Initialize CPU Lists
        self.quantum = quantum
        self.swap = swap
        self.environment = environment

        self.STREAM_THINK = STREAM_THINK
        self.STREAM_SERVICE = 2  # Random-number stream for service times.
        self.mean_service = mean_service
        self.mean_think = mean_think

        self.num_responses_required = num_responses_required

        # define the non-simulation statistical counter.
        self.num_responses = 0

        self.LIST_QUEUE = queue.Queue()
        self.LIST_CPU = queue.Queue(maxsize=1)

        self.stat_monitor = StatisticsMonitor()

    def arrive(self):
        self.LIST_QUEUE.put([self.environment.sim_time, self.environment.rnd.expon(
            self.mean_service, self.STREAM_SERVICE)])
        self.stat_monitor.timest_update(
            float(self.LIST_QUEUE.qsize()), 'LIST_QUEUE', self.environment.sim_time)

        # If the CPU is idle, start a CPU run.
        if self.LIST_CPU.empty():
            self.start_CPU_run()

    def start_CPU_run(self):
        """ Non-event function to start a CPU run of a job. """

        # Remove the first job from the queue.
        transfer = self.LIST_QUEUE.get()
        self.stat_monitor.timest_update(
            float(self.LIST_QUEUE.qsize()), 'LIST_QUEUE', self.environment.sim_time)

        # Determine the CPU time for this pass, including the swap time.
        if (self.quantum < transfer[1]):
            run_time = self.quantum + self.swap
        else:
            run_time = transfer[1] + self.swap

        # Decrement remaining CPU time by a full quantum.  (If less than a full
        # quantum is needed, this attribute becomes negative.  This indicates that
        # the job, after exiting the CPU for the current pass, will be done and is
        # to be sent back to its terminal.)
        transfer[1] -= self.quantum

        # Place the job into the CPU.
        self.LIST_CPU.put(transfer)
        self.stat_monitor.timest_update(
            float(self.LIST_CPU.qsize()), 'LIST_CPU', self.environment.sim_time)

        # Schedule the end of the CPU run.
        self.environment.event_schedule(
            Event(self.environment.sim_time + run_time, EventTypes.EVENT_END_CPU_RUN)
        )

    def end_CPU_run(self):
        """ Event function to end a CPU run of a job. """

        # Remove the job from the CPU.
        transfer = self.LIST_CPU.get()
        self.stat_monitor.timest_update(
            float(self.LIST_CPU.qsize()), 'LIST_CPU', self.environment.sim_time)

        # Check to see whether this job requires more CPU time.
        if (transfer[1] > 0.0):

            # This job requires more CPU time, so place it at the end of the queue
            # and start the first job in the queue.
            self.LIST_QUEUE.put(transfer)
            self.stat_monitor.timest_update(
                float(self.LIST_QUEUE.qsize()), 'LIST_QUEUE', self.environment.sim_time)
            self.start_CPU_run()
        else:

            # This job is finished, so collect response-time statistics and send it
            # back to its terminal, i.e., schedule another arrival from the same
            # terminal.
            self.stat_monitor.sampst_update(self.environment.sim_time - transfer[0])

            self.environment.event_schedule(
                Event(
                    self.environment.sim_time + self.environment.rnd.expon(self.mean_think, self.STREAM_THINK),
                    EventTypes.EVENT_ARRIVAL
                )
            )

            # Increment the number of completed jobs.
            self.num_responses += 1

            # Check to see whether enough jobs are done.
            if self.num_responses >= self.num_responses_required:

                # Enough jobs are done, so schedule the end of the simulation
                # immediately (forcing it to the head of the event list).
                self.environment.event_schedule(
                    Event(self.environment.sim_time, EventTypes.EVENT_END_SIMULATION)
                )

            else:
                # Not enough jobs are done; if the queue is not empty, start
                # another job.
                if not self.LIST_QUEUE.empty():
                    self.start_CPU_run()


class Environment:
    def __init__(self) -> None:
        # Open output files.
        self.outfile = open("tscomp.out", "w")

        # Open input file & read input parameters.
        with open("tscomp.in", "r") as infile:
            min_terms, max_terms, incr_terms, num_responses_required, mean_think, mean_service, quantum, swap \
                = infile.readline().split()
        self.min_terms = int(min_terms)
        self.max_terms = int(max_terms)
        self.incr_terms = int(incr_terms)
        self.num_responses_required = int(num_responses_required)
        self.mean_think = float(mean_think)
        self.mean_service = float(mean_service)
        self.quantum = float(quantum)
        self.swap = float(swap)

        # Write report heading and input parameters.
        self.outfile.write("Time-shared computer model\n\n")
        self.outfile.write(f"Number of terminals{self.min_terms:9d} to{self.max_terms:4d} by {self.incr_terms:4d}\n\n")
        self.outfile.write(f"Mean think time  {self.mean_think:11.3f} seconds\n\n")
        self.outfile.write(f"Mean service time{self.mean_service:11.3f} seconds\n\n")
        self.outfile.write(f"Quantum          {self.quantum:11.3f} seconds\n\n")
        self.outfile.write(f"Swap time        {self.swap:11.3f} seconds\n\n")
        self.outfile.write(f"Number of jobs processed{self.num_responses_required:12d}\n\n\n")
        self.outfile.write("Number of      Average         Average")
        self.outfile.write("       Utilization\n")
        self.outfile.write("terminals   response time  number in queue     of CPU")

        self.LIST_EVENTS: List  # define simulation environment
        self.rnd = RandomNumberGenerator()  # random number generation
        self.STREAM_THINK = 1  # Random-number stream for think times.
        self.sim_time: int

        # keep result of simulation for visualization
        self.result_dict = {}

    def reset_environment(self):
        self.LIST_EVENTS = []
        self.sim_time = 0

    def run(self):
        # Run the simulation varying the number of terminals.
        for num_terms in range(self.min_terms, self.max_terms + 1, self.incr_terms):
            # Initialize simulation environment
            cpu_obj = Computer(
                quantum=self.quantum,
                swap=self.swap,
                environment=self,
                num_responses_required=self.num_responses_required,
                STREAM_THINK=self.STREAM_THINK,
                mean_service=self.mean_service,
                mean_think=self.mean_think,
            )
            self.reset_environment()

            # Schedule the first arrival to the CPU from each terminal.
            for term in range(1, num_terms + 1):
                self.event_schedule(
                    Event(
                        self.rnd.expon(self.mean_think, self.STREAM_THINK),
                        EventTypes.EVENT_ARRIVAL
                    )
                )

            # Run the simulation until it terminates after an end-simulation event
            # (type EVENT_END_SIMULATION) occurs.
            while True:
                # Determine the next event.
                current_event = self.timing()

                # Invoke the appropriate event function.
                if current_event.EVENT_TYPE == EventTypes.EVENT_ARRIVAL:
                    cpu_obj.arrive()
                elif current_event.EVENT_TYPE == EventTypes.EVENT_END_CPU_RUN:
                    cpu_obj.end_CPU_run()
                elif current_event.EVENT_TYPE == EventTypes.EVENT_END_SIMULATION:
                    self.report(num_terms, cpu_obj)
                    break

        self.outfile.close()

    def event_schedule(self, event: Event):
        """
        Schedule an event at time event_time of type event_type.  If attributes
        beyond the first two (reserved for the event time and the event type) are
        being used in the event list, it is the user's responsibility to place
        their values into the transfer array before invoking event_schedule.
        """
        self.LIST_EVENTS.append(event)
        self.LIST_EVENTS.sort(key=lambda e: e.EVENT_TIME)

    def timing(self):
        """
        Remove & return next event from event list
        """

        # Remove the first event from the event list.
        transfer = self.LIST_EVENTS.pop(0)

        # Check for a time reversal.
        if transfer.EVENT_TIME < self.sim_time:
            raise Exception(
                f"\nAttempt to schedule event type {transfer.EVENT_TYPE} for "
                f"time {transfer.EVENT_TIME} at time {self.sim_time}\n"
            )

        self.sim_time = transfer.EVENT_TIME

        return transfer

    def report(self, num_terms, cpu_obj):
        """Report generator function."""

        # Get and write out estimates of desired measures of performance.
        average_response_time = cpu_obj.stat_monitor.sampst_report()
        average_number_in_queue = cpu_obj.stat_monitor.timest_report('LIST_QUEUE', self.sim_time)
        cpu_utilization = cpu_obj.stat_monitor.timest_report('LIST_CPU', self.sim_time)
        self.result_dict[num_terms] = [average_response_time, average_number_in_queue, cpu_utilization]
        self.outfile.write(
            f"\n\n{num_terms:5d}{average_response_time:16.3f}"
            f"{average_number_in_queue:16.3f}"
            f"{cpu_utilization:16.3f}"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--no-plot', '-p', dest="plot", action="store_false", help="don't plot results", default=True
    )
    parsed = parser.parse_args()

    e = Environment()
    e.run()

    if parsed.plot:
    
        import matplotlib.pyplot as plt

        ypoints = []
        x_average_response_times = []
        x_average_number_in_queues = []
        x_cpu_utilizations = []
        for k in e.result_dict:
            ypoints.append(k)
            x_average_response_times.append(e.result_dict[k][0])
            x_average_number_in_queues.append(e.result_dict[k][1])
            x_cpu_utilizations.append(e.result_dict[k][2])

        plt.figure(figsize=(10,4))
        plt.subplot(121)
        plt.xlabel("Number Of Terminals")
        plt.plot(ypoints, x_average_response_times, 'r-', label='Average Response Times')
        plt.plot(ypoints, x_average_number_in_queues, 'b-', label='Average Number In Queues')
        plt.legend(loc='lower right', frameon=False, fontsize="small")
        plt.subplot(122)
        plt.xlabel("Number Of Terminals")
        plt.plot(ypoints, x_cpu_utilizations, 'g-', label='Cpu Utilizations')
        plt.legend(loc='lower right', frameon=False, fontsize="small")
        plt.show()


if __name__ == "__main__":
    main()
