from multiprocessing import Queue
import multiprocessing
from multiprocessing import Process,Pipe
import threading
import random
import time
import os

from os.path import isfile, join

from models.sensitivity_solution import *
from models.HDMR import *
from models.HDMR_out import *
from models.stat_file import *
from lib.sobol_lib import *


SCAN_DIR = "scan/"

solution_name = "misa_test"

class AblatorScan:
        def __init__(self):
                self.configs = self._load_configs()

                self.number_of_simulations = 2
                self.time_limit = 0.5 #cas v sekundach

                self.solutions = []

                self.processing_queues = []
                self.post_process_queues = []

                print(len(self.configs))

                for config in self.configs:
                        self.processing_queues.append(Queue(maxsize=self.number_of_simulations*len(self.configs)))
                        self.post_process_queues.append(Queue(maxsize=self.number_of_simulations*len(self.configs)))

        def _load_configs(self):
                return [ SCAN_DIR + f for f in os.listdir(SCAN_DIR) if isfile(join(SCAN_DIR,f)) ]

        def prepare(self):
                seed = 1
                while(seed <= (self.number_of_simulations)):
                        for i_cfg in range(len(self.configs)):
                                while(os.getcwd().split("\\")[-1] == "temp"):
                                        time.sleep(0.2)

                                solution = SensitivitySolutionScan(solution_name,self.configs[i_cfg],str(seed))
                                solution.load()

                                [sob,seed] = i4_sobol(solution.get_param_count(),seed)
                                #vygeneruj reseni a pridej jej do fronty
                        
                                scaled = solution.scale_matdata(sob)
                                solution.save_matdata()
                                #hdmr_in.add(sob)
                                
                                self.processing_queues[i_cfg].put(solution)
                                seed -= 1
                        seed += 1


        def run(self):
                i = 0
                while(i<self.number_of_simulations):                    
                        for i_cfg in range(len(self.configs)):
                                current_solution = self.processing_queues[i_cfg].get()
                                current_solution.run(i_cfg)
                                #worker = Process(target=current_solution.run())
                                #worker.start()
                                #start_time = time.time()
                        
                                #while(worker.is_alive() and (start_time + self.time_limit) > time.time()):
                                #        time.sleep(1.0)   
                        
                                #worker.terminate()
                                #worker.join()
                                self.post_process_queues[i_cfg].put(current_solution)           
                                time.sleep(0.5)

                        i += 1

        
        def finish(self):
                hdmr_in = HDMR(solution_name)
                hdmr_out = HDMRScanOut(solution_name)
                stat = StatFile(solution_name,hdmr_in, self.number_of_simulations)
                i = 0
                while(i<self.number_of_simulations):
                        temp = []
                        for i_cfg in range(len(self.configs)):
                                current_solution = self.post_process_queues[i_cfg].get()
                                while(os.getcwd().split("\\")[-1] == "temp"):
                                        time.sleep(0.2)

                                #ulozime pouze jednou
                                if i_cfg == 0:
                                        hdmr_in.add(current_solution.get_sobol())

                                try:
                                        temp.append(current_solution.read_depth(i_cfg)) 
                                except(IOError):
                                        temp.append("ERR")
                                        #print('nevytvoril reseni: ' + str(i) + ' pri reseni configu ' + self.configs[i])
                        hdmr_out.add(temp)
                        
                        i += 1
                
                hdmr_in.print_hdmr()
                hdmr_out.print_hdmr()
                stat.save()
                
                
a = AblatorScan()
a.prepare()
a.run()
a.finish()
