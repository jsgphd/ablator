from multiprocessing import Queue
import multiprocessing
from multiprocessing import Process,Pipe, freeze_support
import threading
import random
import time
import os

from models.sensitivity_solution import *
from models.HDMR import *
from models.HDMR_out import *
from models.stat_file import *
from lib.sobol_lib import *

solution_name = "0_CsI_1.0nm_125mJsens_1024"

class AblatorSensitivity:
        def __init__(self, config_src="config_sensitivity.txt"):                
                self.config = config_src

                self.number_of_simulations = 1024

                self.time_limit = 2.5 #cas v sekundach
                self.include_err = False

                self.solutions = []

                self.processing_queue = Queue(maxsize=self.number_of_simulations)
                self.post_process_queue = Queue(maxsize=self.number_of_simulations)

        def prepare(self):
                seed = 1
                while(seed <= self.number_of_simulations):  
                        while(os.getcwd().split("\\")[-1] == "temp"):
                                time.sleep(0.2)
                        solution = SensitivitySolution(solution_name,"config_sensitivity.txt",str(seed))
                        solution.load()

                        [sob,seed] = i4_sobol(solution.get_param_count(),seed)
                        #vygeneruj reseni a pridej jej do fronty
                        
                        scaled = solution.scale_matdata(sob)
                        solution.save_matdata()
                        #hdmr_in.add(sob)

                        self.processing_queue.put(solution)                     


        def run(self):
                i = 0
                while(i<self.number_of_simulations):                    
                        #print('novy process start')
                        current_solution = self.processing_queue.get()        
                        current_solution.run()
                        #worker = Process(target=current_solution.run)
                        #worker.start()
                        #start_time = time.time()
                        
                        #while(worker.is_alive() and (start_time + self.time_limit) > time.time()):
                        #        print('monitoring')
                        #        time.sleep(0.1)   

                        #if((start_time + self.time_limit) < time.time()):
                        #if(worker.is_alive()):
                                #print('process: ' + str(i) + ' terminating')
                                #worker.terminate()
                                
                        #print('joinin')
                        #worker.join()
                        #print('neprosel')
                        #print(worker)
                        #if(worker.is_alive()):
                                #worker.terminate()
                        #        time.sleep(5.0)
                        #        print('process: ' + str(i) + ' terminated')

                        self.post_process_queue.put(current_solution)           
                        time.sleep(0.5)

                        #print('konci jedno reseni')
                        i += 1

        
        def finish(self):
                hdmr_in = HDMR(solution_name)
                hdmr_out = HDMROut(solution_name)
                stat = StatFile(solution_name,hdmr_in, self.number_of_simulations)
                i = 0
                while(i<self.number_of_simulations):
                        current_solution = self.post_process_queue.get()
                        while(os.getcwd().split("\\")[-1] == "temp"):
                                time.sleep(0.2)
                        try:                  
                                hdmr_out.add([current_solution.read_depth()])
                                hdmr_in.add(current_solution.get_sobol())
                        except(IOError):
                                if self.include_err:
                                        hdmr_in.add(current_solution.get_sobol())
                                        hdmr_out.add(["ERR","err2"])
                                
                                print('nevytvoril reseni: ' + str(i))
                        
                        i += 1
                
                hdmr_in.print_hdmr()
                hdmr_out.print_hdmr()
                stat.save()

if __name__ == '__main__':
        freeze_support()
                
        a = AblatorSensitivity()
        a.prepare()
        a.run()
        a.finish()
        
