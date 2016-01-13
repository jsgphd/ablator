from models.solution import Solution
import os,shutil,subprocess

#struktura adresaru pro adresar solution
#name bude slouzit jako vychozi nazev, solutions/name/
#ve slozce name, pak budou pouzity config a soubory pro hdmr

class SensitivitySolution(Solution):
    def __init__(self, name, config_src, job_id ):
        Solution.__init__(self,name=name)
        
        self.config_src = config_src        
        self.config_params = []
        self.sampling = "lin"
        self.kill_timeout = 10 #seconds
        
        self.job_id = job_id
        self.sobol = None
        
    def load(self):
        if(self.name == ""):
            raise "Nezadano reseni"        
        try:
            solution_file = open(self.config_src, "r")
        except IOError as e:
            try:
                solution_file = open("../" + self.config_src, "r")
            except IOError as e:
                print ('Nemohu otevrit soubor: ' + self.config_src)
             
        lines = solution_file.readlines()
        
        for line in lines:            
            line = line.replace("\n","")
            if line[0] == "#":
              continue
            arguments = line.split("=")            
            
            if arguments[0] == "executable":
                self.set_ablator_exe(arguments[1])
            elif arguments[0] == "sampling":
                self.sampling = arguments[1]
            elif arguments[0] == "init_data":
                self.set_init_data(arguments[1])
            elif arguments[0] == "mat_data":
                self.set_mat_data(arguments[1])            
            elif arguments[0] == "opac_data":
                self.set_opac_data(arguments[1])
            else:
                self.config_params.append(arguments)  

        solution_file.close()
        
    def save(self):        
        directory = "solutions/" + self.name + "/" + self.job_id
        if(os.path.exists(directory) == False):
            os.makedirs(directory)

        path = ("solutions/" + self.name + "/sensitivity_solution.txt")
        if(os.path.exists(path)):
            try:
                os.remove(path)
            except:
                print("nemohu odstranit path")
              
        shutil.copy(self.config_src, directory)
        #solution_file = open(path, "w")
        #solution_file.write(self.name + "\n")
        #solution_file.write(self.ablator_exe + "\n")
        #solution_file.write(self.init_data.filepath + "\n")
        #solution_file.write(self.mat_data.filepath + "\n")        
        #solution_file.write(self.opac_data + "\n")        
        #solution_file.close()

    def get_param_count(self):
        i = 0
        for params in self.config_params:
            param_arr = params[1].split(';')

            if param_arr[0] == 'p':
                i += 1
            else:
                indexes = param_arr[0].split(',')                
                for index in indexes:
                    i += 1
        return i

              
    def scale_matdata(self, sobol):
        self.sobol = sobol
        self.mat_data.scale(sobol,self.config_params,self.sampling)

    def save_matdata(self):
        directory = "solutions/" + self.name + "/" + self.job_id
        if(os.path.exists(directory) == False):
            os.makedirs(directory)
        self.mat_data.save(directory + "/matdata")

    def get_sobol(self):
        return self.sobol

    def read_depth(self):
        src = open("solutions\\" + self.name + "\\" + self.job_id + "\\" + "summary", 'r')
        depth_line = src.readlines()[11]
        src.close()
        depth_arr = depth_line.split("=")
        return depth_arr[1].strip()

    def run(self):
        if(self.name == None):
            raise "Nebylo zadano jmeno"
        
        self.save()
        self.copy_solution_to("temp")
        if(self.prepared_to_run()):
            os.chdir("temp")
            #i = os.system(self.ablator_exe)
            try:
                subprocess.call(self.ablator_exe, timeout=self.kill_timeout)
            except:
                print('killed')

            os.system("del " + self.ablator_exe)
            #os.system("del initdata")
            #os.system("del matdata")
            #os.system("del opacdata")
            os.system("copy *.* " + "..\\solutions\\" + self.name + "\\" + self.job_id + "\\")
            os.chdir("..")
            os.system("del /Q temp")

class SensitivitySolutionScan(SensitivitySolution):

    def read_depth(self,i):
        src = open("solutions\\" + self.name + "\\" + self.job_id + "\\" + "summary-" + str(i), 'r')
        depth_line = src.readlines()[11]
        src.close()
        depth_arr = depth_line.split("=")
        return depth_arr[1].strip()


    def run(self,i):
        if(self.name == None):
            raise "Nebylo zadano jmeno"
        
        self.save()
        self.copy_solution_to("temp")
        if(self.prepared_to_run()):
            os.chdir("temp")
            #i = os.system(self.ablator_exe)
            try:
                subprocess.call(self.ablator_exe, timeout=self.kill_timeout)
            except:
                print('killed')

            os.system("del " + self.ablator_exe)
            #os.system("del initdata")
            #os.system("del matdata")
            #os.system("del opacdata")
            os.system("copy summary " + "..\\solutions\\" + self.name + "\\" + self.job_id + "\\summary-" + str(i))
            os.chdir("..")
            os.system("del /Q temp")
