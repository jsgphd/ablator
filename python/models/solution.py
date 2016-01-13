#trida, ktera bude uchovavat solution
#atributy:
#  name - jmeno solutionu, podle nej vznikne adresar ve slozce solutions
#  init_data - prirazeny soubor init_data
#  mat_data - prirazeny soubor mat_data
#  opac_data - prirazeny soubor opac_data
#  ablator_exe - prirazena binarka, ktera provede vypocet

#metody:
#  load - pro nacteni solutionu
#  save - pro ulozeni solutionu do slozky solutions

import os, shutil
from models.init_data import InitData
from models.mat_data import MatData

class Solution:
    def __init__(self, name = None, init_data = None, mat_data = "", opac_data = "", ablator_exe = "Ablator.exe"):
        self.name = name
        #zatim necham jen nacist soubor, v budoucnu se do nej bude zapisovat a bude mit  vlastni tridu
        #self.init_data = InitData(init_data)
        #self.mat_data = MatData(mat_data)
        #self.opac_data = opac_data
        self.init_data = None
        self.mat_data = None
        self.opac_data = None
        self.ablator_exe = ablator_exe        

    def set_init_data(self, value):        
        self.init_data = InitData(value,self.ablator_exe)

    def get_init_data(self):
        return self.init_data

    def set_mat_data(self, value):
        self.mat_data = MatData(value, self.ablator_exe)

    def get_mat_data(self):
        return self.mat_data

    def set_opac_data(self, value):
        self.opac_data = value

    def set_ablator_exe(self, value):        
        self.ablator_exe = value        

    def set_name(self,value):
        self.name = value

    #nacteni solutiony bude probihat jednoduse ze souboru solution.txt:
    #  -prvni radek souboru nazev(cesta) init_data
    #  -druhy radek souboru nazev(cesta) mat_data
    #  -treti radek souboru nazev(cesta) opac_data
    def load(self):
        if(self.name == ""):
            raise "Nezadano reseni"        
        try:
            solution_file = open(self.name, "r")
        except IOError as e:
            raise ('Nemohu otevrit soubor: solutions/' + self.name + "/solution.txt")

        self.name = solution_file.readline().replace("\n","")
        self.ablator_exe = solution_file.readline().replace("\n","")
        self.init_data = InitData(solution_file.readline().replace("\n",""), self.ablator_exe)
        self.mat_data = MatData(solution_file.readline().replace("\n",""), self.ablator_exe)
        self.opac_data = solution_file.readline().replace("\n","")

        solution_file.close()

    #ulozeni struktura solution.txt stejna jako 
    def save(self):        
        directory = "solutions/" + self.name
        if(os.path.exists(directory) == False):
            os.makedirs(directory)

        path = ("solutions/" + self.name + "/solution.txt")
        if(os.path.exists(path)):
            try:
                os.remove(path)
            except:
                print("nemohu odstranit path")

        solution_file = open(path, "w")
        solution_file.write(self.name + "\n")
        solution_file.write(self.ablator_exe + "\n")
        solution_file.write(self.init_data.filepath + "\n")
        solution_file.write(self.mat_data.filepath + "\n")
        solution_file.write(self.opac_data)
        solution_file.close()      

    def copy_solution_to(self,target):
        #shutil.copy("bin/Ablator.exe", target)
        shutil.copy("bin/" + self.ablator_exe, target)
        shutil.copyfile(self.init_data.filepath,target+"/initdata")
        shutil.copyfile(self.mat_data.filepath,target+"/matdata")
        shutil.copyfile(self.opac_data,target+"/opacdata")        
    
    def prepared_to_run(self):
        if(os.path.exists(self.init_data.filepath) == False):
            return False
        if(os.path.exists(self.mat_data.filepath) == False):
            return False
        if(os.path.exists(self.opac_data) == False):
            return False
        return True

    #otevre slozku s resenimi v exploreru
    def open_folder(self):        
        os.system("explorer " + "solutions\\" + self.name)

    def run(self):
        if(self.name == None):
            raise "Nebylo zadano jmeno"
        
        self.save()
        self.copy_solution_to("temp")
        if(self.prepared_to_run()):
            os.chdir("temp")
            i = os.system(self.ablator_exe)
            os.system("del " + self.ablator_exe)
            os.system("del initdata")
            os.system("del matdata")
            os.system("del opacdata")
            os.system("copy *.* " + "..\\solutions\\" + self.name)
            os.chdir("..")
            os.system("del /Q temp")

