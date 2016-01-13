#trida, ktera bude reprezentovat strukturu pro initdata

class InitData:
    def __init__(self, filepath = None, ablator_exe = "Ablator.exe"):
        self.filepath = filepath
        self.ablator_exe = ablator_exe

        if(filepath == None):
            self.set_default()
        else:
            if ablator_exe == "Ablator.exe":
                self.load()
            else:
                self.load2()

    def set_filepath(self, value):
        self.filepath = value
        if ablator_exe == "Ablator.exe":
            self.load()
        else:
            self.load2()

    def set_default(self):
        self.source = 1
        self.absorbtion_coeff = 0.080
        self.first_zone_size = 1
        self.source_type = 1
        self.pulse_duration = 0.001
        self.pulse_energy = 0.6
        self.blackbody_temperature = 0
        self.spot_radius = 5.6e3
        self.tstop = 500

    def set_source(self, value):        
        self.source = value

    def set_absorbtion(self, value):
        self.absorbtion_coeff = value;

    def set_first_zone_size(self, value):
        self.first_zone_size = value;

    def set_source_type(self, value):
        self.source_type = value

    def set_pulse_duration(self, value):
        self.pulse_duration = value

    def set_pulse_energy(self, value):
        self.pulse_energy = value

    def set_blackbody_temperature(self, value):
        self.blackbody_temperature = value

    def set_spot_radius(self, value):
        self.spot_radius = value

    def set_yield_test(self, value):
        self.yield_test = value

    def set_tstop(self, value):
        self.tstop = value

    def get_source_type(self):        
        return self.source_type

    def get_absorbtion(self):
        return self.absorbtion_coeff

    def get_first_zone_size(self):
        return self.first_zone_size

    def get_source(self):
        return self.source

    def get_pulse_duration(self):
        return self.pulse_duration

    def get_pulse_energy(self):
        return self.pulse_energy

    def get_blackbody_temperature(self):
        return self.blackbody_temperature

    def get_spot_radius(self):
        return self.spot_radius

    def get_yield_test(self):
        return self.yield_test

    def get_tstop(self):
        return self.tstop

    #v teto metode budeme ukladat data do slozky odpovidajici prislusnemu reseni
    def save(self,filepath=None):
        try:
            if(filepath == None):
                path = ("solutions/" + self.name + "/initdata")
                if(os.path.exists(path)):
                    os.remove(path)
            else:
                path = filepath
                
            data_file = open(path, "w")
            
            data_file.write("# BB X-ray or laser source\n")            
            data_file.write(str(self.source))
            data_file.write("\n# absorbtion coefficient in microns for laser source\n")
            data_file.write(str(self.absorbtion_coeff))
            data_file.write("\n# the first zone size (nm)\n")
            data_file.write(str(self.first_zone_size))
            data_file.write("\n# source type (1=square pulse;2=Gaussian pulse)\n")
            data_file.write(str(self.source_type))
            data_file.write("\n# pulse duration or pulse FWHM (ns)\n")
            data_file.write(str(self.pulse_duration))
            data_file.write("\n# total energy in the pulse (J)\n")
            data_file.write(str(self.pulse_energy))
            data_file.write("\n# blackbody temperature for X-ray source (keV)\n")
            data_file.write(str(self.blackbody_temperature))
            data_file.write("\n# radius of the irradiated spot (microns) (5.6E3 microns to give the area 1cm2)\n")
            data_file.write(str(self.spot_radius))
            data_file.write("\n# Tstop (ns)\n")
            data_file.write(str(self.tstop))
            
            data_file.close()            

        except(IOError):
            raise "nastala chyba, pri ukladani initdat"

    #pro druhou variantu exe souboru
    def save2(self,filepath=None):
        try:
            if(filepath == None):
                path = ("solutions/" + self.name + "/initdata")
                if(os.path.exists(path)):
                    os.remove(path)
            else:
                path = filepath
                
            data_file = open(path, "w")
            
            data_file.write("# BB X-ray or laser source\n")            
            data_file.write(str(self.source))
            data_file.write("\n# absorbtion coefficient in microns for laser source\n")
            data_file.write(str(self.absorbtion_coeff))
            data_file.write("\n# the first zone size (nm)\n")
            data_file.write(str(self.first_zone_size))
            data_file.write("\n# source type (1=square pulse;2=Gaussian pulse)\n")
            data_file.write(str(self.source_type))
            data_file.write("\n# pulse duration or pulse FWHM (ns)\n")
            data_file.write(str(self.pulse_duration))
            data_file.write("\n# total energy in the pulse (J)\n")
            data_file.write(str(self.pulse_energy))
            data_file.write("\n# blackbody temperature for X-ray source (keV)\n")
            data_file.write(str(self.blackbody_temperature))
            data_file.write("\n# radius of the irradiated spot (microns) (5.6E3 microns to give the area 1cm2)\n")
            data_file.write(str(self.spot_radius))
            data_file.write("\n# whether to use YIELDTEST value for zones with stress bigger than YIELDTEST limit\n")            
            data_file.write(str(self.yield_test))
            data_file.write("\n# Tstop (ns)\n")
            data_file.write(str(self.tstop))
            
            data_file.close()            

        except(IOError):
            raise "nastala chyba, pri ukladani initdat"

    #nacteme data na zaklade solution
    def load(self):
        try:
            data_file = open(self.filepath, "r")
        except(IOError):
            raise "nepodarilo se otevrit soubor: " + self.filename
            
        
        #komentar
        data_file.readline()
        #zdroj
        try:
            self.source = int(data_file.readline())
        except ValueError:
            raise "spatny datovy typ zdroje(1=Laser)"

        #komentar
        data_file.readline()
        #absorbcni koeficient
        try:
            self.absorbtion_coeff = float(data_file.readline())
        except:
            raise "spatny datovy typ absorbcniho koeficientu(realne cislo)"

        #komentar
        data_file.readline()
        #velikost zony
        try:
            self.first_zone_size = float(data_file.readline())
        except ValueError:
            raise "spatny datovy typ velikosti prvni zony(cele cislo znacici nanometry)"

        #komentar
        data_file.readline()
        #typ pulzu
        try:
            self.source_type = int(data_file.readline())
        except ValueError:
            raise "spatny datovy typ zdroje(1=ctvercovy pulz nebo 2=Guasovsky pulz)"

        #komentar
        data_file.readline()
        #delka pulzu
        try:
            self.pulse_duration = float(data_file.readline())
        except:
            raise "spatny datovy typ delky trvani pulzu(realne cislo v ns)"

        #komentar
        data_file.readline()
        #celkova energie pulzu
        try:
            self.pulse_energy = float(data_file.readline())
        except:
            raise "spatny datovy typ energie pulzu(realne cislo v Joulech)"

        #komentar
        data_file.readline()
        #teplota cerneho telesa
        try:
            self.blackbody_temperature = float(data_file.readline())
        except ValueError:
            raise "spatny datovy typ zdroje(cele cislo)"

        #komentar
        data_file.readline()        
        #polomer ozareneho bodu
        try:
            self.spot_radius = float(data_file.readline().replace("D","E"))
        except:
            raise "spatny datovy typ polomeru bodu(realne cislo)"

        #komentar
        data_file.readline()
        #cas konce
        try:
            self.tstop = int(data_file.readline())
        except ValueError:
            raise "spatny datovy casu konce(cele cislo)"

        data_file.close()

    #nacteme data na zaklade solution
    def load2(self):
        try:
            data_file = open(self.filepath, "r")
        except(IOError):
            raise "nepodarilo se otevrit soubor: " + self.filename
            
        
        #komentar
        data_file.readline()
        #zdroj
        try:
            self.source = int(data_file.readline())
        except ValueError:
            raise "spatny datovy typ zdroje(1=Laser)"

        #komentar
        data_file.readline()
        #absorbcni koeficient
        try:
            self.absorbtion_coeff = float(data_file.readline())
        except:
            raise "spatny datovy typ absorbcniho koeficientu(realne cislo)"

        #komentar
        data_file.readline()
        #velikost zony
        try:
            self.first_zone_size = float(data_file.readline())
        except ValueError:
            raise "spatny datovy typ velikosti prvni zony(cele cislo znacici nanometry)"

        #komentar
        data_file.readline()
        #typ pulzu
        try:
            self.source_type = int(data_file.readline())
        except ValueError:
            raise "spatny datovy typ zdroje(1=ctvercovy pulz nebo 2=Guasovsky pulz)"

        #komentar
        data_file.readline()
        #delka pulzu
        try:
            self.pulse_duration = float(data_file.readline())
        except:
            raise "spatny datovy typ delky trvani pulzu(realne cislo v ns)"

        #komentar
        data_file.readline()
        #celkova energie pulzu
        try:
            self.pulse_energy = float(data_file.readline())
        except:
            raise "spatny datovy typ energie pulzu(realne cislo v Joulech)"

        #komentar
        data_file.readline()
        #teplota cerneho telesa
        try:
            self.blackbody_temperature = float(data_file.readline())
        except ValueError:
            raise "spatny datovy typ zdroje(cele cislo)"

        #komentar
        data_file.readline()        
        #polomer ozareneho bodu
        try:
            self.spot_radius = float(data_file.readline().replace("D","E"))
        except:
            raise "spatny datovy typ polomeru bodu(realne cislo)"

        #komentar
        data_file.readline()
        #cas konce
        try:
            self.yield_test = int(data_file.readline())
        except ValueError:
            raise "spatny datovy yield test(cele cislo)"

        #komentar
        data_file.readline()
        #cas konce
        try:
            self.tstop = int(data_file.readline())
        except ValueError:
            raise "spatny datovy casu konce(cele cislo)"

        data_file.close()
