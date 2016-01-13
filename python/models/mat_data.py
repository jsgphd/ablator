#trida, ktera bude reprezentovat objekt pro matdata
from array import array
import math

class MatData:
    def __init__(self, filepath = None, ablator_exe = "Ablator.exe"):
        self.filepath = filepath
        self.ablator_exe = ablator_exe

        if(filepath == None):
            self.set_default()
        else:
            self.load()


    def set_filepath(self, value):
        self.filepath = value
        self.load()

    def set_default():
        self.name = ""
        self.density = 0.0
        self.t_melt = 0
        self.conductivity_solid = array([0.0,0.0,0.0,0.0,0.0,0.0])
        self.conductivity_liquid = array([0.0,0.0,0.0,0.0,0.0,0.0])
        self.conductivity_vapor = array([0.0,0.0,0.0,0.0])
        self.th_solid = array([0.0,0.0,0.0,0.0,0.0,0.0])
        self.th_liquid = array([0.0,0.0])
        self.th_vapor = array([0.0,0.0])
        self.log_coeffs = array([0.0,0.0])
        self.gruneisen_coeffs = array([0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0])
        self.poissons_ratio = 0.0
        self.yield_strength = array([0.0,0.0,0.0])
        self.surface_tension = array([0.0,0.0,0.0,0.0,0.0,0.0,0.0])

    def set_name(self, value):
        self.name = value

    def set_density(self, value):
        self.density = value
        
    def set_t_melt(self, value):
        self.t_melt = value

    def set_conductivity_solid(self, value):
        self.conductivity_solid = value

    def set_conductivity_liquid(self, value):
        self.conductivity_liquid = value

    def set_conductivity_vapor(self, value):
        self.conductivity_vapor = value

    def set_th_solid(self, value):
        self.th_solid = value

    def set_th_liquid(self, value):
        self.th_liquid = value

    def set_th_vapor(self, value):
        self.th_vapor = value

    def set_log_coeffs(self, value):
        self.log_coeffs = value

    def set_gruneisen_coeffs(self, value):
        self.gruneisen_coeffs = value

    def set_poissons_ratio(self, value):
        self.poissons_ratio = value

    def set_yield_strength(self, value):
        self.yield_strength = value

    def set_surface_tension(self, value):
        self.surface_tension = value
				
    def get_name(self):
        return self.name
        
    def get_density(self):
        return self.density
        
    def get_t_melt(self):
        return self.t_melt

    def get_conductivity_solid(self):
        return self.conductivity_solid

    def get_conductivity_liquid(self):
        return self.conductivity_liquid

    def get_conductivity_vapor(self):
        return self.conductivity_vapor

    def get_th_solid(self):
        return self.th_solid

    def get_th_liquid(self):
        return self.th_liquid

    def get_th_vapor(self):
        return self.th_vapor

    def get_log_coeffs(self):
        return self.log_coeffs

    def get_gruneisen_coeffs(self):
        return self.gruneisen_coeffs

    def get_poissons_ratio(self):
        return self.poissons_ratio

    def get_yield_strength(self):
        return self.yield_strength

    def get_surface_tension(self):
        return self.surface_tension

    def save(self, filepath=None):
        try:
            if(filepath == None):
                path = ("solutions/" + self.name + "/matdata")
                if(os.path.exists(path)):
                    os.remove(path)
            else:
                path = filepath
                self.filepath = filepath
                
            data_file = open(path, "w")
            
            data_file.write(self.name)
            data_file.write("# density [kg/m3]\n")
            data_file.write(str(self.density))
            data_file.write("\n#T melt [K]\n")
            data_file.write(str(self.t_melt))

            data_file.write("\n# k (conductivity) vs. T, solid [cal/m.K.s]\n")            
            for val in self.conductivity_solid:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# k vs. T, liquid [cal/m.K.s]\n")
            for val in self.conductivity_liquid:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# k vs. T, vapor [cal/m.K.s]\n")
            for val in self.conductivity_vapor:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# T vs. H, solid [K & cal/g]\n")
            for val in self.th_solid:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# T vs. H, liquid [K & cal/g]\n")
            for val in self.th_liquid:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# T vs. H, vapor [K & cal/g]\n")
            for val in self.th_vapor:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# coeffs for log(P) = A - B/T [micron & K]\n")
            for val in self.log_coeffs:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# Gruneisen EOS coeff & gamma & Rbar in [J/kg.K]\n")
            for val in self.gruneisen_coeffs:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# Poisson's Ratio\n")
            data_file.write(str(self.poissons_ratio))

            data_file.write("\n# Yield strength (flag + coeffs) [Pa & eV]\n")
            for val in self.yield_strength:
                data_file.write(str(val))
                data_file.write("\n")

            data_file.write("# Surface tension (flag + 6 coeff)")
            for val in self.surface_tension:
                data_file.write("\n")
                data_file.write(str(val))
                
            
            data_file.close()
        except(IOError):
            raise "Nastala chyba pri ukladani souboru matdata"

    def load(self):
        try:
            data_file = open(self.filepath, "r")
        except(IOError):
            print("nepodarilo se otevrit soubor: " + self.filepath)
            raise

        #nazev materialu
        try:
            self.name = data_file.readline()
        except:
            raise Exception("nepodarilo se nacist parametr jmeno materialu")

        #comment
        data_file.readline()
        #density
        try:
            self.density = float(data_file.readline().replace("D","E"))
        except:
            raise Exception("spatny datovy typ density(realne cislo)")

        #comment
        data_file.readline()
        #teplota taveni
        try:
            self.t_melt = float(data_file.readline())
        except ValueError:
            raise Exception("spatny datovy teploty taveni(realne cislo)")

        #comment
        data_file.readline()
        #vodivost pevne latky
        try:
            self.conductivity_solid = self.read_to_array(6,data_file)
        except ValueError:
            raise "spatny datovy vodivost pevne latky(realne cislo)"

        #comment
        data_file.readline()
        #vodivost kapalne latky
        try:
            self.conductivity_liquid = self.read_to_array(6,data_file)
        except ValueError:
            raise "spatny datovy vodivost kapalnne latky(realne cislo)"

        #comment
        data_file.readline()
        #vodivost par
        try:
            if self.ablator_exe == "Ablator.exe":
                self.conductivity_vapor = self.read_to_array(4,data_file)
            else:
                self.conductivity_vapor = self.read_to_array(1,data_file)
        except ValueError:
            raise "spatny datovy vodivost par(realne cislo)"

        #comment
        data_file.readline()
        #t vs. h pevne latky
        try:
            self.th_solid = self.read_to_array(6,data_file)
        except ValueError:
            raise "spatny datovy T vs H pevne latky(realne cislo)"

        #comment
        data_file.readline()
        #th liquid
        try:
            self.th_liquid = self.read_to_array(2,data_file)
        except ValueError:
            raise "spatny datovy th liquid(realne cislo)"
        
        #comment
        data_file.readline()
        #vodivost pevne latky
        try:
            self.th_vapor = self.read_to_array(2,data_file)
        except ValueError:
            raise "spatny datovy th vapor(realne cislo)"

        #comment
        data_file.readline()
        #log coeffs
        try:
            self.log_coeffs = self.read_to_array(2,data_file)
        except ValueError:
            raise "spatny datovy log koeffs(realne cislo)"

        #comment
        data_file.readline()
        #gruneisen coeffs
        try:
            self.gruneisen_coeffs = self.read_to_array(10,data_file)
        except ValueError:
            raise "spatny datovy typ gruneisen coeffs(realne cislo)"

        #comment
        data_file.readline()
        #poissons ratio
        try:
            self.poissons_ratio = float(data_file.readline())
        except ValueError:
            raise "spatny datovy poissons ratio(realne cislo)"

        #comment
        data_file.readline()
        #yield strength
        try:
            self.yield_strength = self.read_to_array(3,data_file)
        except ValueError:
            raise "spatny datovy yield strength(realne cislo)"

        #comment
        data_file.readline()
        #povrchove pnuti
        try:
            self.surface_tension = self.read_to_array(7,data_file)
        except ValueError:
            raise "spatny datovy typ povrchove pnuti(realne cislo)"

    def read_to_array(self, number_of_lines, file):
        x = array('f')
        for i in range(0,number_of_lines):
            a = file.readline().replace("D","E")            
            x.append(float(a))

        return x

    def scale(self, sobol, config_params, sampling):
        i = 0
        for params in config_params:
            key = params[0]
            param_arr = params[1].split(';')

            setter = getattr(self, "set_" + key)
            getter = getattr(self, "get_" + key)

            if param_arr[0] == 'p':
                if sampling == "lin":
                    #linearni uprava
                    min_val = getter()*(1-float(param_arr[1]))
                    max_val = getter()*(1+float(param_arr[1]))
                    setter(min_val+sobol[i]*(max_val-min_val))
                elif sampling == "log":
                    #log-puvodni
                    setter(getter()*math.pow(float(param_arr[1]),2*sobol[i]-1))
                i += 1
            else:
                default = getter()
                indexes = param_arr[0].split(',')
                for index in indexes:
                    if sampling == "lin":
                        min_val = default[int(index)]*(1-float(param_arr[1]))
                        max_val = default[int(index)]*(1+float(param_arr[1]))
                        default[int(index)] = min_val+sobol[i]*(max_val-min_val)
                    else:
                        default[int(index)] = default[int(index)]*math.pow(float(param_arr[1]),2*sobol[i]-1)
                    i += 1

                setter(default)

            
