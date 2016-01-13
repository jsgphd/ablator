class HDMROut:
    def __init__(self, dir):
        #self.max = []
        #for i in xrange(0,dim):
        #    self.max.append(0.0)

        #print self.max

        self.dir = dir
        self.list = []
        #self.normalized = False

    def add(self,value):
        self.list.append(value)

    def normalize(self):
        for item in self.list:
            j = 0
            for record in item:
                if record > self.max[j]:
                    self.max[j] = record
                    j+=1

        i = 0
        for item in self.list:
            j = 0
            for record in item:
                self.list[i][j] = record/self.max[j]
                j += 1
            i+=1

    def print_hdmr(self):
        #if self.normalized == False:
        #    self.normalize()
        #    self.normalized = True
        histogram = {}

        in_file = open("solutions/" + self.dir + "/hdmr_out.txt","w")        
        for item in self.list:
            for record in item:
                if record in histogram:
                    histogram[record] = histogram[record] + 1
                else:
                    histogram[record] = 1
                in_file.write(str(record) + "\t")
            in_file.write("\n")

        in_file.close

        histogram_f = open("solutions/" + self.dir + "/histogram.csv", "w")

        for key in histogram:            
            histogram_f.write(str(key) + "\t" + str(histogram[key]) + "\n")
            
        histogram_f.close()


class HDMRScanOut(HDMROut):
    def __init__(self, dir):
        #self.max = []
        #for i in xrange(0,dim):
        #    self.max.append(0.0)

        #print self.max

        self.dir = dir
        self.list = []

        #self.normalized = False

    def add(self,value):
        self.list.append(value)


    def print_hdmr(self):
        #if self.normalized == False:
        #    self.normalize()
        #    self.normalized = True
        histogram = {}

        in_file = open("solutions/" + self.dir + "/hdmr_out.txt","w")        
        for item in self.list:
            for record in item:
                #if record in histogram:
                #    histogram[record] = histogram[record] + 1
                #else:
                #    histogram[record] = 1
                in_file.write(str(record) + "\t")
            in_file.write("\n")

        in_file.close

        #histogram_f = open("solutions/" + self.dir + "/histogram.csv", "w")

        #for key in histogram:            
        #    histogram_f.write(str(key) + "\t" + str(histogram[key]) + "\n")
            
        #histogram_f.close()
        
