class HDMR:
    def __init__(self, dir=""):
        #self.max = []
        #for i in xrange(0,dim):
        #    self.max.append(0.0)

        #print self.max
        self.dir = dir
        self.list = []
        #self.normalized = False

    def add(self,value):
        self.list.append(value)

    def get_size(self):
        return len(self.list)

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

        in_file = open("solutions/" + self.dir + "/hdmr_in.txt","w")        
        for item in self.list:
            for record in item:
                in_file.write(str(record) + "\t")
            in_file.write("\n")

        in_file.close
        
#hdmr = HDMR()
#hdmr.add(10.0)
#hdmr.add(1.0)
#hdmr.print_hdmr()
