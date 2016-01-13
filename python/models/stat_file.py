class StatFile:
    def __init__(self, directory, hdmr, sim_num):
        self.directory = directory
        self.hdmr = hdmr
        self.sim_num = sim_num

    def save(self):
        in_file = open("solutions/" + self.directory + "/stat.csv","w")
        in_file.write(str(self.sim_num) + "\t")
        in_file.write(str(self.hdmr.get_size()))

        in_file.close()
