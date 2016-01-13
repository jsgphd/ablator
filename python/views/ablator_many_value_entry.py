import gtk

class AblatorManyValueEntry (gtk.Window):

    def __init__(self, pocetpolicek, solution, name):

        gtk.Window.__init__(self)
        self.attribute = name
        self.solution = solution
        self.pocetpolicek = pocetpolicek
        self.resize (400 ,250)
        self.set_title(name)
        self.set_position(gtk.WIN_POS_CENTER)
        self.set_resizable(True)
        

        self.tlacitko_ulozit_data = gtk.Button("Ulozit", stock = gtk.STOCK_SAVE)
        self.tlacitko_ulozit_data.connect("clicked", self.ulozit_data_clicked)

        self.tlacitko_zrusit_okno = gtk.Button("Zrusit", stock = gtk.STOCK_CLOSE)
        self.tlacitko_zrusit_okno.connect("clicked", self.zrusit_okno_clicked)

        self.hbox_manyvalueenry_tlacitka = gtk.HBox()
        self.hbox_manyvalueenry_tlacitka.pack_start(self.tlacitko_ulozit_data, False, True,0)
        self.hbox_manyvalueenry_tlacitka.pack_start(self.tlacitko_zrusit_okno, False, True,0)

        self.vbox = gtk.VBox(False,0)
        self.vbox.pack_start(self.hbox_manyvalueenry_tlacitka, False, True, 0)

        self.seznam_poli = []

        for i in range(0, pocetpolicek):
            entry = gtk.Entry()
            self.seznam_poli.append(entry)
            self.vbox.pack_start(entry, False, False,0)

        if(self.solution.get_mat_data() != None):
            self.set_initial_data()

        self.add(self.vbox)

    def set_initial_data(self):
        if self.attribute == "conductivity_solid":
            self.set_data(self.solution.get_mat_data().get_conductivity_solid())
        elif self.attribute == "conductivity_liquid":
            self.set_data(self.solution.get_mat_data().get_conductivity_liquid())            
        elif self.attribute == "conductivity_vapor":
            self.set_data(self.solution.get_mat_data().get_conductivity_vapor())                
        elif self.attribute ==  "th_solid":
            self.set_data(self.solution.get_mat_data().get_th_solid())            
        elif self.attribute == "th_liquid":
            self.set_data(self.solution.get_mat_data().get_th_liquid())            
        elif self.attribute == "th_vapor":
            self.set_data(self.solution.get_mat_data().get_th_vapor())                
        elif self.attribute == "log_coeffs":
            self.set_data(self.solution.get_mat_data().get_log_coeffs())
        elif self.attribute ==  "gruneisen_coeffs":
            self.set_data(self.solution.get_mat_data().get_gruneisen_coeffs())
        elif self.attribute == "yield_strength":
            self.set_data(self.solution.get_mat_data().get_yield_strength())
        elif self.attribute == "surface_tension":
            self.set_data(self.solution.get_mat_data().get_surface_tension())        
        

    def set_data (self, pole):        
        for i in range(0, len(pole)):
            self.seznam_poli[i].set_text(str(pole[i]))

    def ulozit_data_clicked (self, widget):
        values = []
        for i in range(0,len(self.seznam_poli)):
            values.append(float(self.seznam_poli[i].get_text()))

        if self.attribute == "conductivity_solid":
            self.solution.get_mat_data().set_conductivity_solid(values)            
        elif self.attribute == "conductivity_liquid":
            self.solution.get_mat_data().set_conductivity_liquid(values)            
        elif self.attribute == "conductivity_vapor":
            self.solution.get_mat_data().set_conductivity_vapor(values)                
        elif self.attribute ==  "th_solid":
            self.solution.get_mat_data().set_th_solid(values)            
        elif self.attribute == "th_liquid":
            self.solution.get_mat_data().set_th_liquid(values)            
        elif self.attribute == "th_vapor":
            self.solution.get_mat_data().set_th_vapor(values)                
        elif self.attribute == "log_coeffs":
            self.solution.get_mat_data().set_log_coeffs(values)            
        elif self.attribute ==  "gruneisen_coeffs":
            self.solution.get_mat_data().set_gruneisen_coeffs(values)            
        elif self.attribute == "yield_strength":
            self.solution.get_mat_data().set_yield_strength(values)            
        elif self.attribute == "surface_tension":
            self.solution.get_mat_data().set_surface_tension(values)        

    def zrusit_okno_clicked (self, widget):
        self.destroy()
            
        

        

        
        
