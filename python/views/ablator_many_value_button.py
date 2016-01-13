import gtk

from ablator_many_value_entry import AblatorManyValueEntry


class AblatorManyValueButton (gtk.Button):

    def __init__(self, pocetpolicek, solution = None, name = ""):        
        gtk.Button.__init__(self, "nastav")
        self.c_name = name
        self.solution = solution
        self.pocetpolicek = pocetpolicek
        self.connect("clicked", self.otevri_okno_clicked)        

    def otevri_okno_clicked (self, widget):
        self.okno = AblatorManyValueEntry(self.pocetpolicek, self.solution, self.c_name)
        self.okno.show_all()

    #def set_data_to_child (self, pole):
        #self.okno = AblatorManyValueEntry(self.pocetpolicek, self.solution, self.c_name)

        

        #if self.name == "conductivity_solid":
         #   self.okno.set_data(self.solution.get_mat_data().get_conductivity_solid())
        
        #self.policka["conductivity_solid"].set_data_to_child(self.solution.get_mat_data().get_conductivity_solid())
	#self.policka["conductivity_liquid"].set_data_to_child(self.solution.get_mat_data().get_conductivity_liquid())
	#self.policka["conductivity_vapor"].set_data_to_child(self.solution.get_mat_data().get_conductivity_vapor())
	#self.policka["th_solid"].set_data_to_child(self.solution.get_mat_data().get_th_solid())
	#self.policka["th_liquid"].set_data_to_child(self.solution.get_mat_data().get_th_liquid())
        #self.policka["th_vapor"].set_data_to_child(self.solution.get_mat_data().get_th_vapor())
	#self.policka["log_coeffs"].set_data_to_child(self.solution.get_mat_data().get_log_coeffs())
	#self.policka["gruneisen_coeffs"].set_data_to_child(self.solution.get_mat_data().get_gruneisen_coeffs())
	#self.policka["yield_strength"].set_data_to_child(self.solution.get_mat_data().get_yield_strength())
	#self.policka["surface_tension"].set_data_to_child(self.solution.get_mat_data().get_surface_tension())
        
        #self.okno.set_data(pole)
        
