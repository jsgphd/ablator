# -*- coding: utf8 -*-
import gtk

class AblatorSolutionDescription (gtk.VBox):
   
    def __init__(self, settings):

        gtk.VBox.__init__(self)
        
        self.settings = settings
        self.label_opacdata_titulek = gtk.Label("Cesta k opacdatum: ")
        self.label_opacdata_cesta = gtk.Label()
        self.hbox_label_opacdata = gtk.HBox(False,0)
        self.hbox_label_opacdata.pack_start(self.label_opacdata_titulek, False, True,0)
        self.hbox_label_opacdata.pack_start(self.label_opacdata_cesta, False, True,0)

        self.label_matdata_titulek = gtk.Label("Cesta matdatum:     ")
        self.label_matdata_cesta = gtk.Label()
        self.hbox_label_matdata = gtk.HBox(False,0)
        self.hbox_label_matdata.pack_start(self.label_matdata_titulek, False, True,0)
        self.hbox_label_matdata.pack_start(self.label_matdata_cesta, False, True,0)

        self.label_initdata_titulek = gtk.Label("Cesta k initdatum:    ")
        self.label_initdata_cesta = gtk.Label()
        self.hbox_label_initdata = gtk.HBox(False,0)
        self.hbox_label_initdata.pack_start(self.label_initdata_titulek, False, True,0)
        self.hbox_label_initdata.pack_start(self.label_initdata_cesta, False, True,0)

        self.label_solver = gtk.Label("Řešic:    ")
    
        solution = self.settings["solution"].ablator_exe

        print solution
        
        if solution == None:
            self.label_solver_name = gtk.Label()
        else:
            self.label_solver_name = gtk.Label(solution)
        self.hbox_solver = gtk.HBox(False,0)
        self.hbox_solver.pack_start(self.label_solver, False, True,0)
        self.hbox_solver.pack_start(self.label_solver_name, False, True,0)

        #spojeni labelu do vboxu
        self.pack_start(self.hbox_label_opacdata, False, True, 0)
        self.pack_start(self.hbox_label_matdata, False, True, 0)
        self.pack_start(self.hbox_label_initdata, False, True, 0)
        self.pack_start(self.hbox_solver, False, True, 0)
        

    def nastav_init_data (self, filepath):
            
        self.label_initdata_cesta.set_text(filepath)
                

    def nastav_mat_data (self, filepath):


        self.label_matdata_cesta.set_text(filepath)


    def nastav_opac_data (self, filepath):

        self.label_opacdata_cesta.set_text(filepath)

    def nastav_solver(self, name):
        self.label_solver = name


    
