# -*- coding: utf8 -*-
import gtk

from ablator_many_value_button import AblatorManyValueButton

class AblatorMatdataPanel (gtk.VBox):   
    def __init__(self, solution, settings):
			
        gtk.VBox.__init__(self)
        
        self.solution = solution
        self.settings = settings
	self.policka = self.entries()

	self.pack_start(gtk.Label("MatData"),False,False,10)

        for klic, hodnota in self.policka.iteritems():
            hbox = gtk.HBox(True, 0)
            label = gtk.Label(klic)
            hbox.pack_start(label, False, False, 0)
            hbox.pack_start(hodnota, False, False,0)
            self.pack_start(hbox)

        tlacitko_uloz = gtk.Button("Ulozit mat data", stock = gtk.STOCK_SAVE)
        tlacitko_uloz.connect("clicked", self.tlacitko_uloz_clicked)
        hboxx = gtk.HBox(True, 0)
        hboxx.pack_start(tlacitko_uloz, False, False, 0)
        self.pack_start(hboxx, False, False, 0)
        
            

    def entries (self):
        es = {}
        es["name"] = gtk.Entry()
        es["density"] = gtk.Entry()
        es["t_melt"] = gtk.Entry()
        es["conductivity_solid"] = AblatorManyValueButton(6, self.solution, "conductivity_solid")
        es["conductivity_liquid"] = AblatorManyValueButton(6, self.solution, "conductivity_liquid")
        es["conductivity_vapor"] = AblatorManyValueButton(4, self.solution, "conductivity_vapor")
        es["th_solid"] = AblatorManyValueButton(6, self.solution, "th_solid")
        es["th_liquid"] = AblatorManyValueButton(2, self.solution, "th_liquid")
        es["th_vapor"] = AblatorManyValueButton(2, self.solution, "th_vapor")
        es["log_coeffs"] = AblatorManyValueButton(2, self.solution, "log_coeffs")
        es["gruneisen_coeffs"] = AblatorManyValueButton(10, self.solution, "gruneisen_coeffs")
        es["poissons_ratio"] = gtk.Entry()
        es["yield_strength"] = AblatorManyValueButton(3, self.solution, "yield_strength")
        es["surface_tension"] = AblatorManyValueButton(7, self.solution, "surface_tension")

        return es

    def set_data (self):        
	self.policka['name'].set_text(str(self.solution.get_mat_data().get_name()))
	self.policka['density'].set_text(str(self.solution.get_mat_data().get_density()))
	self.policka['t_melt'].set_text(str(self.solution.get_mat_data().get_t_melt()))
        self.policka['poissons_ratio'].set_text(str(self.solution.get_mat_data().get_poissons_ratio()))
        
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


    def tlacitko_uloz_clicked(self, widget):
        file_chooser_dialog = gtk.FileChooserDialog(title="Zvolte soubor, který chcete otevřít", action=gtk.FILE_CHOOSER_ACTION_SAVE ,parent=None, buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK, gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))

        #Vytvoreni a pridani filtru do filechooserdialogu
        def dialog_filtr (dialog, nazev, predloha):
            filtr = gtk.FileFilter()
            filtr.set_name(nazev)
            filtr.add_pattern(predloha)
            dialog.add_filter(filtr)
            dialog.set_filter(filtr)
       
        dialog_filtr(file_chooser_dialog, "Zdrojové kódy v pythonu", "*.py")
        dialog_filtr(file_chooser_dialog, "Textové soubory", "*.txt")
        dialog_filtr(file_chooser_dialog, "Všechny soubory", "*")


        #file_chooser_dialog.set_current_folder(data)
        reakce_dialogu = file_chooser_dialog.run()
        nazev_souboru = file_chooser_dialog.get_filename()
        file_chooser_dialog.destroy()

        if reakce_dialogu == gtk.RESPONSE_OK:

            #try:            
            self.set_data()
            self.settings["description"].nastav_mat_data(nazev_souboru)       
            self.solution.get_mat_data().save(nazev_souboru)
            

            #except:    
                #self.settings["description"].nastav_init_data("Nepovedlo se nacist initdata")
