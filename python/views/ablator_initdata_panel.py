# -*- coding: utf8 -*-
import gtk

class AblatorInitdataPanel (gtk.VBox):   
    def __init__(self,solution,settings):
        gtk.VBox.__init__(self)

        if(solution.ablator_exe == "Ablator.exe"):
            self.policka = self.entries()
        else:            
            self.policka = self.entries2()
            
        self.solution = solution
        self.settings = settings

        self.pack_start(gtk.Label("InitData"),False,False,10)
        
        for klic, hodnota in self.policka.iteritems():
            hbox = gtk.HBox(True, 0)
            label = gtk.Label(klic)
            hbox.pack_start(label, False, False, 0)
            hbox.pack_start(hodnota, False, False,0)
            self.pack_start(hbox)

        self.tlacitko_uloz = gtk.Button("Ulozit init data", stock = gtk.STOCK_SAVE)
        self.tlacitko_uloz.connect("clicked", self.tlacitko_uloz_clicked)
        
        hboxx = gtk.HBox(True, 0)
        vboxx = gtk.VBox(True, 0)
        vboxx.pack_start(self.tlacitko_uloz, False, False, 0)
        hboxx.pack_start(vboxx, False, False, 0)
        self.pack_start(hboxx, False, False, 0)       
            

    def entries (self):
        es = {}
        es["source"] = gtk.Entry()
        es["absorbtion_coeff"] = gtk.Entry()
        es["first_zone_size"] = gtk.Entry()
        es["source_type"] = gtk.Entry()
        es["pulse_duration"] = gtk.Entry()
        es["pulse_energy"] = gtk.Entry()
        es["blackbody_temperature"] = gtk.Entry()
        es["spot_radius"] = gtk.Entry()
        es["tstop"] = gtk.Entry()

        return es

    def entries2 (self):
        es = {}
        es["source"] = gtk.Entry()
        es["absorbtion_coeff"] = gtk.Entry()
        es["first_zone_size"] = gtk.Entry()
        es["source_type"] = gtk.Entry()
        es["pulse_duration"] = gtk.Entry()
        es["pulse_energy"] = gtk.Entry()
        es["blackbody_temperature"] = gtk.Entry()
        es["spot_radius"] = gtk.Entry()
        es["yield_test"] = gtk.Entry()
        es["tstop"] = gtk.Entry()

        return es

    def set_data (self,init_data):       
        self.policka["source"].set_text(str(init_data.get_source()))
        self.policka["absorbtion_coeff"].set_text(str(init_data.get_absorbtion()))
        self.policka["first_zone_size"].set_text(str(init_data.get_first_zone_size()))
        self.policka["source_type"].set_text(str(init_data.get_source_type()))
        self.policka["pulse_energy"].set_text(str(init_data.get_pulse_energy()))
        self.policka["pulse_duration"].set_text(str(init_data.get_pulse_duration()))        
        self.policka["blackbody_temperature"].set_text(str(init_data.get_blackbody_temperature()))
        self.policka["spot_radius"].set_text(str(init_data.get_spot_radius()))        
        self.policka["tstop"].set_text(str(init_data.get_tstop()))

    def set_data2 (self,init_data):       
        self.policka["source"].set_text(str(init_data.get_source()))
        self.policka["absorbtion_coeff"].set_text(str(init_data.get_absorbtion()))
        self.policka["first_zone_size"].set_text(str(init_data.get_first_zone_size()))
        self.policka["source_type"].set_text(str(init_data.get_source_type()))
        self.policka["pulse_energy"].set_text(str(init_data.get_pulse_energy()))
        self.policka["pulse_duration"].set_text(str(init_data.get_pulse_duration()))        
        self.policka["blackbody_temperature"].set_text(str(init_data.get_blackbody_temperature()))
        self.policka["spot_radius"].set_text(str(init_data.get_spot_radius()))
        self.policka["yield_test"].set_text(str(init_data.get_yield_test()))
        self.policka["tstop"].set_text(str(init_data.get_tstop()))

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
            self.set_data(self.solution.get_init_data())
            self.settings["description"].nastav_init_data(nazev_souboru)
            #self.settings["initdatapanel"].set_data(self.init_data)
            #self.solution.set_init_data(nazev_souboru)
            if self.solution.ablator_exe == "Ablator.exe":
                self.solution.get_init_data().save(nazev_souboru)
            else:
                self.solution.get_init_data().save2(nazev_souboru)
            

            #except:    
                #self.settings["description"].nastav_init_data("Nepovedlo se nacist initdata")
        
        
        
        


        
