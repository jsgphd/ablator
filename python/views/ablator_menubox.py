# -*- coding: utf8 -*-

import gtk

from models.solution import Solution
from ablator_save_window import AblatorNewSaveWindow

class AblatorMenuBox (gtk.HBox):

    def __init__(self, settings):

        
        gtk.HBox.__init__(self, False, 0)

        self.settings = settings

        self.tlacitko_novadata = gtk.Button("Nova Data..", stock = gtk.STOCK_NEW)
        self.tlacitko_novadata.connect("clicked", self.pridat_stranku_clicked)

        self.tlacitko_nahratdata = gtk.Button("Nahrat data..", stock = gtk.STOCK_OPEN)  
        self.tlacitko_nahratdata.connect("clicked", self.nahrat_data_clicked)
               
        self.pack_start(self.tlacitko_novadata, False, True, 0)
        self.pack_start(self.tlacitko_nahratdata, False, True, 0)        


#____________________________________________________________________________#

        
    def pridat_stranku_clicked (self, widget):
        solution = Solution()
        AblatorNewSaveWindow(solution, self.settings["notebook"])
        #self.settings["notebook"].set_current_page(self.settings["notebook"].get_n_pages()-1)        
        #self.settings["notebook"].set_current_page(-1)
        #self.settings["notebook"].pridat_stranku()



    def nahrat_data_clicked (self, widget):
        
        file_chooser_dialog = gtk.FileChooserDialog(title="Zvolte soubor, ktery chcete otevrit", parent=None, buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK, gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))


     #Vytvoreni a pridani filtru do filechooserdialogu
        def dialog_filtr (dialog, nazev, predloha):
            filtr = gtk.FileFilter()
            filtr.set_name(nazev)
            filtr.add_pattern(predloha)
            dialog.add_filter(filtr)
            dialog.set_filter(filtr)
       
        dialog_filtr(file_chooser_dialog, "Zdrojove kody v pythonu", "*.py")
        dialog_filtr(file_chooser_dialog, "Vsechny soubory", "*")
        dialog_filtr(file_chooser_dialog, "Textove soubory", "*.txt")


        #file_chooser_dialog.set_current_folder(solution)
        reakce_dialogu = file_chooser_dialog.run()
        nazev_souboru = file_chooser_dialog.get_filename()
        file_chooser_dialog.destroy()

        if reakce_dialogu == gtk.RESPONSE_OK:
            solution = Solution(name=nazev_souboru)
            solution.load()
            page = self.settings["notebook"].pridat_stranku_from_load(solution)
            page.settings["initdatapanel"].set_data(solution.get_init_data())
            page.settings["matdatapanel"].set_data()

            page.settings["description"].nastav_init_data(solution.init_data.filepath)
            page.settings["description"].nastav_mat_data(solution.mat_data.filepath)
            page.settings["description"].nastav_opac_data(solution.opac_data)

            self.settings["notebook"].set_current_page(self.settings["notebook"].get_n_pages()-1)
            
        



        

  

        
    

            
