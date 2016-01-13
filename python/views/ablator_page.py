# -*- coding: utf8 -*-

import gtk

from models.solution import Solution
from models.init_data import InitData
from models.mat_data import MatData

from ablator_solution_description import AblatorSolutionDescription
from ablator_initdata_panel import AblatorInitdataPanel
from ablator_matadata_panel import AblatorMatdataPanel
from ablator_save_window import AblatorSaveWindow
#from ablator_opacdata_panel import AblatorOpacdataPanel



class AblatorPage (gtk.VBox):

    def __init__(self, parent_notebook, solution = None):
        
        gtk.VBox.__init__(self, False, 0)

        self.notebook = parent_notebook

        if(solution == None):
            self.solution = Solution()
        else:
            self.solution = solution

        self.settings = {}
        self.settings["initdatapanel"]=AblatorInitdataPanel(self.solution,self.settings)
        self.settings["matdatapanel"]=AblatorMatdataPanel(self.solution, self.settings)
        self.settings["opacdatacesta"]= ""
        self.settings["solution"]=self.solution
        self.settings["description"]=AblatorSolutionDescription(self.settings)
        

        self.Init_tlacitko = gtk.Button("Init Data..")
        self.Init_tlacitko.connect("clicked", self.cteni_initdat_clicked)

        self.tlacitko_smazat_page = gtk.Button ("Zavřít záložku", stock = gtk.STOCK_CLOSE)
        self.tlacitko_smazat_page.connect("clicked", self.smazat_page_clicked)
        
        self.Mat_tlacitko = gtk.Button("Mat Data..")
        self.Mat_tlacitko.connect("clicked", self.cteni_matdata_clicked)
        
        self.Opac_tlacitko = gtk.Button("Opac Data..")
        self.Opac_tlacitko.connect("clicked", self.cteni_opacdata_clicked)

        self.tlacitko_vypocet = gtk.Button("Vypocet", stock = gtk.STOCK_EXECUTE)
        self.tlacitko_vypocet.connect("clicked", self.spustit_vypocet_clicked)

        self.tlacitko_otevrit_slozku = gtk.Button("Otevřít složku s řešením")
        self.tlacitko_otevrit_slozku.connect("clicked", self.otevrit_slozku_clicked)

        self.tlacitko_ulozit_data = gtk.Button("Ulozit data..", stock = gtk.STOCK_SAVE)
        self.tlacitko_ulozit_data.connect("clicked", self.ulozit_data_clicked)

        self.combobox = gtk.combo_box_new_text()
        self.combobox.append_text("Ablator - 1.0")
        self.combobox.append_text("Ablator - 2.0")
        self.combobox.set_active(0)


        #spojeni tlacitek do hboxu
        self.hbox_zalozka_tlacitka = gtk.HBox()
        self.hbox_zalozka_tlacitka.pack_start(self.Init_tlacitko, False, True,0)
        self.hbox_zalozka_tlacitka.pack_start(self.Mat_tlacitko, False, True,0)
        self.hbox_zalozka_tlacitka.pack_start(self.Opac_tlacitko, False, True,0)
        self.hbox_zalozka_tlacitka.pack_start(self.tlacitko_otevrit_slozku, False, True,0)
        self.hbox_zalozka_tlacitka.pack_start(self.tlacitko_ulozit_data, False, True, 0)
        self.hbox_zalozka_tlacitka.pack_start(self.tlacitko_smazat_page, False, True, 0)
        self.hbox_zalozka_tlacitka.pack_start(self.tlacitko_vypocet, False, True,0)
        #self.hbox_zalozka_tlacitka.pack_start(self.combobox, False, True,0)
        
        #spojeni initdata a matdata do hboxu

        self.hbox_initdata_matdata_opacdata = gtk.HBox()
        self.hbox_initdata_matdata_opacdata.pack_start(self.settings["initdatapanel"], True, True,0)
        self.hbox_initdata_matdata_opacdata.pack_start(self.settings["matdatapanel"], True, True,0)
        #self.hbox_initdata_matdata_opacdata.pack_start(self.settings["opacdatapanel"], True, True,0)
        

        self.pack_start(self.hbox_zalozka_tlacitka, False, False, 0)
        self.pack_start(self.hbox_initdata_matdata_opacdata, False, False, 0)      
        self.pack_end(self.settings["description"], False, True, 0)




    def cteni_initdat_clicked (self, widget):

        file_chooser_dialog = gtk.FileChooserDialog(title="Zvolte soubor, který chcete otevřít", parent=None, buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK, gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))

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
            #self.init_data = InitData(nazev_souboru)
            self.solution.set_init_data(nazev_souboru)
            self.settings["description"].nastav_init_data(nazev_souboru)
            
            if self.solution.ablator_exe == "Ablator.exe":
                self.settings["initdatapanel"].set_data(self.solution.init_data)                
            else:
                self.settings["initdatapanel"].set_data2(self.solution.init_data)
                
            self.solution.set_init_data(nazev_souboru)
            

            #except:    
                #self.settings["description"].nastav_init_data("Nepovedlo se nacist initdata")
        

    def cteni_matdata_clicked (self, widget):

        file_chooser_dialog = gtk.FileChooserDialog(title="Zvolte soubor, který chcete otevřít", parent=None, buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK, gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))

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
            self.solution.set_mat_data(nazev_souboru)            
            self.settings["description"].nastav_mat_data(nazev_souboru)            
            self.settings["matdatapanel"].set_data()
            

            #except:    
                #self.settings["description"].nastav_mat_data("Nepovedlo se nacist matdata")        


    def cteni_opacdata_clicked (self, widget):
        file_chooser_dialog = gtk.FileChooserDialog(title="Zvolte soubor, ktery chcete otevrit", parent=None, buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK, gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))
    
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
            self.solution.set_opac_data(nazev_souboru)  
            self.settings["description"].nastav_opac_data(nazev_souboru)

            #except:
            #    self.settings["description"].nastav_opac_data("Nepovedlo se nacist opacdata")
                

    def spustit_vypocet_clicked (self, widget):        
        #s = Solution("Kremik","data/initdata","data/matdata-Si","data/opacdata-PMMA","Ablator2.exe")
        
        #model = self.combobox.get_model()
        #active = self.combobox.get_active()
        
        #item = model[active][0]
        #if(item == "Ablator - 1.0"):
        #    self.solution.set_ablator_exe("Ablator.exe")
        #else:
        #    self.solution.set_ablator_exe("Ablator2.exe")
        try:
            self.solution.run()
            self.solution.open_folder()
        except Exception, e:            
            md = gtk.MessageDialog(None, 
            gtk.DIALOG_DESTROY_WITH_PARENT, gtk.MESSAGE_INFO, 
            gtk.BUTTONS_CLOSE, "Nejdříve řešení uložte.")
            md.run()
            md.destroy()
        

    def otevrit_slozku_clicked (self, widget):
        self.solution.open_folder()

    def smazat_page_clicked (self,widget):        
        strana = self.notebook.get_current_page()
        self.notebook.remove_page(strana)
        self.notebook.queue_draw_area(0,0,-1,-1)

    def ulozit_data_clicked (self, widget):        
        AblatorSaveWindow(self.solution, self.notebook)
        

        
        
        
