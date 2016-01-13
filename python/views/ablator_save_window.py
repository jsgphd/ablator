# -*- coding: utf8 -*-
import gtk

class AblatorSaveWindow (gtk.Window):

    def __init__(self, solution, notebook):

        gtk.Window.__init__(self)        
        self.solution = solution
        self.notebook = notebook
        self.resize (400 ,250)
        self.set_title("Název řešení")
        self.set_position(gtk.WIN_POS_CENTER)
        self.set_resizable(True)

        self.vbox = gtk.VBox(False,0)

        self.tlacitko_ulozit_data = gtk.Button(label="Vytvořit", stock = gtk.STOCK_SAVE)
        self.tlacitko_ulozit_data.set_label("Vytvořit")
        self.tlacitko_ulozit_data.connect("clicked", self.ulozit_data_clicked)

        self.tlacitko_zrusit_okno = gtk.Button("Zrušit", stock = gtk.STOCK_CLOSE)
        self.tlacitko_zrusit_okno.connect("clicked", self.zrusit_okno_clicked)

        self.hbox = gtk.HBox(True,0)
        self.hbox.pack_start(self.tlacitko_ulozit_data,False,False,0)
        self.hbox.pack_start(self.tlacitko_zrusit_okno,False,False,0)

        self.label = gtk.Label("Název řešení")

        if(self.solution.name == None):
            self.name_entry = gtk.Entry()
        else:
            self.name_entry = gtk.Entry()
            self.name_entry.set_text(self.solution.name)

        self.vbox.pack_start(self.label,False,False,10)
        self.vbox.pack_start(self.name_entry,False,False,10)
        self.vbox.pack_end(self.hbox, False,False,0)

        self.add(self.vbox)
        self.show_all()
        
    def zrusit_okno_clicked (self, widget):
        self.destroy()

    def ulozit_data_clicked (self, widget):
        self.solution.set_name(self.name_entry.get_text())
        self.notebook.set_tab_label_text(self.notebook.get_nth_page(self.notebook.get_current_page()), self.name_entry.get_text())
        self.destroy()
            

class AblatorNewSaveWindow (gtk.Window):

    def __init__(self, solution, notebook):

        gtk.Window.__init__(self)        
        self.solution = solution
        self.notebook = notebook
        self.resize (400 ,250)
        self.set_title("Název řešení")
        self.set_position(gtk.WIN_POS_CENTER)
        self.set_resizable(True)

        self.vbox = gtk.VBox(False,0)

        self.tlacitko_ulozit_data = gtk.Button("Uložit", stock = gtk.STOCK_SAVE)
        self.tlacitko_ulozit_data.connect("clicked", self.ulozit_data_clicked)

        self.tlacitko_zrusit_okno = gtk.Button("Zrušit", stock = gtk.STOCK_CLOSE)
        self.tlacitko_zrusit_okno.connect("clicked", self.zrusit_okno_clicked)

        self.hbox = gtk.HBox(True,0)
        self.hbox.pack_start(self.tlacitko_ulozit_data,False,False,0)
        self.hbox.pack_start(self.tlacitko_zrusit_okno,False,False,0)

        self.label = gtk.Label("Název řešení")

        if(self.solution.name == None):
            self.name_entry = gtk.Entry()
        else:
            self.name_entry = gtk.Entry()
            self.name_entry.set_text(self.solution.name)

        self.combobox = gtk.combo_box_new_text()
        self.combobox.append_text("Ablator")
        self.combobox.append_text("Ablator - XUV")
        self.combobox.set_active(0)

        self.vbox.pack_start(self.label,False,False,10)
        self.vbox.pack_start(self.name_entry,False,False,10)
        self.vbox.pack_start(self.combobox, False, False, 10)
        self.vbox.pack_end(self.hbox, False,False,0)

        self.add(self.vbox)
        self.show_all()
        
    def zrusit_okno_clicked (self, widget):
        self.destroy()

    def ulozit_data_clicked (self, widget):
        model = self.combobox.get_model()
        active = self.combobox.get_active()
        
        item = model[active][0]
        if(item == "Ablator"):
            self.solution.set_ablator_exe("Ablator.exe")            
        else:            
            self.solution.set_ablator_exe("AblatorXUV.exe")        

        self.solution.set_name(self.name_entry.get_text())
        self.notebook.pridat_stranku_name(name=self.name_entry.get_text(),solution=self.solution)
        self.notebook.set_current_page(self.notebook.get_n_pages()-1)        
        #self.notebook.set_tab_label_text(self.notebook.get_nth_page(self.notebook.get_current_page()), self.name_entry.get_text())
        self.destroy()
        

        

        
        
