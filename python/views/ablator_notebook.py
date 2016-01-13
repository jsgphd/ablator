# -*- coding: utf8 -*-

import gtk
from ablator_page import AblatorPage

class AblatorNotebook (gtk.Notebook):
   
    def __init__(self):

        gtk.Notebook.__init__(self)
        self.pridat_stranku()
        

    def pridat_stranku (self,widget=None):
        self.cislo_strany = "Nové řešení"
        self.vbox_zalozka = AblatorPage(self)            

        self.popisek_cislo_strany = gtk.Label(self.cislo_strany)           
        self.insert_page(self.vbox_zalozka, self.popisek_cislo_strany, -1)
        self.show_all()

    def pridat_stranku_name (self,widget=None, name="", solution=None):
        self.cislo_strany = name
        self.vbox_zalozka = AblatorPage(self,solution)            

        self.popisek_cislo_strany = gtk.Label(self.cislo_strany)           
        self.insert_page(self.vbox_zalozka, self.popisek_cislo_strany, -1)
        self.show_all()

    def pridat_stranku_from_load (self,solution):        
        self.vbox_zalozka = AblatorPage(self, solution)            

        self.popisek_cislo_strany = gtk.Label(solution.name)           
        self.insert_page(self.vbox_zalozka, self.popisek_cislo_strany, -1)
        self.show_all()

        return self.vbox_zalozka


        
