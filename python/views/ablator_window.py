import gtk

#coding:utf8

#from models.solution import Solution
#from models.init_data import InitData
#from models.mat_data import MatData
from ablator_menubox import AblatorMenuBox
from ablator_notebook import AblatorNotebook

class AblatorWindow (gtk.Window):

    def __init__(self):

        gtk.Window.__init__(self)

        self.settings = {}
        self.settings["notebook"]=AblatorNotebook()
        self.settings["menubox"]=AblatorMenuBox(self.settings)
        self.set_title("Ablator")
        self.connect ("destroy", gtk.main_quit)
        self.set_position(gtk.WIN_POS_CENTER)
        self.set_border_width(1)
        self.set_resizable(True)
        self.resize (800 ,600)
        self.vbox = gtk.VBox(False, 0)
        self.vbox.pack_start(self.settings["menubox"], False, False, 0)
        self.vbox.pack_start(self.settings["notebook"], False, False, 0)
        
        self.add(self.vbox)
        

        self.show_all()
