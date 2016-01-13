from views.ablator_window import *

from models.mat_data import *
from models.init_data import *

def spust():
        settings = gtk.settings_get_default()
        settings.props.gtk_button_images = True
        gtk.main()    

if __name__ == "__main__" :
    aplikace = AblatorWindow()
    spust()
