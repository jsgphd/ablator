#slozka ve ktere se vytvori adresar se vstupy
OUT_DIR = "scan_configs"
#slozka v OUT_DIR, mela by obsahovat nazev reseni
SOL_NAME = "lif"

import os

def prepare(min_val, max_val, step):
    work_dir,init_dir = prepare_work_dir()
    
    val_range = max_val - min_val
    steps_no = int(val_range/step) + 1

    for i in range(steps_no):
        #vygeneruj initdata
        s_i = prepare_init_data(min_val + i*step)
        s_i_path = init_dir + "/" + str(i) + "_" + SOL_NAME + "_initdata"
        f = open(s_i_path,"w")
        f.write(s_i)
        f.close()
                 
        #vygeneruj config sensitivity
        s_c = prepare_config(s_i_path)
        f = open(work_dir + "/" + str(i) + "_config_sensitivity.txt","w")
        f.write(s_c)
        f.close()
        
    print("Done!")


def prepare_work_dir():
    if(not os.path.isdir(OUT_DIR)):
        os.mkdir(OUT_DIR)

    work_dir = OUT_DIR + "/" + SOL_NAME
    if(not os.path.isdir(work_dir)):
        os.mkdir(work_dir)

    init_work_dir = work_dir + "/init" 
    if(not os.path.isdir(init_work_dir)):
        os.mkdir(init_work_dir)

    return [work_dir,init_work_dir]
        
    
def prepare_init_data(energy):
    s = """\
# BB X-ray or laser source (1=Laser)
1
# absorbtion coefficient in microns for laser source
0.080
# the first zone size (nm)
1
# source type (1=square pulse;2=Gaussian pulse)
1
# pulse duration or pulse FWHM (ns)
0.001
# total energy in the pulse (J)
"""

    s += str(energy) + "\n"

    s += """\
# blackbody temperature for X-ray source (keV)
0
# radius of the irradiated spot (microns) (5.6E3 microns to give the area 1cm2)
5.6D+3
# Tstop (ns)
500
"""
    return s

def prepare_config(init_path):
    s = """\
#atribut nemenit, pokud je dano atribut=p, pak se bere primo
#hodnota(ma se za to, ze hodnota je popsana pouze jednim parametrem)
#pokud ma pole vice hodnot, pak se jeho indexem(prvni prvek je oznacen 0, dalsi 1,2,3,...)
#zapisuji porade prvky, oddelene carkou, bez mezer, ktere se budou vybirat, viceprvkove atributy 
executable=Ablator.exe
#pro vzorkovani v logaritm. prostoru-> sampling=lin
sampling=lin
init_data="""

    s += init_path + "\n"

    s += """\
mat_data=data/sio2/sio2_matdata
opac_data=data/sio2/sio2_opacdata
density=p;0.1
t_melt=p;0.1
conductivity_solid=0;0.1
conductivity_liquid=0;0.1
conductivity_vapor=0;0.1
th_solid=0;0.1
th_liquid=0;0.1
th_vapor=0;0.1
log_coeffs=0,1;0.05
gruneisen_coeffs=1,2,5,6,8,9;0.1
poissons_ratio=p;0.1
surface_tension=2,3;0.1
    """

    return s

min_value = 0.02
max_value = 0.1
step_size = 0.02
prepare(min_value, max_value, step_size)
