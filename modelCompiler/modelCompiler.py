


from pymodelica import compile_fmux,compile_fmu
from casadi import SymbolicOCP
from pyfmi import load_fmu

import os
import zipfile



def compileModelica(class_name,models_list):
    
    fmux_name = compile_fmux(class_name=class_name,
                             file_name=models_list,
                             compiler='optimica',
                             compiler_options={'inline_functions':'all',
                                               'generate_html_diagnostics':True,
                                               'eliminate_alias_variables':False},
                             compile_to='.',
                             compiler_log_level='w')

    class_name = class_name.replace('.','_')
    sfile = zipfile.ZipFile('./'+class_name+'.fmux','r')
    mfile = sfile.extract('modelDescription.xml','.')

    if os.path.isfile(class_name+'.xml'):
        os.remove(class_name+'.xml')       
    os.rename('modelDescription.xml',class_name+'.xml')       
    


def getOCP(class_name,models_list):
    
    #try:
    compileModelica(class_name,models_list) 
    #except:
    #     print "Oops!  compilation was not possible"   

    ocp = SymbolicOCP()
    ocp.parseFMI(class_name+'.xml')
    return ocp;


def getModelList(path='.',withSubDirs=False):

    models_list = [];
    
    if withSubDirs:
    
        for path, subdirs, files in os.walk(path):
           for filename in files:
             f = os.path.join(path, filename)
             if f.endswith(".mop") or f.endswith(".mo"):
                 models_list.append(f)    
        
    else:    
        
        files = [f for f in os.listdir(path) if os.path.isfile(f)]
        for f in files:
             if f.endswith(".mop") or f.endswith(".mo"):
                 models_list.append(f)
    
             
    return models_list;
    
    
def compileForSimulation(class_name,models_list):
    
    init_sim_fmu = compile_fmu(class_name=class_name,
                               file_name=models_list,
                               compiler_options={'inline_functions':'all',
                                               'generate_html_diagnostics':True,
                                               'eliminate_alias_variables':False},
                               compile_to='.',
                               compiler_log_level='w')
                               
    init_sim_model = load_fmu(init_sim_fmu)

                               
    return init_sim_model
