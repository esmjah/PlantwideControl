###OLGA_SIM
import OpenOPC
import time
import datetime
import numpy as NP

model = 'OlgaNet.ServerDemo.'
minimumPauseSec = 0.01
maxPauseSec = 1



def OLGA_restart_time(da,block=True):
    da.write((model+'PAUSE',1))
    tinit = datetime.datetime(2000,01,13,00,00,00)
    initTime = tinit.strftime("%m/%d/%y %H:%M:%S")
    da.write((model+'ExternalClock',initTime))
    da.write((model+'INITTIME',initTime))
    da.write((model+'SIMTIME',initTime))
    da.write((model+'PAUSE',0))
    ok = False;
    time.sleep(1)    
    if block:
        for k in range(10):
            m = da.read(model+"PAUSE")[0] # can we do a better check?
            if m is False:
                ok = True;                
                break
            else:
                time.sleep(1)
    else:
        ok = True

    return ok

def OLGA_connect(serverIP='127.0.0.1',modelPrefix='Sim.ServerDemo.'):  #Create Data Access Object from OPC server

    model = modelPrefix;                   # turn his into a class
    da = OpenOPC.client()
    id = 'SPT.OLGAOPCServer.7'
    da.connect(id)
    time.sleep(1)
    #OLGA_restart_time(da,True)
    return da


def OLGA_load_state(da,fileName='initial_simulation_conditions',block=True):
    da.write((model+'PAUSE',1))
    da.write((model+'ExternalClock',da.read(model+'INITTIME')[0]))
    da.write((model+'LoadSnap.File',fileName))
    da.write((model+'LoadSnap',1))
    ok = False;
    time.sleep(1)    
    if block:
        for k in range(10):
            m = da.read(model+"LastMessage")[0]
            if m == ('INFO Snap file '+"'"+fileName+"'"+', successfully loaded.'):
                ok = True;                
                break
            else:
                time.sleep(1)
    else:
        ok = True

    da.write((model+'PAUSE',0))
    return ok                
        

    
def OLGA_save_state(da,fileName='initial_simulation_conditions',block=True):
    da.write((model+'PAUSE',1))
    da.write((model+'SaveSnap.File',fileName))
    da.write((model+'SaveSnap',1))
    da.write((model+'PAUSE',0))
    ok = False;
    time.sleep(1)    
    if block:
        for k in range(10):
            m = da.read(model+"PAUSE")[0] # can we do a better check?
            if m == False:
                ok = True;                
                break
            else:
                time.sleep(1)
    else:
        ok = True
    
    return ok
    

def OLGA_simulate_step(da,stepsize=1,block=True):

    #Read External Clock On Olga Server
    ex_time_read = da.read(model+'ExternalClock')[0]

    #Convert read variables to datetime datatype
    ex_time = datetime.datetime.strptime(ex_time_read,"%m/%d/%y %H:%M:%S")
    
    DT = datetime.timedelta(seconds=stepsize)

    step = ex_time + DT
    
    
    if (step-ex_time).seconds is not DT.seconds:
        print 'Something went wrong with the time calculation.  Investigate what is going on with the OPC server'
        print 'Time read from olga', ex_time_read
        print 'Time to be sent for simulation', step
        print 'Dt in seconds', stepsize
        assert((step-ex_time).seconds is DT.seconds)

    #Convert back to string
    ex_time_write = step.strftime("%m/%d/%y %H:%M:%S")

    #Update External Clock

    da.write((model+'ExternalClock',ex_time_write))

    while block:
        sim_time_read = da.read(model+'SIMTIME')[0]
        sim_time = datetime.datetime.strptime(sim_time_read,"%m/%d/%y %H:%M:%S")
        deltaTime = ex_time - sim_time
        if (deltaTime.total_seconds() <= 0):
            if da.read(model+"State")[0]=='STATE_RUNNING_READY':
                break

        if (deltaTime.total_seconds() > 1.1*stepsize):
            print 'Something went wrong with the time calculation.  Investigate what is going on with the OPC server'
            print 'Time read from olga', ex_time_read
            print 'Time to be sent for simulation', step
            print 'current time being simulated', sim_time_read
            print 'Dt in seconds', stepsize
            assert(False)
                  
        speed_read = da.read((model+'SPEED'))[0]
        
        if speed_read == 0:
            time.sleep(minimumPauseSec)
        else:
            remaingWallClockSec = deltaTime.seconds/speed_read
            time.sleep(min(maxPauseSec,max(0.9*remaingWallClockSec,minimumPauseSec)))
            
        
    


def OLGA_get_init(da):
    t = da.read(model+"TIME")[0]
    m_ga0 = da.read(model+"Annulus.GASMASS")[0]
    m_gt0 = da.read(model+"PRODUCTION Tubing.GASMASS")[0]
    m_ot0 = da.read(model+"PRODUCTION Tubing.OILMASS")[0]
    m_wt0 = da.read(model+"PRODUCTION Tubing.WATMASS")[0]
    m_lt0 = m_ot0 + m_wt0
    T_a = 273.15 + da.read(model+"Annulus.MEANTMBRCT")[0]
    T_t = 273.15 + da.read(model+"PRODUCTION Tubing.MEANTMBRCT")[0]
    p_ai = da.read(model+"Annulus.PT")[0][142]
    p_ti = da.read(model+"PRODUCTION Tubing.PT")[0][4]

    return [t,m_ga0,m_gt0,m_ot0,m_wt0,m_lt0,T_a,T_t,p_ai,p_ti]

def OLGA_get_pipeline_measurements_not_working(da):

    results = {}
    results['t'] = da.read(model+"TIME")[0]
    results['p_inlet'] = da.read(model+'BRAN-1.PT')[0][0]
    results['p_lp'] = da.read(model+'BRAN-1.PT')[0][460]
    results['p_tp'] = da.read(model+'BRAN-1.PT')[0][520]
    results['p_outlet'] = da.read(model+'BRAN-1.PT')[0][540]
    results['w_g_out'] =  da.read(model+'BRAN-1.GG')[0][540]
    results['w_o_out'] =   da.read(model+'BRAN-1.GLTHL')[0][540]
    results['w_w_out'] =  da.read(model+'BRAN-1.GLTWT')[0][540]

    return results

def OLGA_get_pipeline_measurements(da):

    t = da.read(model+"TIME")[0]
    p_inlet = da.read(model+'BRAN-1.PT')[0][0]
    p_lp = da.read(model+'BRAN-1.PT')[0][460]
    p_tp = da.read(model+'BRAN-1.PT')[0][520]
    p_outlet = da.read(model+'BRAN-1.PT')[0][540]
    w_g_out =  da.read(model+'BRAN-1.GG')[0][540]
    w_o_out =   da.read(model+'BRAN-1.GLTHL')[0][540]
    w_w_out =  da.read(model+'BRAN-1.GLTWT')[0][540]

    return [t,p_inlet,p_lp,p_tp,p_outlet,w_g_out,w_o_out,w_w_out,w_w_out+w_o_out]

    
def OLGA_get_pipeline_measurements2(da):
    N = 16
    results = NP.zeros(N)
    results[0] = da.read(model+"TIME")[0]
    results[1] = da.read(model+'INLET.PT')[0]
    results[2] = da.read(model+'LOW POINT.PT')[0]
    results[3] = da.read(model+'OUTLET.PT')[0]
    results[4] = da.read(model+'OUTLET.PT')[0]
    results[5] =  da.read(model+'OUTLET.GG')[0]
    results[6] =  da.read(model+'OUTLET.GLTHL')[0]
    results[7] =  da.read(model+'OUTLET.GLTWT')[0]
    results[8] =  da.read(model+'OUTLET.GLT')[0]
    results[9] =  da.read(model+'BRAN-1.GASMASS')[0]
    results[10] = da.read(model+'BRAN-1.WATMASS')[0]
    results[11] = da.read(model+'BRAN-1.LIQMASS')[0]
    results[12] = da.read(model+'SOUR-1-1.MASSFLOW')[0]
    results[13] = da.read(model+'SOUR-1-1.GASFRACTION')[0]
    results[14] = da.read(model+'SOUR-1-1.WATERFRACTION')[0]
    results[15] = da.read(model+'CHOKE-1-1.OPENING')[0] 
    
    return results



def OLGA_read_tags(da,tags):
    
        
    
    N = len(tags)
    results = NP.zeros(N)
    
    i = 0
    for i in range(N):
        results[i] = da.read(model+tags[i])[0]
        
    return results


def OLGA_write_tags(da,newValues):
    ok = False
    for k in newValues.keys():
        ok = da.write((model+k,float(newValues[k])))
        if ok is False:
            return ok
    
    return ok




def OLGA_read_pipeline(da):
    
    tags = ["TIME",
    'INLET.PT',
    'LOW POINT.PT',
    'OUTLET.PT',
    'OUTLET.GG',
    'OUTLET.GLTHL',
    'OUTLET.GLTWT',
    'OUTLET.GLT',
    'BRAN-1.GASMASS',
    'BRAN-1.WATMASS',
    'BRAN-1.LIQMASS',
    'SOUR-1-1.MASSFLOW',
    'SOUR-1-1.GASFRACTION',
    'SOUR-1-1.WATERFRACTION',
    'CHOKE-1-1.OPENING'] 
        
    
    results = OLGA_read_tags(da,tags)

    data = {}        
    for k in range(len(tags)):
        data[tags[k]] = results[k]
           
    return data



def OLGA_get_well_measurements(da):
    N = 18
    results = NP.zeros(N)

    results[0]  = da.read(model+"TIME")[0]
    results[1] = da.read(model+"Annulus.GASMASS")[0] #m_ga
    results[2] = da.read(model+"PRODUCTION Tubing.GASMASS")[0] #m_gt
    results[3] = da.read(model+"PRODUCTION Tubing.OILMASS")[0] #m_ot
    results[4] = da.read(model+"PRODUCTION Tubing.WATMASS")[0] #m_wt
    results[5] = results[3]+results[4] #m_lt
    results[6] = da.read(model+"Casinghead.MASSFLOW")[0] #w_ga
    results[7] = -da.read(model+"Operation valve.GTLEAK")[0] #w_gi
    results[8] = 0 #w_gp
    results[9] = da.read(model+"PRODUCTION Tubing.GLT")[0][1] #w_lr
    results[10] = 0 #w_lp
    results[11] = da.read(model+"PRODUCTION Tubing.PT")[0][0] #p_bh
    results[12] = da.read(model+"Operation valve.PTLKUP")[0] #p_ai
    results[13] = da.read(model+"Operation valve.PTLEAK")[0] #p_ti
    results[14] = da.read(model+"Prod_Choke.GVALVE")[0] #w_p
    results[15] = da.read(model+"Prod_Choke.PVALVE")[0] #p_p ish
    results[16] = da.read(model+"Prod_Choke.VALVDP")[0] #dp_valve
    results[17] = da.read(model+"PRODUCTION Tubing.PT")[0][96] #p_cm
    

    return results
    
    



def OLGA_simulate_pipeline(da):
    n = 500
    results = []
    for i in range(n):
        results.append(OLGA_get_pipeline_measurements(da))

        while(da.read(model+"State")[0]!='STATE_RUNNING_READY'):
            pass
                    
        OLGA_simulate_step(da,10)

           
    return results


    
    

def OLGA_finished(da):
    da.write((model+'Stop',1))
    da.close()
    print('Disconnected from OLGA')



def OLGA_simulate(da):
    n = 4000;
    results = NP.zeros([n,len(OLGA_get_init(da))])
    for i in range(n):
        results[i] = OLGA_get_init(da)
        OLGA_simulate_step(da)
    return results


    