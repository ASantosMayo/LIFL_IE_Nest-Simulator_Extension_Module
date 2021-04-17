#
#  A bio-inspired model of the early visual system, based on parallel spike-sequence detection, showing orientation selectivity

#  Simulation script; Results plots; Comparison simulation - MEG signal

import matplotlib.pyplot as plt
import matplotlib
import pickle
import gzip
import numpy
import nest
import nest.raster_plot
import numpy as np
import scipy.io
import pickle


if not 'lifl_psc_exp_ie' in nest.Models(): # Load our customized NEST Simulator module for MNSD.
    nest.Install('LIFL_IEmodule')

from OrientedColumnV1 import column # Load the cortical column properties and building


totalsimulations = 1  # In this file we will only simulate 1 trial. We could create an bucle FOR and repeat Visual stimulations


FiringRate_0 = np.zeros([totalsimulations, 37]); FiringRate_45 = np.zeros([totalsimulations, 37]);
FiringRate_90 = np.zeros([totalsimulations, 37]); FiringRate_135 = np.zeros([totalsimulations, 37])
SynLFP = np.zeros([totalsimulations, 4999])
OSI = np.zeros([totalsimulations, 3])

simulations = 0 # Current simulation iteration. In this file only will be 1.

nest.ResetKernel()
nest.SetKernelStatus({'resolution': 0.1})  # Time resolution set at 0.1 ms
plt.close('all')

GCells = 324  # Number of Ganglionar Cells (Same number for #Retina = #LGN = #SS4 )
spikingGanglion = range(0, 324, 1)
inputs = nest.Create("spike_generator", GCells)
InputsDetector = nest.Create('spike_detector')
nest.Connect(inputs, InputsDetector)

LGN = nest.Create('parrot_neuron', GCells) # Parrot neurons will fire at same time than Ganglionar Cells (retina)
nest.Connect(inputs, LGN, 'one_to_one')

# We create the V1 columns defined on the function
Detector0, Spikes0, Multimeter0, SomaMultimeter0, Pyr230, SS40, Pyr50, Pyr60, In230, In40, In50, In60 = column(0, LGN)
Detector45, Spikes45, Multimeter45, SomaMultimeter45, Pyr2345, SS445, Pyr545, Pyr645, In2345, In445, In545, In645 = column(45, LGN)
Detector90, Spikes90, Multimeter90, SomaMultimeter90, Pyr2390, SS490, Pyr590, Pyr690, In2390, In490, In590, In690 = column(90, LGN)
Detector135, Spikes135, Multimeter135, SomaMultimeter135, Pyr23135, SS4135, Pyr5135, Pyr6135, In23135, In4135, In5135, In6135 = column(135, LGN)


print(simulations) # In this file wil be only 1
nest.Simulate(400 + np.round(np.random.rand(1)*200)) # We simulate about half a second to initialize the network and assure randomness on results


# file 19 corresponds with angle 90º ;   We give some randomness to spike times. (Avoiding spikes to be exactly in same time step)
exec('file = "./files/spikes_reponse_gabor_randn02_19.pckl"')
with open(file, 'rb') as f:
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    resultados = u.load()
    senders = resultados['senders'];
    times = resultados['times']

currtime = nest.GetKernelStatus('time')
for j in range(0, len(spikingGanglion)):
    spike_time = times[senders == min(senders)+j] + currtime
    nest.SetStatus([inputs[j]], {'spike_times': spike_time})

nest.Simulate(200)


###  Orientation Selectivity Index (OSI)
# Now we Calculate the OSI index by extracting the firing rate from layer 2/3 in both 0º an 90ª
senders = nest.GetStatus(Spikes0)[0]['events']['senders'][nest.GetStatus(Spikes0)[0]['events']['times'] > currtime]
times = nest.GetStatus(Spikes0)[0]['events']['times'][nest.GetStatus(Spikes0)[0]['events']['times'] > currtime]
rate230 = len(times[np.isin(senders, Pyr230)])
senders = nest.GetStatus(Spikes90)[0]['events']['senders'][nest.GetStatus(Spikes90)[0]['events']['times'] > currtime]
times = nest.GetStatus(Spikes90)[0]['events']['times'][nest.GetStatus(Spikes90)[0]['events']['times'] > currtime]
rate2390 = len(times[np.isin(senders, Pyr2390)])
OSI[simulations, 0] = (rate2390 - rate230) / (rate2390 + rate230)

#  ####### PREFERRED COLUMN:  V_m membrane potential
idneuron = 15; # Any id can be selected to get membrane potential data for plot.

events = nest.GetStatus(Multimeter90)[0]['events']; eventos = events['V_m']
events_V_mPref = events
t = range(0, (2000-1),)
j = 0; v_1 = np.zeros([96,2000-1]); v_2 = np.zeros([96,2000-1]); v_3 = np.zeros([96,2000-1]); v_4 = np.zeros([96,2000-1]); v_5 = np.zeros([96,2000-1]); v_6 = np.zeros([96,2000-1]); v_7 = np.zeros([96,2000-1])
k = 0
for j in range(0,324,36):
    for i in range(0,18,2):
        #Extract Membrane Potential (V) of neurons to plot
        tmp = events['V_m'][events['senders'] == SS490[i+j]];
        v_1[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == SS490[i+j+1]];
        v_2[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == SS490[i+j+18]];
        v_3[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == SS490[i + j + 19]];
        v_4[k,:] = tmp[len(tmp)-1999:]
        tmp= events['V_m'][events['senders'] == Pyr2390[k]];
        v_5[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == Pyr590[k]];
        v_6[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == Pyr690[k]];
        v_7[k,:] = tmp[len(tmp)-1999:]
        print('V_m calculating: ', k)
        k +=1

plt.figure(10); ax = plt.subplot(4,1,1); plt.plot(t, v_5[idneuron,:], color = [0.4, 0.2, 0.3]); plt.xticks([])
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none'); plt.gca().yaxis.tick_right();  plt.ylabel('Pyr II/III'); plt.ylim([-70, 10])
plt.subplot(4,1,2);plt.plot(t, v_1[idneuron,:]); plt.plot(t, v_2[idneuron,:]); plt.plot(t, v_3[idneuron,:]); plt.plot(t, v_4[idneuron,:]); plt.xticks([])
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel('SS IV'); plt.ylim([-70, 10])
plt.subplot(4,1,3); plt.plot(t, v_6[idneuron,:], color = [0.2, 0.2, 0.2]); plt.ylabel('L'); plt.xticks([])
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel('Pyr V'); plt.ylim([-70, 10])
plt.subplot(4, 1, 4); plt.plot(t, v_7[idneuron, :], color = [0.2, 0.2, 0.2])
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel('Pyr VI'); plt.ylim([-70, 10])
#plt.savefig('figures/V_m_prefered.pdf'); plt.savefig('figures/V_m_prefered.png')


#  ####### NON PREFERRED COLUMN:  V_m membrane potential
####### :  V_m membrane potential
events = nest.GetStatus(Multimeter0)[0]['events']; eventos = events['V_m']
events_V_mNonPref = events
t = range(0, (2000-1),)
j = 0; v_1 = np.zeros([96,2000-1]); v_2 = np.zeros([96,2000-1]); v_3 = np.zeros([96,2000-1]); v_4 = np.zeros([96,2000-1]); v_5 = np.zeros([96,2000-1]); v_6 = np.zeros([96,2000-1]); v_7 = np.zeros([96,2000-1])
k = 0
for j in range(0,324,36):
    for i in range(0,18,2):
        tmp = events['V_m'][events['senders'] == SS40[i+j]]
        v_1[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == SS40[i+j+1]]
        v_2[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == SS40[i+j+18]]
        v_3[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == SS40[i + j + 19]]
        v_4[k,:] = tmp[len(tmp)-1999:]
        tmp= events['V_m'][events['senders'] == Pyr230[k]]
        v_5[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == Pyr50[k]]
        v_6[k,:] = tmp[len(tmp)-1999:]
        tmp = events['V_m'][events['senders'] == Pyr60[k]]
        v_7[k,:] = tmp[len(tmp)-1999:]
        print('V_m calculating: ', k)
        k +=1

plt.figure(); ax = plt.subplot(4,1,1); plt.plot(t, v_5[idneuron,:], color = [0.4, 0.2, 0.3]); plt.xticks([])
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none'); plt.gca().yaxis.tick_right();  plt.ylabel('Pyr II/III'); plt.ylim([-70, 10])
plt.subplot(4,1,2);plt.plot(t, v_1[idneuron,:]); plt.plot(t, v_2[idneuron,:]); plt.plot(t, v_3[idneuron,:]); plt.plot(t, v_4[idneuron,:]); plt.xticks([])
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel('SS IV'); plt.ylim([-70, 10])
plt.subplot(4,1,3); plt.plot(t, v_6[idneuron,:], color = [0.2, 0.2, 0.2]); plt.ylabel('L'); plt.xticks([])
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel('Pyr V'); plt.ylim([-70, 10])
plt.subplot(4, 1, 4); plt.plot(t, v_7[idneuron, :], color = [0.2, 0.2, 0.2])
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel('Pyr VI'); plt.ylim([-70, 10])
#plt.savefig('figures/V_m_nonprefered.pdf'); plt.savefig('figures/V_m_nonprefered.png')

del v_1, v_2, v_3, v_4, v_5, v_6, v_7

##  Here SYNAPTIC INPUT CURRENT is extracted and computed in relation with Dendrite lenght. Thus we obtain signal y pA (picoAmpere). See Method in the article.
t = np.linspace(0,498, 4999)
LFP_pinwheel = np.zeros([1,4999]); plt.figure(); s = 1
for d in [0, 45, 90, 135]:
    exec('eventsSyn = nest.GetStatus(Multimeter' + str(d) + ')[0]["events"];')
    exec('Pyr23 = Pyr23' + str(d) + '; Pyr5 = Pyr5' + str(d) + '; Pyr6 = Pyr6' + str(d) + ';')
    syn_23 = np.zeros([96,5000-1]); syn_5 = np.zeros([96,5000-1]); syn_6 = np.zeros([96,5000-1]); k = 0-1
    for j in range(0, 324, 36):
        for i in range(0, 18, 2):
            syn_23tmp = eventsSyn['I_syn_ex'][eventsSyn['senders'] == Pyr23[k]]
            syn_23[k,:] = syn_23tmp[len(syn_23tmp)-4999:]* 0.000482
            syn_5tmp = eventsSyn['I_syn_ex'][eventsSyn['senders'] == Pyr5[k]]
            syn_5[k,:] = syn_5tmp[len(syn_5tmp)-4999:]  * 0.001342
            syn_6tmp = eventsSyn['I_syn_ex'][eventsSyn['senders'] == Pyr6[k]]
            syn_6[k,:] = syn_6tmp[len(syn_6tmp)-4999:]  * 0.000963
            k += 1
    exec('syn_exc_' + str(d) + ' = (np.mean(syn_23, axis=0) + np.mean(syn_5, axis=0) + np.mean(syn_6, axis=0))')
    plt.subplot(4,1,s)
    eval('plt.plot(t, syn_exc_' + str(d) + ', "k")')
    plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel(str(d) + "º")
    exec('LFP_pinwheel = LFP_pinwheel + syn_exc_' + str(d))
    del syn_23, syn_23tmp, syn_5, syn_5tmp, syn_6, syn_6tmp
    s += 1
SynLFP[simulations, :] = syn_exc_0 + syn_exc_45 + syn_exc_90 + syn_exc_135


plt.figure(); plt.plot(syn_exc_0, "k")
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel(str(d) + "º")
plt.ylim()[1]
plt.figure(); plt.plot(syn_exc_90, "k")
plt.gca().spines['left'].set_color('none'); plt.xticks([]); plt.gca().spines['top'].set_color('none'); plt.xticks([]); plt.gca().spines['bottom'].set_color('none');plt.gca().yaxis.tick_right(); plt.ylabel(str(d) + "º")
plt.ylim()[1]
print(currtime)

plt.figure(); plt.plot(SynLFP[simulations, :]); plt.legend(['LFP of Pinwheel'])

import scipy.io
LFP_V1 = SynLFP
scipy.io.savemat('LFP_V1.mat', dict(x = t, y = LFP_V1))




###  RASTER PLOT of 0degree COLUMN

events = nest.GetStatus(Detector0)[0]['events']
times = events['times']*10; times = times.astype(int)
senders = events['senders']
raster23 = []; raster4 = []; raster5 = []; raster6 = []; raster23i = []; raster4i = []; raster5i = []; raster6i = []; rasterdata = []

for i in In60: rasterdata.append(times[senders == i])
raster6i = range(0, len(rasterdata));tmp = len(rasterdata)
for i in Pyr60: rasterdata.append(times[senders == i])
raster6 = range(tmp, len(rasterdata));tmp = len(rasterdata)

for i in In50: rasterdata.append(times[senders == i])
raster5i = range(tmp, len(rasterdata)); tmp = len(rasterdata)
for i in Pyr50: rasterdata.append(times[senders == i])
raster5 = range(tmp, len(rasterdata)); tmp = len(rasterdata)

for i in In40: rasterdata.append(times[senders == i])
raster4i = range(tmp, len(rasterdata)); tmp = len(rasterdata)
for i in SS40: rasterdata.append(times[senders == i])
raster4 = range(tmp, len(rasterdata)); tmp = len(rasterdata)

for i in In230: rasterdata.append(times[senders == i])
raster23i = range(tmp, len(rasterdata)); tmp = len(rasterdata)
for i in Pyr230: rasterdata.append(times[senders == i])
raster23 = range(tmp, len(rasterdata)); tmp = len(rasterdata)

rasterdata.append([0])

keycolors = np.array([[0.0, 0.0, 0.0, 0.0],
[0.0,0.2, 0.2, 0.2],
[0.0,0.5, 0.5, 0.5],
[0.0,0.2, 0.2, 0.2],
[0.,0.5, 0.5, 0.5],
[0.0,0.2, 0.2, 0.2],
[0.0,0.5, 0.5, 0.5],
[0.0,0.2, 0.2, 0.2],
[0.0,0.5, 0.5, 0.5]])


color_layers = np.zeros([tmp+1, 3]); color_layers[tmp, :] = [1,1,1]
color_layers[min(raster23):max(raster23), 0] = keycolors[1,1];color_layers[min(raster23):max(raster23), 1] = keycolors[1,2];color_layers[min(raster23):max(raster23), 2] = keycolors[1,3];
color_layers[min(raster4):max(raster4), 0] = keycolors[3,1];color_layers[min(raster4):max(raster4), 1] = keycolors[3,2];color_layers[min(raster4):max(raster4), 2] = keycolors[3,3];
color_layers[min(raster5):max(raster5), 0] = keycolors[5,1];color_layers[min(raster5):max(raster5), 1] = keycolors[5,2];color_layers[min(raster5):max(raster5), 2] = keycolors[5,3];
color_layers[min(raster6):max(raster6), 0] = keycolors[7,1];color_layers[min(raster6):max(raster6), 1] = keycolors[7,2];color_layers[min(raster6):max(raster6), 2] = keycolors[7,3];

color_layers[min(raster23i):max(raster23i), 0] = keycolors[2,1];color_layers[min(raster23i):max(raster23i), 1] = keycolors[2,2];color_layers[min(raster23i):max(raster23i), 2] = keycolors[2,3];
color_layers[min(raster4i):max(raster4i), 0] = keycolors[4,1];color_layers[min(raster4i):max(raster4i), 1] = keycolors[4,2];color_layers[min(raster4i):max(raster4i), 2] = keycolors[4,3];
color_layers[min(raster5i):max(raster5i), 0] = keycolors[6,1];color_layers[min(raster5i):max(raster5i), 1] = keycolors[6,2];color_layers[min(raster5i):max(raster5i), 2] = keycolors[6,3];
color_layers[min(raster6i):max(raster6i), 0] = keycolors[8,1];color_layers[min(raster6i):max(raster6i), 1] = keycolors[8,2];color_layers[min(raster6i):max(raster6i), 2] = keycolors[8,3];


plt.figure()
plt.eventplot(rasterdata, color = color_layers, linewidths = 3, linelengths=7)
plt.xlabel('Time (ms)'); plt.ylabel('Neuron')
plt.xticks([currtime*10,  currtime*10+250, currtime*10+500,  currtime*10+750, currtime*10+1000,   currtime*10+1250, currtime*10+1500,  currtime*10+1750,  currtime*10+2000 ], [0, 25, 50, 75, 100, 125, 150, 175, 200])
plt.yticks([1055, 855, 617, 455, 370, 313, 181, 30], ["Pyr23", "Inh23", "SS4", "Inh4", "Pyr5", "Inh5", "Pyr6", "Inh6"])
plt.gca().spines['right'].set_color('none'); plt.gca().spines['top'].set_color('none');
#plt.savefig('figures/Raster_nonpreferedCentered.pdf'); plt.savefig('figures/Raster_nonpreferedCentered.png')
rasterdata1 = rasterdata

###  RASTER PLOT of 90degree COLUMN

events = nest.GetStatus(Detector90)[0]['events']
times = events['times']*10; times = times.astype(int)
senders = events['senders']
raster23 = []; raster4 = []; raster5 = []; raster6 = []; raster23i = []; raster4i = []; raster5i = []; raster6i = []; rasterdata = []


for i in In690: rasterdata.append(times[senders == i])
raster6i = range(0, len(rasterdata));tmp = len(rasterdata)
for i in Pyr690: rasterdata.append(times[senders == i])
raster6 = range(tmp, len(rasterdata));tmp = len(rasterdata)

for i in In590: rasterdata.append(times[senders == i])
raster5i = range(tmp, len(rasterdata)); tmp = len(rasterdata)
for i in Pyr590: rasterdata.append(times[senders == i])
raster5 = range(tmp, len(rasterdata)); tmp = len(rasterdata)

for i in In490: rasterdata.append(times[senders == i])
raster4i = range(tmp, len(rasterdata)); tmp = len(rasterdata)
for i in SS490: rasterdata.append(times[senders == i])
raster4 = range(tmp, len(rasterdata)); tmp = len(rasterdata)

for i in In2390: rasterdata.append(times[senders == i])
raster23i = range(tmp, len(rasterdata)); tmp = len(rasterdata)
for i in Pyr2390: rasterdata.append(times[senders == i])
raster23 = range(tmp, len(rasterdata)); tmp = len(rasterdata)

rasterdata.append([0])

keycolors = np.array([[0.0, 0.0, 0.0, 0.0],
[0.0,0.2, 0.2, 0.2],
[0.0,0.5, 0.5, 0.5],
[0.0,0.2, 0.2, 0.2],
[0.,0.5, 0.5, 0.5],
[0.0,0.2, 0.2, 0.2],
[0.0,0.5, 0.5, 0.5],
[0.0,0.2, 0.2, 0.2],
[0.0,0.5, 0.5, 0.5]])


color_layers = np.zeros([tmp+1, 3]); color_layers[tmp, :] = [1,1,1]
color_layers[min(raster23):max(raster23), 0] = keycolors[1,1];color_layers[min(raster23):max(raster23), 1] = keycolors[1,2];color_layers[min(raster23):max(raster23), 2] = keycolors[1,3];
color_layers[min(raster4):max(raster4), 0] = keycolors[3,1];color_layers[min(raster4):max(raster4), 1] = keycolors[3,2];color_layers[min(raster4):max(raster4), 2] = keycolors[3,3];
color_layers[min(raster5):max(raster5), 0] = keycolors[5,1];color_layers[min(raster5):max(raster5), 1] = keycolors[5,2];color_layers[min(raster5):max(raster5), 2] = keycolors[5,3];
color_layers[min(raster6):max(raster6), 0] = keycolors[7,1];color_layers[min(raster6):max(raster6), 1] = keycolors[7,2];color_layers[min(raster6):max(raster6), 2] = keycolors[7,3];

color_layers[min(raster23i):max(raster23i), 0] = keycolors[2,1];color_layers[min(raster23i):max(raster23i), 1] = keycolors[2,2];color_layers[min(raster23i):max(raster23i), 2] = keycolors[2,3];
color_layers[min(raster4i):max(raster4i), 0] = keycolors[4,1];color_layers[min(raster4i):max(raster4i), 1] = keycolors[4,2];color_layers[min(raster4i):max(raster4i), 2] = keycolors[4,3];
color_layers[min(raster5i):max(raster5i), 0] = keycolors[6,1];color_layers[min(raster5i):max(raster5i), 1] = keycolors[6,2];color_layers[min(raster5i):max(raster5i), 2] = keycolors[6,3];
color_layers[min(raster6i):max(raster6i), 0] = keycolors[8,1];color_layers[min(raster6i):max(raster6i), 1] = keycolors[8,2];color_layers[min(raster6i):max(raster6i), 2] = keycolors[8,3];


plt.figure()
plt.eventplot(rasterdata, color = color_layers, linewidths = 3, linelengths=7)
plt.xlabel('Time (ms)'); plt.ylabel('Neuron')
plt.xticks([currtime*10,  currtime*10+250, currtime*10+500,  currtime*10+750, currtime*10+1000,   currtime*10+1250, currtime*10+1500,  currtime*10+1750,  currtime*10+2000 ], [0, 25, 50, 75, 100, 125, 150, 175, 200])
plt.yticks([1055, 855, 617, 455, 370, 313, 181, 30], ["Pyr23", "Inh23", "SS4", "Inh4", "Pyr5", "Inh5", "Pyr6", "Inh6"])
plt.gca().spines['right'].set_color('none'); plt.gca().spines['top'].set_color('none');

plt.savefig('figures/Raster_preferedCentered.pdf'); plt.savefig('figures/Raster_preferedCentered.png')

rasterdata2 = rasterdata



#######  Load MEG data. Filter and prepare both signals for
#######         COMPARISON SIMULATION vs MEG SIGNAL

from scipy.io import loadmat
MEG = loadmat("./files/MEGGaborStimScoutTSeries.mat")
MEGdata = MEG['Value']; MEGtime = MEG['Time']

plt.figure(); plt.plot(MEGtime.T, MEGdata.T)

#  CORRELATION

from scipy import signal
datameg = MEGdata[:,420:720].T;
datasim = SynLFP[simulations, :];
datasim = datasim - np.mean(datasim)
datasim = scipy.signal.resample(datasim, 300);
datasim = datasim - np.mean(datasim)
datameg = np.squeeze(datameg - np.mean(datameg))
datatime = MEGtime[:,420:720].T;

from FUNCIONES import filtro

datasim = filtro(np.flip(datasim), 0.5, 60, 600); datasim = filtro(np.flip(datasim), 0.1, 60, 600);

datameg = filtro(np.flip(datameg), 0.5, 60, 600); datameg = filtro(np.flip(datameg), 0.1, 60, 600)


# Now we have preprocessed timeseries

datameg = datameg[180:252]; datasim = datasim[180:252]; datatime = datatime[180:252]

datacorr = numpy.array([datameg, datasim])
CORR = numpy.corrcoef(datacorr)
import scipy.stats
scipy.stats.pearsonr(datameg, datasim)

datameg01 = (datameg - min(datameg))/(max(datameg) - min(datameg))
datasim01 = (datasim - min(datasim))/(max(datasim) - min(datasim))
colors = ['#08F7FE',  # teal/cyan
    '#FE53BB',  # pink
          ]
plt.figure(); plt.plot(datatime, datameg01, 'b', color= 'olive', linewidth=2); plt.plot(datatime, datasim01, color='orchid', linewidth=2)
plt.legend(['Real data', 'Simulated data']);
plt.gca().spines['right'].set_color('none'); plt.gca().spines['top'].set_color('none');
plt.xlabel('Time (ms)'); plt.ylabel('Current (pA·m)')
plt.yticks([0, 1], [0, 1]); plt.xticks([0, 0.025, 0.05, 0.075, 0.1, 0.125], [0, 25, 50, 75, 100, 125])

