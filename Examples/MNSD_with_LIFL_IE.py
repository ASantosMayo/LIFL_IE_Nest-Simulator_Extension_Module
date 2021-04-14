#%% 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
     ----- MNSD with LIF-SL-I -----
Multi-Neuron Spike Detector with Leaky and Integrate-Fire
with Spike Latency and Intrinsic Plasticity

@author: Alejandro Santos-Mayo // Complutense University of Madrid (UCM)
"""


import nest
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
if not 'lifl_psc_exp_ie' in nest.Models():
    nest.Install('LIFL_IEmodule')

nest.ResetKernel()
nest.SetKernelStatus({'resolution' : 0.1})  # Time resolution set at 0.1 ms

# First of all, we create the dc generator for ordered inputs
E1 = nest.Create("dc_generator"); E2 = nest.Create("dc_generator"); E3 = nest.Create("dc_generator"); E4 = nest.Create("dc_generator")
# Neurons which send the initial supra-threshold input. Connections between Generator-Neuron.
N1 = nest.Create('iaf_psc_alpha'); N2 = nest.Create('iaf_psc_alpha'); N3 = nest.Create('iaf_psc_alpha'); N4 = nest.Create('iaf_psc_alpha')
nest.Connect(E1,N1, {"rule": "one_to_one"}, {"weight": 725.5});
nest.Connect(E2,N2, {"rule": "one_to_one"}, {"weight": 725.5});
nest.Connect(E3,N3, {"rule": "one_to_one"}, {"weight": 725.5});
nest.Connect(E4,N4, {"rule": "one_to_one"}, {"weight": 725.5});


# Main leaky integrate and fire with spike latency and intrinsic plasticity. 
D1 = nest.Create('lifl_psc_exp_ie', 1, { 'lambda': 0.0005 }) # plasticity parameter lambda
D2 = nest.Create('lifl_psc_exp_ie', 1, { 'lambda': 0.0005 }) # Tau parameter is 12.5 as default
D3 = nest.Create('lifl_psc_exp_ie', 1, { 'lambda': 0.0005 })
D4 = nest.Create('lifl_psc_exp_ie', 1, { 'lambda': 0.0005 })
nest.Connect(N1, D1, {"rule": "one_to_one"}, {"weight": 3700.0}) # Connections Source Neuron <-> Detector (LIF_SL_IP)
nest.Connect(N2, D2, {"rule": "one_to_one"}, {"weight": 3700.0})
nest.Connect(N3, D3, {"rule": "one_to_one"}, {"weight": 3700.0})
nest.Connect(N4, D4, {"rule": "one_to_one"}, {"weight": 3700.0})

# We set the GID of neurons which modulate the Intrinsic Plasticity of the neuron.
nest.SetStatus(D1, {'stimulator': [D2[0]]});
nest.SetStatus(D2, {'stimulator': [D3[0], D1[0]]});
nest.SetStatus(D3, {'stimulator': [D2[0], D4[0]]});
nest.SetStatus(D4, {'stimulator': [D3[0]]});

# Target neuron. Connections are set in order to produce a target spike only in pattern detection.
Target = nest.Create('lifl_psc_exp_ie', 1, {'I_e' : 0.0, "tau_minus": 10.0})
nest.Connect(D1,Target, {"rule": "one_to_one"}, {"weight": 410.0});
nest.Connect(D2,Target, {"rule": "one_to_one"}, {"weight": 410.0});
nest.Connect(D3,Target, {"rule": "one_to_one"}, {"weight": 410.0});
nest.Connect(D4,Target, {"rule": "one_to_one"}, {"weight": 410.0});

# Connections between 
nest.Connect(D1,D2, {"rule": "one_to_one"}, { "model": "stdp_synapse", 'delay': 0.1});
nest.Connect(D2,D1, {"rule": "one_to_one"}, { "model": "stdp_synapse", 'delay': 0.1});
nest.Connect(D2,D3, {"rule": "one_to_one"}, { "model": "stdp_synapse", 'delay': 0.1});
nest.Connect(D3,D2, {"rule": "one_to_one"}, { "model": "stdp_synapse", 'delay': 0.1});
nest.Connect(D3,D4, {"rule": "one_to_one"}, { "model": "stdp_synapse", 'delay': 0.1});
nest.Connect(D4,D3, {"rule": "one_to_one"}, { "model": "stdp_synapse", 'delay': 0.1});
        
# We create a Detector so that we can get spike times and raster plot
detector = nest.Create('spike_detector'); 
nest.Connect(D1,detector); nest.Connect(D2, detector); nest.Connect(D3, detector); nest.Connect(D4, detector); nest.Connect(Target, detector);
nest.Connect(N1, detector); nest.Connect(N2, detector); nest.Connect(N3, detector); nest.Connect(N4, detector);
# We create a multim to record the Voltage potential of the membrane (V_m)
# and the excitability parameter of Intrinsic Plasticity (soma_exc)
multim = nest.Create('multimeter', params = {'withtime': True, 'record_from': ['V_m', 'soma_exc'], 'interval': 0.1})
nest.Connect(multim, Target); nest.Connect(multim, D1); nest.Connect(multim, D2); nest.Connect(multim, D3); nest.Connect(multim, D4) #Target

# We create a matrix to save the spike times of each trial
iterations = 300
times = np.zeros((iterations+1, 6)); times[:,4]=range(0,iterations+1)

for i in range(1,iterations): #450
      # First, we set the amplitude to the generator and create a pattern (order) of spikes
    nest.SetStatus(E1, {"amplitude": 0.6575, 'start': (i*1000)-1000+30.0, 'stop': (i*1000)-1000+55.0}) #30 55   0.5635
    nest.SetStatus(E2, {"amplitude": 0.6575, 'start': (i*1000)-1000+33.3, 'stop': (i*1000)-1000+58.3}) #40 65
    nest.SetStatus(E3, {"amplitude": 0.6575, 'start': (i*1000)-1000+36.6, 'stop': (i*1000)-1000+61.6}) #50 75
    nest.SetStatus(E4, {"amplitude": 0.6575, 'start': (i * 1000) - 1000 + 40.0, 'stop': (i * 1000) - 1000 + 65.0})  # 50 75
    nest.SetStatus(D1, {"V_m": -70.0}); nest.SetStatus(D2, {"V_m": -70.0}); nest.SetStatus(D3, {"V_m": -70.0}); nest.SetStatus(D4, {"V_m": -70.0}); nest.SetStatus(Target, {"V_m": -70.0})
    
    nest.Simulate(1000) # Simulate 1000 ms per trial
    
    # We get the status of the detector in order to get the spike times and save them
    spikes = nest.GetStatus(detector)[0]['events']
    t_spikes = spikes['times'][spikes['times']>(1000*i)-1000]; sends_spikes = spikes['senders'][len(spikes['senders'])-len(t_spikes):];
    times1 = t_spikes[sends_spikes == 9]; times2 = t_spikes[sends_spikes == 10]; times3 = t_spikes[sends_spikes == 11]; times4 = t_spikes[sends_spikes == 12]
    if len(times1)==1: times[i,0] = times1-(1000*i)+1000; times[i,4] = times[i,4]+1
    if len(times2)==1: times[i,1] = times2-(1000*i)+1000; times[i,4] = times[i,4]+1
    if len(times3)==1: times[i,2] = times3-(1000*i)+1000; times[i,4] = times[i,4]+1
    if len(times4) == 1: times[i, 3] = times4 - (1000 * i) + 1000; times[i, 4] = times[i, 4] + 1

# Raster plot
nest.raster_plot.from_device(detector)

# We get events of multimeter (V_m and Soma_exc)
events = nest.GetStatus(multim)[0]['events']
t = np.linspace(int(min(events['times'])), int(max(events['times']))+1, (int(max(events['times']))+1)*10-1);
v_1 = events['V_m'][events['senders']==9];v_2 = events['V_m'][events['senders']==10];v_3 = events['V_m'][events['senders']==11];v_4 = events['V_m'][events['senders']==12]; v_5 = events['V_m'][events['senders']==13];
se1 = events['soma_exc'][events['senders']==9]-1;se2 = events['soma_exc'][events['senders']==10]-1;se3 = events['soma_exc'][events['senders']==11]-1; se4 = events['soma_exc'][events['senders']==12]-1;
plt.figure(); plt.plot(t, v_1); plt.plot(t, v_2); plt.plot(t, v_3); plt.plot(t, v_4); plt.plot(t, v_5); plt.ylabel('Membrane potential [mV]')

from mpl_toolkits.mplot3d import Axes3D


delta1 = times[:,0]- times[:,1]
delta2 = times[:,0]- times[:,2]
delta3 = times[:,0]- times[:,3]
delta4 = times[:,1]- times[:,2]
delta5 = times[:,1]- times[:,3]
delta6 = times[:,2]- times[:,3]



fig = plt.figure();
ax = Axes3D(fig);
ax.plot(delta1, delta2, delta3)
ax.scatter(delta1, delta2, delta3)
ax.set_xlabel('Time between D1-D2'); ax.set_ylabel('Time between D2-D3'); ax.set_zlabel('Time between D1-D3')

# 3D plot of Intrinsic Plasticity  excitability parameter for each neuron
import matplotlib.pyplot as plt
fig = plt.figure()
#ax = Axes3D(fig)
ax = fig.add_subplot(111, projection='3d')
ax.plot(se1, se2, se3)
ax.scatter(se1, se2, se3, c = se4, cmap = plt.hot(), s= 0.5)
ax.set_xlabel('Intrinsic Excitability D1'); ax.set_ylabel('Intrinsic Excitability D2'); ax.set_zlabel('Intrinsic Excitability D3')
ax.locator_params(axis='y', nbins=6)
ax.locator_params(axis='x', nbins=6)
ax.locator_params(axis='z', nbins=6)

plt.figure();    #  PLOT IE in time
plt.plot(se1)
plt.plot(se2)
plt.plot(se3)
plt.plot(se4)
plt.xticks([0, 1000000, 2000000, 3000000], [0, 100, 200, 300]); plt.xlabel('Training trials'); plt.ylabel('Intrinsic Excitability (IE)')

plt.figure();    #  PLOT DELAY deltas
plt.plot(delta1)
plt.plot(delta2)
plt.plot(delta3)
plt.plot(delta4)
plt.plot(delta5)
plt.plot(delta6)
plt.xlim([1,300]); plt.xticks([1, 100, 200, 300], [0, 100, 200, 300]); plt.xlabel('Training trials'); plt.ylabel('Spike time difference between neurons')

################################################################################
################   TEST new Patterns   #########################################

# In order to test other possible patterns, we swich off the Intrinsic Plasticity mode (std_mod: synaptic time dependent modificator)
nest.SetStatus(D1, {'std_mod': False})
nest.SetStatus(D2, {'std_mod': False})
nest.SetStatus(D3, {'std_mod': False})
nest.SetStatus(D4, {'std_mod': False})

i = iterations + 2; #453
nest.SetStatus(E1, {"amplitude": 0.6575, 'start': (i*1000)-1000+30.0, 'stop': (i*1000)-1000+55.0}) #30 55
nest.SetStatus(E2, {"amplitude": 0.6575, 'start': (i*1000)-1000+33.3, 'stop': (i*1000)-1000+58.3}) #40 65
nest.SetStatus(E3, {"amplitude": 0.6575, 'start': (i*1000)-1000+40.0, 'stop': (i*1000)-1000+65.0}) #50 75
nest.SetStatus(E4, {"amplitude": 0.6575, 'start': (i*1000)-1000+36.6, 'stop': (i*1000)-1000+61.6}) #50 75
nest.SetStatus(D1, {"V_m": -70.0}); nest.SetStatus(D2, {"V_m": -70.0}); nest.SetStatus(D3, {"V_m": -70.0}); nest.SetStatus(Target, {"V_m": -70.0})
nest.Simulate(3000)

i = i+1; #453
nest.SetStatus(E1, {"amplitude": 0.6575, 'start': (i*1000)-1000+40.0, 'stop': (i*1000)-1000+65.0}) #30 55
nest.SetStatus(E2, {"amplitude": 0.6575, 'start': (i*1000)-1000+36.6, 'stop': (i*1000)-1000+61.6}) #40 65
nest.SetStatus(E3, {"amplitude": 0.6575, 'start': (i*1000)-1000+33.3, 'stop': (i*1000)-1000+58.3}) #50 75
nest.SetStatus(E4, {"amplitude": 0.6575, 'start': (i*1000)-1000+30.0, 'stop': (i*1000)-1000+55.0}) #50 75
nest.SetStatus(D1, {"V_m": -70.0}); nest.SetStatus(D2, {"V_m": -70.0}); nest.SetStatus(D3, {"V_m": -70.0}); nest.SetStatus(Target, {"V_m": -70.0})
nest.Simulate(1000)
i = i+1; #453
nest.SetStatus(E3, {"amplitude": 0.6575, 'start': (i*1000)-1000+35.0, 'stop': (i*1000)-1000+60.0}) #30 55 
nest.SetStatus(E2, {"amplitude": 0.6575, 'start': (i*1000)-1000+30.0, 'stop': (i*1000)-1000+55.0}) #40 65
nest.SetStatus(E1, {"amplitude": 0.6575, 'start': (i*1000)-1000+40.0, 'stop': (i*1000)-1000+65.0}) #50 75
nest.SetStatus(D1, {"V_m": -70.0}); nest.SetStatus(D2, {"V_m": -70.0}); nest.SetStatus(D3, {"V_m": -70.0}); nest.SetStatus(Target, {"V_m": -70.0})
nest.Simulate(1000)
i = i+1; #453
nest.SetStatus(E1, {"amplitude": 0.6575, 'start': (i*1000)-1000+35.0, 'stop': (i*1000)-1000+60.0}) #30 55 
nest.SetStatus(E3, {"amplitude": 0.6575, 'start': (i*1000)-1000+30.0, 'stop': (i*1000)-1000+55.0}) #40 65
nest.SetStatus(E2, {"amplitude": 0.6575, 'start': (i*1000)-1000+40.0, 'stop': (i*1000)-1000+65.0}) #50 75
nest.SetStatus(D1, {"V_m": -70.0}); nest.SetStatus(D2, {"V_m": -70.0}); nest.SetStatus(D3, {"V_m": -70.0}); nest.SetStatus(Target, {"V_m": -70.0})
nest.Simulate(1000)
i = i+1; #453
nest.SetStatus(E2, {"amplitude": 0.6575, 'start': (i*1000)-1000+35.0, 'stop': (i*1000)-1000+60.0}) #30 55 
nest.SetStatus(E3, {"amplitude": 0.6575, 'start': (i*1000)-1000+30.0, 'stop': (i*1000)-1000+55.0}) #40 65
nest.SetStatus(E1, {"amplitude": 0.6575, 'start': (i*1000)-1000+40.0, 'stop': (i*1000)-1000+65.0}) #50 75
nest.SetStatus(D1, {"V_m": -70.0}); nest.SetStatus(D2, {"V_m": -70.0}); nest.SetStatus(D3, {"V_m": -70.0}); nest.SetStatus(Target, {"V_m": -70.0})
nest.Simulate(1000)




# We get events of multimeter (V_m and Soma_exc)
events = nest.GetStatus(multim)[0]['events']
t = np.linspace(int(min(events['times'])), int(max(events['times']))+1, (int(max(events['times']))+1)*10-1);
v_1 = events['V_m'][events['senders']==9];v_2 = events['V_m'][events['senders']==10];v_3 = events['V_m'][events['senders']==11];v_4 = events['V_m'][events['senders']==12]; v_5 = events['V_m'][events['senders']==13];
se1 = events['soma_exc'][events['senders']==9]-1;se2 = events['soma_exc'][events['senders']==10]-1;se3 = events['soma_exc'][events['senders']==11]-1; se4 = events['soma_exc'][events['senders']==12]-1;
plt.figure(); plt.plot(t, v_1); plt.plot(t, v_2); plt.plot(t, v_3); plt.plot(t, v_4); plt.plot(t, v_5); plt.ylabel('Membrane potential [mV]')

plt.figure()
ax = plt.subplot2grid((3, 2), (2, 0), )
ax.plot(t, v_1); ax.plot(t, v_2); ax.plot(t, v_3); ax.plot(t, v_4); ax.plot(t, v_5);
ax.set_xlim([303040,303100])
plt.xticks(range(305000), np.repeat(range(1000), 305))

ax1 = plt.subplot2grid((3, 2), (2, 1), sharey  = ax)
ax1.plot(t, v_1); ax1.plot(t, v_2); ax1.plot(t, v_3); ax1.plot(t, v_4); ax1.plot(t, v_5);
ax1.set_xlim([304040,304100])
plt.setp(ax1.get_yticklabels(), visible=False)
plt.xticks(range(305000), np.repeat(range(1000), 305))

ax2 = plt.subplot2grid((3, 2), (1, 0))
ax2.plot(t, v_1); ax2.plot(t, v_2); ax2.plot(t, v_3); ax2.plot(t, v_4); ax2.plot(t, v_5); plt.ylabel('Membrane potential [mV]')
ax2.set_xlim([302040,302100])
plt.setp(ax2.get_xticklabels(), visible=False)

ax3 = plt.subplot2grid((3, 2), (1, 1), sharey = ax2)
ax3.plot(t, v_1); ax3.plot(t, v_2); ax3.plot(t, v_3); ax3.plot(t, v_4); ax3.plot(t, v_5);
ax3.set_xlim([305040,305100])
plt.setp(ax3.get_xticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

ax4 = plt.subplot2grid((3, 2), (0, 0))
ax4.plot(t, v_1); ax4.plot(t, v_2); ax4.plot(t, v_3); ax4.plot(t, v_4); ax4.plot(t, v_5);
ax4.set_xlim([298040,298100])
plt.setp(ax4.get_xticklabels(), visible=False)

ax5 = plt.subplot2grid((3, 2), (0, 1), sharey = ax4)
ax5.plot(t, v_1); ax5.plot(t, v_2); ax5.plot(t, v_3); ax5.plot(t, v_4);ax5.plot(t, v_5);
ax5.set_xlim([301040,301100])
plt.setp(ax5.get_xticklabels(), visible=False)
plt.setp(ax5.get_yticklabels(), visible=False)

plt.show()