def column(setdegree, LGN):
    """"
    Visual Cortex (V1) Oriented Column -

    This function creates a cortical column trained with a specific orientation (pre-trained 'soma_exc' IE values)

    The column is composed by:
        - Pyramidal (Pyr)               layer 2/3
        - Inhibitory interneurons(inh)  layer 2/3
        - Spiny Stellate Cells(SS)      layer 4
        - Inhibitory interneurons       layer 4
        - Pyramidal                     layer 5
        - Inhibitory interneurons       layer 5
        - Pyramidal                     layer 6
        - Inhibitory interneurons       layer 6

    """
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
    if not 'lifl_psc_exp_ie' in nest.Models():
        nest.Install('LIFL_IEmodule')


    SS4 = nest.Create('lifl_psc_exp_ie', 324, {'I_e': 0.0,  # 122.1
                           'V_m': -70.0,
                           'E_L': -65.0,
                           'V_th': -50.0,
                           'V_reset': -65.0,
                           'C_m': 250.0,
                           'tau_m': 10.0,
                           'tau_syn_ex': 2.0,
                           'tau_syn_in': 2.0,
                           't_ref': 2.0,
                           'std_mod': False,
                           'lambda': 0.0005,
                           'tau': 12.5, })


    lgn2v1_delay = 1.0 # Delay from LGN to Cortex

    nest.Connect(LGN, SS4, {'rule': 'one_to_one'}, {
        "weight": 15000.0,
        "delay": lgn2v1_delay})

    # Target neuron. Connections are set in order to produce a target spike only in pattern detection.
    Pyr23 = nest.Create('aeif_psc_exp_peak', 324, {
                            'I_e': 0.0,
                            'V_m': -70.0,
                            'E_L': -70.0,
                            'V_th': -50.0,
                            'V_reset': -55.0,
                            'C_m': 250.0,
                            'tau_syn_ex': 2.0,
                            'tau_syn_in': 2.0,
                            't_ref': 2.0,
                            'g_L': 980.0
    })
    for idn in Pyr23:
        nest.SetStatus([idn], {'V_m': (-65.0 + np.random.rand()*10.0)})

    Pyr5 = nest.Create('aeif_psc_exp_peak', 81,
                       {   'I_e': 0.0,
                           'V_m': -70.0,
                           'E_L': -70.0,
                           'V_th': -50.0,
                           'V_reset': -55.0,
                           'C_m': 250.0,
                           'tau_syn_ex': 2.0,
                           'tau_syn_in': 2.0,
                           't_ref': 2.0,
                           'g_L': 980.0
                       })
    for idn in Pyr5:
        nest.SetStatus([idn], {'V_m': (-65.0 + np.random.rand()*10.0)})
    Pyr6 = nest.Create('aeif_psc_exp_peak', 243,
                       {  'I_e': 0.0,
                           'V_m': -70.0,
                           'E_L': -70.0,
                           'V_th': -50.0,
                           'V_reset': -55.0,
                           'C_m': 250.0,
                           'tau_syn_ex': 2.0,
                           'tau_syn_in': 2.0,
                           't_ref': 2.0,
                           'g_L': 980.0
                       })
    for idn in Pyr6:
        nest.SetStatus([idn], {'V_m': (-65.0 + np.random.rand()*10.0)})

    # Poisson Noise Generators
    poisson_activityL23 = nest.Create('poisson_generator', 1)
    nest.SetStatus(poisson_activityL23, {'rate': 1721500.0})
    poisson_activityL5 = nest.Create('poisson_generator', 1)
    nest.SetStatus(poisson_activityL5, {'rate': 1740000.0})
    poisson_activityL6 = nest.Create('poisson_generator', 1)
    nest.SetStatus(poisson_activityL6, {'rate': 1700000.0})
    poisson_activityInh = nest.Create('poisson_generator', 1)
    nest.SetStatus(poisson_activityInh, {'rate': 1750000.0})
    nest.Connect(poisson_activityL23, Pyr23, {'rule': 'all_to_all'}, {"weight": 5.0})
    nest.Connect(poisson_activityL5, Pyr5, {'rule': 'all_to_all'}, {"weight": 5.0})
    nest.Connect(poisson_activityL6, Pyr6, {'rule': 'all_to_all'}, {"weight": 5.0})


    # FeedForward
    nest.Connect(Pyr23, Pyr5, {'rule': 'fixed_indegree', 'indegree': 15}, {"weight": 100.0, "delay": 1.0})
    nest.Connect(Pyr5, Pyr6, {'rule': 'fixed_indegree', 'indegree': 20}, {"weight": 100.0, "delay": 1.0})

    ## Connections between layers
    nest.Connect(Pyr23, Pyr23, {'rule': 'fixed_indegree', 'indegree': 36}, {"weight": 100.0, "delay": 1.0})
    nest.Connect(Pyr5, Pyr5, {'rule': 'fixed_indegree', 'indegree': 10}, {"weight": 100.0, "delay": 1.0})
    nest.Connect(Pyr6, Pyr6, {'rule': 'fixed_indegree', 'indegree': 20}, {"weight": 100.0, "delay": 1.0})

    In4 = nest.Create('aeif_psc_exp_peak', 65,
                      {
                          'I_e': 0.0,
                          'V_m': -70.0,
                          'E_L': -70.0,
                          'V_th': -50.0,
                          'V_reset': -55.0,
                          'C_m': 250.0,
                          'tau_syn_ex': 2.0,
                          'tau_syn_in': 2.0,
                          't_ref': 1.0,
                          'g_L': 980.0
                      })
    nest.Connect(poisson_activityInh, In4, {'rule': 'all_to_all'}, {"weight": 4.9})
    nest.Connect(SS4, In4, {'rule': 'fixed_indegree', 'indegree': 32}, {"weight": 100.0, "delay": 1.0})
    nest.Connect(In4, SS4, {'rule': 'fixed_indegree', 'indegree': 6}, {"weight": -100.0, "delay": 1.0})
    nest.Connect(In4, In4, {'rule': 'fixed_indegree', 'indegree': 6}, {"weight": -100.0, "delay": 1.0})

    poisson_activity_inh = nest.Create('poisson_generator', 1)
    nest.SetStatus(poisson_activity_inh, {'rate': 340000.0})
    In23 = nest.Create('aeif_psc_exp_peak', 65,
                       {
                           'I_e': 0.0,
                           'V_m': -70.0,
                           'E_L': -70.0,
                           'V_th': -50.0,
                           'V_reset': -55.0,
                           'C_m': 250.0,
                           'tau_syn_ex': 2.0,
                           'tau_syn_in': 2.0,
                           't_ref': 1.0,
                           'g_L': 980.0
                       })

    nest.Connect(poisson_activityInh, In23, {'rule': 'all_to_all'}, {"weight": 5.0});
    nest.Connect(Pyr23, In23, {'rule': 'fixed_indegree', 'indegree': 35}, {"weight": 100.0, "delay": 1.0})
    nest.Connect(In23, Pyr23, {'rule': 'fixed_indegree', 'indegree': 8}, {"weight": -100.0, "delay": 1.0})
    nest.Connect(In23, In23, {'rule': 'fixed_indegree', 'indegree': 8}, {"weight": -100.0, "delay": 1.0})

    In5 = nest.Create('aeif_psc_exp_peak', 16,
                      {
                          'I_e': 0.0,
                          'V_m': -70.0,
                          'E_L': -70.0,
                          'V_th': -50.0,
                          'V_reset': -55.0,
                          'C_m': 250.0,
                          'tau_syn_ex': 2.0,
                          'tau_syn_in': 2.0,
                          't_ref': 1.0,
                          'g_L': 980.0
                      })

    nest.Connect(poisson_activityInh, In5, {'rule': 'all_to_all'}, {"weight": 5.0});
    nest.Connect(Pyr5, In5, {'rule': 'fixed_indegree', 'indegree': 30}, {"weight": 100.0, "delay": 1.0})
    nest.Connect(In5, Pyr5, {'rule': 'fixed_indegree', 'indegree': 8}, {"weight": -100.0, "delay": 1.0})
    nest.Connect(In5, In5, {'rule': 'fixed_indegree', 'indegree': 8}, {"weight": -100.0, "delay": 1.0})

    In6 = nest.Create('aeif_psc_exp_peak', 49,
                      {
                          'I_e': 0.0,
                          'V_m': -70.0,
                          'E_L': -70.0,
                          'V_th': -50.0,
                          'V_reset': -55.0,
                          'C_m': 250.0,
                          'tau_syn_ex': 2.0,
                          'tau_syn_in': 2.0,
                          't_ref': 1.0,
                          'g_L': 980.0
                      })
    nest.Connect(poisson_activityInh, In6, {'rule': 'all_to_all'}, {"weight": 5.0});
    nest.Connect(Pyr6, In6, {'rule': 'fixed_indegree', 'indegree': 32}, {"weight": 100.0, "delay": 1.0})
    nest.Connect(In6, Pyr6, {'rule': 'fixed_indegree', 'indegree': 6}, {"weight": -100.0, "delay": 1.0})
    nest.Connect(In6, In6, {'rule': 'fixed_indegree', 'indegree': 6}, {"weight": -100.0, "delay": 1.0})

    # Here we load the Soma_exc (IE value) trained before for each preferred angle.
    exec(
        'file = "./files/soma_exc_15_' + str(
            setdegree) + '.pckl"', None, globals())
    with open(file, 'rb') as f:
        u = pickle._Unpickler(f)
        u.encoding = 'latin1'
        SS4_soma_exc_raw = u.load()
        SS4_soma_exc = SS4_soma_exc_raw[:]
        del SS4_soma_exc_raw, setdegree

    for i in range(0, 324):
        nest.SetStatus([SS4[i]], {'soma_exc': SS4_soma_exc[i]})

    k = 0
    for j in range(0, 324, 36):
        for i in range(0, 18, 2):

            # Set stimulator (i.e. Neuromodulator) of each SS cell
            nest.SetStatus([SS4[i + j]],
                           {'stimulator': [SS4[i + j + 1], SS4[i + j + 18]]})
            nest.SetStatus([SS4[i + j + 1]],
                           {'stimulator': [SS4[i + j], SS4[i + j + 19]]})
            nest.SetStatus([SS4[i + j + 18]],
                           {'stimulator': [SS4[i + j], SS4[i + j + 19]]})
            nest.SetStatus([SS4[i + j + 19]],
                           {'stimulator': [SS4[i + j + 18], SS4[i + j + 1]]})

            # Connect betwen neuromodulators of SS cell (groups of 4 SS)
            nest.Connect([SS4[i + j]], [SS4[i + j + 1]], {"rule": "one_to_one"},
                         {"model": "stdp_synapse", 'delay': 0.1})
            nest.Connect([SS4[i + j]], [SS4[i + j + 18]], {"rule": "one_to_one"},
                         {"model": "stdp_synapse", 'delay': 0.1})
            nest.Connect([SS4[i + j + 1]], [SS4[i + j]], {"rule": "one_to_one"},
                         {"model": "stdp_synapse", 'delay': 0.1})
            nest.Connect([SS4[i + j + 1]], [SS4[i + j + 19]], {"rule": "one_to_one"},
                         {"model": "stdp_synapse", 'delay': 0.1})
            nest.Connect([SS4[i + j + 18]], [SS4[i + j]], {"rule": "one_to_one"},
                         {"model": "stdp_synapse", 'delay': 0.1})
            nest.Connect([SS4[i + j + 18]], [SS4[i + j + 19]], {"rule": "one_to_one"},
                         {"model": "stdp_synapse", 'delay': 0.1})
            nest.Connect([SS4[i + j + 19]], [SS4[i + j + 1]], {"rule": "one_to_one"},
                         {"model": "stdp_synapse", 'delay': 0.1})
            nest.Connect([SS4[i + j + 19]], [SS4[i + j + 18]], {"rule": "one_to_one"},
                         {"model": "stdp_synapse", 'delay': 0.1})

            # Connect each group of 4 SS to 4 layer 2/3 Pyramidal Cells so that each pyramidal only fire on polychrony arrival of SS fires
            nest.Connect([SS4[i + j], SS4[i + j + 1], SS4[i + j + 18], SS4[i + j + 19]],
                         [Pyr23[i+j], Pyr23[i+j+1], Pyr23[i+j+2], Pyr23[i+j+3]], {"rule": "all_to_all"},
                         {"weight": 100.0, "delay": 1.0})


            k += 1

    # We create a Detector so that we can get spike times and raster plot
    Detector = nest.Create('spike_detector')
    nest.Connect(Pyr5, Detector)
    nest.Connect(Pyr6, Detector)
    nest.Connect(SS4, Detector)
    nest.Connect(Pyr23, Detector)
    nest.Connect(In23, Detector)
    nest.Connect(In4, Detector)
    nest.Connect(In5, Detector)
    nest.Connect(In6, Detector)

    Spikes = nest.Create('spike_detector')
    nest.Connect(Pyr23, Spikes)
    nest.Connect(Pyr5, Spikes)
    nest.Connect(Pyr6, Spikes)

    Multimeter = nest.Create('multimeter',
                             params={'withtime': True, 'record_from': ['V_m', 'I_syn_ex'], 'interval': 0.1})
    nest.Connect(Multimeter, Pyr23)
    nest.Connect(Multimeter, Pyr5)
    nest.Connect(Multimeter, SS4)
    nest.Connect(Multimeter, Pyr6)

    SomaMultimeter = nest.Create('multimeter', params={'withtime': True, 'record_from': ['soma_exc'], 'interval': 0.1})
    nest.Connect(SomaMultimeter, SS4)

    return Detector, Spikes, Multimeter, SomaMultimeter, Pyr23, SS4, Pyr5, Pyr6, In23, In4, In5, In6
