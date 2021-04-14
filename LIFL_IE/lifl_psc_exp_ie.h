/*
 *  lifl_psc_exp_ie.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef LIFL_PSC_EXP_IE_H
#define LIFL_PSC_EXP_IE_H

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "recordables_map.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"


// Includes from sli:
#include "dictdatum.h"


namespace mynest
{
/* BeginDocumentation
   Name: lifl_psc_exp_ie - Leaky integrate-and-fire with spike latency neuron model with exponential
                       PSCs integrated in Exponential function and Intrinsic excitability plasticity.

   Description:
   lifl_psc_exp_iep is an implementation of a leaky integrate-and-fire model
   with exponential shaped postsynaptic currents (PSCs) according to [1].
   Thus, postsynaptic currents have an infinitely short rise time.

   The threshold crossing is followed by an absolute refractory period (t_ref)
   during which the membrane potential is clamped to the resting potential
   and spiking is prohibited.

   The linear subthresold dynamics is integrated by the Exact
   Integration scheme [2]. The neuron dynamics is solved on the time
   grid given by the computation step size. Incoming as well as emitted
   spikes are forced to that grid.

   An additional state variable and the corresponding differential
   equation represents a piecewise constant external current.

   The general framework for the consistent formulation of systems with
   neuron like dynamics interacting by point events is described in
   [2]. A flow chart can be found in [3].

   ADDITIONAL FEATURES:

   This neuron model shows SPIKE LATENCY and INTRINSIC EXCITABILITY.
   
   SPIKE LATENCY implies an over-threshold state until the spike generation, and allows the time between threshold crossing and spike generation to vary depending on the membrane potential reached during the summation process. Thereby, SL enables a neuron to adjust its spike time by modulating its internal state. The SL mechanism present in the LIFL has been extracted from the biological Hodgkin-Huxley (HH) model (see [5]), representing a relevant biologicalresource to implement temporal coding mechanism by many types of neurons. 

    INTRINSIC EXCITABILITY plasticity is a biological concept that has emerged during the last years. It can be envisaged as an internal modulator of the synaptic input strength that is called intrinsic neuronal excitability (IE) [6]. This effect occurs in parallel to the induction of long-term synaptic modifications and produces changes in the neuron’s conductivity, by modulating the h-channel activity). It is defined as a non-synaptic plasticity process because changes are produced within the postsynaptic cell and not at the synapse. Thus, in presence of IE, the temporal difference between presynaptic and postsynaptic spike not only produces synaptic long term potentiation/depression (LTP-s/LTD-s), but also an intrinsic excitability  increase/decrease (LTP-IE/LTD-IE). This mechanism causes a change in the neuron’s conductance, thus strengthens the modulation of  the membrane potential, and finally controls the Spike Latency (i.e. recreating the Multi Neural Spike-sequence Detector MNSD [7]).
   

   Remarks:
   The present implementation uses individual variables for the
   components of the state vector and the non-zero matrix elements of
   the propagator.  Because the propagator is a lower triangular matrix
   no full matrix multiplication needs to be carried out and the
   computation can be done "in place" i.e. no temporary state vector
   object is required.

   The template support of recent C++ compilers enables a more succinct
   formulation without loss of runtime performance already at minimal
   optimization levels. A future version of lifl_psc_exp_ie will probably
   address the problem of efficient usage of appropriate vector and
   matrix objects.

   Parameters:
   The following parameters can be set in the status dictionary.

   E_L          double - Resting membrane potential in mV.
   C_m          double - Capacity of the membrane in pF
   tau_m        double - Membrane time constant in ms.
   tau_syn_ex   double - Time constant of postsynaptic excitatory currents in ms
   tau_syn_in   double - Time constant of postsynaptic inhibitory currents in ms
   t_ref        double - Duration of refractory period (V_m = V_reset) in ms.
   V_m          double - Membrane potential in mV
   V_th         double - Spike threshold in mV.
   V_reset      double - Reset membrane potential after a spike in mV.
   I_e          double - Constant input current in pA.
   t_spike      double - Point in time of last spike in ms.
 
   Specific for spike latency and intrinsic excitability

   nueromodulator array - array with GID of neurons which modulate the Intrinsic plasticity.
   lambda   double . Value of change per intrinsic plasticity effect.
   tau      double . value of window of the intrinsic plasticity effect.
   std_mod  bool   . Swhich ON (true) or OFF (false) the intrinsic plasticity effect.

Remarks:

   If tau_m is very close to tau_syn_ex or tau_syn_in, the model
   will numerically behave as if tau_m is equal to tau_syn_ex or
   tau_syn_in, respectively, to avoid numerical instabilities.
   For details, please see IAF_neurons_singularity.ipynb in the
   NEST source code (docs/model_details).

   lifl_psc_exp_ie can handle current input in two ways: Current input
   through receptor_type 0 are handled as stepwise constant current
   input as in other iaf models, i.e., this current directly enters
   the membrane potential equation. Current input through
   receptor_type 1, in contrast, is filtered through an exponential
   kernel with the time constant of the excitatory synapse,
   tau_syn_ex. For an example application, see [4].

   References:
   [1] Misha Tsodyks, Asher Uziel, and Henry Markram (2000) Synchrony Generation
   in Recurrent Networks with Frequency-Dependent Synapses, The Journal of
   Neuroscience, 2000, Vol. 20 RC50 p. 1-5
   [2] Rotter S & Diesmann M (1999) Exact simulation of time-invariant linear
   systems with applications to neuronal modeling. Biologial Cybernetics
   81:381-402.
   [3] Diesmann M, Gewaltig M-O, Rotter S, & Aertsen A (2001) State space
   analysis of synchronous spiking in cortical neural networks.
   Neurocomputing 38-40:565-571.
   [4] Schuecker J, Diesmann M, Helias M (2015) Modulated escape from a
   metastable state driven by colored noise.
   Physical Review E 92:052119
   [5] Salerno M, Susi G, Cristini A. Accurate latency characterization for very
   large asynchronous spiking neural networks. International Conference on Bioinformatics
    Models, Methods andAlgorithms (BIOINFORMATICS 2011). SciTePress; 2011. pp. 116–124.
   [6] Debanne D. Spike-timing dependent plasticity beyond synapse - pre- and post-synaptic plasticity of intrinsic neuronal excitability. Frontiers in Synaptic Neuroscience. 2010. doi:10.3389/fnsyn.2010.00021
   [7] Susi G, Antón Toro L, Canuet L, López ME, Maestú F, Mirasso CR, et al. A Neuro-Inspired System for Online Learning and Recognition of Parallel Spike Trains, Based on Spike Latency, and Heterosynaptic STDP. Front Neurosci. 2018;12: 780.

   Sends: SpikeEvent

   Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

   SeeAlso: lifl_psc_exp_ie_ps

   FirstVersion: 2019-2020
   Author: Alejandro Santos-Mayo, based on iaf_psc_exp
*/

/**
 * Leaky integrate-and-fire neuron with exponential PSCs.
 */
class lifl_psc_exp_ie : public nest::Archiving_Node
{

public:
  lifl_psc_exp_ie();
  lifl_psc_exp_ie( const lifl_psc_exp_ie& );

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and
   * Hiding
   */
  using nest::Node::handle;
  using nest::Node::handles_test_event;

  nest::port send_test_event( nest::Node&, nest::rport, nest::synindex, bool );

  void handle( nest::SpikeEvent& );
  void handle( nest::CurrentEvent& );
  void handle( nest::DataLoggingRequest& );

  nest::port handles_test_event( nest::SpikeEvent&, nest::rport );
  nest::port handles_test_event( nest::CurrentEvent&, nest::rport );
  nest::port handles_test_event( nest::DataLoggingRequest&, nest::rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();

  void update( const nest::Time&, const long, const long );

  // The next two classes need to be friends to access the State_ class/member
  friend class nest::RecordablesMap< lifl_psc_exp_ie >;
  friend class nest::UniversalDataLogger< lifl_psc_exp_ie >;

  // ----------------------------------------------------------------

  /**
   * Independent parameters of the model.
   */
  struct Parameters_
  {

    /** Membrane time constant in ms. */
    double Tau_;

    /** Membrane capacitance in pF. */
    double C_;

    /** Refractory period in ms. */
    double t_ref_;

    /** Resting potential in mV. */
    double E_L_;

    /** External current in pA */
    double I_e_;

    /** Threshold, RELATIVE TO RESTING POTENTAIL(!).
        I.e. the real threshold is (E_L_+Theta_). */
    double Theta_;

    /** reset value of the membrane potential */
    double V_reset_;

    /** Time constant of excitatory synaptic current in ms. */
    double tau_ex_;

    /** Time constant of inhibitory synaptic current in ms. */
    double tau_in_;

	// Spike latency and Intrinsic Excitability parameters.

    /** Auxiliar temporal variable to calculate the Spike Latency */
    long stepms;

    /** Time step resolution to take care of in operations */
    double dt;

    /** Global ID of neuro-modulator neurons */
    std::vector< long > stimulator_;

    /** Lambda value indicates the change amplitude in IE plasticity */
    double lambda;

    /** Tau value indicates the time window for IE plasticity */
    double tau;

    /** std_mod can swich on / off the IE plasticity mechanism */
    bool std_mod;







    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary

    /** Set values from dictionary.
     * @returns Change in reversal potential E_L, to be passed to State_::set()
     */
    double set( const DictionaryDatum& );
  };

  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   */
  struct State_
  {
    // state variables
    //! synaptic stepwise constant input current, variable 0
    double i_0_;
    double i_1_;      //!< presynaptic stepwise constant input current
    double i_syn_ex_; //!< postsynaptic current for exc. inputs, variable 1
    double i_syn_in_; //!< postsynaptic current for inh. inputs, variable 1
    double V_m_;      //!< membrane potential, variable 2
   
       // Spike latency and Intrinsic Excitability States

    double Vpositive; //!< Auxiliar value used to correctly calculate spike latency
    double enhancement; //!< Intrinsic Excitability value modulator of incoming current.
        
    std::vector< long > hist_;

    std::vector< double > t_lastspike_;



    //! absolute refractory counter (no membrane potential propagation)
    int r_ref_;

    State_(); //!< Default initialization

    void get( DictionaryDatum&, const Parameters_& ) const;

    /** Set values from dictionary.
     * @param dictionary to take data from
     * @param current parameters
     * @param Change in reversal potential E_L specified by this dict
     */
    void set( const DictionaryDatum&, const Parameters_&, const double );
  };

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( lifl_psc_exp_ie& );
    Buffers_( const Buffers_&, lifl_psc_exp_ie& );

    /** buffers and sums up incoming spikes/currents */
    nest::RingBuffer spikes_ex_;
    nest::RingBuffer spikes_in_;
    std::vector< nest::RingBuffer > currents_;

    //! Logger for all analog data
    nest::UniversalDataLogger< lifl_psc_exp_ie > logger_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    /** Amplitude of the synaptic current.
        This value is chosen such that a post-synaptic potential with
        weight one has an amplitude of 1 mV.
        @note mog - I assume this, not checked.
    */
    //    double PSCInitialValue_;

    // time evolution operator
    double P20_;
    double P11ex_;
    double P11in_;
    double P21ex_;
    double P21in_;
    double P22_;

    double weighted_spikes_ex_;
    double weighted_spikes_in_;

    int RefractoryCounts_;
  };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out the real membrane potential
  inline double
  get_V_m_() const
  {
    return S_.V_m_ + P_.E_L_;
  }

  // INTRINSIC EXCTIABILITY value
  double
  get_soma_exc_() const
  {
    return S_.enhancement;
  }

  inline double
  get_weighted_spikes_ex_() const
  {
    return V_.weighted_spikes_ex_;
  }

  inline double
  get_weighted_spikes_in_() const
  {
    return V_.weighted_spikes_in_;
  }

  inline double
  get_I_syn_ex_() const
  {
    return S_.i_syn_ex_;
  }

  inline double
  get_I_syn_in_() const
  {
    return S_.i_syn_in_;
  }

  // ----------------------------------------------------------------

  /**
   * @defgroup lifl_psc_exp_ie_data
   * Instances of private data structures for the different types
   * of data pertaining to the model.
   * @note The order of definitions is important for speed.
   * @{
   */
  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;
  /** @} */

  //! Mapping of recordables names to access functions
  static nest::RecordablesMap< lifl_psc_exp_ie > recordablesMap_;
};


inline nest::port
mynest::lifl_psc_exp_ie::send_test_event( nest::Node& target,
  nest::rport receptor_type,
  nest::synindex,
  bool )
{
  nest::SpikeEvent e;
  e.set_sender( *this );
  return target.handles_test_event( e, receptor_type );
}

inline nest::port
mynest::lifl_psc_exp_ie::handles_test_event( nest::SpikeEvent&, nest::rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw nest::UnknownReceptorType( receptor_type, get_name() );
  }
  return 0;
}

inline nest::port
mynest::lifl_psc_exp_ie::handles_test_event( nest::CurrentEvent&, nest::rport receptor_type )
{
  if ( receptor_type == 0 )
  {
    return 0;
  }
  else if ( receptor_type == 1 )
  {
    return 1;
  }
  else
  {
    throw nest::UnknownReceptorType( receptor_type, get_name() );
  }
}

inline nest::port
mynest::lifl_psc_exp_ie::handles_test_event( nest::DataLoggingRequest& dlr, nest::rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw nest::UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
lifl_psc_exp_ie::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d, P_ );
  Archiving_Node::get_status( d );

  ( *d )[ nest::names::recordables ] = recordablesMap_.get_list();
}

inline void
lifl_psc_exp_ie::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_;                 // temporary copy in case of errors
  const double delta_EL = ptmp.set( d ); // throws if BadProperty
  State_ stmp = S_;                      // temporary copy in case of errors
  stmp.set( d, ptmp, delta_EL );         // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace

#endif // lifl_psc_exp_ie_H
