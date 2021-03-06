/*
 *  lifl_psc_exp_ie.cpp
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

#include "lifl_psc_exp_ie.h"

// C++ includes:
#include <limits>

// Includes from libnestutil:
#include "numerics.h"
#include "propagator_stability.h"

#include <algorithm>
#include <iterator>

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

#include "histentry.h"
#include "node.h"

#include "common_synapse_properties.h"//*_
#include "connection.h"
#include "connector_base.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< mynest::lifl_psc_exp_ie > mynest::lifl_psc_exp_ie::recordablesMap_;

namespace nest
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< mynest::lifl_psc_exp_ie >::create()
{
  // use standard names whereever you can for consistency!
  insert_( names::V_m, &mynest::lifl_psc_exp_ie::get_V_m_ );
  insert_( names::weighted_spikes_ex, &mynest::lifl_psc_exp_ie::get_weighted_spikes_ex_ );
  insert_( names::weighted_spikes_in, &mynest::lifl_psc_exp_ie::get_weighted_spikes_in_ );
  insert_( names::I_syn_ex, &mynest::lifl_psc_exp_ie::get_I_syn_ex_ );
  insert_( names::I_syn_in, &mynest::lifl_psc_exp_ie::get_I_syn_in_ );
  insert_( names::soma_exc, &mynest::lifl_psc_exp_ie::get_soma_exc_ );
}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

mynest::lifl_psc_exp_ie::Parameters_::Parameters_()
  : Tau_( 10.0 )             // in ms
  , C_( 250.0 )              // in pF
  , t_ref_( 2.0 )            // in ms
  , E_L_( -70.0 )            // in mV
  , I_e_( 0.0 )              // in pA
  , Theta_( -55.0 - E_L_ )   // relative E_L_
  , V_reset_( -70.0 - E_L_ ) // in mV
  , tau_ex_( 2.0 )           // in ms
  , tau_in_( 2.0 )           // in ms

  // Latency and Intrinsic excitability 
  , dt( 0.0 )    // ms Parameter of time resolution
  , lambda(0.0001) // Lambda - Intrinsic Plasticity enhance
  , tau(12.5) // Tau - Intrinsic Plasticity window
  , std_mod(true) // ON/OFF of the spike time dependent modification
  , stimulator_()

{
}

mynest::lifl_psc_exp_ie::State_::State_()
  : i_0_( 0.0 )
  , i_syn_ex_( 0.0 )
  , i_syn_in_( 0.0 )
  , V_m_( 0.0 )
  , r_ref_( 0 )

  // Latency and Intrinsic excitability 
  //, refr_count( 0 )
  , enhancement(1.0)
  , t_lastspike_()
  , hist_(0.0)
{
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
mynest::lifl_psc_exp_ie::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, nest::names::E_L, E_L_ ); // resting potential
  def< double >( d, nest::names::I_e, I_e_ );
  def< double >( d, nest::names::V_th, Theta_ + E_L_ ); // threshold value
  def< double >( d, nest::names::V_reset, V_reset_ + E_L_ );
  def< double >( d, nest::names::C_m, C_ );
  def< double >( d, nest::names::tau_m, Tau_ );
  def< double >( d, nest::names::tau_syn_ex, tau_ex_ );
  def< double >( d, nest::names::tau_syn_in, tau_in_ );
  def< double >( d, nest::names::t_ref, t_ref_ );

  // Spike latency and Intrinsc Excitability
  def< double >( d, nest::names::lambda, lambda );
  def< double >( d, nest::names::tau, tau );
  def< bool >(d, nest::names::std_mod, std_mod );

const size_t n_stims = stimulator_.size();
std::vector< long >* stims = new std::vector< long >();
stims -> reserve(  n_stims );
for ( size_t n = 0; n < n_stims; ++n)
{
    stims->push_back( stimulator_[ n ]);
}

(*d )[ nest::names::stimulator] = IntVectorDatum( stims );

}

double
mynest::lifl_psc_exp_ie::Parameters_::set( const DictionaryDatum& d )
{
  // if E_L_ is changed, we need to adjust all variables defined relative to
  // E_L_
  const double ELold = E_L_;
  updateValue< double >( d, nest::names::E_L, E_L_ );
  const double delta_EL = E_L_ - ELold;

  if ( updateValue< double >( d, nest::names::V_reset, V_reset_ ) )
  {
    V_reset_ -= E_L_;
  }
  else
  {
    V_reset_ -= delta_EL;
  }

  if ( updateValue< double >( d, nest::names::V_th, Theta_ ) )
  {
    Theta_ -= E_L_;
  }
  else
  {
    Theta_ -= delta_EL;
  }

  updateValue< double >( d, nest::names::I_e, I_e_ );
  updateValue< double >( d, nest::names::C_m, C_ );
  updateValue< double >( d, nest::names::tau_m, Tau_ );
  updateValue< double >( d, nest::names::tau_syn_ex, tau_ex_ );
  updateValue< double >( d, nest::names::tau_syn_in, tau_in_ );
  updateValue< double >( d, nest::names::t_ref, t_ref_ );

  // SL and IE
  updateValue< double >( d, nest::names::lambda, lambda );
  updateValue< double >( d, nest::names::tau, tau );
  updateValue< std::vector< long > >( d, nest::names::stimulator, stimulator_ );
  updateValue< bool >(d,nest::names::std_mod, std_mod );


  if ( V_reset_ >= Theta_ )
  {
    throw nest::BadProperty( "Reset potential must be smaller than threshold." );
  }
  if ( C_ <= 0 )
  {
    throw nest::BadProperty( "Capacitance must be strictly positive." );
  }
  if ( Tau_ <= 0 || tau_ex_ <= 0 || tau_in_ <= 0 )
  {
    throw nest::BadProperty(
      "Membrane and synapse time constants must be strictly positive." );
  }
  if ( t_ref_ < 0 )
  {
    throw nest::BadProperty( "Refractory time must not be negative." );
  }

  return delta_EL;
}



void
mynest::lifl_psc_exp_ie::State_::get( DictionaryDatum& d, const Parameters_& p ) const
{
  def< double >( d, nest::names::V_m, V_m_ + p.E_L_ ); // Membrane potential
  ( *d )[ nest::names::soma_exc ] = enhancement;
}

void
mynest::lifl_psc_exp_ie::State_::set( const DictionaryDatum& d,
  const Parameters_& p,
  double delta_EL )
{
  if ( updateValue< double >( d, nest::names::V_m, V_m_ ) )
  {
    V_m_ -= p.E_L_;
  }
  else
  {
    V_m_ -= delta_EL;
  }
  updateValue< double >( d, nest::names::soma_exc, enhancement );
}

mynest::lifl_psc_exp_ie::Buffers_::Buffers_( lifl_psc_exp_ie& n )
  : logger_( n )
{
}

mynest::lifl_psc_exp_ie::Buffers_::Buffers_( const Buffers_&, lifl_psc_exp_ie& n )
  : logger_( n )
{
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

mynest::lifl_psc_exp_ie::lifl_psc_exp_ie()
  : Archiving_Node()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

mynest::lifl_psc_exp_ie::lifl_psc_exp_ie( const lifl_psc_exp_ie& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
mynest::lifl_psc_exp_ie::init_state_( const Node& proto )
{
  const lifl_psc_exp_ie& pr = downcast< lifl_psc_exp_ie >( proto );
  S_ = pr.S_;
}

void
mynest::lifl_psc_exp_ie::init_buffers_()
{
  B_.spikes_ex_.clear(); // includes resize
  B_.spikes_in_.clear(); // includes resize
  B_.currents_.clear();  // includes resize
  B_.logger_.reset();
 /// Archiving_Node::clear_history_();
}

void
mynest::lifl_psc_exp_ie::calibrate()
{
  B_.currents_.resize( 2 );
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  const double h = nest::Time::get_resolution().get_ms();

  P_.dt = h;      // We save resolution for further calculations.

  // numbering of state vaiables: i_0 = 0, i_syn_ = 1, V_m_ = 2

  // commented out propagators: forward Euler
  // needed to exactly reproduce Tsodyks network

  // these P are independent
  V_.P11ex_ = std::exp( -h / P_.tau_ex_ );
  // P11ex_ = 1.0-h/tau_ex_;

  V_.P11in_ = std::exp( -h / P_.tau_in_ );
  // P11in_ = 1.0-h/tau_in_;

  V_.P22_ = std::exp( -h / P_.Tau_ );
  // P22_ = 1.0-h/Tau_;

  // these are determined according to a numeric stability criterion
  V_.P21ex_ = propagator_32( P_.tau_ex_, P_.Tau_, P_.C_, h );
  V_.P21in_ = propagator_32( P_.tau_in_, P_.Tau_, P_.C_, h );

  // P21ex_ = h/C_;
  // P21in_ = h/C_;

  V_.P20_ = P_.Tau_ / P_.C_ * ( 1.0 - V_.P22_ );
  // P20_ = h/C_;

  // t_ref_ specifies the length of the absolute refractory period as
  // a double in ms. The grid based lifl_psc_exp_ie can only handle refractory
  // periods that are integer multiples of the computation step size (h).
  // To ensure consistency with the overall simulation scheme such conversion
  // should be carried out via objects of class nest::Time. The conversion
  // requires 2 steps:
  //     1. A time object r is constructed, defining representation of
  //        t_ref_ in tics. This representation is then converted to computation
  //        time steps again by a strategy defined by class nest::Time.
  //     2. The refractory time in units of steps is read out get_steps(), a
  //        member function of class nest::Time.
  //
  // Choosing a t_ref_ that is not an integer multiple of the computation time
  // step h will lead to accurate (up to the resolution h) and self-consistent
  // results. However, a neuron model capable of operating with real valued
  // spike time may exhibit a different effective refractory time.

  V_.RefractoryCounts_ = nest::Time( nest::Time::ms( P_.t_ref_ ) ).get_steps();
  // since t_ref_ >= 0, this can only fail in error
  assert( V_.RefractoryCounts_ >= 0 );

  // INITIALIZATIONS
  const size_t n_stims = P_.stimulator_.size();
  for(size_t i = 0; i < n_stims; i++)
  {
   S_.t_lastspike_.push_back(0.0);
  }

}

void
mynest::lifl_psc_exp_ie::update( const nest::Time& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( nest::delay ) from < nest::kernel().connection_manager.get_min_delay() );
  assert( from < to );

  // evolve from timestep 'from' to timestep 'to' with steps of h each
  for ( long lag = from; lag < to; ++lag )
  {
     if (S_.V_m_ >= 105.0)
     {
     S_.V_m_ = 105.0;
     // send spike, and set spike time in archive.
     set_spiketime( nest::Time::step( origin.get_steps() + lag + 1 ) );
     nest::SpikeEvent se;
     nest::kernel().event_delivery_manager.send( *this, se, lag );
     S_.r_ref_ = V_.RefractoryCounts_;

     S_.hist_.push_back(nest::Archiving_Node::get_spiketime_ms());

     }

    if ( S_.r_ref_ == 0 ) // neuron not refractory, so evolve V
    {
      // Implementing SPIKE LATENCY feature

      if (S_.V_m_ > 15.6) // 15.6 is the value calculated for this specific SpikeLatency
      {
        S_.Vpositive = S_.V_m_ / 15;	
      	S_.V_m_ = S_.V_m_ + (pow((S_.Vpositive-1),2)*P_.dt)/(1-(S_.Vpositive - 1)*P_.dt) * 15;

      if (S_.V_m_ >= 105.0)
      {
      S_.V_m_ = 105.0;  // Put peak of spike and send it, setting spike time in archive.
      S_.V_m_ = 105.0;
      set_spiketime( nest::Time::step( origin.get_steps() + lag + 1 ) );
      nest::SpikeEvent se;
      nest::kernel().event_delivery_manager.send( *this, se, lag );
      S_.r_ref_ = V_.RefractoryCounts_;

      S_.hist_.push_back(lag);
      }

      }
      else
      {
        S_.V_m_ = S_.V_m_ * V_.P22_ + S_.i_syn_in_ * V_.P21in_ 
          + (S_.i_syn_ex_ * V_.P21ex_ + ( P_.I_e_ + S_.i_0_ ) * V_.P20_) * S_.enhancement; // Compute V_m of neuron
      }
    }
    else
    {
      --S_.r_ref_;
      S_.V_m_ = 0.1;

    } // neuron is absolute refractory

    // exponential decaying PSCs
    S_.i_syn_ex_ *= V_.P11ex_;
    S_.i_syn_in_ *= V_.P11in_;

    // add evolution of presynaptic input current
    S_.i_syn_ex_ += ( 1. - V_.P11ex_ ) * S_.i_1_;

    // the spikes arriving at T+1 have an immediate effect on the state of the
    // neuron

    V_.weighted_spikes_ex_ = B_.spikes_ex_.get_value( lag );
    V_.weighted_spikes_in_ = B_.spikes_in_.get_value( lag );

    S_.i_syn_ex_ += V_.weighted_spikes_ex_;
    S_.i_syn_in_ += V_.weighted_spikes_in_;

    // set new input current
    S_.i_0_ = B_.currents_[ 0 ].get_value( lag );
    S_.i_1_ = B_.currents_[ 1 ].get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
mynest::lifl_psc_exp_ie::handle( nest::SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  if ( e.get_weight() >= 0.0 )
  {
    B_.spikes_ex_.add_value( e.get_rel_delivery_steps(
                               nest::kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
  else
  {
    B_.spikes_in_.add_value( e.get_rel_delivery_steps(
                               nest::kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }






  if (P_.std_mod)   // Implementing INTRINSIC EXCITABILITY (IE) Plasticity
  {
    size_t origSize = P_.stimulator_.size();
    
    for (size_t i=0; i < origSize; i++)
    {
      long modulator = P_.stimulator_[ i ];
      long source_gid = e.get_sender_gid();

      if (source_gid == modulator) // If gID of input current is from an Stimulator (IE modulator)
      {

  // For a new synapse, t_lastspike_ contains the point in time of the last
  // spike. So we initially read the
  // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of
  // history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
  // details.
  //target->get_history( 0.0, t_spike, &start, &finish );


	double t_spike = e.get_stamp().get_ms(); 

	size_t histsize = S_.hist_.size();


 	// Take last spikes (history) and compute the LTP-IE or LTD-IE Plasticity changes

	long lstspk = S_.t_lastspike_[i];
	auto r = std::find_if(std::begin(S_.hist_), std::end(S_.hist_), [lstspk](int its){return its > lstspk;});
	while (r != std::end(S_.hist_)){ 

		S_.enhancement = S_.enhancement + std::exp(((S_.t_lastspike_[i]) - S_.hist_[std::distance(std::begin(S_.hist_), r)] )/P_.tau)*P_.lambda;

		S_.enhancement = S_.enhancement - std::exp((S_.hist_[std::distance(std::begin(S_.hist_), r)]  - t_spike)/P_.tau)*P_.lambda;

		r = std::find_if(std::next(r), std::end(S_.hist_), [lstspk](int its){return its > lstspk;});

	}

	S_.t_lastspike_[i] = t_spike; // Save the last spike of this neuron for next occasion
      }    
    }  
  }
}


void
mynest::lifl_psc_exp_ie::handle( nest::CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  if ( 0 == e.get_rport() )
  {
    B_.currents_[ 0 ].add_value(
      e.get_rel_delivery_steps(
        nest::kernel().simulation_manager.get_slice_origin() ),
      w * c );
  }
  if ( 1 == e.get_rport() )
  {
    B_.currents_[ 1 ].add_value(
      e.get_rel_delivery_steps(
        nest::kernel().simulation_manager.get_slice_origin() ),
      w * c );
  }
}

void
mynest::lifl_psc_exp_ie::handle( nest::DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
