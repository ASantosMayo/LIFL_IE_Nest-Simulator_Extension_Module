/*
 *  aeif_psc_exp_peak.cpp
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

#include "aeif_psc_exp_peak.h"

#ifdef HAVE_GSL

// C++ includes:
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "nest_names.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< mynest::aeif_psc_exp_peak > mynest::aeif_psc_exp_peak::recordablesMap_;

namespace nest
{
/*
 * template specialization must be placed in namespace
 *
 * Override the create() method with one call to RecordablesMap::insert_()
 * for each quantity to be recorded.
 */
template <>
void
RecordablesMap< mynest::aeif_psc_exp_peak >::create()
{
  // use standard names whereever you can for consistency!
  insert_(
    nest::names::V_m, &mynest::aeif_psc_exp_peak::get_y_elem_< mynest::aeif_psc_exp_peak::State_::V_M > );
  insert_( nest::names::I_syn_ex,
    &mynest::aeif_psc_exp_peak::get_y_elem_< mynest::aeif_psc_exp_peak::State_::I_EXC > );
  insert_( nest::names::I_syn_in,
    &mynest::aeif_psc_exp_peak::get_y_elem_< mynest::aeif_psc_exp_peak::State_::I_INH > );
  insert_( nest::names::w, &mynest::aeif_psc_exp_peak::get_y_elem_< mynest::aeif_psc_exp_peak::State_::W > );
}
}


extern "C" int
mynest::aeif_psc_exp_peak_dynamics( double, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef mynest::aeif_psc_exp_peak::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const mynest::aeif_psc_exp_peak& node =
    *( reinterpret_cast< mynest::aeif_psc_exp_peak* >( pnode ) );

  const bool is_refractory = node.S_.r_ > 0;

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away...

  // Clamp membrane potential to V_reset while refractory, otherwise bound
  // it to V_peak. Do not use V_.V_peak_ here, since that is set to V_th if
  // Delta_T == 0.
  const double& V =
    is_refractory ? node.P_.V_reset_ : std::min( y[ S::V_M ], node.P_.V_peak_ );
  // shorthand for the other state variables
  const double& I_syn_ex = y[ S::I_EXC ];
  const double& I_syn_in = y[ S::I_INH ];
  const double& w = y[ S::W ];

  const double I_spike = node.P_.Delta_T == 0.
    ? 0.
    : ( node.P_.g_L * node.P_.Delta_T
        * std::exp( ( V - node.P_.V_th ) / node.P_.Delta_T ) );

  // dv/dt
  f[ S::V_M ] = is_refractory
    ? 0.
    : ( -node.P_.g_L * ( V - node.P_.E_L ) + I_spike + I_syn_ex - I_syn_in - w
        + node.P_.I_e + node.B_.I_stim_ ) / node.P_.C_m;

  f[ S::I_EXC ] = -I_syn_ex / node.P_.tau_syn_ex; // Exc. synaptic current (pA)

  f[ S::I_INH ] = -I_syn_in / node.P_.tau_syn_in; // Inh. synaptic current (pA)

  // Adaptation current w.
  f[ S::W ] = ( node.P_.a * ( V - node.P_.E_L ) - w ) / node.P_.tau_w;

  return GSL_SUCCESS;
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

mynest::aeif_psc_exp_peak::Parameters_::Parameters_()
  : V_peak_( 0.0 )    // mV
  , V_reset_( -60.0 ) // mV
  , t_ref_( 0.0 )     // ms
  , g_L( 30.0 )       // nS
  , C_m( 281.0 )      // pF
  , E_L( -70.6 )      // mV
  , Delta_T( 2.0 )    // mV
  , tau_w( 144.0 )    // ms
  , a( 4.0 )          // nS
  , b( 80.5 )         // pA
  , V_th( -50.4 )     // mV
  , tau_syn_ex( 0.2 ) // ms
  , tau_syn_in( 2.0 ) // ms
  , I_e( 0.0 )        // pA
  , gsl_error_tol( 1e-6 )
{
}

mynest::aeif_psc_exp_peak::State_::State_( const Parameters_& p )
  : r_( 0 )
{
  y_[ 0 ] = p.E_L;
  for ( size_t i = 1; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = 0;
  }
}

mynest::aeif_psc_exp_peak::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
}

mynest::aeif_psc_exp_peak::State_& mynest::aeif_psc_exp_peak::State_::operator=(
  const State_& s )
{
  assert( this != &s ); // would be bad logical error in program
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
  r_ = s.r_;
  return *this;
}

/* ----------------------------------------------------------------
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
mynest::aeif_psc_exp_peak::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, nest::names::C_m, C_m );
  def< double >( d, nest::names::V_th, V_th );
  def< double >( d, nest::names::t_ref, t_ref_ );
  def< double >( d, nest::names::g_L, g_L );
  def< double >( d, nest::names::E_L, E_L );
  def< double >( d, nest::names::V_reset, V_reset_ );
  def< double >( d, nest::names::tau_syn_ex, tau_syn_ex );
  def< double >( d, nest::names::tau_syn_in, tau_syn_in );
  def< double >( d, nest::names::a, a );
  def< double >( d, nest::names::b, b );
  def< double >( d, nest::names::Delta_T, Delta_T );
  def< double >( d, nest::names::tau_w, tau_w );
  def< double >( d, nest::names::I_e, I_e );
  def< double >( d, nest::names::V_peak, V_peak_ );
  def< double >( d, nest::names::gsl_error_tol, gsl_error_tol );
}

void
mynest::aeif_psc_exp_peak::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, nest::names::V_th, V_th );
  updateValue< double >( d, nest::names::V_peak, V_peak_ );
  updateValue< double >( d, nest::names::t_ref, t_ref_ );
  updateValue< double >( d, nest::names::E_L, E_L );
  updateValue< double >( d, nest::names::V_reset, V_reset_ );

  updateValue< double >( d, nest::names::C_m, C_m );
  updateValue< double >( d, nest::names::g_L, g_L );

  updateValue< double >( d, nest::names::tau_syn_ex, tau_syn_ex );
  updateValue< double >( d, nest::names::tau_syn_in, tau_syn_in );

  updateValue< double >( d, nest::names::a, a );
  updateValue< double >( d, nest::names::b, b );
  updateValue< double >( d, nest::names::Delta_T, Delta_T );
  updateValue< double >( d, nest::names::tau_w, tau_w );

  updateValue< double >( d, nest::names::I_e, I_e );

  updateValue< double >( d, nest::names::gsl_error_tol, gsl_error_tol );

  if ( V_reset_ >= V_peak_ )
  {
    throw nest::BadProperty( "Ensure that V_reset < V_peak ." );
  }

  if ( Delta_T < 0. )
  {
    throw nest::BadProperty( "Delta_T must be positive." );
  }
  else if ( Delta_T > 0. )
  {
    // check for possible numerical overflow with the exponential divergence at
    // spike time, keep a 1e20 margin for the subsequent calculations
    const double max_exp_arg =
      std::log( std::numeric_limits< double >::max() / 1e20 );
    if ( ( V_peak_ - V_th ) / Delta_T >= max_exp_arg )
    {
      throw nest::BadProperty(
        "The current combination of V_peak, V_th and Delta_T"
        "will lead to numerical overflow at spike time; try"
        "for instance to increase Delta_T or to reduce V_peak"
        "to avoid this problem." );
    }
  }

  if ( V_peak_ < V_th )
  {
    throw nest::BadProperty( "V_peak >= V_th required." );
  }

  if ( C_m <= 0 )
  {
    throw nest::BadProperty( "Ensure that C_m > 0" );
  }

  if ( t_ref_ < 0 )
  {
    throw nest::BadProperty( "Ensure that t_ref >= 0" );
  }

  if ( tau_syn_ex <= 0 || tau_syn_in <= 0 || tau_w <= 0 )
  {
    throw nest::BadProperty( "All time constants must be strictly positive." );
  }

  if ( gsl_error_tol <= 0. )
  {
    throw nest::BadProperty( "The gsl_error_tol must be strictly positive." );
  }
}

void
mynest::aeif_psc_exp_peak::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, nest::names::V_m, y_[ V_M ] );
  def< double >( d, nest::names::I_syn_ex, y_[ I_EXC ] );
  def< double >( d, nest::names::I_syn_in, y_[ I_INH ] );
  def< double >( d, nest::names::w, y_[ W ] );
}

void
mynest::aeif_psc_exp_peak::State_::set( const DictionaryDatum& d, const Parameters_& )
{
  updateValue< double >( d, nest::names::V_m, y_[ V_M ] );
  updateValue< double >( d, nest::names::I_syn_ex, y_[ I_EXC ] );
  updateValue< double >( d, nest::names::I_syn_in, y_[ I_INH ] );
  updateValue< double >( d, nest::names::w, y_[ W ] );
  if ( y_[ I_EXC ] < 0 || y_[ I_INH ] < 0 )
  {
    throw nest::BadProperty( "Conductances must not be negative." );
  }
}

mynest::aeif_psc_exp_peak::Buffers_::Buffers_( mynest::aeif_psc_exp_peak& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

mynest::aeif_psc_exp_peak::Buffers_::Buffers_( const Buffers_&, mynest::aeif_psc_exp_peak& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

mynest::aeif_psc_exp_peak::aeif_psc_exp_peak()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

mynest::aeif_psc_exp_peak::aeif_psc_exp_peak( const aeif_psc_exp_peak& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

mynest::aeif_psc_exp_peak::~aeif_psc_exp_peak()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ )
  {
    gsl_odeiv_step_free( B_.s_ );
  }
  if ( B_.c_ )
  {
    gsl_odeiv_control_free( B_.c_ );
  }
  if ( B_.e_ )
  {
    gsl_odeiv_evolve_free( B_.e_ );
  }
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
mynest::aeif_psc_exp_peak::init_state_( const Node& proto )
{
  const aeif_psc_exp_peak& pr = downcast< aeif_psc_exp_peak >( proto );
  S_ = pr.S_;
}

void
mynest::aeif_psc_exp_peak::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize
  B_.currents_.clear();  // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = nest::Time::get_resolution().get_ms();

  // We must integrate this model with high-precision to obtain decent results
  B_.IntegrationStep_ = std::min( 0.01, B_.step_ );

  if ( B_.s_ == 0 )
  {
    B_.s_ =
      gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_step_reset( B_.s_ );
  }

  if ( B_.c_ == 0 )
  {
    B_.c_ = gsl_odeiv_control_yp_new( P_.gsl_error_tol, P_.gsl_error_tol );
  }
  else
  {
    gsl_odeiv_control_init(
      B_.c_, P_.gsl_error_tol, P_.gsl_error_tol, 0.0, 1.0 );
  }

  if ( B_.e_ == 0 )
  {
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.e_ );
  }

  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );
  B_.sys_.function = aeif_psc_exp_peak_dynamics;

  B_.I_stim_ = 0.0;
}

void
mynest::aeif_psc_exp_peak::calibrate()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  // set the right threshold and GSL function depending on Delta_T
  if ( P_.Delta_T > 0. )
  {
    V_.V_peak = P_.V_peak_;
  }
  else
  {
    V_.V_peak = P_.V_th; // same as IAF dynamics for spikes if Delta_T == 0.
  }

  V_.refractory_counts_ = nest::Time( nest::Time::ms( P_.t_ref_ ) ).get_steps();
  // since t_ref_ >= 0, this can only fail in error
  assert( V_.refractory_counts_ >= 0 );
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
mynest::aeif_psc_exp_peak::update( const nest::Time& origin, const long from, const long to )
{
  assert(
    to >= 0 && ( nest::delay ) from < nest::kernel().connection_manager.get_min_delay() );
  assert( from < to );
  assert( State_::V_M == 0 );

  for ( long lag = from; lag < to; ++lag )
  {

	//if ( S_.y_[ State_::V_M ] == V_.refractory_counts_ )
	//if ( S_.y_[ State_::V_M ] == 10.0 )
	//{
	//	S_.y_[ State_::V_M ] = P_.V_reset_;
	//}


    double t = 0.0;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply( B_.e_,
        B_.c_,
        B_.s_,
        &B_.sys_,             // system of ODE
        &t,                   // from t
        B_.step_,             // to t <= step
        &B_.IntegrationStep_, // integration step size
        S_.y_ );              // neuronal state
      if ( status != GSL_SUCCESS )
      {
        throw nest::GSLSolverFailure( get_name(), status );
      }

      // check for unreasonable values; we allow V_M to explode
      if ( S_.y_[ State_::V_M ] < -1e3 || S_.y_[ State_::W ] < -1e6
        || S_.y_[ State_::W ] > 1e6 )
      {
        throw nest::NumericalInstability( get_name() );
      }

      // spikes are handled inside the while-loop
      // due to spike-driven adaptation
      if ( S_.r_ > 0 )
      {
        S_.y_[ State_::V_M ] = P_.V_reset_;
      }
      else if ( S_.y_[ State_::V_M ] >= V_.V_peak ) // If V_m is on peak, reset V_m.
      {
        S_.y_[ State_::W ] += P_.b; // spike-driven adaptation

        /* Initialize refractory step counter.
         * - We need to add 1 to compensate for count-down immediately after
         *   while loop.
         * - If neuron has no refractory time, set to 0 to avoid refractory
         *   artifact inside while loop.
         */
        S_.r_ = V_.refractory_counts_ > 0 ? V_.refractory_counts_ + 1 : 0;

        set_spiketime( nest::Time::step( origin.get_steps() + lag + 1 ) );
        nest::SpikeEvent se;
        nest::kernel().event_delivery_manager.send( *this, se, lag );
	
	if ( S_.r_ ==  V_.refractory_counts_ )
	{
		S_.y_[ State_::V_M ] = 10.0;
	}

      }
    }

    // decrement refractory count
    if ( S_.r_ > 0 )
    {
	if ( S_.r_ ==  V_.refractory_counts_ )  // V_m is on peak at first refractory count 
	{
		S_.y_[ State_::V_M ] = 10.0;
	}
	else
	{
      		S_.y_[ State_::V_M ] = P_.V_reset_;
	}

      --S_.r_;
    }

    S_.y_[ State_::I_EXC ] += B_.spike_exc_.get_value( lag );
    S_.y_[ State_::I_INH ] += B_.spike_inh_.get_value( lag );

    // set new input current
    B_.I_stim_ = B_.currents_.get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
mynest::aeif_psc_exp_peak::handle( nest::SpikeEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  if ( e.get_weight() > 0.0 )
  {
    B_.spike_exc_.add_value( e.get_rel_delivery_steps(
                               nest::kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
  else
  {
    B_.spike_inh_.add_value( e.get_rel_delivery_steps(
                               nest::kernel().simulation_manager.get_slice_origin() ),
      -e.get_weight() * e.get_multiplicity() );
  } // keep conductances positive
}

void
mynest::aeif_psc_exp_peak::handle( nest::CurrentEvent& e )
{
  assert( e.get_delay_steps() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(
    e.get_rel_delivery_steps( nest::kernel().simulation_manager.get_slice_origin() ),
    w * c );
}

void
mynest::aeif_psc_exp_peak::handle( nest::DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL
