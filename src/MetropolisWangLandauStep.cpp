/*
 * MetropolisWangLandauStep.cpp
 *
 *  Created on: 09.09.2016
 *      Author: dschmidt
 */

#include "MetropolisWangLandauStep.h"

namespace FermiOwn {

MetropolisWangLandauStep::MetropolisWangLandauStep( Field<int>& spinField, const IsingHamiltonian& hamilton, std::ranlux48* rndGen ) :
			MetropolisStep( rndGen ),
			spin( spinField ),
			H( hamilton ),
			x_dist( 0, spin.getSize()-1 ),
			new_x( 0 ),
			newE( 0. ),
			oldE( -H.calculateEnergy() )
{}

MetropolisWangLandauStep::~MetropolisWangLandauStep() {
	// TODO Auto-generated destructor stub
}

void MetropolisWangLandauStep::propose() {
	new_x = x_dist( *rnd ); //change spin at a random point
	spin( new_x ) *= -1;
	newE = -H.calculateEnergy();
}

double MetropolisWangLandauStep::change() {
	return std::exp( getDos( oldE ) - getDos( newE ) );
}

void MetropolisWangLandauStep::accept() {
	oldE = newE;
}

void MetropolisWangLandauStep::reject() {
	spin( new_x ) *= -1;
}

double MetropolisWangLandauStep::getDos( const double energy ) const {
	auto it = dos.find( energy );
	if( it == dos.end() ) {
		return 0;
	} else {
		return it->second;
	}
}

} /* namespace FermiOwn */
