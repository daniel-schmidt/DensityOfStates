/*
 * MetropolisWangLandauStep.cpp
 *
 *  Created on: 09.09.2016
 *      Author: dschmidt
 */

#include "MetropolisWangLandauStep.h"

namespace FermiOwn {

MetropolisWangLandauStep::MetropolisWangLandauStep( Field<int>& spinField, const IsingHamiltonian& hamilton, double initialF, const double finalFTol, const double flatness, const size_t checkHistEvery, std::ranlux48* rndGen ) :
			MetropolisStep( rndGen ),
			spin( spinField ),
			H( hamilton ),
			x_dist( 0, spin.getSize()-1 ),
			new_x( 0 ),
			newE( 0. ),
			oldE( -H.calculateEnergy() ),
			flat( flatness ),
			checkHist( checkHistEvery ),
			f( initialF ),
			Ftol( finalFTol )
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

bool MetropolisWangLandauStep::onConfig( size_t confNum ) {
	dos[oldE] += f;
	hist[oldE] ++;

	if( (confNum+1)%checkHist == 0 && histFlat() ) {
		f *= 0.5;
		if( f < Ftol ) return false;
		hist.clear();
	}
	return true;
}

double MetropolisWangLandauStep::getDos( const double energy ) const {
	auto it = dos.find( energy );
	if( it == dos.end() ) {
		return 0;
	} else {
		return it->second;
	}
}

bool MetropolisWangLandauStep::histFlat() {
	double meanHist = 0.;
	double minHist = hist.begin()->second;
	for( auto histPair : hist ) {
		meanHist += histPair.second;
		if( histPair.second < minHist ) {
			minHist = histPair.second;
		}
	}
	meanHist /= hist.size();
	std::cout << "meanHist: " << meanHist << " minHist: " << minHist << " > " << meanHist*flat << "?";
	if( meanHist*flat < minHist ) {
		std::cout << " Histogram flat!" << std::endl;
		return true;
	} else {
		std::cout << " Not flat!" << std::endl;
		return false;
	}
}

} /* namespace FermiOwn */
