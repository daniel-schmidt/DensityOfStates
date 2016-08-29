/*
 * IsingHamiltonian.cpp
 *
 *  Created on: 29.08.2016
 *      Author: dschmidt
 */

#include "IsingHamiltonian.h"

namespace FermiOwn {

IsingHamiltonian::IsingHamiltonian( const Lattice& lattice, Field<int>& spinField, double coupling ) :
 lat( lattice ),
 spin( spinField ),
 J( coupling )
{}

IsingHamiltonian::~IsingHamiltonian() {
	// TODO Auto-generated destructor stub
}

double IsingHamiltonian::calculateEnergy() const {
	double sum = 0;
	for( size_t x = 0; x < lat.getVol(); x++ ) {
		const std::vector< size_t > nn = lat.getNeighbours( x );
		for( size_t y : nn ) {
			sum += spin(y)*spin(x);
		}
	}
	sum *= J;
	return sum;
}

} /* namespace FermiOwn */
