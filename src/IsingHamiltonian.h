/*
 * IsingHamiltonian.h
 *
 *  Created on: 29.08.2016
 *      Author: dschmidt
 */

#ifndef SRC_ISINGHAMILTONIAN_H_
#define SRC_ISINGHAMILTONIAN_H_

#include "Lattice.h"
#include "Field.h"

namespace FermiOwn {

class IsingHamiltonian {
public:
	IsingHamiltonian( const Lattice& lattice, Field<int>& spinField, const double coupling );
	virtual ~IsingHamiltonian();

	/**
	 * @brief Returns the negative value of the Hamiltonian
	 */
	double calculateEnergy() const;

private:
	const Lattice& lat;
	Field<int> & spin;
	const double J;
};

} /* namespace FermiOwn */

#endif /* SRC_ISINGHAMILTONIAN_H_ */
