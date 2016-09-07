/*
 * TrunkatedMetropolis_test.cpp
 *
 *  Created on: 07.09.2016
 *      Author: dschmidt
 */

#include "Lattice.h"
#include "IsingHamiltonian.h"
#include "MetropolisStep.h"
#include "ConfigGenerator.h"

//#include "TrunkatedMetropolis.h"

int main( int argc, char** argv ) {
	using namespace FermiOwn;

	size_t Nt = 16;
	size_t Ns = 16;
	size_t dim = 2;

	Lattice lat( Nt, Ns, dim );

	size_t dofPerPoint = 1;
	double J = 1.;
	double beta = atof( argv[1] );

	size_t numThermal = 1;
	size_t numConfs = 10000;
//	size_t numUpPerConf = lat.getVol() * 20;
	size_t numUpPerConf = 1;

	std::ranlux48 rndGen;

	std::uniform_int_distribution<int> spin_dist(0,1);
	std::uniform_int_distribution<int> x_dist(0, lat.getVol()-1 );

	int new_x = 0;

	// initialization of physical parameters

	Field<int> spin( lat.getVol(), dofPerPoint, &rndGen, oneInit);

	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( spin_dist( rndGen) ) {
			spin( x ) *= -1;
		}
	}
	IsingHamiltonian H( lat, spin, J );

	double E0 = H.calculateEnergy();
	double currE = E0;
	double dE = atof( argv[2] );

	std::cout << "#Running algorithm for E0 = " << E0 << " and dE = " << dE << std::endl;
	double Echange = 0.;
	// metropolis initialization
	auto propose = [&] () {
//		std::cout << "Propose" << std::endl;
		do {
			new_x = x_dist( rndGen );
			Echange = beta*H.changeAt( new_x );
//			std::cout << " change in E: " << Echange << std::endl;
		} while( currE+Echange > E0+dE || currE+Echange < E0-dE );
	};
	auto change = [&]() {
		return std::exp( -Echange );
	};
	auto accept = [&]() {
		spin( new_x ) *= -1;	// flip spin to accept config.
		currE += Echange;
		std::cout << currE << std::endl;
	};
	auto reject = [&]() {
//		std::cout << "rejected" << std::endl;
	};
	MetropolisStep met( propose, change, accept, reject, &rndGen );

	auto step = [&](){ met.step(); };
	auto onConfig = [&](){
		// average spin on lattice
//		int averageSpinOnConf = 0.;
//		for(size_t x = 0; x < lat.getVol(); x++ )
//			averageSpinOnConf += spin(x);
//
//		averageSpin += averageSpinOnConf/double(lat.getVol());
//		avSpinAbs += abs( averageSpinOnConf )/double(lat.getVol());
//		avSpinOnConfig << averageSpinOnConf/double(lat.getVol()) << std::endl;
	};
	ConfigGenerator confGen( numThermal, numConfs, numUpPerConf, step, onConfig );
	confGen.run();

	return 0;
}
