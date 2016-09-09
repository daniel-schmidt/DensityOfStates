//============================================================================
// Name        : WangLandauIsing.cpp
// Author      : Daniel Schmidt
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <random>
#include "Lattice.h"
#include "Field.h"
#include "IsingHamiltonian.h"
#include "MetropolisWangLandauStep.h"
#include "ConfigGenerator.h"

int main() {
	using namespace FermiOwn;

	size_t Nt = 4;
	size_t Ns = 4;
	size_t dim = 2;
	size_t dofPerPoint = 1;
	double J = 0.3;

	double f = 2.;
	double flatness = 0.90;
	double finalTol = 1e-10;

	size_t numUpdates = 10000000;
	size_t histCheckEvery = 1000;
	Lattice lat( Nt, Ns, dim );

	std::ranlux48 rndGen;
	std::uniform_int_distribution<int> spin_dist(0,1);

	Field<int> spin( lat.getVol(), dofPerPoint, &rndGen, oneInit);

	IsingHamiltonian H( lat, spin, J );

	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( spin_dist(rndGen) == 1 )
			spin(x) -= 2;
	}

	spin.Print();

	MetropolisWangLandauStep wlStep( spin, H, f, finalTol, flatness, histCheckEvery, &rndGen );

	size_t numThermal = 1;
	size_t numUpPerConf = 1;

	auto measure = [](){};
	ConfigGenerator confGen( numThermal, numUpdates, numUpPerConf, &wlStep , measure );

	confGen.run();

	std::cout << std::endl;
	std::cout << "Density of states: " << std::endl;
	std::cerr.precision( 17 );
	auto dos = wlStep.getDos();
	for( auto dosPair : dos ) {
		std::cerr <<  dosPair.first << "\t" << dosPair.second << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Histogram:" << std::endl;
	for( auto histPair : wlStep.getHist() ) {
		std::cout << histPair.first << "\t" << histPair.second << std::endl;
	}

	std::cout << "Final f: " << wlStep.getF() << std::endl;

	return 0;
}
