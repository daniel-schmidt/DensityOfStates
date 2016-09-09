/*
 * TrunkatedMetropolis_test.cpp
 *
 *  Created on: 07.09.2016
 *      Author: dschmidt
 */

#include "Lattice.h"
#include "IsingHamiltonian.h"
#include "MetropolisWangLandauStep.h"
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
	double a = 0.8;

	size_t numThermal = 1;
	size_t numConfs = 1000;
	size_t numUpPerConf = lat.getVol();
	//	size_t numUpPerConf = 1;

	std::ranlux48 rndGen;

	std::uniform_int_distribution<int> spin_dist(0,1);
	std::uniform_int_distribution<int> x_dist(0, lat.getVol()-1 );

	int new_x = 0;

	// initialization of physical parameters

	Field<int> spin( lat.getVol(), dofPerPoint, &rndGen, oneInit);

	for( size_t x = 0; x < 0.5*lat.getVol(); x+=2 ) {
//		if( spin_dist( rndGen) ) {
			spin( x ) *= -1;
//		}
	}
	IsingHamiltonian H( lat, spin, J );

	double E0 = -H.calculateEnergy();
	double currE = E0;
	double dE = atof( argv[2] );

	std::cout << "#Running algorithm for E0 = " << E0 << " and dE = " << dE << std::endl;
//	double Echange = 0.;
//	// metropolis initialization
//	auto propose = [&] () {
//		//		std::cout << "Propose" << std::endl;
//		do {
//			new_x = x_dist( rndGen );
//			Echange = beta*H.changeAt( new_x );
//			//			std::cout << " change in E: " << Echange << std::endl;
//		} while( currE+Echange > E0+dE/2. || currE+Echange < E0-dE/2. );
//	};
//	auto change = [&]() {
//		return std::exp( a*Echange );
//	};
//	auto accept = [&]() {
//		spin( new_x ) *= -1;	// flip spin to accept config.
//		currE += Echange;
//	};
//	auto reject = [&]() {
//		//		std::cout << "rejected" << std::endl;
//	};

	double f =2.0;
	double fFinal = 1e-8;
	double flatness = 0.8;
	size_t checkHist = 100;
	MetropolisWangLandauStep met( spin, H, f, fFinal, flatness, checkHist, &rndGen );

	double EdiffFromE0 = 0.;
//	auto step = [&](){ met.step(); };
	auto measure = [&]( ){
		EdiffFromE0 += currE-E0;
		//		std::cout << currE-E0 << std::endl;
		// average spin on lattice
		//		int averageSpinOnConf = 0.;
		//		for(size_t x = 0; x < lat.getVol(); x++ )
		//			averageSpinOnConf += spin(x);
		//
		//		averageSpin += averageSpinOnConf/double(lat.getVol());
		//		avSpinAbs += abs( averageSpinOnConf )/double(lat.getVol());
		//		avSpinOnConfig << averageSpinOnConf/double(lat.getVol()) << std::endl;
		return true;
	};

	ConfigGenerator confGen( numThermal, numConfs, numUpPerConf, &met, measure );
	for( size_t n = 0; n < 500; n++ ){
		confGen.run();

		a += 12./(4.*dE + dE*dE)*EdiffFromE0/numConfs;
//		a += 12./(dE*dE)*EdiffFromE0/numConfs;

		std::cout << EdiffFromE0/numConfs << "\t" << a << std::endl;
	}

	return 0;
}
