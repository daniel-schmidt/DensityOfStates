//============================================================================
// Name        : WangLandauIsing.cpp
// Author      : Daniel Schmidt
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <random>
#include <map>
#include "Lattice.h"
#include "Field.h"
#include "IsingHamiltonian.h"
#include "MetropolisStep.h"
#include "ConfigGenerator.h"

namespace FermiOwn {

double getDos( const std::map<double,double>& dos, const double energy ) {
	auto it = dos.find( energy );
	if( it == dos.end() ) {
		return 0;
	} else {
		return it->second;
	}
}

bool histFlat( const std::map<double, size_t>& hist, const double flatness ) {
	double meanHist = 0.;
	double minHist = hist.begin()->second;
	for( auto histPair : hist ) {
		meanHist += histPair.second;
		if( histPair.second < minHist ) {
			minHist = histPair.second;
		}
	}
	meanHist /= hist.size();
	std::cout << "meanHist: " << meanHist << " minHist: " << minHist << " > " << meanHist*flatness << "?";
	if( meanHist*flatness < minHist ) {
		std::cout << " Histogram flat!" << std::endl;
		return true;
	} else {
		std::cout << " Not flat!" << std::endl;
		return false;
	}
}

}

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
	std::uniform_int_distribution<int> x_dist( 0, lat.getVol()-1 );
	std::uniform_real_distribution<double> real_dist;

	Field<int> spin( lat.getVol(), dofPerPoint, &rndGen, oneInit);

	IsingHamiltonian H( lat, spin, J );

	std::map<double, double> dos;
	std::map<double, size_t> hist;

	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( spin_dist(rndGen) == 1 )
			spin(x) -= 2;
	}

	spin.Print();

	double E = -H.calculateEnergy();

	double newE;
	size_t x;

	auto propose = [&]() {
		x = x_dist( rndGen ); //change spin at a random point
		spin( x ) *= -1;
		newE = -H.calculateEnergy();
	};

	auto change = [&]() {
		return std::exp( getDos(dos, E) - getDos( dos, newE ) );
	};

	auto accept = [&]() {
		E=newE;
	};

	auto reject = [&]() {
		spin( x ) *= -1;
	};

	auto onConfig = [&](int confNum) {
		dos[E] += f;
		hist[E]++;

		if( (confNum+1)%histCheckEvery == 0 && histFlat( hist, flatness ) ) {
			f *= 0.5;
			if( f < finalTol ) return false;
			hist.clear();
		}
		return true;
	};

	MetropolisStep metStep( propose, change, accept, reject, &rndGen );
	auto step = [&]() {
		metStep.step();
	};
	size_t numThermal = 1;
	size_t numUpPerConf = 1;

	ConfigGenerator confGen( numThermal, numUpdates, numUpPerConf, step , onConfig );

	confGen.run();

	std::cout << std::endl;
	std::cout << "Density of states: " << std::endl;
	std::cerr.precision( 17 );
	for( auto dosPair : dos ) {
		std::cerr <<  dosPair.first << "\t" << dosPair.second << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Histogram:" << std::endl;
	for( auto histPair : hist ) {
		std::cout << histPair.first << "\t" << histPair.second << std::endl;
	}

	std::cout << "Final f: " << f << std::endl;

	return 0;
}
