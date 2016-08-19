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

namespace FermiOwn {

double calculateEnergy( Lattice& lat, Field<int>& spin, double J ) {
	//returns -Hamiltonian
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

	size_t Nt = 16;
	size_t Ns = 16;
	size_t dim = 2;
	size_t dofPerPoint = 1;
	double J = 0.3;

	double f = 1.;
	double flatness = 0.9;

	size_t numUpdates = 10000000;
	size_t histCheckEvery = 1000;
	Lattice lat( Nt, Ns, dim );

	std::ranlux48 rndGen;
	std::uniform_int_distribution<int> spin_dist(0,1);
	std::uniform_int_distribution<int> x_dist( 0, lat.getVol()-1 );
	std::uniform_real_distribution<double> real_dist;

	Field<int> spin( lat.getVol(), dofPerPoint, &rndGen, oneInit);
	std::map<double, double> dos;
	std::map<double, size_t> hist;

	for( size_t x = 0; x < lat.getVol(); x++ ) {
		if( spin_dist(rndGen) == 1 )
			spin(x) -= 2;
	}

	spin.Print();

	double E = calculateEnergy( lat, spin, J );

//	std::cout << "Energy: " << E << " DoS: " << getDos( dos, E );


//	std::cout << " DoS again: " << getDos( dos, E );

	for( size_t i = 0; i < numUpdates; i++ ) {
		size_t x = x_dist( rndGen ); //change spin at a random point
		spin( x ) *= -1;

		double newE = calculateEnergy( lat, spin, J );
//		std::cout << " new Energy: " << newE;

		double prop = std::exp( getDos(dos, E) - getDos( dos, newE ) );
		double r = real_dist( rndGen );
//		std::cout << " r=" << r << " < " << prop;
		if( r < prop ) {
			E=newE;
//			std::cout << " accepted." << std::endl;
		} else {
			// we reset the spin flip
			spin( x ) *= -1;
//			std::cout << " rejected." << std::endl;
		}
		dos[E] += f;
		hist[E]++;

		if( (i+1)%histCheckEvery == 0 && histFlat( hist, flatness ) ) {
			f *= 0.5;
			hist.clear();
//			for( auto histPair : hist ) {
//				histPair.second = 0;
//			}
		}

	}

	std::cout << std::endl;
	std::cout << "Density of states: " << std::endl;
	for( auto dosPair : dos ) {
		std::cerr << dosPair.first << "\t" << dosPair.second << std::endl;
	}

	std::cout << std::endl;
	std::cout << "Histogram:" << std::endl;
	for( auto histPair : hist ) {
		std::cout << histPair.first << "\t" << histPair.second << std::endl;
	}

	std::cout << "Final f: " << f << std::endl;

	return 0;
}
