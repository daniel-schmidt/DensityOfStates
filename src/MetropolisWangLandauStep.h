/*
 * MetropolisWangLandauStep.h
 *
 *  Created on: 09.09.2016
 *      Author: dschmidt
 */

#ifndef SRC_METROPOLISWANGLANDAUSTEP_H_
#define SRC_METROPOLISWANGLANDAUSTEP_H_

#include "Field.h"
#include "IsingHamiltonian.h"
#include "MetropolisStep.h"

namespace FermiOwn {

class MetropolisWangLandauStep: public MetropolisStep {
public:
	MetropolisWangLandauStep( Field<int>& spinField, const IsingHamiltonian& hamilton, double initialF, const double finalFTol, const double flatness, const size_t checkHistEvery, std::ranlux48* rndGen );
	virtual ~MetropolisWangLandauStep();

	virtual bool onConfig( size_t confNum );

	inline std::map< double, size_t > & getHist();
	inline std::map< double, double > & getDos();

protected:
	virtual void propose();
	virtual double change();
	virtual void accept();
	virtual void reject();

	double getDos( double E ) const;
	bool histFlat();

	Field<int>& spin;
	const IsingHamiltonian& H;
	std::uniform_int_distribution<int> x_dist;
	size_t new_x;
	double newE;
	double oldE;

	std::map<double,double> dos;
	std::map<double, size_t> hist;

	const double flat;
	const size_t checkHist;
	double f;
	double Ftol;
};

std::map< double, size_t > & MetropolisWangLandauStep::getHist() {
	return hist;
}

std::map< double, double > & MetropolisWangLandauStep::getDos() {
	return dos;
}

} /* namespace FermiOwn */

#endif /* SRC_METROPOLISWANGLANDAUSTEP_H_ */
