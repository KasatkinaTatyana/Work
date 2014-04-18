#ifndef BRACKETFUNCTIONS_H
#define BRACKETFUNCTIONS_H

#include "bracket.h"

#include <math.h>

#include <vector>
#include <cstdlib>
#include <ctime>

struct GainsPowers {
	double g;
    int p1;
    int p2;
    int p3;
    int p4;
	GainsPowers(double g, int p1, int p2, int p3, int p4) :
		g(g),
		p1(p1),
		p2(p2),
		p3(p3),
		p4(p4)
	{
	}
private:
};

//typedef struct _GainsPowers GainsPowers_t;

bool comparefun(GainsPowers x, GainsPowers y);

bool cond(GainsPowers x);

void VectBracketValue (std::vector<Bracket>& br, std::vector<double>& vect, double ksi1, double ksi2, double ksi3);
void LocalPowersChange(Power_t& local_powers, unsigned ind,unsigned value);
void SimplifyBracket(Bracket& br);

#endif // BRACKETFUNCTIONS_H