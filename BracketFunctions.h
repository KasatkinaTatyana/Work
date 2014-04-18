#ifndef BRACKETFUNCTIONS_H
#define BRACKETFUNCTIONS_H

#include "bracket.h"

#include <math.h>

#include <vector>
#include <cstdlib>
#include <ctime>

struct _GainsPowers {
	double g;
    int p1;
    int p2;
    int p3;
    int p4;
};

typedef struct _GainsPowers GainsPowers_t;

bool comparefun(GainsPowers_t x, GainsPowers_t y);

bool cond(GainsPowers_t x);

void VectBracketValue (std::vector<Bracket>& br, std::vector<double>& vect, double ksi1, double ksi2, double ksi3);
void LocalPowersChange(Power_t& local_powers, unsigned ind,unsigned value);
void SimplifyBracket(Bracket& br);

#endif // BRACKETFUNCTIONS_H