#ifndef BRACKETFUNCTIONS_H
#define BRACKETFUNCTIONS_H

#include "bracket.h"

#include <math.h>

#include <vector>
#include <cstdlib>
#include <ctime>


bool comparefun(GainPower_t x, GainPower_t y);

bool cond(GainPower_t x);

bool equality(GainPower_t x, GainPower_t y);

void VectBracketValue (std::vector<Bracket>& br, std::vector<double>& vect, double ksi1, double ksi2, double ksi3);
void LocalTermsChange(GainPower_t& local_terms, unsigned ind,unsigned value);
void SimplifyBracket(Bracket& br);

#endif // BRACKETFUNCTIONS_H