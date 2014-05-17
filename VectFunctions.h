#ifndef VECTFUNCTIONS_H
#define VECTFUNCTIONS_H

#include "bracket.h"

#include <vector>
#include <cstdlib>
#include <ctime>

std::vector<double> VectProduct(std::vector<double> a, std::vector<double> b);
void VectBracketProduct(std::vector<Bracket>& a, std::vector<Bracket>& b, std::vector<Bracket>& result);
double ScalarProduct(std::vector<double> a, std::vector<double> b);
Bracket GeneralVectorTensorVectorProduct(std::vector<Bracket>& vect1, std::vector<Bracket>& vect2, double** M);

void MultNumber(std::vector<double>& vect, double number);
void SumVector(std::vector<double>& sum,const std::vector<double> added);

double NumericalVectorTensorVectorProduct(std::vector<double>& vect1,
                                          std::vector<double>& vect2, double** M);

double NumericalVectorTensorVectorProduct(double* arr1,
                                          double* arr2, double** M);
#endif // VECTFUNCTIONS_H
