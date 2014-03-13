#ifndef FINITEELEMENTMATRIX_H
#define FINITEELEMENTMATRIX_H

#endif // FINITEELEMENTMATRIX_H

#include <vector>
#include <list>
#include "bracket.h"
#include <stdlib.h>
using namespace std;

typedef class Bracket Bracket_t;

class FiniteElementMatrix
{
public:
    static const unsigned m_CountPeaks=4;
    static const unsigned m_Dim=3;
    static const unsigned m_MaxOrderFact=25;

    FiniteElementMatrix(unsigned p, double simplex_peaks[m_CountPeaks][m_Dim], double Eps[m_Dim][m_Dim],
                        double Mu[m_Dim][m_Dim]);
    virtual ~FiniteElementMatrix();
    void MatrixInit();

    std::vector<Bracket>& GetVectBracket() {return m_VectBracket;}

    void AddToVectBracket(std::vector<double>& gains, std::vector<Power_t>& powers, std::vector<Bracket>& vect_bracket);

    void ShowVectBracket(std::vector<Bracket>& vect_bracket);

    void AddRotMatrixProduct(Bracket vect1[m_Dim], Bracket vect2[m_Dim], double M[][m_Dim]);

private:

        unsigned m_P;

        Power_t DefPowers(unsigned n,unsigned m);
        std::vector<Bracket> m_VectBracket;       //��� ������ ������������� ���� ���� ���������
        std::vector<Bracket> m_CurVectBracket;    //��� ������ ������������ �������� �������� � ����� ���������
                                                  //��� ������� �� ���� ��������� ���������

        void AddVTVProduct(std::vector<double> a, std::vector<double> b, std::vector<double> c,
                           std::vector<double> d, double M[][m_Dim], unsigned n1, unsigned m1, unsigned n2, unsigned m2,
                           std::vector<Bracket>& vect_bracket);
        std::vector<double> DefVector(unsigned ind);
        std::vector<unsigned> Def_nm(unsigned gamma, unsigned beta);
        double m_Peaks[m_CountPeaks][m_Dim];       //������ ��������� (4) ������ ���������

        std::vector<double> VectProduct(std::vector<double> a, std::vector<double> b);
        double ScalarProduct(std::vector<double> a, std::vector<double> b);

        double m_MatrixMu[m_Dim][m_Dim];    //������� ������� ��������� �������������
        double m_MatrixEps[m_Dim][m_Dim];    //������� ������� ��������������� �������������

        void AddSilvester(unsigned gamma, unsigned beta, unsigned numb, unsigned ind, std::vector<Bracket>& vect_bracket);

        void LocalPowersChange(Power_t& local_powers, unsigned ind,unsigned value);

        double Integrate (Bracket_t& br);

        double** m_MetrMatrix;                    //����������� �������
        double** m_EulerMatrix;                   //������� ������
        unsigned m_MatrixSize;    //����������� ����������� ������� � ������� ������

        //���������� �����������--------------------------------------------------------------------------------------------
        double Fact(unsigned N);

        double m_ArrayFact[m_MaxOrderFact+1];

        inline double CalcFact(unsigned N);
        //------------------------------------------------------------------------------------------------------------------
};
