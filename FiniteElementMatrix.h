#ifndef FINITEELEMENTMATRIX_H
#define FINITEELEMENTMATRIX_H

#include "bracket.h"

#include <vector>
#include <list>
#include <stdlib.h>

class FiniteElementMatrix
{
public:
    static const unsigned m_CountPeaks=4;
    static const unsigned m_Dim=3;
    static const unsigned m_MaxOrderFact=25;

    FiniteElementMatrix(unsigned p, double** simplex_peaks, double** Eps, double** Mu);
    virtual ~FiniteElementMatrix();

	void FormArrayAnalyt();     //����������� ������ ������������ �������� ����������� �������

	void FormArrayNum();        //���������� ������ �������� �������� ����������� ������� �� ���� �����

    void MatrixInit();        //������������ ����������� ������� � ������� ������ � ������������� ������

	void NumMatrixInit();    //������������ ����������� ������� � ������� ������ � �������������� ���������� ��������������

	// ���������� ������
    void AddToVectBracket(std::vector<GainPower_t>& terms, std::vector<Bracket>& vect_bracket);

    void ShowVectBracket(std::vector<Bracket>& vect_bracket);

    std::vector<Bracket> FormVectEigFunc(Bracket& br, unsigned n, unsigned m);

	std::vector<Bracket> RotorCalc(Bracket& br, unsigned n, unsigned m);

	void NumFormVectEigFunc(std::vector<Bracket>& vect_br, unsigned n, unsigned m, 
							double ksi1, double ksi2, double ksi3, 
							std::vector<double>& result);

	double Integrate(Bracket* br);

	// ����������� �����������
	void ShowMatrixes();

	void ExportFuncValues(std::string s1, std::string s2, std::string s3);

	//----------&----------&----------&----------&----------&----------&----------&----------&
	// ���� �������, ������� �����, ����� ��������� ����� � �������������� �������� ����������
	void FormArrayAnalyt_LinComb();
	void FormArrayNum_LinComb();
	//----------&----------&----------&----------&----------&----------&----------&----------&
	void FormEdjeFunc(Bracket* bracket,unsigned* index_array,
					  std::vector<Bracket>* eig_func,
					  std::vector<Bracket>* rot_eigfunc);

	void FormFaceFunc(Bracket* bracket,unsigned* index_array,
					  unsigned* non_zero_array, std::vector<Bracket>* eig_func,
					  std::vector<Bracket>* rot_eigfunc, unsigned ind);

	void FormInsideFunc(Bracket* bracket,unsigned* node_array,
					    std::vector<Bracket>* eig_func, std::vector<Bracket>* rot_eigfunc, unsigned ind);

	void CalcUnitRot(Bracket* bracket, unsigned numb, std::vector<Bracket>* rot);

private:
	    unsigned m_QuadOrder;            //������� ������������ ������
		int Q2;                          // m_QuadOrder ^ 2
		int Q3;                          // m_QuadOrder ^ 3 

        unsigned m_P;

        std::vector<double> DefVector(unsigned ind);

		double** m_Peaks;       //������ ��������� (4) ������ ��������� 4(=m_CountPeaks) �� 3(=m_Dim)

        double** m_MatrixMu;    //������� ������� ��������� ������������� 3 �� 3
        double** m_MatrixEps;    //������� ������� ��������������� �������������

        void AddSilvester(unsigned gamma, unsigned beta, unsigned numb, unsigned ind, std::vector<Bracket>& vect_bracket);


        void findIndex(unsigned gamma, unsigned beta, unsigned order, unsigned& idxLow, unsigned& idxHigh);

        double *m_MetrMatrix;                    //����������� �������
        double *m_EulerMatrix;                   //������� ������
		double *m_NumEulerMatrix;                //������� ������ ��� ��������� �������������� 
		double *m_NumMetrMatrix;                //����������� ������� ��� ��������� �������������� 
        int m_MatrixSize;    //����������� ����������� ������� � ������� ������

        //���������� �����������--------------------------------------------------------------------------------------------
        double Fact(unsigned N);

        double m_ArrayFact[m_MaxOrderFact+1];

        inline double CalcFact(unsigned N);
        //------------------------------------------------------------------------------------------------------------------
        //��������� ���������������� ���������
        std::vector<double> m_GradKsi_1;
        std::vector<double> m_GradKsi_2;
        std::vector<double> m_GradKsi_3;
        std::vector<double> m_GradKsi_4;

        //�������� ��������� ��������, ������������ ����� ����� ���������
        std::vector<double> m_UnVect_1;
        std::vector<double> m_UnVect_2;
        std::vector<double> m_UnVect_3;

        //������� ������ ���������
        double m_Icob;
		//-----------------------------------------------------------------------------------------------------------------
		//������� ������ ��������� �������� � ����� ��� ���������� ��������������
		std::vector<double> m_Weights;
		std::vector<double> m_Roots;

		void CompareMatrixs();
		//-------------------------------------------------------------------------------------------------------------
		double** m_Arr_AllNodes;     //������ ���� ����������� ������� � �� ������� �� ���� ����� 
		double** m_Arr_RotAllNodes;

		bpp_t m_ArrAnalyt_EigFunc;   //������ ������������� ����������� ������� � �� �������
		bpp_t m_ArrAnalyt_RotEigFunc;

		bp_t* m_ArrAnalyt_Shapes;
		bp_t* m_ArrAnalyt_CurlShapes;

		void display(unsigned rows, unsigned columns, double** arr);
};



#endif // FINITEELEMENTMATRIX_H
