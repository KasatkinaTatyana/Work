#ifndef BRACKET_H
#define BRACKET_H

#include <vector>


struct _Power {
    int p1;
    int p2;
    int p3;
    int p4;
};

typedef struct _Power Power_t;

class Bracket
{
public:
    Bracket(unsigned N = 0);
    Bracket(std::vector<double> gains, std::vector<Power_t> powers);
    virtual ~Bracket();
    void BracketInit(unsigned N);
    void BracketInit(std::vector<double> gains, std::vector<Power_t> powers);
    void BracketCleanUp();
    unsigned BracketSize() {return m_N;}

    std::vector<double> GetGains() {return m_Gains;}
    void SetGains(std::vector<double> gains) {this->m_Gains=gains;}

    //std::vector<Power_t> GetGains() {return m_Powers;}    ????????
    void SetPowers(std::vector<Power_t> powers){this->m_Powers=powers;}


private:
    std::vector<double> m_Gains;
    std::vector<Power_t> m_Powers;
    unsigned m_N;

    Bracket operator*(Bracket& b)
    {
        unsigned M = b.BracketSize()*BracketSize();
        unsigned L = b.BracketSize();
        unsigned N = BracketSize();

        Bracket result(M);

        for(unsigned i = 0; i < N; i++)
        {

            unsigned offset = i*L;

            float Ai = m_Gains.at(i);
            Power_t q1 = m_Powers.at(i);

            for(unsigned j = 0; j < L; j++)
             {

                float Bj = b.m_Gains.at(j);

                float Cij = Ai * Bj;

                result.m_Gains[offset + j] = Cij;

                Power_t q2 = m_Powers.at(j);

                Power_t q3;

                q3.p1 = q1.p1+q2.p1;
                q3.p2 = q1.p2+q2.p2;
                q3.p3 = q1.p3+q2.p3;
                q3.p4 = q1.p4+q2.p4;

                result.m_Powers[offset+j]=q3;

            }
        }
        return result;
    }
};

#endif // BRACKET_H
