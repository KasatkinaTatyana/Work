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
    Bracket(std::vector<double>& gains, std::vector<Power_t>& powers);
    Bracket(void);

    virtual ~Bracket();
    void BracketInit(unsigned N);
    void BracketInit(std::vector<double>& gains, std::vector<Power_t>& powers);
    void BracketCleanUp();
    unsigned BracketSize() {return m_N;}

    std::vector<double>& GetGains() {return m_Gains;}
    void SetGains(std::vector<double>& gains);

    std::vector<Power_t>& GetPowers() {return m_Powers;}
    void SetPowers(std::vector<Power_t>& powers);

    void ShowElements();

    Bracket operator*(Bracket& b);

    Bracket operator*(double number);

    Bracket operator*=(Bracket& b);

    Bracket operator+(Bracket& b);

    //Bracket operator+=(Bracket& b);

    Bracket operator-(Bracket& b);

private:
    std::vector<double> m_Gains;
    std::vector<Power_t> m_Powers;
    unsigned m_N;
};

#endif // BRACKET_H
