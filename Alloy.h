#pragma once
#include <string>
#include <vector>


class Alloy //alloy class
{

	public:
		Alloy();
		Alloy(double _a, double _base_a, double _Eg, double _VBO, double _av, double _ac, double _b, double _C11, double _C12);//constructor when tempretaure doesn't change
		Alloy(double _a, double _base_a, double _Eg, double _VBO, double _av, double _ac, double _b, double _C11, double _C12, double _E0, double _alpha, double _beta, double _a_lin_a, double _a_lin_b); //constructor when temperature changes
		Alloy(double _a, double _base_a, double _Eg, double _VBO, double _av, double _ac, double _b, double _C11, double _C12, double _E0, double _alpha, double _beta, double _a_lin_a, double _a_lin_b, double _gamma1, double _gamma2, double _gamma3); //extended constructor using constants used for effective mass aproximation - not used 
		Alloy(const Alloy&);
		void Mod_base(double); //substrate lattice constant modification
		void UpdateTemp(double,bool);  // updating all values dependent on temperature, when temperature changes
		double get_a(); // gives lattice constant
		~Alloy();
		friend class TwoAlloy;

	private:
		double a; //lattice constatnt
		double base_a;//substrate lattice constant
		double Eg; // band gap
		double VBO; //experimentally aquired constants used for bands calculation
		double av;
		double ac;
		double b;
		double C11;
		double C12;
		double E0; // band gap in 0K - used for Varschni equation
		double alpha; // alpha in Varschni equation
		double beta; // beta in Varschni equation
		double Temp; 
		double EpsilonII; //stress variable
		double a_lin_a; //linear coeficient in linear dependence of lattice constant from temperature
		double a_lin_b; //constant coeficient in linear dependence of lattice constant from temperature
		double gamma1; //coefficient for effective mass aproximation - not used
		double gamma2;
		double gamma3;
		 
};

class TwoAlloy
{
public:
	TwoAlloy();
	TwoAlloy(Alloy AC_, Alloy BC_, double _base_a, double b_a, double b_Eg, double b_VBO, double b_av, double b_ac, double b_b, double b_C11, double b_C12); //only one constructor 
	//bowings don't depent on temperature so there is no nedd to make separate constructors. There is also no effective mass aproximation constructor as it was abandoned
	TwoAlloy(const TwoAlloy&);
	~TwoAlloy();
	void Mod_base(double newbase); //substrate lattice constant modification
	void UpdateTemp(double newtemp); // updating all values dependent on temperature, when temperature changes
	void UpdateEpsilon(double new_a); // updating stress variable when lattice constant changes
	double alloy_var_bow(double, int); //calculate effective variable for trinary alloy using binary alloy variables and bowings values. Goes from 1 to 10. List is in definition.
	double alloy_var_bow(double, int, double); // same as previous but with specified temperature
	std::vector<double> calc_bands(double x); //calculate bands
	std::vector<double> calc_bands(double x,double temp); //same as previous one but with temperature
	friend class TriAlloy;

private:
	Alloy AC; //first alloy
	Alloy BC; //second alloy
	double base_a; //substrate lattice constant
	double bowing_a; //bowing of trinary alloy for all variables 
	double bowing_Eg;
	double bowing_VBO;
	double bowing_av;
	double bowing_ac;
	double bowing_b;
	double bowing_C11;
	double bowing_C12;
	double Temp; //temperature
	double EpsilonII; // stress varaible
};

class TriAlloy
{
public:
	TriAlloy();
	TriAlloy(const TriAlloy&);
	TriAlloy(TwoAlloy ABC_, TwoAlloy ABD_, TwoAlloy ACD_, TwoAlloy BCD_, double _base_a); //main constructor, need all trinary combinations and  substrate lattice constant
	double alloy_var_bow(double x, double y, int var); //same as in trinary aloys , calculate selected variable using appropraite equation
	double alloy_var_bow(double x, double y, int var, double T); // same but with temperature
	void Mod_base(double new_base); //substrate lattice constant modification
	void UpdateTemp(double new_temp); // updating all values dependent on temperature, when temperature changes
	~TriAlloy();

private:
	TwoAlloy ABC; //all trialloy combinations
	TwoAlloy ABD;
	TwoAlloy ACD;
	TwoAlloy BCD;
	double base_a; //substrate lattice constant
	double Temp; //temperature
	double EpsilonII; // stress varaible
};


