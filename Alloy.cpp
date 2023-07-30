#include "Alloy.h"
double two_alloy_bowing_abc(double x, double all_ac, double all_bc, double bow_abc) //second order bowing formula for trinary alloy
{
	return all_bc * x + all_ac * (1 - x) + x * (1 - x) * bow_abc;
}

double four_alloy_bow_abcd(double x, double y, double xG_abc, double yG_acd, double yG_bcd, double xG_abd) // bowing formula for quadrary alloy of type 
{

	if (x*(1-x) == 0) //case when x part in denominator is zero and y parts cancel out.
	{
		return (x * yG_acd + (1 - x) * yG_bcd);
	}
	else if (y * (1 - y) == 0) //case when y part in denominator is zero and x parts cancel out.
	{
		return ((1 - y) * xG_abd + y * xG_abc) ;
	}
	else //full formula when nothing cancels out. Above cases implemented to avoid divading by zero
	{
		return (x * (1 - x) * ((1 - y) * xG_abd + y * xG_abc) + y * (1 - y) * (x * yG_acd + (1 - x) * yG_bcd)) / (x * (1 - x) + y * (1 - y));
	}
	
}

double Varschni(double E0, double alpha, double beta, double  T) // Varschni equation
{
	return E0 - alpha * T * T / (T + beta);
}

double Ev_EpsilonII(double apod, double a0) //equation for stress variable
{
	return (apod - a0) / a0;
}

double linear_a(double a, double b, double T) // linear dependence of lattice constant from temperature
{
	return b + a * (T - 300.0);
}


Alloy::Alloy()
{
	a = 0;
	base_a = 0;
	Eg = 0;
	VBO = 0;
	av = 0;
	ac = 0;
	b = 0;
	C11 = 0;
	C12 = 0;
	E0 = 0;
	alpha = 0;
	beta = 0;
	a_lin_a = 0;
	a_lin_b = 0;
	Temp = 300;
	EpsilonII = 0;
	gamma1 = 0;
	gamma2 = 0;
	gamma3 = 0;
}

Alloy::Alloy(double _a , double _base_a, double _Eg , double _VBO, double _av, double _ac, double _b,double _C11, double _C12)
{
	a = _a;
	base_a = _base_a;
	Eg = _Eg;
	VBO = _VBO;
	av = _av;
	ac = _ac;
	b = _b;
	C11 = _C11;
	C12 = _C12;
	E0 = Eg;
	alpha = 0;
	beta = 0;
	a_lin_a = 0;
	a_lin_b = 0;
	Temp = 300;
	EpsilonII = Ev_EpsilonII(base_a, a);
	gamma1 = 0;
	gamma2 = 0;
	gamma3 = 0;
}

Alloy::Alloy(double _a, double _base_a, double _Eg, double _VBO, double _av, double _ac, double _b, double _C11, double _C12,double _E0, double _alpha,double _beta, double _a_lin_a, double _a_lin_b)
{
	a = _a;
	base_a = _base_a;
	Eg = _Eg;
	VBO = _VBO;
	av = _av;
	ac = _ac;
	b = _b;
	C11 = _C11;
	C12 = _C12;
	E0 = _E0;
	alpha = _alpha;
	beta = _beta;
	Temp = 300;
	EpsilonII = Ev_EpsilonII(base_a, a);
	a_lin_a = _a_lin_a;
	a_lin_b = _a_lin_b;
	gamma1 = 0;
	gamma2 = 0;
	gamma3 = 0;
}

Alloy::Alloy(double _a, double _base_a, double _Eg, double _VBO, double _av, double _ac, double _b, double _C11, double _C12, double _E0, double _alpha, double _beta, double _a_lin_a, double _a_lin_b, double _gamma1, double _gamma2, double _gamma3)
{
	a = _a;
	base_a = _base_a;
	Eg = _Eg;
	VBO = _VBO;
	av = _av;
	ac = _ac;
	b = _b;
	C11 = _C11;
	C12 = _C12;
	E0 = _E0;
	alpha = _alpha;
	beta = _beta;
	Temp = 300;
	EpsilonII = Ev_EpsilonII(base_a, a);
	a_lin_a = _a_lin_a;
	a_lin_b = _a_lin_b;
	gamma1 = _gamma1;
	gamma2 = _gamma2;
	gamma3 = _gamma3;
}

Alloy::Alloy(const Alloy &A)
{
	a = A.a;
	base_a = A.base_a;
	Eg = A.Eg;
	VBO = A.VBO;
	av = A.av;
	ac = A.ac;
	b = A.b;
	C11 = A.C11;
	C12 = A.C12;
	E0 = A.E0;
	alpha = A.alpha;
	beta = A.beta;
	a_lin_a = A.a_lin_a;
	a_lin_b = A.a_lin_b;
	Temp = A.Temp;
	EpsilonII = A.EpsilonII;
	gamma1 = A.gamma1;
	gamma2 = A.gamma2;
	gamma3 = A.gamma3;
}

void Alloy::Mod_base(double newbase)
{
	base_a = newbase;
}

void Alloy::UpdateTemp(double T, bool final)
{
	if (T != Temp)
	{
		a = linear_a(a_lin_a, a_lin_b, T);
		if (final) EpsilonII = Ev_EpsilonII(base_a, a);
		Eg = Varschni(E0, alpha, beta, T);
		Temp = T;
	}
	
}

double Alloy::get_a()
{
	return a;
}

Alloy::~Alloy()
{

}

TwoAlloy::TwoAlloy()
{
	AC = Alloy();
	BC = Alloy();
	base_a = 0;
	bowing_a = 0;
	bowing_Eg = 0;
	bowing_VBO = 0;
	bowing_av = 0;
	bowing_ac = 0;
	bowing_b = 0;
	bowing_C11 = 0;
	bowing_C12 = 0;
	Temp = 300;
	EpsilonII = 0;
}

TwoAlloy::TwoAlloy(const TwoAlloy& C)
{
	AC = C.AC;
	BC = C.BC;
	base_a = C.base_a;
	bowing_a = C.bowing_a;
	bowing_Eg = C.bowing_Eg;
	bowing_VBO = C.bowing_VBO;
	bowing_av = C.bowing_av;
	bowing_ac = C.bowing_ac;
	bowing_b = C.bowing_b;
	bowing_C11 = C.bowing_C11;
	bowing_C12 = C.bowing_C12;
	Temp = C.Temp;
	EpsilonII = C.EpsilonII;
}

TwoAlloy::TwoAlloy(Alloy AC_, Alloy BC_,double _base_a, double b_a, double b_Eg, double b_VBO, double b_av, double b_ac, double b_b, double b_C11, double b_C12)
{
	AC = AC_;
	BC = BC_;
	base_a = _base_a;
	bowing_a = b_a;
	bowing_Eg = b_Eg;
	bowing_VBO = b_VBO;
	bowing_av = b_av;
	bowing_ac = b_ac;
	bowing_b = b_b;
	bowing_C11 = b_C11;
	bowing_C12 = b_C12;
	Temp = 300;
	EpsilonII = 0;
}

void TwoAlloy::Mod_base(double newbase)
{
	base_a = newbase;
}

void TwoAlloy::UpdateTemp(double T)
{
	AC.UpdateTemp(T, false);
	BC.UpdateTemp(T, false);
}

void TwoAlloy::UpdateEpsilon(double new_a)
{
	EpsilonII = Ev_EpsilonII(base_a, new_a);
}

double TwoAlloy::alloy_var_bow(double x, int var)
{
	double result;
	switch (var)//calculate given variable using bowing formula
	{
	case 0:
		result = two_alloy_bowing_abc(x,AC.a,BC.a,bowing_a);
		break;
	case 1:
		result = two_alloy_bowing_abc(x, AC.Eg, BC.Eg, bowing_Eg);
		break;
	case 2:
		result = two_alloy_bowing_abc(x, AC.VBO, BC.VBO, bowing_VBO);
		break;
	case 3:
		result = two_alloy_bowing_abc(x, AC.av, BC.av, bowing_av);
		break;
	case 4:
		result = two_alloy_bowing_abc(x, AC.ac, BC.ac, bowing_ac);
		break;
	case 5:
		result = two_alloy_bowing_abc(x, AC.b, BC.b, bowing_b);
		break;
	case 6:
		result = two_alloy_bowing_abc(x, AC.C11, BC.C11, bowing_C11);
		break;
	case 7:
		result = two_alloy_bowing_abc(x, AC.C12, BC.C12, bowing_C12);
		break;
	case 8:
		result = two_alloy_bowing_abc(x, AC.gamma1, BC.gamma1, 0);
		break;
	case 9:
		result = two_alloy_bowing_abc(x, AC.gamma2, BC.gamma2, 0);
		break;
	case 10:
		result = two_alloy_bowing_abc(x, AC.gamma3, BC.gamma3, 0);
		break;
	default:
		result =  0;
		break;
	}
	return result;
}

double TwoAlloy::alloy_var_bow(double x, int var, double T)
{
	UpdateTemp(T);
	return alloy_var_bow(x, var);
}

std::vector<double> TwoAlloy::calc_bands(double x)
{
	double a0 = alloy_var_bow(x, 0);//calculate all needed variables with bowing
	UpdateEpsilon(a0);
	double Eg = alloy_var_bow(x, 1);
	double VBO = alloy_var_bow(x, 2);
	double av = alloy_var_bow(x, 3);
	double ac = alloy_var_bow(x, 4);
	double b = alloy_var_bow(x, 5);
	double C11 = alloy_var_bow(x, 6);
	double C12 = alloy_var_bow(x, 7);
	double deltaE__ck = 2 * ac * (1 - C12 / C11) * EpsilonII; //calculate deltas
	double deltaE__vh = 2 * av * (1 - C12 / C11) * EpsilonII;
	double deltaE__s = b * (1 - 2 * C12 / C11) * EpsilonII;
	double Ec = VBO + Eg + deltaE__ck; //calculate condutction band and ligh and heavy holes valence bands
	double Ehh = VBO + deltaE__vh + deltaE__s;
	double Elh = VBO + deltaE__vh - deltaE__s;
	std::vector<double> result;
	result.push_back(Ec); result.push_back(Ehh); result.push_back(Elh);
	return result;
	//the formulas used to calculates bands are probably avaible in common sources, however due to the fact that those formulas were just given to us during project classes, I don't have source for it at the moment
}

std::vector<double> TwoAlloy::calc_bands(double x, double temp)
{
	UpdateTemp(temp);
	return calc_bands(x);
}

TwoAlloy::~TwoAlloy()
{
	AC.~Alloy();
	BC.~Alloy();
}

TriAlloy::TriAlloy()
{
	ABC = TwoAlloy();
	ABD = TwoAlloy();
	ACD = TwoAlloy();
	BCD = TwoAlloy();
	base_a = 0;
	Temp = 0;
	EpsilonII = 0;
}

TriAlloy::TriAlloy(const TriAlloy& ABCD)
{
	ABC = ABCD.ABC;
	ABD = ABCD.ABD;
	ACD = ABCD.ACD;
	BCD = ABCD.BCD;
	base_a = ABCD.base_a;
	Temp = ABCD.Temp;
	EpsilonII = ABCD.EpsilonII;
}

TriAlloy::TriAlloy(TwoAlloy ABC_, TwoAlloy ABD_, TwoAlloy ACD_, TwoAlloy BCD_, double _base_a)
{
	ABC = ABC_;
	ABD = ABD_;
	ACD = ACD_;
	BCD = BCD_;
	base_a = _base_a;
	Temp = 0;
	EpsilonII = 0;
}

TriAlloy::~TriAlloy()
{
	ABC.~TwoAlloy();
	ABD.~TwoAlloy();
	ACD.~TwoAlloy();
	BCD.~TwoAlloy();
}

double TriAlloy::alloy_var_bow(double x,double y, int var)
{
	double result;
	result = four_alloy_bow_abcd(x, y, ABC.alloy_var_bow(x, var), ACD.alloy_var_bow(y, var),  BCD.alloy_var_bow(y, var), ABD.alloy_var_bow(x, var));//the list of varaibles is take care of in trinary function, so here there is only one function
	return result;
}

double TriAlloy::alloy_var_bow(double x, double y, int var,double T)
{
	double result;
	result = four_alloy_bow_abcd(
		x, y, ABC.alloy_var_bow(x, var, T), ACD.alloy_var_bow(y, var, T),
		BCD.alloy_var_bow(y, var, T), ABD.alloy_var_bow(x, var, T)
	);
	return result;
}

void TriAlloy::Mod_base(double newbase)
{
	base_a = newbase;
}

void TriAlloy::UpdateTemp(double T)
{
	ABC.UpdateTemp(T);
	ABD.UpdateTemp(T);
	ACD.UpdateTemp(T);
	BCD.UpdateTemp(T);
}
