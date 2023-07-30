#include <iostream>
#include "Alloy.h"
#include <fstream>

using namespace std;

int main()
{
	string punkt;
	fstream plik;
	cout << "zadanie?" << endl;
	cin >> punkt;
	if (punkt == "1")
	{
		//1.task
		//For the selected ternary alloy, determine the relationship between the energy gap and the lattice constant
		// as a function of the composition. For x = 0.5 plot the temperature dependencies of these parameters. (*) Determine the dependency
		// of the above parameters as a function of the composition for a four-component alloy, the results are presented
		// on the chart (e.g. in the form as it was discussed in class).

		//declaring binary alloys
		Alloy InAs = Alloy(6.0583 + 2.74e-5 * (-300), 6.0583 + 2.74e-5 * (-300), 0.41, -0.59, -1.00, -5.08, -1.8, 832.9, 452.6, 0.417, 2.76e-4, 93, 2.74e-5, 6.0583);
		Alloy GaAs = Alloy(5.65325 + 3.88e-5 * (-300), 5.65325 + 3.88e-5 * (-300),1.519, -0.84, -1.16, -7.17, -2.0, 1221, 566, 1.519, 5.405e-4, 204, 3.88e-5, 5.65325);
		//declaring trinary alloy and it's bowing parameters
		TwoAlloy GaInAs = TwoAlloy(GaAs, InAs, 0, 0, -0.477, 0.38, 0, -2.61, 0, 0, 0);
		cout << "sub point" << endl;
		cin >> punkt;
		if (punkt == "a")
		{
			vector<double> x, a, Eg;
			GaInAs.UpdateTemp(300);
			GaInAs.Mod_base(GaInAs.alloy_var_bow(0, 0));
			for (int i = 0; i <= 1000000; i++) //loop going from 0 to 1 with very big precision for alloy doping
			{
				x.push_back(i * 1e-6);
				a.push_back(GaInAs.alloy_var_bow(x.back(), 0));
				Eg.push_back(GaInAs.alloy_var_bow(x.back(), 1));
			}
			plik.open("zad1/a/result.txt",ios::out);
			if (plik.is_open()) //simple loop for saving vector contents in file
			{
				cout << "saving to file..." << endl;
				for (int i = 0; i < x.size(); i++)
				{
					plik << x.at(i) << " " << a.at(i) << " " << Eg.at(i) << endl;
				}
				plik.close();
				cout << "saving finished" << endl;
			}
			else
			{
				cout << "File not opened" << endl;
			}
		}
		else if (punkt == "b") //this sub point is similiar as last one but here the temperature changes instead of doping
		{
			vector<double> T, a, Eg;
			GaInAs.UpdateTemp(0);
			GaInAs.Mod_base(GaInAs.alloy_var_bow(0.5, 0));
			for (int t = 0; t <= 2000; t++)
			{
				T.push_back(t);
				GaInAs.UpdateTemp(t); //after each change of tempreture, program need to update temperature for correct values of constants
				a.push_back(GaInAs.alloy_var_bow(0.5, 0));
				Eg.push_back(GaInAs.alloy_var_bow(0.5, 1));
			}
			plik.open("zad1/b/result.txt", ios::out);
			if (plik.is_open())//simple loop for saving vector contents in file
			{
				cout << "saving to file..." << endl;
				for (int i = 0; i < T.size(); i++)
				{
					plik << T.at(i) << " " << a.at(i) << " " << Eg.at(i) << endl;
				}
				plik.close();
				cout << "saving finished" << endl;
			}
			else
			{
				cout << "File not opened" << endl;
			}
		}
		else if (punkt == "c")
		{
			//InSb GaSb InAsSb GaAsSb GaInSb
			Alloy InSb = Alloy(6.4794 + 3.48e-5 * (-300), 6.4794 + 3.48e-5 * (-300), 0.235, -0, -0.36, -6.94, -2.0, 8684.7, 373.5, 0.235, 3.2e-4, 170, 3.48e-5, 6.4794);
			Alloy GaSb = Alloy(6.0959 + 4.72e-5 * (-300), 6.0959 + 4.72e-5 * (-300), 0.812, -0.03, -0.8, -7.5, -2.0, 884.2, 402.6, 0.812, 4.17e-4, 140, 4.72e-5, 6.0959);
			//declaring additional binary alloys
			TwoAlloy GaInAs = TwoAlloy(InAs, GaAs, 0, 0, -0.477, 0.38, 0, -2.61, 0, 0, 0);
			TwoAlloy InAsSb = TwoAlloy(InSb, InAs, 0, 0, -0.67, 0, 0, 0, 0, 0, 0);
			TwoAlloy GaAsSb = TwoAlloy(GaSb, GaAs, 0, 0, -1.43, 1.06, 0, 0, 0, 0, 0);
			TwoAlloy GaInSb = TwoAlloy(InSb, GaSb, 0, 0, -0.415, 0, 0, 0, 0, 0, 0);
			//declaring rest of combinations for trinary alloys
			TriAlloy GaInAsSb = TriAlloy(GaInAs, GaInSb, GaAsSb, InAsSb,0);
			//declaring quaternary alloy
			cout << GaInAsSb.alloy_var_bow(0, 0, 1) << " " << GaInAsSb.alloy_var_bow(0, 1, 1) << " " << GaInAsSb.alloy_var_bow(1, 0, 1) << " " << GaInAsSb.alloy_var_bow(1, 1, 1)<<endl;
			//checking if vertex points behave correctly
			vector<double> x, y, a, Eg;
			GaInAsSb.UpdateTemp(0);
			GaInAsSb.Mod_base(GaInAsSb.alloy_var_bow(0,0,0));
			
			for (int i = 0; i < 100000; i++)
			{
				x.push_back(0);
				y.push_back(i * 1e-5);
				a.push_back(GaInAsSb.alloy_var_bow(x.back(),y.back(),0));
				Eg.push_back(GaInAsSb.alloy_var_bow(x.back(), y.back(), 1));
			}
			for (int i = 0; i < 100000; i++)
			{
				x.push_back(i*0.5e-5);
				y.push_back(1.0);
				a.push_back(GaInAsSb.alloy_var_bow(x.back(), y.back(), 0));
				Eg.push_back(GaInAsSb.alloy_var_bow(x.back(), y.back(), 1));
			}
			for (int i = 0; i < 100000; i++)
			{
				x.push_back(0.5);
				y.push_back(1-i * 1e-5);
				a.push_back(GaInAsSb.alloy_var_bow(x.back(), y.back(), 0));
				Eg.push_back(GaInAsSb.alloy_var_bow(x.back(), y.back(), 1));
			}
			for (int i = 0; i <= 100000; i++)
			{
				x.push_back(0.5 - i * 0.5e-5);
				y.push_back(0);
				a.push_back(GaInAsSb.alloy_var_bow(x.back(), y.back(), 0));
				Eg.push_back(GaInAsSb.alloy_var_bow(x.back(), y.back(), 1));
			}
			//running 4 loops changing only one variable at time to allow the data to be graphed in 2d. It follows a rectangle from (x,y)=(0,0) to (0,1), (0.5,1), (0.5,0), back to (0,0)
			plik.open("zad1/c/result1.txt", ios::out);
			if (plik.is_open())//simple loop for saving vector contents in file
			{
				cout << "saving to file..." << endl;
				for (int i = 0; i < x.size(); i++)
				{
					plik << x.at(i) << " " << y.at(i) << " " << a.at(i) << " " << Eg.at(i) << endl;
				}
				plik.close();
				cout << "saving finished" << endl;
			}
			else
			{
				cout << "File not opened" << endl;
			}
		}
	}
	else if (punkt == "2")
	{
		//2.task
		//For the previously selected ternary material on a GaSb substrate, compare the energy relationships of
		// conduction bands and heavy and light hole at point Γ as a function of the composition taking into account
		// stresses and without them.

		Alloy InAs = Alloy(6.0583 + 2.74e-5 * (-300), 6.0583 + 2.74e-5 * (-300), 0.41, -0.59, -1.00, -5.08, -1.8, 832.9, 452.6, 0.417, 2.76e-4, 93, 2.74e-5, 6.0583);
		Alloy GaAs = Alloy(5.65325 + 3.88e-5 * (-300), 5.65325 + 3.88e-5 * (-300), 1.519, -0.84, -1.16, -7.17, -2.0, 1221, 566, 1.519, 5.405e-4, 204, 3.88e-5, 5.65325);
		TwoAlloy GaInAs = TwoAlloy(GaAs, InAs, 0, 0, -0.477, 0.38, 0, -2.61, 0, 0, 0);
		cout << "sub point" << endl;
		cin >> punkt;
		if (punkt == "a")
		{
			Alloy GaSb = Alloy(6.0959+4.72e-5*(0.0-300.0), 0, 0.812, -.03, -2, -7.5, -2.0, 884.2, 402.6);
			GaInAs.Mod_base(GaSb.get_a());
			GaInAs.UpdateTemp(0);			
			//GaInAs.UpdateEpsilon(GaInAs)
			vector<double> x,tempor;
			vector<vector<double>> bands;
			for (int i = 0; i <= 1000000; i++)//calculating bands vaule in gamma for each doping value
			{
				x.push_back(i * 1e-6);
				tempor=GaInAs.calc_bands(x.back());
				bands.push_back(tempor);
				
			}
			plik.open("zad2/a/result.txt", ios::out);
			if (plik.is_open())// saving to file
			{
				cout << "saving to file..." << endl;
				for (int i = 0; i < x.size(); i++)
				{
					plik << x.at(i) << " " << bands.at(i).at(0) << " " << bands.at(i).at(1) << " " << bands.at(i).at(2) << endl;
				}
				plik.close();
				cout << "saving finished" << endl;
			}
			else
			{
				cout << "File not opened" << endl;
			}
		}
		else if(punkt == "b") //this is same as previous one, but I set base's lattice constant same as alloy that calculations run for, to set stress to 0
		{
			//Alloy GaSb = Alloy(6.0959 + 4.72e-5 * (0.0 - 300.0), 0, 0.812, -.03, -2, -7.5, -2.0, 884.2, 402.6);
			GaInAs.Mod_base(GaInAs.alloy_var_bow(0, 0));
			GaInAs.UpdateTemp(0);
			//GaInAs.UpdateEpsilon(GaInAs)
			vector<double> x, tempor;
			vector<vector<double>> bands;
			for (int i = 0; i <= 1000000; i++)
			{
				x.push_back(i * 1e-6);
				GaInAs.Mod_base(GaInAs.alloy_var_bow(x.back(), 0));
				tempor = GaInAs.calc_bands(x.back());
				bands.push_back(tempor);

			}
			plik.open("zad2/b/result.txt", ios::out);
			if (plik.is_open())
			{
				cout << "saving to file..." << endl;
				for (int i = 0; i < x.size(); i++)
				{
					plik << x.at(i) << " " << bands.at(i).at(0) << " " << bands.at(i).at(1) << " " << bands.at(i).at(2) << endl;
				}
				plik.close();
				cout << "saving finished" << endl;
			}
			else
			{
				cout << "File not opened" << endl;
			}
		}
	}
	else if (punkt == "3")
	{
		//3.task
		//Design a type I quantum well with an AB barrier material whose energy width
		//at 300K corresponds to the wavelength λ. Assume the influence of the dimension constraint(quantum
		//confinement) for an 8nm wide well as equal to 100meV. Include stresses and try to minimize it.
		// (*) Solve the Schrödinger equation for the well and find the actual impact quantum confinement.
		// the (*) is morepart for more ambitious that was not worked on in this program

		// Parameters for task 3:
		// InP substrate
		//Wavelength 1550nm
		Alloy InP = Alloy(5.8697, 5.8697, 1.4236, -0.94, -0.6, -6.0, -2.0, 1011, 561, 1.4236, 3.63e-4, 162, 2.79e-5, 5.8697);
		InP.UpdateTemp(300,false);
		double lambda_nm = 1550;
		double DeltE = 6.626e-34 * 299792000 / (lambda_nm * 1e-9) / 1.602e-19;
		double needed_DeltE = DeltE - 0.1;

		Alloy InAs = Alloy(6.0583 + 2.74e-5 * (-300), InP.get_a(), 0.41, -0.59, -1.00, -5.08, -1.8, 832.9, 452.6, 0.417, 2.76e-4, 93, 2.74e-5, 6.0583);
		Alloy GaAs = Alloy(5.65325 + 3.88e-5 * (-300), InP.get_a(), 1.519, -0.84, -1.16, -7.17, -2.0, 1221, 566, 1.519, 5.405e-4, 204, 3.88e-5, 5.65325);
		TwoAlloy GaInAs = TwoAlloy(GaAs, InAs, InP.get_a(), 0, -0.477, 0.38, 0, -2.61, 0, 0, 0);
		GaInAs.UpdateTemp(300);
		double konc;
		vector<double> x,diffhh,difflh,pom;
		for (int i = 0; i < 100000; i++)//this loop was performed over commented ranges to find dpoing giving smallest diffrence in energy. For automating this one should use  bisection algorithm
		{
			//konc = i * 1e-5;
			//konc = 0.65 + i * 5e-7;
			//konc = 0.686 + i * 1e-8;
			//konc = 0.6864 + i * 1e-9;
			//konc = 0.68643 + i * 2e-10;
			//konc = 0.6875 + i * 5e-9;
			konc = 0.6876 + i * 1e-9;
			pom = GaInAs.calc_bands(konc);
			x.push_back(konc);
			diffhh.push_back(abs(pom.at(0) - pom.at(1) - needed_DeltE));
			difflh.push_back(abs(pom.at(0) - pom.at(2) - needed_DeltE));
		}
		plik.open("zad3/result.txt", ios::out);
		if (plik.is_open())
		{
			cout << "saving to file..." << endl;
			for (int i = 0; i < x.size(); i++)
			{
				plik << x.at(i) << " " << diffhh.at(i) << " " << difflh.at(i) << endl;
			}
			plik.close();
			cout << "saving finished" << endl;
		}
		else
		{
			cout << "File not opened" << endl;
		}

	}
	return 0;
}

