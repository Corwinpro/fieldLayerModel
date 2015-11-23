#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <complex>      // std::complex, std::abs
#include <cstdlib>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

const double PI = 3.14159265359;
const double cspeed = 3.0e8;

const double ampl = 3.7e-5;
const double n0 = 1.0+ampl*10.0;


const double bulkLength = 2.0e2;
const int discretSize = 1000;
const int size = int(bulkLength)*discretSize;


const double latticeSpeed = 0.0;

const double lambda0 = 1.0*double(discretSize);
const double k0 = 2.0*PI/lambda0;

const double lambda1 = lambda0*(1.0 + 2.0*latticeSpeed/cspeed*n0);
const double k1 = 2.0*PI/lambda1;

const double alpha = 2.88e-40;
const double Dens0 = 1.98e25;
const double eps0 = 8.85e-12;
const double mu0 = 1.26e-6;
const double kb = 1.38e-23;
const double T = 300.0;
const double Eo2 = 1.0e17 / sqrt(eps0/mu0 * (1.0 + alpha*Dens0/eps0) );

template <typename T>
    std::string to_string(T value)
    {
      //create an output string stream
      std::ostringstream os ;

      //throw the value into the string stream
      os << value ;

      //convert the string stream into a string and return
      return os.str() ;
    }

void SetNCentralized(double * n){
	int BulkCenter = size / 2;
	
	for (int i = 0; i <= BulkCenter; i++)
	{
		*(n+BulkCenter+i) = n0 + ampl*cos((k0+k1)*n0*double(i));
		*(n+BulkCenter-i) = n0 + ampl*cos((k0+k1)*n0*double(i));
	}
	return;
}

void set_n (double * n,
			double phase,
			bool random)
{

	std::srand(std::time(0)); // use current time as seed for random generator
    
    cout << "bulk Length is: " << bulkLength*532.0e-7 << " cm\n" << endl;
    cout << "lattice speed: " << cspeed*(lambda1 - lambda0) / (lambda1 + lambda0) / n0 << " m/s" << endl; 
	
	SetNCentralized(n);

	/*while ( abs(cos((k0+k1)*n0*double(size)) - 1.0) > 1.0e-2)
		size++;
	if ( abs(cos((k0+k1)*n0*double(size)) - cos((k0+k1)*n0*double(size+1))) < 1.0e-2)
		size++;
	cout << "added " << size - bulkLength*discretSize << " points." << endl;*/

	/*for(int i = 0; i <= size; i++){

		double random_variable = double(std::rand()) / double(RAND_MAX);

		*(n+i) = n0 + ampl*(cos((k0+k1)*n0*double(i) + phase) )*( exp(-4.0*pow(double(i-size/2)/(double(size)/2.0),2.0) ) - exp(-4.0) );
		//*(n+i) = n0 + ampl*(cos(n0*(k0+k1)*double(i)) + 0.0*(-2.0*random_variable+1.0));
		//(exp(-4.0*log(10.0)*pow((double(i)-double(size)/2)/double(size),2.0)) - 0.1)*
	}*/
	return;
}

void NormalizeFields (complex<double> * EL_plus, 
					  complex<double> * EL_minus, 
					  complex<double> * ER_plus, 
					  complex<double> * ER_minus)
{	
	for (int i = size; i >= 0; i--)
	{	
		*(EL_minus+i) = *(EL_minus+i) / *(EL_plus);
		*(EL_plus+i) = *(EL_plus+i) / *(EL_plus);
	}
		for (int i = 0; i <= size; i++)
	{
		*(ER_plus+i) = *(ER_plus+i) / *(ER_minus+size);
		*(ER_minus+i) = *(ER_minus+i) / *(ER_minus+size);
	}
	return;
}

void CalcFullField (  double * n,
					  complex<double> * EL_plus, 
					  complex<double> * EL_minus, 
					  complex<double> * ER_plus, 
					  complex<double> * ER_minus,
					  double * E_full)
{
	
	for (int i = 0; i <= size; i++){
		complex<double> phase0 (0.0, *(n+i)*k0*(double(i)-double(size)/2.0));
		complex<double> phase1 (0.0, *(n+i)*k1*(double(i)-double(size)/2.0));

		*(E_full+i) = pow(abs( *(EL_plus+i)*exp(phase0) + *(ER_plus+i) * exp(phase1) + *(EL_minus+i)*exp(-phase0)  + *(ER_minus+i)*exp(-phase1) ) , 2.0);
	}

	return;
}

void CalcFields (double * n, 
				complex<double> * EL_plus, 
				complex<double> * EL_minus, 
				complex<double> * ER_plus, 
				complex<double> * ER_minus,
				double * E_full)
{
	complex<double> Eleft_plus;
	complex<double> Eleft_minus;
	complex<double> Eright_plus;
	complex<double> Eright_minus;

	*(EL_plus+size) = Eright_plus = complex<double>(1.0,0.0);		//Волна, которая вылетает справа из среды
	*(EL_minus+size) = Eright_minus = complex<double>(0.0,0.0);		//Волна, которая влетает в среду справа (нулевая)

	for (int i = size; i >= 1; i--)	//волна справа налево
	{
		double nr = *(n+i);	// Текущие показатели преломления
		double nl = *(n+i-1);

		complex<double> phaseR (0.0, k0*nr*(double(i)-double(size+1)/2.0)); // i * kR * z
		complex<double> phaseL (0.0, k0*nl*(double(i)-double(size+1)/2.0)); // i * kL * z

		//Расчет полей
		Eleft_plus = (Eright_plus*exp(phaseR)*(nr+nl)/(2.0*nl) + Eright_minus*exp(-phaseR)*(nl-nr)/(2.0*nl)) * exp(-phaseL); 
		Eleft_minus = (Eright_plus*exp(phaseR)*(nl-nr)/(2.0*nl) + Eright_minus*exp(-phaseR)*(nr+nl)/(2.0*nl)) * exp(phaseL);
	
		//Заполнение и сохранение массивов
		*(EL_plus+i-1) = Eleft_plus;
		*(EL_minus+i-1) = Eleft_minus;

		// Переход на следующую границу
		Eright_plus = Eleft_plus;
		Eright_minus = Eleft_minus;
	}

	*(ER_plus) = Eleft_plus = complex<double>(0.0,0.0);	//Волна, которая вылетает справа из среды 
	*(ER_minus) = Eleft_minus = complex<double>(1.0,0.0);	//Волна, которая влетает в среду справа (нулевая)

	for (int i = 1; i<= size; i++) //волна слева направо
	{
		double nr = *(n+i);	// Текущие показатели преломления
		double nl = *(n+i-1);

		complex<double> phaseR (0.0, k1*nr*(double(i)-double(size+1)/2.0)); // i * kR * z !!! -1.0
		complex<double> phaseL (0.0, k1*nl*(double(i)-double(size+1)/2.0)); // i * kL * z

		//Расчет полей
		Eright_plus = (Eleft_plus*exp(phaseL)*(nr+nl)/(2.0*nr) + Eleft_minus*exp(-phaseL)*(nr-nl)/(2.0*nr)) * exp(-phaseR);
		Eright_minus = (Eleft_plus*exp(phaseL)*(nr-nl)/(2.0*nr) + Eleft_minus*exp(-phaseL)*(nr+nl)/(2.0*nr)) * exp(phaseR);
	
		*(ER_plus+i) = Eright_plus;
		*(ER_minus+i) = Eright_minus;

		Eleft_plus = Eright_plus;
		Eleft_minus = Eright_minus;
	}

	NormalizeFields(EL_plus,EL_minus,ER_plus,ER_minus);

	CalcFullField(n,EL_plus,EL_minus,ER_plus,ER_minus,E_full);

	
	return;
}

void PrintResults ( int j,
					double * n,
					complex<double> * EL_plus, 
					complex<double> * EL_minus, 
					complex<double> * ER_plus, 
					complex<double> * ER_minus,
					double * E_full)
{
	std::string file;
	file = "file_n_" + to_string(j) + ".txt";
	const char *cstr = file.c_str();
	FILE * file_n = fopen(cstr,"w");


	FILE * file_E_plus = fopen("file_E_plus.txt","w");
	FILE * file_E_minus = fopen("file_E_minus.txt","w");
	FILE * file_Energy = fopen("file_Energy.txt","w");

	std::string file2;
	file2 = "file_Efull_" + to_string(j) + ".txt";
	const char *cstr2 = file2.c_str();
	FILE * file_Efull = fopen(cstr2,"w");



	/*fprintf(file_E_plus,  "%e %e %e %e %e\n", double(i)/lambda0,  real(*(EL_plus+i)), imag(*(EL_plus+i)),  real(*(ER_plus+i)), imag(*(ER_plus+i))  );
	//fprintf(file_E_minus, "%e %e %e %e %e\n", double(i)/lambda0,  real(*(EL_minus+i)),imag(*(EL_minus+i)), real(*(ER_minus+i)), imag(*(ER_minus+i)) );
	//fprintf(file_diff, "%e %e %e\n", double(i)/lambda0, abs(*(EL_plus+i)) - abs(*(ER_minus+size-i)), abs(*(ER_plus+i)) - abs(*(EL_minus+size-i)));

	// file_Energy хранит вектор Пойнтинга	
	//fprintf(file_Energy, "%e %e %e\n", double(i)/lambda0, ((pow(abs(*(EL_plus+i)),2.0)) - pow(abs(*(EL_minus+i)),2.0))*(*(n+i)), ((pow(abs(*(ER_plus+i)),2.0)) - pow(abs(*(ER_minus+i)),2.0))*(*(n+i)));

	//complex<double> phase (0.0, *(n+i)*k0*(double(i)-double(size)/2.0));
	//Efield = pow(abs( exp(phase) * (*(EL_plus+i) + *(ER_plus+i)) + exp(-phase)*( *(EL_minus+i) + *(ER_minus+i)) ) , 2.0); */

	for (int i = 1; i <= size; i++)
		{
			//if ( *(E_full+i) > *(E_full+i+1) && *(E_full+i) > *(E_full+i-1))
				fprintf(file_Efull, "%e %e %e %e\n", double(i)/lambda0, abs (*(EL_plus+i) + (*(ER_plus+i)) ), abs( *(EL_minus+i) + (*(ER_minus+i))), *(E_full+i) );
				//fprintf(file_Efull, "%e %e %e\n", double(i)/lambda0, abs(*(EL_plus+i) ), abs( *(EL_minus+i) ));	
			//if ( *(n+i) > *(n+i+1) && *(n+i) > *(n+i-1))
				fprintf(file_n, "%e %e\n", double(i)/lambda0, *(n+i) );
		}
		fprintf(file_n, "\n");
		cout << "Written!" << endl;

		fclose(file_n);
		fclose(file_E_plus);
		fclose(file_E_minus);
		fclose(file_Energy);
		fclose(file_Efull);
}

double CalcAveragePotent(double * E_full)
{
	double Average = 0.0;
	for (int i = 0; i <= size; i++)
		Average += exp( *(E_full+i) * alpha*Eo2/(4.0*kb*T) );

	return Average / double(size+1);
}

void PrintFullDensity(int j, double * Density)
{
	std::string fileName;
	fileName = "dens_" + to_string(j) + ".txt";
	const char *cstr = fileName.c_str();
	FILE * file = fopen(cstr,"w");

	for (int i = 0; i <= size; i++)
		fprintf(file, "%e %e\n", double(i)/lambda0, *(Density+i));
	fclose(file);
}

void ReadDensitFromFile(double * Density)
{
	ifstream fin("dens_14800.txt");
	double d;
	for (int i = 0; i <= size; i++)
	{
   		fin >> d;
    	fin >> d;
    	*(Density+i) = d;
	}

	/*FILE * file = fopen("compDens.txt", "w");
	for (int i = 0; i <= size ; i++)
		fprintf(file, "%e %e\n", double(i)/lambda0, *(Density+i));

	fprintf(file, "\n");

	cout << setprecision(15) << *(Density+size) << endl;*/

	fin.close();
	cout << "Reading from file complete!" << endl;
}
int main()
{
	double * n = new double [size + 1];	
	
	complex<double> * EL_plus = new complex<double> [size+1]; 	// Массив амплитуд волн, которые распространяются слева направо от левого лазера
	complex<double> * EL_minus = new complex<double> [size+1];	// Массив амплитуд волн, которые распространяются справа налево
	complex<double> * ER_plus = new complex<double> [size+1]; 	// Массив амплитуд волн, которые распространяются слева направо от правого лазера
	complex<double> * ER_minus = new complex<double> [size+1];	// Массив амплитуд волн, которые распространяются справа налево
	double * E_full = new double [size+1];	// Массив полного поля

	/*int j = 0;
	bool random = false;
	set_n(n,0.0,random);
	CalcFields(n,EL_plus,EL_minus,ER_plus,ER_minus,E_full);
	PrintResults(j,n,EL_plus,EL_minus,ER_plus,ER_minus,E_full);*/

	//Setting 0 iteration medium - no perturbations
	for (int i = 0; i <= size; i++)
		*(n+i) = sqrt(1.0 + alpha*Dens0/eps0);
	

	double DensAverage;
	double * Density = new double [size+1];
	for (int i = 0; i <= size; i++)
		*(Density+i) = Dens0;

	//ReadDensitFromFile(Density);

	//iterations Field->Medium->...
	for (int j = -1; j <= 4000000; j++)
	{
		if (j%100 == 0)
			cout << "Iteration #" << j << endl;
	
		CalcFields(n,EL_plus,EL_minus,ER_plus,ER_minus,E_full);

		double PotentAverage = CalcAveragePotent(E_full);

		for (int i = 0; i<=size; i++)
		{
			*(Density+i) =  *(Density+i)*0.999 + 0.001 * Dens0 * exp( *(E_full+i) * alpha*Eo2/(4.0*kb*T))/ PotentAverage;
			*(n+i) = sqrt( 1.0 + *(Density+i)*alpha/eps0) ;
		}

		if(j%1000==0)
			PrintResults(j,n,EL_plus,EL_minus,ER_plus,ER_minus,E_full);
		if (j%1000==0)
			PrintFullDensity(j,Density);

	}

	delete [] EL_plus;
	delete [] EL_minus;
	delete [] ER_plus;
	delete [] ER_minus;
	delete [] E_full;

	//delete [] n;
	//delete [] Density;

	return 0;
}
