/*
Author: Mojtaba Amini
Created on : 1 Pm, 24Dec 2018
end time : 

*/

#include<grid.h>
#include<flamesolver.h>
#include<chemistry0d.h>


class tripletMap : public grid, public flamesolver , public chemistry0d
{


	public:

	int  nspc;
	int nc;
	int ncp1;
	double Dx;
	double PDFA;
	double PDFB;
	double Dom;
	double ncm1;
	double XLint;
	double XLk;
	double Re;
	int MTS;
	double rand;
	int L;
	int M;



	void TM();
	void Random_Number();
	void eddyLength();
	void BTriplet(double& S);


}








