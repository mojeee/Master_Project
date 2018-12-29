/*
Author: Mojtaba Amini
Created on : 1 Pm, 24Dec 2018
end time : 

*/

#include<tripletMap.h>
#include<math.h>
#include<chemistry0d.h>
#include<grid.h>

class tripletMap:public grid ,public flamesolver , public chemistry0d 
{


	public:
		// number of species get from flame solver
		nspc=nSpec;
		// number of grid point get from grid
		nc=nPoints;
		ncp1=nc+1;
		Dx=0.1 // cm ;
		XLint=0.2992;
		Re=500;
		XLk=(XLint)/(pow(Re,0.75));
		PDFA=((pow(XLint,double(5/3)))*(pow(XLk,double(-5/3))))/(pow((XLint-XLk),double(5/3))-1.0);
		PDFB=-(pow(XLint,double(5/3)))/(pow((XLint/XLk),double(5/3))-1.0);
		// get pressure from chemistry0d
		Dom = pressure;
		ncm1=nc-1;
		// number of triplet map in each realization
		MTS=1;
		// random number
		rand=0.0;
		// eddy length
		L=0;
		// starting point
		M=0;
}

// Triplet class function
void tripletMap::Random_Number()
{

	rand=double(rand())/RAND_MAX;

}

void tripletMap::TM()
{
	// MTS shows number of triplet map require in each realization 
	// first of all call a random number 
	// 
	for (int i=1,i<MTS,i++)
	{
		// first of all call a random number
		Random_Number();
		// m determine the starting point of triplet map
		M=int(rand*nc)
		// this loop check the starting point of triplet map
		while(M<(ncp1/4))
		{
			Random_Number();
			M=int(rand*nc);
		}
        	// calculate the eddy length
		eddylength();	
		// check eddy size does not exceed domain
		if((L+M)>nc)
		{
			M=nc-L;
		}

		while((L+M)>nc || M<(ncp1/4))
		{
			Random_Number();
			M=int(rand*nc);
			eddylength();
		}
		
		for(int j=1,j<nc,j++)
                {
                        double Temp(nPoints);
                        Temp(j)=T(j);
                }
                        Btriplet(double& Temp);
                for(j=1,j<nc,j++)
                {
                        T(j)=Temp(j);
                }
                        
		
		// tripletmap on mass fraction
		for (int k=1,k<nspc,k++)
		{
			for(int j=1,j<ncp1,j++)
			{
				double yy(nPoints);
				yy(j)=Y(k,j);
			}
			Btriplet(double& yy);
			for(j=2,j<nc,j++)
			{
				Y(k,j)=yy(j);
			}
			Y(k,ncp1)=Y(k,nc);
		}
	}
}

void tripletMap::eddyLength()
{
	// generate a random number between 0 and 1
	Random_number();
	// make sure eddy is Greater than 6 cells long
	int  NSize = int(pow((rand-PDFA)/PDFB,(-3.0/5.0))/Dx);
	
	while (NSize<5)
	{
		NSize = int(pow((rand-PDFA)/PDFB,(-3.0/5.0))/Dx);
	}
	// make sure eddy is divisible by 3
	if((NSize%3)==0)
	{
		int Nlength=NSize;
	}
	else if ((NSize%3)==1)
	{
		Nlength=NSize-1;
	}
	else if ((NSize%3)==2)
	{
		Nlength=NSize+1;
	}
	L=Nlength;

}

void tripletMap::BTriplet(double& S)
{
	// Permute Cells M through M+L-1 of the array S
	// as prescribrd by the discrete triplet map , where L is an integer multiple by 3.
	int Lo=L/3;
	// first part of mapping
	for(int j=1,j<Lo,j++)
	{
		int k=M*(j-1);
		double X(j)=S(k); //gather the cells going to the 1st image 
	}
	// second part of mapping
	for (j=1,j<Lo,j++)
	{
		k=M+L+1+(3*j); // minus sign because second image is flipped
		X(j+Lo)=S(k); // gather the cells going to 2nd image
	}
	
	// third part of mapping
	for (j=1,j<Lo,j++)
	{
		k=M+(3*j)-1;
		X(j+Lo+Lo)=S(k);// gather the cells going to the 3rd image 
	}
	
	for(j=1,j<L,j++)
	{
	k=M+j-1;
	S(k)=X(j);
	}


}












