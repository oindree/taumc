// Monte Carlo for tau propagation and energy loss in rock (later, air) 
// 2-24-2016 - initial attempt
// 4-18-2016 - converting from matlab
// Oindree Banerjee (banerjee.104@osu.edu) 

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include "TCanvas.h"
#include "TGraph.h" 
#include "TTreeIndex.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TText.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include <unistd.h>
#include "TVector3.h"
#include "TRotation.h"
#include "TSpline.h"
#include "Math/InterpolationTypes.h"
#include "TGaxis.h"
#include "TPaveStats.h"


using namespace std;

static const int middle_num_rows = 26;
static const int middle_num_cols = 2;
static const double tau_lifetime = 2.906e-13; // in seconds 
static const double tau_mass = 1.7e9; // in eV/c^2
static const double c = 3.0e8; // m/s
double curve_middle[middle_num_rows][middle_num_cols];
double x_axis_middle[middle_num_rows];
double y_axis_middle[middle_num_rows]; 



int main()
{ 
// define constants and declare variables
//
  double x_max = 50000.0;                                                // max total dist in meters that tau can travel
  double N = 6.023e23;                                                   // Avogadro number
  double A = 22.0;                                                       // "mass number" of rock as given in paper, are the units in grams? 
  double unit = 1.0e-6;                                                  // g^-1 cm^2
  double proton_mass = 1.67e-24;                                         // units of grams
  double rho = 4.43;                                                     // for titanium (rock) whose mass number is 22
  double energy = 1.0e18;                                                // initial energy of tau in eV 
  double sigma_middle = 0.0; 

  int num_interpl = 3000;                                                // number of points to interpolate
  double x_axis_middle_interpl[num_interpl];                             // new x axis with interpolated points 
  double y_axis_middle_interpl[num_interpl];                             // new y axis with interpolated points 
  
  double step_sum_middle[25];	
  double y_cdf[25]; 

  double x = 0.0;                                                            // distance traveled by tau along its path
  double y = 0.0;                                                            // inelasticity: fractional energy loss of tau at each "POWWW!"
  double ycdf2 = 0.0; 
  double ycdf1 = 0.0;
  double y2 = 0.0; 
  double y1 = 0.0; 
  int top_index = 0; 
  double dE; 
  int num_int = 0;                                                           // number of interactions
 
  std::vector <double> dx_vector;                                             // will hold dx for each interaction
  std::vector <double> y_vector;                                              // will hold y for each interaction
  std::vector <double> dE_vector;                                             // will hold change in energy for each interaction
  std::vector <double> E_vector;                                              // will hold energy after each interaction
                          
 
  cout << "Welcome to oindreeMC" << endl; 
  cout << "The initial energy of tau is " << energy << " eV" << endl; 
  cout << "x_max is " << x_max << " m" << endl;

std::ifstream file;
    file.open("ALLM_curve_middle.txt");

    if(file.is_open())
    {
            std::cout << "File Opened successfully!!! Reading data from file into array." << std::endl;
            while(!file.eof())
            {
                    for(int i = 0; i < middle_num_rows; i++)
                    {
                            for(int j = 0; j < middle_num_cols; j++)
                            {
                                  
                                    file >> curve_middle[i][j];
                                    
                            }
                            
                    }
            }

    }
    file.close();

    for(int i = 0; i < middle_num_rows; i++)
      {
	for(int j = 0; j < middle_num_cols; j++)
	  {
	    //cout << curve_middle[i][j] << endl;
	    x_axis_middle[i] = curve_middle[i][0];
	    y_axis_middle[i] = curve_middle[i][1];
	    
	  }
      }

//to check that I got the correct arrays having numbers from x-axis and y-axis of paper
cout << "x_axis_middle " << endl; 
cout << x_axis_middle[0] << endl;
cout << x_axis_middle[1] << endl;
cout << x_axis_middle[2] << endl; 

cout << "y_axis_middle " << endl; 
cout << y_axis_middle[0] << endl;
cout << y_axis_middle[1] << endl;
cout << y_axis_middle[2] << endl; 
 


// TAKE LOG10 of X-AXIS OF PAPER AND CALL THAT X-AXIS FROM NOW ON
 for(int i = 0; i < middle_num_rows; i++)
   {
     x_axis_middle[i] = log10(x_axis_middle[i]);
     y_axis_middle[i] = y_axis_middle[i] * (A / N) * unit * log(10); // CORRECTING Y-AXIS TO HAVE UNITS OF CM^2
   }

 double xrange = x_axis_middle[25] - x_axis_middle[0]; 

 cout << "log10 of x_axis_middle " << endl;
 cout << x_axis_middle[0] << endl;
 cout << x_axis_middle[1] << endl;
 cout << x_axis_middle[2] << endl; 
 cout << x_axis_middle[25] << endl; 
 cout << "x range " << endl; 
 cout << xrange << endl; 


{
   TCanvas *c1 = new TCanvas("c1","yaxis vs log10 xaxis",200,10,700,500);
   TGraph *gr = new TGraph(middle_num_rows,x_axis_middle,y_axis_middle);
   gr -> SetTitle("y dsig/dy vs log10 y");
   gr -> GetXaxis()->SetTitle("log10 y");
   gr -> GetYaxis()->SetTitle("y dsig/dy (cm^2)");
   gr -> Draw("AC*");
   c1 -> SaveAs("plots/fig8paper_usinglog10y.jpg");
   c1 -> Update();
   c1 -> Modified();
}


// INTERPOLATE

//ROOT::Math::Interpolator inter(middle_num_rows, ROOT::Math::Interpolation::kLINEAR);

 for ( int i = 0; i < num_interpl; ++i )
   {
     x_axis_middle_interpl[i]  = ((double) i * (xrange)/(num_interpl - 1)) + x_axis_middle[0];
     //cout << "x axis interpolated " << endl; 
     //cout << x_axis_middle_interpl[i] << endl;      
 //iy[i] = inter.Eval(x[i]);
   }


 cout << "length of x-axis, y-axis, x-axis-interpl" << endl; 
 cout << sizeof(x_axis_middle) / sizeof(x_axis_middle[0]) << endl; 
 cout << sizeof(y_axis_middle) / sizeof(y_axis_middle[0]) << endl; 
 cout << sizeof(x_axis_middle_interpl) / sizeof(x_axis_middle_interpl[0]) << endl; 
 


 for ( int i = 0; i < int((sizeof(x_axis_middle)/sizeof(x_axis_middle[0]))) - 1; i++ )
   {
     step_sum_middle[i] = (x_axis_middle[i+1] - x_axis_middle[i]) * (((y_axis_middle[i+1] + y_axis_middle[i]))/2);
     sigma_middle = sigma_middle + step_sum_middle[i]; 

   } 

 cout << "sigma_middle" << endl;
 cout << sigma_middle << endl; 
 
 
 for ( int i = 0; i < int((sizeof(step_sum_middle)/sizeof(step_sum_middle[0]))); i++ )
   {
    for ( int j = 0; j <= i; j++ )
      { 
    	y_cdf[i] += step_sum_middle[j]; 
      }
   }	

 // find max of y_cdf

 double ycdfmax = TMath::MaxElement(25, y_cdf);  // ycdfmax should equal sigma and it does 
 cout << "ycdfmax " << endl; 
 cout << ycdfmax << endl; 

 // divide out ycdfmax from each element of ycdf 

 for ( int i = 0; i < int((sizeof(y_cdf)/sizeof(y_cdf[0]))); i++ )
   {
     y_cdf[i] = y_cdf[i] / ycdfmax; 
   }

 // find interaction length

 double int_length_middle = proton_mass / (rho * sigma_middle); 
 double int_length_m_middle = int_length_middle / (1e2); 
 cout << "interaction length in meter " << endl; 
 cout << int_length_m_middle << endl; 


 // WHILE LOOP x < x_max

 TRandom r; // for picking dx, generate a number in interval [0,1] (0 is excluded) 
 TRandom2 s; // for picking y, need to check that these are good to use
 TRandom3 t; // for determining if tau will decay

 while (x <= x_max)
   {
     // RANDOMLY PICK THE dx's TRAVELED BY THE TAU ALONG ITS PATH


     // pick random number between 0 and 1 and save it

     double rand_x = r.Rndm();
     cout << "rand_x is " << endl; 
     cout << rand_x << endl; 
     double dx = -int_length_m_middle * log(1 - rand_x); 
     x = x + dx;  
     dx_vector.push_back( dx );
     
     //cout << " dx " << endl; 
     //cout << dx << endl; 

     // DETERMINE IF TAU WILL DECAY FOR THIS dx 
     
     double gamma = energy / tau_mass; 
     double rand_decay = t.Rndm();
     cout << "rand_decay is " << endl;
     cout << rand_decay << endl; 
     double dx_decay = -(c * gamma * tau_lifetime) * log(1 - rand_decay); 
     dx_decay = dx; // for testing
     if (dx_decay == dx)
       {
	 cout << "TAU DECAYS FOR THIS dx" << endl;
	 break; 
       }
     

     

     
   
     // RANDOMLY PICK THE y (fractional energy lost by the tau) for this dx 
     // first find random number 
     double rand_y = s.Rndm(); 
     cout << "rand_y is " << endl; 
     cout << rand_y << endl; 

     // use cdf of y to randomly pick y 
     
     
     if (rand_y < y_cdf[0]) // case where random number less than min ycdf 
       {
	 ycdf2 = y_cdf[0];
	 ycdf1 = 0; 
	 cout << "rand_y less than min y_cdf" << endl; 
	 //cout << ycdf2 << endl; 
	 //cout << ycdf1 << endl; 
	 y2 = x_axis_middle[0]; 
	 y1 = 0.0; 
	 double slope = (ycdf2 - ycdf1) / (y2 - y1); 
	 double intercept = ycdf2 - (slope * y2); 
	 y = (rand_y - intercept) / slope;
       }
     else if (rand_y == y_cdf[0])  // case where random number equal to min ycdf 
       { 
	 cout << "rand_y equal to min y_cdf" << endl; 
	 y = x_axis_middle[0]; 
       }
     else     // case where random number greater than min ycdf
       {
	  for (int i = 1; i < int(sizeof(y_cdf) / sizeof(y_cdf[0])); i++)
	   {
	     if (y_cdf[i] > rand_y)
	       {
		 ycdf2 = y_cdf[i]; 
		 ycdf1 = y_cdf[i-1];
		 top_index = i; 
		 break; 
	       }
	   }
	  y2 = x_axis_middle[top_index];
	  y1 = x_axis_middle[top_index - 1]; 
	  double topminusrand = ycdf2 - rand_y; 
	  double randminusbot = rand_y - ycdf1; 
	  if (topminusrand == randminusbot)                   // case where random number fell right in the middle between top and bot ycdf points 
	    {
	      y = (1./2.) * (y2 + y1); 
	    }
	  else if (topminusrand > randminusbot)              // case where random number closer to bot ycdf point
	    { 
	      y = y1; 
	    }
	  else                                               // case where random number closer to top ycdf point 
	    { 
	      y = y2;   
	    } 
       }

     // push back to y_array
     y_vector.push_back( y ); 
     
   
	       
     // have a y, unlog it to get the actual fractional value 

     double y_notlog = pow(10,y);
     cout << "y now is" << endl;
     cout << y_notlog << endl; 
	 
     num_int = num_int + 1;   // to get number of interactions


     cout << "energy before minusing this y:" << endl; 
     cout << energy << endl; 

     
     // find change in energy
     dE = -y_notlog * energy; 
     
     // push back to dE array
     dE_vector.push_back( dE ); 

     // update energy
     energy = energy + dE;

     // push back to energy array
     E_vector.push_back( energy ); 
     


     cout << "energy after minusing this y:" << endl; 
     cout << energy << endl; 

      
   }

 cout << "number of interactions" << endl; 
 cout << num_int << endl; 

 //declare array for dEdx and energy for TGraph plotting later 
 double dEdx_array[num_int]; 
 double E_array[num_int]; 


 // Fill array dEdx_array with log10 of -dE/dx and E_array with log10 of E_vector

 for (int i = 0; i < num_int; i++)
   {
     dEdx_array[i] = ( log10(-dE_vector[i] / dx_vector[i]) );
     E_array[i] = log10(E_vector[i]); // turning into array for plotting
   }


 //cout << "y_array [0]" << endl; 
 //cout << y_array[0] << endl; 


 // Make hist of dx's 
 {
   TCanvas *c3 = new TCanvas("c3", "hist_dx", 200, 10, 700, 500);
   TH1D *hdx_array = new TH1D("hist_dx", "; dx; counts", 100, 0, 100); 

   for (int i = 0; i < int(dx_vector.size()); i++) 
     { 
       hdx_array->Fill(dx_vector[i]);
     } 
   hdx_array -> Draw(); 
   hdx_array -> SetTitle("Histogram of dx");
   hdx_array -> GetXaxis() -> SetTitle("dx (meter)");
   hdx_array -> GetYaxis() -> SetTitle("counts"); 
   c3 -> SaveAs("plots/dx_hist.jpg"); 
   c3 -> Update();
   c3 -> Modified(); 
 }

 // Make hist of y's 
 {
   TCanvas *c4 = new TCanvas("c4", "hist_y", 200, 10, 700, 500);
   TH1D *hy_array = new TH1D("hist_y", "; y; counts",10,-4,0); 

   for (int i = 0; i < int(y_vector.size()); i++) 
     { 
       hy_array->Fill(y_vector[i]);
     } 
   hy_array -> Draw(); 
   hy_array -> SetTitle("Histogram of y");
   hy_array -> GetXaxis() -> SetTitle("log10 y");
   hy_array -> GetYaxis() -> SetTitle("counts"); 
   c4 -> SaveAs("plots/y_hist.jpg"); 
   c4 -> Update();
   c4 -> Modified(); 
 }

 
 // Make plot of -dE/dx vs energy 


{
   TCanvas *c5 = new TCanvas("c5","dEdx vs E", 200, 10, 700, 500);
   TGraph *dEdxvsE = new TGraph(num_int,E_array,dEdx_array);
   dEdxvsE -> SetTitle("log10(-dE/dx) vs log10(E)");
   dEdxvsE -> GetXaxis()->SetTitle("log10( E / eV )");
   dEdxvsE -> GetYaxis()->SetTitle("log10( (-dE/dx) / (eV/m) )");
   dEdxvsE -> Draw("AC*");
   c5 -> SaveAs("plots/dEdxvsE.jpg");
   c5 -> Update();
   c5 -> Modified();
}


 // Make plot of y cdf vs log10 of y

{
   TCanvas *c2 = new TCanvas("c2","y_cdf", 200, 10, 700, 500);
   TGraph *ycdf = new TGraph(middle_num_rows-1,x_axis_middle,y_cdf);
   ycdf -> SetTitle("cdf of y");
   ycdf -> GetXaxis()->SetTitle("log10 y");
   ycdf -> GetYaxis()->SetTitle("cdf of y (cm^2)");
   ycdf -> Draw("AC*");
   c2 -> SaveAs("plots/y_cdf.jpg");
   c2 -> Update();
   c2 -> Modified();
}




  return 0;
}


 











