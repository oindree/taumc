// Monte Carlo program for tau propagation and energy loss in rock (later, air) 
// 2-24-2016 - initial attempt
// 4-18-2016 - converting from matlab
// 7-12-2016 - adding more comments, preparing readme, putting on github
// Oindree Banerjee (banerjee.104@osu.edu) 
//
//
//
// to run: first compile by typing "make" and then once compiled type "./taumc" without the quotes 


// includes 
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm> 
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <cmath> 
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
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

//static const double norm = 1.0; 
static const int middle_num_rows = 26;
static const int bottom_num_rows = 36; 
static const int top_num_rows = 28; 
static const int num_cols = 2; 
static const double tau_lifetime = 2.906e-13; // in seconds 
static const double tau_mass = 1.7e9; // in eV/c^2
static const double c = 3.0e8; // m/s

double curve_middle[middle_num_rows][num_cols];
double curve_bottom[bottom_num_rows][num_cols]; 
double curve_top[top_num_rows][num_cols];
 
double x_axis_middle[middle_num_rows];
double y_axis_middle[middle_num_rows]; 

double x_axis_bottom[bottom_num_rows];
double y_axis_bottom[bottom_num_rows];

double x_axis_top[top_num_rows];
double y_axis_top[top_num_rows]; 


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
  double sigma_middle = 0.0; 
  double sigma_top = 0.0; 
  double sigma_bottom = 0.0; 
  double sigma = 0.0; 
  int energy_index = 0; 
  double int_length = 0.0; 
  double int_length_m = 0.0; 

  int num_interpl = 5000;                                                // number of points to interpolate
  double x_axis_middle_interpl[num_interpl]; 
  double y_axis_middle_interpl[num_interpl]; 
  double x_axis_top_interpl[num_interpl];
  double y_axis_top_interpl[num_interpl];
  double x_axis_bottom_interpl[num_interpl];
  double y_axis_bottom_interpl[num_interpl]; 



  double paper_energy[3] = {15.0, 18.0, 21.0};     
  double y_axis[3]; 
  double x_axis[3]; 
  int num_energy_interpl = 12; 
  double energy_interpl[num_energy_interpl]; 
  double y_axis_interpl[num_energy_interpl][num_interpl]; 
  double x_axis_interpl[num_energy_interpl][num_interpl]; 

  double energydiff[num_energy_interpl]; 
  
  double step_sum_middle[num_interpl-1];
  double step_sum_top[num_interpl-1];
  double step_sum_bottom[num_interpl-1];
  double step_sum[num_interpl-1];  
  double y_cdf[num_interpl-1]; 

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
  std::vector <double> x_vector;                           
 
  cout << "Welcome to taumc!! This program is written by Oindree Banerjee - Department of Physics - The Ohio State University - banerjee.104@osu.edu" << endl; 
  cout << "Maximum distance x_max the taus are allowed to travel is currently set at " << x_max << " meter" << endl;

  double initial_energy;
  cout << "Enter initial energy in eV that you want the taus to have. Choose between 10^15 and 10^21 eV. For 10^21 eV type 1e21." << endl; 
  cin >> initial_energy; 

  int num_tau;
  cout << "Enter number of taus you want to run." << endl; 
  cin >> num_tau; 


std::ifstream file_middle;
    file_middle.open("ALLM_curve_middle.txt");

    if(file_middle.is_open())
    {
            std::cout << "Middle file Opened successfully!!! Reading data from file into array." << std::endl;
            while(!file_middle.eof())
            {
                    for(int i = 0; i < middle_num_rows; i++)
                    {
                            for(int j = 0; j < num_cols; j++)
                            {
                                  
                                    file_middle >> curve_middle[i][j];
                                    
                            }
                            
                    }
            }

    }
    file_middle.close();

    for(int i = 0; i < middle_num_rows; i++)
      {
	for(int j = 0; j < num_cols; j++)
	  {
	    x_axis_middle[i] = curve_middle[i][0];
	    y_axis_middle[i] = curve_middle[i][1];
	    
	  }
      }


std::ifstream file_bottom;
    file_bottom.open("ALLM_curve_bottom.txt");

    if(file_bottom.is_open())
    {
            std::cout << "Bottom file Opened successfully!!! Reading data from file into array." << std::endl;
            while(!file_bottom.eof())
            {
                    for(int i = 0; i < bottom_num_rows; i++)
                    {
                            for(int j = 0; j < num_cols; j++)
                            {
                                  
                                    file_bottom >> curve_bottom[i][j];
                                    
                            }
                            
                    }
            }

    }
    file_bottom.close();

    for(int i = 0; i < bottom_num_rows; i++)
      {
	for(int j = 0; j < num_cols; j++)
	  {
	    x_axis_bottom[i] = curve_bottom[i][0];
	    y_axis_bottom[i] = curve_bottom[i][1];
	    
	  }
      }


std::ifstream file_top;
    file_top.open("ALLM_curve_top.txt");

    if(file_top.is_open())
    {
            std::cout << "Top file Opened successfully!!! Reading data from file into array." << std::endl;
            while(!file_top.eof())
            {
                    for(int i = 0; i < top_num_rows; i++)
                    {
                            for(int j = 0; j < num_cols; j++)
                            {
                                  
                                    file_top >> curve_top[i][j];
                                    
                            }
                            
                    }
            }

    }
    file_top.close();

    for(int i = 0; i < top_num_rows; i++)
      {
	for(int j = 0; j < num_cols; j++)
	  {
	    x_axis_top[i] = curve_top[i][0];
	    y_axis_top[i] = curve_top[i][1];
	    
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

 double middle_xrange = TMath::MaxElement(sizeof(x_axis_middle)/sizeof(x_axis_middle[0]),x_axis_middle) - TMath::MinElement(sizeof(x_axis_middle)/sizeof(x_axis_middle[0]),x_axis_middle);

 for(int i = 0; i < bottom_num_rows; i++)
   {
     x_axis_bottom[i] = log10(x_axis_bottom[i]);
     y_axis_bottom[i] = y_axis_bottom[i] * (A / N) * unit * log(10); // CORRECTING Y-AXIS TO HAVE UNITS OF CM^2
   }

 double bottom_xrange = TMath::MaxElement(sizeof(x_axis_bottom)/sizeof(x_axis_bottom[0]),x_axis_bottom) - TMath::MinElement(sizeof(x_axis_bottom)/sizeof(x_axis_bottom[0]),x_axis_bottom); 

 for(int i = 0; i < top_num_rows; i++)
   {
     x_axis_top[i] = log10(x_axis_top[i]);
     y_axis_top[i] = y_axis_top[i] * (A / N) * unit * log(10); // CORRECTING Y-AXIS TO HAVE UNITS OF CM^2
   }

 double top_xrange = TMath::MaxElement(sizeof(x_axis_top)/sizeof(x_axis_top[0]),x_axis_top) - TMath::MinElement(sizeof(x_axis_top)/sizeof(x_axis_top[0]),x_axis_top);

 cout << "log10 of x_axis_middle " << endl;
 cout << x_axis_middle[0] << endl;
 cout << x_axis_middle[1] << endl;
 cout << x_axis_middle[2] << endl; 
 cout << x_axis_middle[25] << endl; 

// INTERPOLATE 

// with ROOT 


 // Plot fig8 from paper 


   TCanvas *c1 = new TCanvas("c1","yaxis vs log10 xaxis",200,10,700,500);
   TGraph *gr_middle = new TGraph(middle_num_rows,x_axis_middle,y_axis_middle);
   TGraph *gr_top = new TGraph(top_num_rows,x_axis_top,y_axis_top);
   TGraph *gr_bot = new TGraph(bottom_num_rows,x_axis_bottom,y_axis_bottom); 

   gr_top -> SetTitle("y dsig/dy vs log10 y");
   gr_top -> GetXaxis()->SetTitle("log10 y");
   gr_top -> GetYaxis()->SetTitle("y dsig/dy (cm^2)");

   gr_top -> Draw("AP*"); 
   gr_top -> SetMarkerColor(3);

   gr_middle -> Draw("same");
   gr_middle -> SetMarkerColor(2); 
  

   gr_bot -> Draw("same"); 
   gr_bot -> SetMarkerColor(4); 

   TSpline3 *sp_middle = new TSpline3("sp_middle",gr_middle);
   sp_middle->Draw("same");

   TSpline3 *sp_top = new TSpline3("sp_top",gr_top);
   sp_top->Draw("same"); 

   TSpline3 *sp_bot = new TSpline3("sp_bot",gr_bot);
   sp_bot->Draw("same"); 
   
   cout << "test evaluating spline at specific x values " << endl; 
   cout << sp_middle->Eval(-4.0) << "   " << sp_middle->Eval(-3.2) << endl;  

   for (int i = 0; i < num_interpl; ++i) 
     {
      x_axis_middle_interpl[i] = ((double) i * (middle_xrange)/(num_interpl -1)) + x_axis_middle[0];
      y_axis_middle_interpl[i] = sp_middle -> Eval(x_axis_middle_interpl[i]); 
     }

   for (int i = 0; i < num_interpl; ++i) 
     {
      x_axis_top_interpl[i] = ((double) i * (top_xrange)/(num_interpl -1)) + x_axis_top[0];
      y_axis_top_interpl[i] = sp_top -> Eval(x_axis_top_interpl[i]); 
     }

   for (int i = 0; i < num_interpl; ++i) 
     {
      x_axis_bottom_interpl[i] = ((double) i * (bottom_xrange)/(num_interpl -1)) + x_axis_bottom[0];
      y_axis_bottom_interpl[i] = sp_bot -> Eval(x_axis_bottom_interpl[i]); 
     }
   
   // get a TGraph to verify interpolated arrays

   TGraph *interpl_middle = new TGraph(num_interpl, x_axis_middle_interpl, y_axis_middle_interpl);
   interpl_middle -> Draw("same");
   interpl_middle -> SetMarkerColor(6); 

   

   c1 -> SaveAs("plots/fig8paper_usinglog10y.jpg");
   c1 -> SaveAs("plots/fig8paper_usinglog10y.root");
   c1 -> Update();
   c1 -> Modified();
   gPad -> Update();
   gPad -> Modified(); 


 cout << "length of x-axis, y-axis, x-axis-interpl" << endl; 
 cout << sizeof(x_axis_middle) / sizeof(x_axis_middle[0]) << endl; 
 cout << sizeof(y_axis_middle) / sizeof(y_axis_middle[0]) << endl; 
 cout << sizeof(x_axis_middle_interpl) / sizeof(x_axis_middle_interpl[0]) << endl; 



 // INTERPOLATE TO GET CURVES IN BETWEEN 



   TCanvas *c11 = new TCanvas("c11", "y_axis vs log of paper energies", 200, 10, 700, 500); 

   for (int i = 0; i < num_interpl; i++) 
     {
      y_axis[0] = y_axis_bottom_interpl[i];
      y_axis[1] = y_axis_middle_interpl[i]; 
      y_axis[2] = y_axis_top_interpl[i];
      
      x_axis[0] = x_axis_bottom_interpl[i]; 
      x_axis[1] = x_axis_middle_interpl[i]; 
      x_axis[2] = x_axis_top_interpl[i];
  
      TGraph *gyaxis = new TGraph(3, paper_energy, y_axis);
      gyaxis -> Draw("AP*");
      gyaxis -> SetMarkerColor(6);
      gyaxis -> GetXaxis() -> SetTitle("log10(Energy / eV)");
      gyaxis -> GetYaxis() -> SetTitle("y dsig/dy (cm^2) "); 
      
      TSpline3 *sp_yaxis = new TSpline3("sp_yaxis", gyaxis); 
      sp_yaxis -> Draw("same"); 
      
      TGraph *gxaxis = new TGraph(3, paper_energy, x_axis); 
      
      TSpline3 *sp_xaxis = new TSpline3("sp_xaxis", gxaxis);
       
      for (int j = 0; j < num_energy_interpl; ++j)
        {
         energy_interpl[j] = ((double) j * (21 - 15) / (num_energy_interpl - 1)) + 15; 
         y_axis_interpl[j][i] = sp_yaxis -> Eval(energy_interpl[j]);
         x_axis_interpl[j][i] = sp_xaxis -> Eval(energy_interpl[j]);    
        }

     }

    c11 -> SaveAs("plots/y_axis_vs_energy.jpg");
    c11 -> Update(); 
    c11 -> Modified(); 



   // ******************************************************************************************************************************************************************


 // get sigma in cm^2 for top, middle and bottom curves (as a check)



 for ( int i = 0; i < int(sizeof(x_axis_middle_interpl)/sizeof(x_axis_middle_interpl[0])) - 1; i++ )
   {
     step_sum_middle[i] = (x_axis_middle_interpl[i+1] - x_axis_middle_interpl[i]) * (((y_axis_middle_interpl[i+1] + y_axis_middle_interpl[i]))/2);
     sigma_middle = sigma_middle + step_sum_middle[i]; 

   } 

 cout << "sigma_middle" << endl;
 cout << sigma_middle << endl;


 for ( int i = 0; i < int(sizeof(x_axis_top_interpl)/sizeof(x_axis_top_interpl[0])) - 1; i++ )
   {
     step_sum_top[i] = (x_axis_top_interpl[i+1] - x_axis_top_interpl[i]) * (((y_axis_top_interpl[i+1] + y_axis_top_interpl[i]))/2);
     sigma_top = sigma_top + step_sum_top[i]; 

   } 

 cout << "sigma_top" << endl;
 cout << sigma_top << endl;


 for ( int i = 0; i < int(sizeof(x_axis_bottom_interpl)/sizeof(x_axis_bottom_interpl[0])) - 1; i++ )
   {
     step_sum_bottom[i] = (x_axis_bottom_interpl[i+1] - x_axis_bottom_interpl[i]) * (((y_axis_bottom_interpl[i+1] + y_axis_bottom_interpl[i]))/2);
     sigma_bottom = sigma_bottom + step_sum_bottom[i]; 

   }
 

 cout << "sigma_bottom" << endl;
 cout << sigma_bottom << endl;


 // ***********************************************************************************************************************************************************************


 TRandom r; // for picking dx, generate a number in interval [0,1] (0 is excluded) 
 TRandom2 s; // for picking y, need to check that these are good to use
 TRandom3 t; // for determining if tau will decay


// LOOP OVER TAUS

for (int tau = 0; tau < num_tau; tau++) 

{

 double energy = initial_energy; 

 // WHILE LOOP x < x_max 

 while (x <= x_max)
   {

     double log10energy = log10(energy);

     if (log10energy < 15 || log10energy > 21)
       {
	 cout << "DONE WITH THIS TAU. ENERGY OF TAU OUT OF THE RANGE 10^15 TO 10^21 EV." << endl; 
	 break;
       }
     else
       {
	 for (int i = 0; i < num_energy_interpl; i++)
	   {
	     energydiff[i] = abs(log10energy - energy_interpl[i]);
	   }
	 double energydiffmin = TMath::MinElement(num_energy_interpl, energydiff); 
  
	 for (int k = 0; k < num_energy_interpl; k++)
	   {
	     if (energydiff[k] == energydiffmin)
	       { 
		 energy_index = k; 
	       }
	   }
       }


     // got energy_index now use it to get correct x and y axes for integration to get updated sigma



     for ( int i = 0; i < num_interpl - 1; i++ )
       {
	 step_sum[i] = (x_axis_interpl[energy_index][i+1] - x_axis_interpl[energy_index][i]) * (((y_axis_interpl[energy_index][i+1] + y_axis_interpl[energy_index][i]))/2);
	 sigma = sigma + step_sum[i]; 
       } 


     // FIND THE CDF OF y

 
     for ( int i = 0; i < int(sizeof(step_sum)/sizeof(step_sum[0])); i++ )
       {
	 for ( int j = 0; j <= i; j++ )
	   { 
	     y_cdf[i] += step_sum[j]; 
	   }
       }

	

     // find max of y_cdf

     double ycdfmax = TMath::MaxElement(sizeof(y_cdf)/sizeof(y_cdf[0]), y_cdf);  // ycdfmax should equal sigma and it does 
     //cout << "ycdfmax " << endl; 
     //cout << ycdfmax << endl; 



     // divide out ycdfmax from each element of ycdf 

     for ( int i = 0; i < int(sizeof(y_cdf)/sizeof(y_cdf[0])); i++ )
       {
	 y_cdf[i] = y_cdf[i] / ycdfmax; 
       }



     // FIND INTERACTION LENGTH which will update after each interaction

     int_length = proton_mass / (rho * sigma); 
     int_length_m = int_length / (1e2); 



     // RANDOMLY PICK THE dx's TRAVELED BY THE TAU ALONG ITS PATH


     // pick random number between 0 and 1 and save it

     double rand_x = r.Rndm();
     //cout << "rand_x is " << endl; 
     //cout << rand_x << endl; 
     double dx = -int_length_m * log(1 - rand_x); 
     x = x + dx;  
     dx_vector.push_back( dx );
     
     //cout << " dx " << endl; 
     //cout << dx << endl; 

     // DETERMINE IF TAU WILL DECAY FOR THIS dx 
     
     double gamma = energy / tau_mass; 
     double rand_decay = t.Rndm();
     //cout << "rand_decay is " << endl;
     //cout << rand_decay << endl; 
     double decay_pdf = exp(-dx / (c * gamma * tau_lifetime));  
     if (rand_decay >= decay_pdf)
      {
     	 cout << "TAU DECAYS FOR THIS dx" << endl;
     	 break; 
      }
     
   
     // RANDOMLY PICK THE y (fractional energy lost by the tau) for this dx 
     // first find random number 
     double rand_y = s.Rndm(); 
     //cout << "rand_y is " << endl; 
     //cout << rand_y << endl; 

     // use cdf of y to randomly pick y 
     
     
     if (rand_y < y_cdf[0]) // case where random number less than min ycdf 
       {
	 ycdf2 = y_cdf[0];
	 ycdf1 = 0; 
	 cout << "rand_y less than min y_cdf" << endl; 
	 //cout << ycdf2 << endl; 
	 //cout << ycdf1 << endl; 
	 y2 = x_axis_interpl[energy_index][0]; 
	 y1 = 0.0; 
	 double slope = (ycdf2 - ycdf1) / (y2 - y1); 
	 double intercept = ycdf2 - (slope * y2); 
	 y = (rand_y - intercept) / slope;
       }
     else if (rand_y == y_cdf[0])  // case where random number equal to min ycdf 
       { 
	 cout << "rand_y equal to min y_cdf" << endl; 
	 y = x_axis_interpl[energy_index][0]; 
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
	  y2 = x_axis_interpl[energy_index][top_index];
	  y1 = x_axis_interpl[energy_index][top_index - 1]; 
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
     //cout << "y now is" << endl;
     //cout << y_notlog << endl; 
	 
     num_int = num_int + 1;   // to get number of interactions


     //cout << "energy before minusing this y:" << endl; 
     //cout << energy << endl; 

     
     // find change in energy
     dE = -y_notlog * energy; 
     
     // push back to dE array
     dE_vector.push_back( dE ); 

     // update energy
     energy = energy + dE;

     // push back to energy array
     E_vector.push_back( energy ); 
     


     //cout << "energy after minusing this y:" << endl; 
     //cout << energy << endl; 

      
   }  // end of while loop 


// fill vector containing x's 

x_vector.push_back( x ); 


cout << "number of interactions of this tau:" << endl; 
cout << num_int << endl; 
cout << "total distance traveled by this tau: " << endl; 
cout << x << " meter" << endl; 


}  // end of for loop over taus 


 //declare array for dEdx and energy for TGraph plotting later 
 double dEdx_array[num_int]; 
 double E_array[num_int]; 
 double dx_array[num_int]; 

 // Fill array dEdx_array with log10 of -dE/dx and E_array with log10 of E_vector

 for (int i = 0; i < num_int; i++)
   {
     dEdx_array[i] = ( log10(-dE_vector[i] / dx_vector[i]) );
     E_array[i] = log10(E_vector[i]); // turning into array for plotting
     dx_array[i] = dx_vector[i];  // turning into array for plotting 
 }



 // FIND min and max of dx_array so that I can make a new array of dx's that increments regularly from the min to the max 

 double dx_max = TMath::MaxElement(num_int, dx_array);
 double dx_min = TMath::MinElement(num_int, dx_array); 
 double dx_range = dx_max - dx_min; 
 int num_forpdf = 500;  

// cout << "min and max of dx_array" << endl; 
// cout << dx_min << "     " << dx_max << endl; 


 // make the new array called dx_forpdf where store dx's from min to max regularly spaced 

 double dx_forpdf[num_forpdf]; 
 double pdf_dx[num_forpdf]; 
 double step_pdf[num_forpdf]; 
 double pdf_integral = 0; 

 for (int i = 0; i < num_forpdf; i++)
   {
     dx_forpdf[i] = ((double) i * (dx_range) / (num_forpdf - 1)) + dx_min; 
     pdf_dx[i] = exp(-dx_forpdf[i] / int_length_m);
   }
   
 // integrate the pdf of dx so can divide out later for normalization 

 for (int i = 0; i < num_forpdf-1; i++)
   {
     step_pdf[i] = ( dx_forpdf[i+1] - dx_forpdf[i] ) * (( pdf_dx[i+1] + pdf_dx[i] ) / 2 );
     pdf_integral = pdf_integral + step_pdf[i]; 
   }

 // divide out integral of pdf of dx to get normalized pdf 

 for (int i = 0; i < num_forpdf; i++)
   {
     pdf_dx[i] = pdf_dx[i] / pdf_integral; 
   }

 
 // Make hist of dx's 
 {
   TCanvas *c3 = new TCanvas("c3", "hist_dx", 200, 10, 700, 500);

   TGraph *gpdf_dx = new TGraph(num_forpdf, dx_forpdf, pdf_dx); // plotting normalized pdf of dx here 
   gpdf_dx -> Draw("AP*");
   gpdf_dx -> SetMarkerColor(6); 

   TH1D *hdx_array = new TH1D("hist_dx", "; dx; counts", 1000, 0, 100); 

   for (int i = 0; i < int(num_int); i++) 
     { 
       hdx_array->Fill(dx_array[i]);
     } 
   double dx_integral = hdx_array -> Integral("width"); 
   //cout << "dx_array integral is" << endl; 
   //cout << dx_integral << endl; 
   hdx_array -> Scale(1.0 * (1.0 / dx_integral));  // plotting normalized hist of picked dx here 
   hdx_array -> Draw("same"); 
   gpdf_dx  -> SetTitle("Histogram of dx");
   gpdf_dx -> GetXaxis() -> SetTitle("dx (meter)");
   gpdf_dx -> GetYaxis() -> SetTitle("counts"); 
   hdx_array -> SetFillColor(kBlack); 
 
   c3 -> SaveAs("plots/dx_hist.jpg"); 
   c3 -> Update();
   c3 -> Modified(); 
 }

 // Make hist of y's 
 {

   for (int i = 0; i < top_num_rows; i++)
    {
     y_axis_top[i] =  y_axis_top[i] / sigma_top; 
    } 
       
   TCanvas *c4 = new TCanvas("c4", "hist_y", 200, 10, 700, 500);

   TGraph *norm_fig8 = new TGraph(top_num_rows, x_axis_top, y_axis_top); 
   norm_fig8 -> Draw("AC*"); 
   norm_fig8 -> SetMarkerColor(6); 
 
   TH1D *hy_vector = new TH1D("hist_y", "; y; counts",100,-4,0); 
 
   for (int i = 0; i < int(y_vector.size()); i++) 
     { 
       hy_vector -> Fill(y_vector[i]);
     }
 
   double y_integral = hy_vector -> Integral("width"); 
   hy_vector -> Scale(1.0 * (1.0 / y_integral));  
   hy_vector -> Draw("same"); 
   norm_fig8 -> SetTitle("Histogram of y");
   norm_fig8 -> GetXaxis() -> SetTitle("log10 y");
   norm_fig8 -> GetYaxis() -> SetTitle("Number of Events"); 
   hy_vector -> SetFillColor(kBlack); 
   
   c4 -> SaveAs("plots/y_hist.jpg"); 
   c4 -> Update();
   c4 -> Modified(); 
 }


// Make a hist of x's that is the total distance traveled by each tau

{
 TCanvas *c12 = new TCanvas("c12", "hist_x", 200, 10, 700, 500);
 TH1D *hx_vector = new TH1D("hist_x", "; x; counts", 10, 0, 50); 
 
 for (int i = 0; i < int(x_vector.size()); i++)
   {
    hx_vector -> Fill(x_vector[i]);
   }
 hx_vector -> Draw(""); 
 hx_vector -> SetTitle("Histogram of x");
 hx_vector -> GetXaxis() -> SetTitle("Total distance traveled (m)");
 hx_vector -> GetYaxis() -> SetTitle("Number of taus"); 

 c12 -> SaveAs("plots/x_hist.jpg");
 c12 -> Update();
 c12 -> Modified(); 
}
 
 // Make plot of -dE/dx vs energy 


{
   TCanvas *c5 = new TCanvas("c5","dEdx vs E", 200, 10, 700, 500);
   TGraph *dEdxvsE = new TGraph(num_int,E_array,dEdx_array);
   dEdxvsE -> SetTitle("log10(-dE/dx) vs log10(E)");
   dEdxvsE -> GetXaxis()->SetTitle("log10( Energy / eV )");
   dEdxvsE -> GetYaxis()->SetTitle("log10( (-dE/dx) / (eV/m) )");
   dEdxvsE -> Draw("A*");
   c5 -> SaveAs("plots/dEdxvsE.jpg");
   c5 -> Update();
   c5 -> Modified();
}


 // Make plot of y cdf vs log10 of y

{
   TCanvas *c2 = new TCanvas("c2","y_cdf", 200, 10, 700, 500);
   TGraph *ycdf = new TGraph((num_interpl-1),x_axis_top_interpl,y_cdf);
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


 
