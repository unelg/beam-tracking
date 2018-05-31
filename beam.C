#include "Math/GSLRndmEngines.h"
#include <fstream>
#include <iostream>
#include "TMath.h"
using namespace TMath;

void beam(){
     TH2F *h2=new TH2F( " h2 " , "x-px" , 1000 , -50 , 50 , 1000 , -30 , 30 ) ;
     TH2F *h3=new TH2F( " h3 " , "y-py" , 1000 , -50 , 50 , 1000 , -30 , 30 ) ;
     TH1F *h4=new TH1F( " h4 " , "x-number of particles" , 100 , -100 , 100 ) ;
     TH1F *h5=new TH1F( " h5 " , "y-number of particles" , 100 , -100 , 100 ) ;
     TH2F *h6=new TH2F( " h6 " , "driftx-px" , 1000 , -400 , 400 , 1000 , -8 , 8 ) ;
     TH2F *h7=new TH2F( " h2 " , "quad x-px" , 1000 , -50 , 50 , 1000 , -30 , 30 ) ;
     TH2F *h8=new TH2F( " h3 " , "quad y-py" , 1000 , -50 , 50 , 1000 , -30 , 30 ) ;
     TH1F *h9=new TH1F( " h4 " , "quad x-number of particles" , 100 , -100 , 100 ) ;
     TH1F *h10=new TH1F( " h5 " , "quad y-number of particles" , 100 , -100 , 100 ) ;
     TH2F *h11=new TH2F( " h6 " , "solenoid x-px" , 1000 , -400 , 400 , 1000 , -8 , 8 ) ;
     TGraph *gr = new TGraph();
     TGraph *gr2 = new TGraph();
     TGraph *gr3 = new TGraph();
     TGraph *gr4 = new TGraph();
     TGraph *gr5 = new TGraph();
     TGraph *gr6 = new TGraph();
     TGraph *gr7 = new TGraph();
     TGraph *gr8 = new TGraph();
     TGraph *gr9 = new TGraph();
     TGraph *gr10 = new TGraph();
     TGraph *gr11 = new TGraph();
     TGraph *gr12 = new TGraph();

//-----initial conditions   
    double restmass= 938*1000; //KeV
    double kenergy= 20; //KeV
    double x_max = 10; //milimeter
    double px = 0.6; //miliradyan
    double y_max = 10;
    double py = 0.6;
    double nemitx= 0.01;//(mm*mrad*pi)
    double nemity= 0.01;//(mm*mrad*pi)

//----- derived quantities
    double rgamma= 1+(kenergy/restmass);
    double rbeta= sqrt(1-(1/(rgamma*rgamma)));
    double gemitx= nemitx/(rgamma*rbeta);
    double gemity= nemity/(rgamma*rbeta);;
    double tbetax= (x_max*x_max)/gemitx;
    double tbetay= (y_max*y_max)/gemity;
    double talphax = px/sqrt(gemitx/tbetax);
    double talphay = py/sqrt(gemity/tbetay);
    double tgammax = ((1+(talphax*talphax))/tbetax);
    double tgammay = ((1+(talphay*talphay))/tbetay);
    double emitx ;
    double emity ;
    double Nemitx ;
    double Nemity ;
    double rhox = 0.5*atan((2*talphax)/(tgammax-tbetax));
    double rhoy = 0.5*atan((2*talphay)/(tgammay-tbetay));
    //double rhox = -0.9;
    //double rhoy = -0.9;
    cout<<"alpha: "<<talphax<<"beta: "<<tbetax<<"gamma: "<<tgammax<<"gemit: "<<gemitx<<"rgamma: "<<rgamma<<"rbeta "<<rbeta<<endl;

    double px_max= sqrt(tgammax*gemitx);
    double py_max= sqrt(tgammay*gemity);
    cout<<"pxmax: "<<px_max<<endl;
    double rx, rpx;
    double ry, rpy;
    double nrx, nrpx;
    double nry, nrpy;
   
    const int N=5000;
    double particle [N][4];
    cout<<rhox<<" "<<rhoy<<endl;
    TRandom2 *r2=new TRandom2();
    ofstream myfile2;
    myfile2.open ("beambegin.txt");
    ROOT::Math::GSLRandomEngine *myengine = new ROOT::Math::GSLRandomEngine();
    myengine->Initialize();
    cout<<"gammax="<<tgammax<<" betax="<<tbetax<<" alphax="<< talphax <<" geometricemmitance= "<<gemitx<<" rhox:"<<rhox<<" rhoy: "<<rhoy << endl;
    for (int i = 0; i < N; ++i){
        
        myengine->Gaussian2D(x_max/2, px_max, rhox, rx, rpx);
        myengine->Gaussian2D(y_max/2, py_max, rhoy, ry, rpy);

        //cout<<rx<<" "<<rpx<<endl;
        emitx = (tgammax*rx*rx) + (2*talphax*rx*rpx) + (tbetax*rpx*rpx);
        emity = (tgammay*ry*ry) + (2*talphay*ry*rpy) + (tbetay*rpy*rpy);
        
        particle[i][0]= rx;
        particle[i][1]= rpx;
        particle[i][2]= ry;
        particle[i][3]= rpy;
//cout<<particle[i][0]<<" "<<particle [i][1]<<endl;
        myfile2 <<rx<<" "<<rpx<<" "<<ry<<" "<<rpy<<endl;
        gr2-> SetPoint(i,rx,rpx);
        gr3-> SetPoint(i,ry,rpy);
        gr4-> SetPoint(i,rx,ry);
            h2->Fill(rx,rpx);
            h3->Fill(ry,rpy);
            h4->Fill(rx);
            h5->Fill(ry);
    }
     /*
    double beamRad;
    for (int i = 0; i < N; ++i){
     double rx = r2->Gaus(0,x_max);
     double rpx = r2->Gaus(0,px_max);
     double ry = r2->Gaus(0,y_max);
     double rpy = r2->Gaus(0,py_max);
     
     //cout<<rx<<" "<<ry<<endl;
     emitx = (tgammax*rx*rx) + (2*talphax*rx*rpx) + (tbetax*rpx*rpx);
     emity = (tgammay*ry*ry) + (2*talphay*ry*rpy) + (tbetay*rpy*rpy);
     beamRad=sqrt(rx*rx + ry*ry);
     
     if (gemitx >= emitx  && gemity >= emity && beamRad<=x_max){
     //cout<<emitx<<" "<<emity<<endl;
     
     particle[i][0] = rx;
     particle[i][1] = rpx;
     particle[i][2] = ry;
     particle[i][3] = rpy;
     //cout<<"x: "<<particle [i][0]<<"px: "<<particle [i][1]<<"y: "<<particle [i][2]<<"py: "<<particle [i][3]<<endl;
     myfile2 <<rx<<" "<<rpx<<" "<<ry<<" "<<rpy<<endl;
     //cout<<"firstrx: "<<particle [i][0] <<" firstrpx: "<<particle [i][1]<<endl;
     gr2-> SetPoint(i,rx,rpx);
     gr3-> SetPoint(i,ry,rpy);
     gr4-> SetPoint(i,rx,ry);
     h2->Fill(rx,rpx);
     h3->Fill(ry,rpy);
     h4->Fill(rx);
     h5->Fill(ry);} else --i;
     }*/
    myfile2.close();
//----------------------beam creation ends here

    double l=312;// this is drift distance in mm
    ofstream myfile;
    myfile.open ("beam.txt");
    double ntbetax;
    double ntbetay;
    double ntgammay;
    double ntgammax;
    double ntalphay;
    double ntalphax;
    
    for (int n = 0; n<N; n++){
        //cout<<"old rx: "<<particle [n][2]<<"old rpx: "<<particle [n][2]<<endl;
         rx = particle[n][0];
        rpx = particle[n][1];
         ry = particle[n][2];
        rpy = particle[n][3];
        
        nrx = rx + (l*rpx);
        ntbetax = tbetax -2*l*talphax+l*l*tgammax;
        ntalphax = talphax- l*tgammax;
        ntgammax=tgammax;
        nry = ry + (l*rpy);
        ntbetay = tbetay -2*l*talphay+l*l*tgammay;
        ntalphay = talphay- l*tgammay;
        ntgammay=tgammay;
      
        particle[n][0] = nrx;
        particle[n][1] = rpx;
        particle[n][2] = nry;
        particle[n][3] = rpy;
        //cout<<"newry: "<<particle [n][2]<<"newrpy: "<<particle [n][3]<<endl;
        Nemitx = (tgammax*nrx*nrx) + (2*talphax*nrx*rpx) + (tbetax*rpx*rpx);
        //Nemity = (tgammay*nry*nry) + (2*talphay*ry*rpy) + (tbetay*rpy*rpy);
        //cout<<" "<<Nemitx<<endl;
        myfile <<nrx<<" "<<rpx<<" "<<nry<<" "<<rpy<<endl;
        gr-> SetPoint(n,nrx,rpx);
        gr5-> SetPoint(n,nry,rpy);
        gr6-> SetPoint(n,nrx,nry);
        h6->Fill(nrx,rpx);
    }
//------------------------- drift in free space ends here

    myfile.close();
    double L=0.5;//(meter) ----------QUADRUPOLE
    double g=5 ;//quadrupole gradient(T/m)
    double K = (0.299*g)/(rbeta*kenergy); // quadrupole strength(mm^-2)
    // Quadrupole transfer focusing in x defocusing in y
    double nqtbetax;
    double nqtalphax;
    double nqtgammax;
    ofstream myfile3;
    myfile3.open ("beamquad.txt");
    for (int n = 0; n<N; n++){
        //cout<<"old rx: "<<particle [n][2]<<"old rpx: "<<particle [n][2]<<endl;
        rx = particle[n][0];
       rpx = particle[n][1];
        ry = particle[n][2];
       rpy = particle[n][3];
        
        nrx = rx*cos(sqrt(K)*L) + rpx*(1/sqrt(K))*sin(sqrt(K)*L);
        
        nqtbetax = cos(sqrt(K)*L)*cos(sqrt(K)*L)*ntbetax -2*(1/sqrt(K))*sin(sqrt(K)*L)*cos(sqrt(K)*L)*ntalphax+ (1/K)*sin(sqrt(K)*L)*sin(sqrt(K)*L)*ntgammax;
        
        nqtalphax = -sqrt(K)*sin(sqrt(K)*L)*cos(sqrt(K)*L)*ntbetax-sin(sqrt(K)*L)*sin(sqrt(K)*L)*ntalphax-cos(sqrt(K)*L)*cos(sqrt(K)*L)*ntalphax-(1/sqrt(K))*sin(sqrt(K)*L)*cos(sqrt(K)*L)*ntgammax;
        nqtgammax = K*sin(sqrt(K)*L)*sin(sqrt(K)*L)*ntbetax+2*sin(sqrt(K)*L)*cos(sqrt(K)*L)*ntalphax+cos(sqrt(K)*L)*cos(sqrt(K)*L)*ntgammax;
        
        nrpx = rx*(-1*sqrt(K))*sin(sqrt(K)*L) + rpx*cos(sqrt(K)*L);
        
        nry = ry*cosh(sqrt(K)*L) + rpy*(1/sqrt(K))*sinh(sqrt(K)*L);
        
        nrpy = ry*(sqrt(K))*sinh(sqrt(K)*L) + rpy*cosh(sqrt(K)*L);
        cout<<nrx<<"ry: "<<rx<<" rpy: "<<rpx<<"  "<<nrpx<<endl;
        
        particle[n][0] = nrx;
        particle[n][1] = nrpx;
        particle[n][2] = nry;
        particle[n][3] = nrpy;
        //cout<<"newry: "<<particle [n][2]<<"newrpy: "<<particle [n][3]<<endl;
        Nemitx = (tgammax*nrx*nrx) + (2*talphax*nrx*rpx) + (tbetax*rpx*rpx);
        //Nemity = (tgammay*nry*nry) + (2*talphay*ry*rpy) + (tbetay*rpy*rpy);
        //cout<<" "<<Nemitx<<endl;
        myfile3 <<nrx<<"  "<<nrpx<<"  "<<nry<<"  "<<nrpy<<endl;
        gr7-> SetPoint(n,nrx,nrpx);
        gr8-> SetPoint(n,nry,nrpy);
        gr9-> SetPoint(n,nrx,nry);
        h7->Fill(nrx,nrpx);
        h8->Fill(nry,nrpy);
        h9->Fill(nrx,n);
        h10->Fill(nry,n);
    }
//------------------------QUAD ends here
    myfile3.close();


//--------------------------------- solenoid transfer entry,exit and body all together
    double Ls=0.5; // solenoid length in M
    double k=0.9; //solenoid strength in
    ofstream myfile4;
    myfile4.open ("beamsolenoid.txt");
    for (int n = 0; n<N; n++){
        //cout<<"old rx: "<<particle [n][2]<<"old rpx: "<<particle [n][2]<<endl;
        rx = particle[n][0];
       rpx = particle[n][1];
        ry = particle[n][2];
       rpy = particle[n][3];
        
        nrx =  rx*cos(k*Ls)*cos(k*Ls)   +rpx*(1/k)*sin(k*Ls)*cos(k*Ls)+ry*sin(k*Ls)*cos(k*Ls) +rpy*(1/k)*sin(k*Ls)*sin(k*Ls);
       nrpx = -rx*sin(k*Ls)*cos(k*Ls)*k +rpx*cos(k*Ls)*cos(k*Ls)-ry*k*sin(k*Ls)*sin(k*Ls)     +rpy*sin(k*Ls)*cos(k*Ls);
        nry = -rx*sin(k*Ls)*cos(k*Ls)   -rpx*(1/k)*sin(k*Ls)*sin(k*Ls)+ry*cos(k*Ls)*cos(k*Ls) +rpy*(1/k)*sin(k*Ls)*cos(k*Ls);
        nrpy = rx*sin(k*Ls)*sin(k*Ls)*k -rpx*cos(k*Ls)*sin(k*Ls)-ry*k*sin(k*Ls)*cos(k*Ls)     +rpy*cos(k*Ls)*cos(k*Ls);
        
        particle[n][0] = nrx;
        particle[n][1] = nrpx;
        particle[n][2] = nry;
        particle[n][3] = nrpy;
        //cout<<"newry: "<<particle [n][2]<<"newrpy: "<<particle [n][3]<<endl;
        Nemitx = (tgammax*nrx*nrx) + (2*talphax*nrx*rpx) + (tbetax*rpx*rpx);
        //Nemity = (tgammay*nry*nry) + (2*talphay*ry*rpy) + (tbetay*rpy*rpy);
        //cout<<" "<<Nemitx<<endl;
        myfile4 <<nrx<<" "<<nrpx<<""<<nry<<""<<nrpy<<endl;
        gr10-> SetPoint(n,nrx,nrpx);
        gr11-> SetPoint(n,nry,nrpy);
        gr12-> SetPoint(n,nrx,nry);
        h10->Fill(nrx,nrpx);
    }
//---------------------------solenoid ends.    
    myfile4.close();
    
    
//-----------------------------NOW DRAW ALL    
    TCanvas *c1=new TCanvas();
    c1->Divide(4,4);
 
    //c1->Update();
    c1->cd(1);
    gr2->GetHistogram()->GetXaxis()->SetTitle("px");
    gr2->GetHistogram()->GetYaxis()->SetTitle("x");
    gr2->Draw("AP");
    //h2->Draw();
    c1->cd(2);
    gStyle->SetOptFit();
    gr4->GetHistogram()->GetXaxis()->SetTitle("x");
    gr4->GetHistogram()->GetYaxis()->SetTitle("y");
    gr4->Draw("AP");
    //h3->Draw();
    c1->cd(3);
    h4->Draw();
    c1->cd(4);
    h5->Draw();
    c1->cd(5);
    h2->Draw();
    //gr6->Draw("AP");
    c1->cd(6);
    h3->Draw();
    cout<<h2->GetCorrelationFactor()<<endl;
    cout<<h3->GetCorrelationFactor()<<endl;
    c1->cd(7);
    h7->Draw();
    c1->cd(8);
    h8->Draw();
    c1->cd(9);
    h9->Draw();
    c1->cd(10);
    h10->Draw();
    c1->cd(11);
    gStyle->SetOptFit();
    gr7->GetHistogram()->GetXaxis()->SetTitle("x");
    gr7->GetHistogram()->GetYaxis()->SetTitle("px");
    gr7->Draw("AP");
    c1->cd(12);
    gStyle->SetOptFit();
    gr8->GetHistogram()->GetXaxis()->SetTitle("y");
    gr8->GetHistogram()->GetYaxis()->SetTitle("py");
    gr8->Draw("AP");
    c1->cd(13);
    gStyle->SetOptFit();
    gr9->GetHistogram()->GetXaxis()->SetTitle("x");
    gr9->GetHistogram()->GetYaxis()->SetTitle("y");
    gr9->Draw("AP");
    c1->cd(14);
    gStyle->SetOptFit();
    gr10->GetHistogram()->GetXaxis()->SetTitle("x");
    gr10->GetHistogram()->GetYaxis()->SetTitle("px");
    gr10->Draw("AP");
    c1->cd(15);
    gStyle->SetOptFit();
    gr11->GetHistogram()->GetXaxis()->SetTitle("y");
    gr11->GetHistogram()->GetYaxis()->SetTitle("py");
    gr11->Draw("AP");
    c1->cd(16);
    gStyle->SetOptFit();
    gr12->GetHistogram()->GetXaxis()->SetTitle("x");
    gr12->GetHistogram()->GetYaxis()->SetTitle("y");
    gr12->Draw("AP");
    
    }
