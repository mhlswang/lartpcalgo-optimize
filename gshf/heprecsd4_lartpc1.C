void heprecsd4_lartpc1()
{
    gROOT->Reset();

    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    //gStyle->SetPadLeftMargin(.1);
    //gStyle->SetPadRightMargin(.009);
    gStyle->SetPadTopMargin(.05);
    //gStyle->SetPadBottomMargin(.1);

    gStyle->SetOptTitle(0);
    //gStyle->SetTitleW(0.5);
    //gStyle->SetTitleH(0.08);

    // ... SetOptStat(int,ovr,undr,rms,mean,ent,name),default:0001111
    gStyle->SetOptStat(111110);
    //gStyle->SetOptFit(0101);
    gStyle->SetStatW(.3);
    gStyle->SetStatH(0.18);

    gStyle->SetLabelSize(0.045);
    gStyle->SetTitleSize(0.045);

    // ... define histograms
    TH1F* h1   = new TH1F("dsim1" ,"Default ROOT fitter", 100,-5,5);
    TH1F* h2   = new TH1F("dsim2" ,"New Custom fitter", 100,-5,5);
    TH1F* h3   = new TH1F("drec" ,"myrecTick-recTick", 100,-5,5);

    int ievt,im,ipk;
    double simtck,rectck,rms,mytck2,mysigma;

    ifstream fin("result.txt");
    for(int n=0;n<113386;n++){
      fin>>ievt>>im>>ipk>>simtck>>rectck>>rms>>mytck2>>mysigma;
      //if(ipk==0 && im==0){
      h1->Fill(rectck-simtck);
      h2->Fill(mytck2-simtck);
      h3->Fill(mytck2-rectck);
    }

    TCanvas* canv1 = new TCanvas("c1","plot",0,0,800,400);
    canv1->Divide(2,1);

    h1->GetXaxis()->SetNdivisions(505,kTRUE);
    h1->GetXaxis()->SetTitle("Reconstructed time - True time [#time bins]");
    h1->SetMaximum(h1->GetMaximum()*1.2);
    //h1->SetMaximum(1e6);
    h1->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h1->SetMarkerStyle(4);
    h1->SetMarkerColor(kRed);
    h1->SetMarkerSize(.5);
    h1->SetStats(0);
    canv1->cd(1);h1->Draw("P");
    //gPad->SetLogy();
    h2->GetXaxis()->SetNdivisions(505,kTRUE);
    h2->GetXaxis()->SetTitle("Reconstructed time - True time [#time bins]");
    h2->SetMaximum(h2->GetMaximum()*1.2);
    //h2->SetMaximum(1e6);
    h2->SetLineColor(kBlue);
    h2->SetLineWidth(2);
    h2->SetMarkerStyle(5);
    h2->SetMarkerColor(kBlue);
    h2->SetMarkerSize(1.1);
    //canv1->cd(2);h2->Draw();
    h2->Draw("P,Same");
    //gPad->SetLogy();

   TLegend* leg = new TLegend(0.11,0.75,0.52,0.9);
   //leg->SetHeader("The Legend Title");
   leg->AddEntry(h1,"","p");
   leg->AddEntry(h2,"","p");
   leg->Draw();



    h3->GetXaxis()->SetNdivisions(505,kTRUE);
    //h3->SetMaximum(h3->GetMaximum()*1.15);
    h3->SetMaximum(h3->GetMaximum()*10);
    h3->GetXaxis()->SetTitle("ROOT - Custom recon time [#time bins]");
    h3->SetLineColor(kBlack);
    h3->SetLineWidth(2);
    canv1->cd(2);h3->Draw();
    gPad->SetLogy();

    canv1->cd(0);
    TPad *txtpad = new TPad("txtpad","txtpad",0,0,1,1);
    txtpad->SetFillStyle(4000);
    txtpad->Draw();
    txtpad->cd();
    TLatex ltx;
    ltx.SetTextAlign(12);
    ltx.SetTextSize(0.048);
    //ltx.DrawLatex(.046,.88,"Default ROOT-based fitter");
    //ltx.DrawLatex(.385,.88,"New Custom fitter");
    ltx.DrawLatex(.57,.88,"Default vs. New fitter");


    canv1->Print("heprecosd4-lartpc1.eps");
}
        //fprintf(fout,"%d %d %d %lf %lf %lf %lf %lf\n",n,i,j,simtck,rectck,rms,mytck2,mysigma);
