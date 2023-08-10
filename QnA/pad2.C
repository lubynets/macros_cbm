void CanvasPartition(const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05,
                     Float_t vSpacing = 0., Float_t hSpacing = 0.);
void pad2()
{
   gStyle->SetOptStat(0);
   gROOT->Macro( "/home/oleksii/cbmdir/flow_drawing_tools/example/style_multipad.cc" );
   TCanvas *C = (TCanvas*) gROOT->FindObject("C");
   if (C) delete C;
   C = new TCanvas("C","canvas",4500,3000);
   C->SetFillStyle(4000);
   // Number of PADS
   const Int_t Nx = 3;
   const Int_t Ny = 2;
   // Margins
   Float_t lMargin = 0.05;
   Float_t rMargin = 0.01;
   Float_t bMargin = 0.05;
   Float_t tMargin = 0.01;
   // Canvas setup
   CanvasPartition(Nx,Ny,lMargin,rMargin,bMargin,tMargin);
   // Dummy histogram.

   // TICKS X Axis

   TFile* fileIn = TFile::Open("/home/oleksii/cbmdir/flow_drawing_tools/build/res.psd.root");
   HeapPicture* hp = fileIn->Get<HeapPicture>("abcd");

   TPad *pad[Nx][Ny];
   for (Int_t i=0;i<Nx;i++) {
      for (Int_t j=0;j<Ny;j++) {
         C->cd(0);
         // Get the pads previously created.
         char pname[16];
         sprintf(pname,"pad_%i_%i",i,j);
         pad[i][j] = (TPad*) gROOT->FindObject(pname);
         pad[i][j]->Draw();
         pad[i][j]->SetFillStyle(4000);
         pad[i][j]->SetFrameFillStyle(4000);
//          pad[i][j]->cd();
         // Size factors
         Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
         Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();

//          std::cout << i << "\t" << j << "\n" << pad[i][j]->GetAbsWNDC() << "\t" << pad[i][j]->GetAbsHNDC() << "\n\n";

         char hname[16];
         sprintf(hname,"h_%i_%i",i,j);
//          TH1F *hFrame = (TH1F*) h->Clone(hname);
//
//          // TICKS Y Axis
//          hFrame->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
//
//          hFrame->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
//
//          hFrame->Draw();

         HeapPicture* hpclone = (HeapPicture*)hp->Clone();
         hpclone->DrawPad(pad[i][j]);
      }
   }
   C->SaveAs("C.png", "png");
   C->Print("C.pdf", "pdf");
}
void CanvasPartition(const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin,
                     Float_t vSpacing, Float_t hSpacing)
{
   // Setup Pad layout:
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;
   for (Int_t i=0;i<Nx;i++) {
      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
      for (Int_t j=0;j<Ny;j++) {
         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
         pad->Draw();
      }
   }
}
