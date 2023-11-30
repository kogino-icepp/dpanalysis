//FFTをやるために作ったコードをこちらに退避


/*int startbin = 2970;
            double jfmin = min(Freq1[0],Freq1[nbin-1]);
            double jfmax = max(Freq1[0],Freq1[nbin-1]);
            TH1* thist = new TH1D("thist",";Freq[GHz];Spec",nbin,jfmin,jfmax);
            TH1* onhist = new TH1D("onhist",";Freq[GHz];Spec",nbin,0,nbin);
            TH1* delhist = new TH1D("delhist","test;Freq[GHz];Spec",nbin,jfmin,jfmax);
            TH1* fhist = new TH1D("fhist","MAG;1/Freq[ns];",nbin,0,nbin/2.5);//棘抜く前のFFT
            TH1* tfhist = new TH1D("tfhist","PH;1/Freq[ns];",nbin,0,nbin/2.5);//棘抜いた後のFFT
            TH1* sighist = new TH1D("sighist","sig;Freq[GHz];Spec",nbin,jfmin,jfmax);
            TH1* sighist30 = new TH1D("sighist30","sig;Freq[GHz];Sig",30,Freq1[startbin],Freq1[startbin+30]);
            TH1* sinhist = new TH1D("sinhist","sin;t;f",50,0,50);
            TH1* sinfft = new TH1D("sinfft","sinfft;;",50,0,50);
            TH1* sinrfft = new TH1D("sinrfft","sinr;t;f",50,0,50);
            TH1* parhist = new TH1D("parhsit","partition;Freq[GHz];Spec",30,Freq1[startbin],Freq1[startbin+30]);
            TH1* fparhist = new TH1D("fparhist","fft;1/Freq[ns];Spec",30,0,nbin/2.5);
            TH1* fsighist30 = new TH1D("fsighist30","fft(signal);1/Freq[ns];Spec",30,0,nbin/2.5);
            TH1* rawhist = new TH1D("rawhist","",nbin,0,nbin);
            //30binごとに区切ってFFTを行うための場所
            rep(bin,30){
                parhist -> SetBinContent(bin,Cold[bin+startbin]/pow(10,13));
                sighist30 -> SetBinContent(bin,F_sig2(Freq1[bin],Freq1[1],10,0.5));
            }
            rep(bin,nbin)rawhist -> SetBinContent(bin,Cold[bin]);
            st.Hist(rawhist);
            rawhist -> Draw();
            st.Hist(parhist);
            parhist -> Draw();
            //30binだけのFFTデータをそれぞれに格納する
            double re_full[30],im_full[30],sre_full[30],sim_full[30];
            
            fparhist = parhist -> FFT(fparhist,"MAG");
            TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
            fft -> GetPointsComplex(re_full,im_full);
            fsighist30 = sighist30 -> FFT(fsighist30,"MAG");
            fft = TVirtualFFT::GetCurrentTransform();
            fft -> GetPointsComplex(sre_full,sim_full);
            st.Hist(fsighist30);
            fsighist30 -> Draw();
            double cre = 0;
            double cim = 0;//カットオフ周波数での値、後から調節する
            double afre[30],afim[30];
            TH1* wiehist = new TH1D("wiehist","wiener filter(no cutoff);1/Freq[ns];Spec",30,0,30);
            rep(bin,16){
                double deno = (sre_full[bin]*sre_full[bin]+sim_full[bin]*sim_full[bin]);
                double deno2 = (sre_full[bin]*sre_full[bin]+sim_full[bin]*sim_full[bin])+(cre*cre+cim*cim);
                double sre = sre_full[bin]/deno;
                double sim = -sim_full[bin]/deno;
                double sre2 = sre_full[bin]/deno2;
                double sim2 = -sim_full[bin]/deno2;
                double re = re_full[bin]*sre-im_full[bin]*sim;
                double im = re_full[bin]*sim+im_full[bin]*sre;
                double re2 = re_full[bin]*sre2-im_full[bin]*sim2;
                double im2 = re_full[bin]*sim2+im_full[bin]*sre2;
                wiehist -> SetBinContent(bin,sqrt(sre*sre+sim*sim));
            }
            st.Hist(wiehist);
            wiehist -> Draw();
            //TH1* testrev = nullptr;
            //testrev = wiehist -> FFT(testrev,"MAG");
            //st.Hist(testrev);
            //testrev -> Draw();
            /*double sre_full[50],sim_full[50];
            rep(bin,50)sinhist -> SetBinContent(bin,sin(2*TMath::Pi()*bin/10));
            sinfft = sinhist -> FFT(sinfft,"MAG");
            TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
            fft -> GetPointsComplex(sre_full,sim_full);
            TH1* backhist = new TH1D("","",50,0,50);
            int N = 50;
            TVirtualFFT *fft_back = TVirtualFFT::FFT(1,&N,"C2R M K");
            fft_back->SetPointsComplex(sre_full,sim_full);
            rep(bin,50)cout << bin << ": " <<  sre_full[bin] << " " << sim_full[bin] << endl;
            fft_back -> Transform();
            backhist = TH1::TransformHisto(fft_back,backhist,"Re");
            st.Hist(backhist);
            backhist -> Draw();*/
            //フィルタは前半だけにかける、後半は放置しても問題なく復元できそう(あと絶対値だけ問題にする)

            //TVirtualFFT::SetTransform(nullptr);
            /*
            rep(k,30)parhist -> SetBinContent(k,Cold[25940+k]/pow(10,13));
            st.Hist(parhist);
            parhist -> Draw();
            fhist = parhist -> FFT(fhist,"MAG");
            st.Hist(fhist);
            rep(bin,nbin){
                thist -> SetBinContent(bin,Cold[bin]/pow(10,13));
                sighist -> SetBinContent(bin,F_sig2(Freq1[bin],Freq1[1],10,0.5));
            }
            st.Hist(sighist);
            sighist -> Draw();

            //filter作り、カットオフ有無で二通り
            TH1* filhist = new TH1D("filhist","wiener filter[not cutoff];1/Freq[ns];",nbin,0,nbin/2.5);
            TH1* filhist2 = new TH1D("filhist2","wiener filter[cutoff];1/Freq[ns];",nbin,0,nbin/2.5);
            TH1* tfilhist = new TH1D("tfilhist","after filtered[not cutoff];",nbin,0,nbin/2.5);
            TH1* tfilhist2 = new TH1D("tfilhist2","after filtered[cutoff];",nbin,0,nbin/2.5);

            Double_t *re_full = new Double_t[nbin];//シグナルのFFTデータ
            Double_t *im_full = new Double_t[nbin];//上の虚部
            Double_t *tre_full = new Double_t[nbin];//元データのFFTデータ
            Double_t *tim_full = new Double_t[nbin];//上の虚部
            //シグナルと生データをそれぞれFFTにかける
            fhist = sighist -> FFT(fhist,"MAG");
            TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();//直近でFFTしたデータのあれこれを引っ張り出すポインタ
            fft -> GetPointsComplex(re_full,im_full);
            double cre = re_full[10];
            double cim = im_full[10];
            tfhist = thist -> FFT(tfhist,"MAG");
            TVirtualFFT *tfft = TVirtualFFT::GetCurrentTransform();
            tfft -> GetPointsComplex(tre_full,tim_full);
            

            rep(bin,nbin/2){
                double deno = (re_full[bin]*re_full[bin]+im_full[bin]*im_full[bin]);
                double deno2 = (re_full[bin]*re_full[bin]+im_full[bin]*im_full[bin])+(cre*cre+cim*cim);
                double re = re_full[bin]/deno;
                double im = -im_full[bin]/deno;
                double re2 = re_full[bin]/deno2;
                double im2 = -im_full[bin]/deno2;
                double tre = tre_full[bin]*re-tim_full[bin]*im;
                double tim = tre_full[bin]*im+tim_full[bin]*re;
                double tre2 = tre_full[bin]*re2-tim_full[bin]*im2;
                double tim2 = tre_full[bin]*im2+tim_full[bin]*re2;
                
                tfilhist -> SetBinContent(bin,sqrt(tre*tre+tim*tim));
                tfilhist2 -> SetBinContent(bin,sqrt(tre2*tre2+tim2*tim2));
                filhist -> SetBinContent(bin,sqrt(re*re+im*im));
                filhist2 -> SetBinContent(bin,sqrt(re2*re2+im2*im2));
                tre_full[bin] = tre2;
                tim_full[bin] = tim2;
            }
            st.Hist(filhist2);
            filhist2 -> Draw();
            st.Hist(tfilhist);
            tfilhist -> Draw();
            Int_t Nbin = 32767;
            TH1* backhist = new TH1D("backhist","back;Freq[GHz];Spec",nbin,jfmin,jfmax);
            TVirtualFFT *fft_back = TVirtualFFT::FFT(1,&Nbin,"C2R M K");
            fft_back->SetPointsComplex(tre_full,tim_full);
            fft_back -> Transform();
            backhist = TH1::TransformHisto(fft_back,backhist,"Re");
            st.Hist(backhist);
            backhist -> Draw();

            
            st.Hist(fhist);
            fhist -> Draw();
            
            //これを踏まえた時にフィルタがどう機能しているのかを確認する
            fhist2 = sighist -> FFT(fhist2,"MAG");
            TVirtualFFT *fft2 = TVirtualFFT::GetCurrentTransform();//直近でFFTしたデータのあれこれを引っ張り出すポインタ
            fft2 -> GetPointsComplex(mre_full,mim_full);//各点のDCとナイキスト周波数を取得
            double mcr = mre_full[100];
            double mci = mim_full[100];
            TH1* reshist = new TH1D("reshist","res;1/Freq[ns];Spec",nbin,0,nbin/2.5);
            rep(bin,nbin){
                double tr = tre_full[bin];
                double ti = tim_full[bin];
                double mr = mre_full[bin];
                double mi = mim_full[bin];
                double K = (mcr*mcr+mci*mci)+(mr*mr+mi*mi);
                tre_full[bin] = (tr*mr+ti*mi)/K;
                tim_full[bin] = (ti*mr-tr*mi)/K;
                reshist -> SetBinContent(bin,sqrt(tre_full[bin]*tre_full[bin]+tim_full[bin]*tim_full[bin]));
            }
            st.Hist(reshist);
            reshist -> Draw();*/