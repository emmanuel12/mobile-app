/*
   You can email me at "emmanuel.nsanga@gmail.com".

   Usage: ./freq <filename path>; e.g './freq ~/audio.wav'.

   Copyright [2013] [Emmanuel Nsanga]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License. **/

#include<iostream>
#include<complex>
#include<fstream>
#include<algorithm>
#include<vector>
#include<map>

using namespace std;

class Frequency{
    public:
        // Structure for reading the wave file.
        typedef struct snd_file{
                char chunk_id[4];
                int chunk_size;
                char format[4];
                char subchunk1_id[4];
                int subchunk1_size;
                short int audioformat;
                short int channels;
                int samplerate;
                int byte_rate;
                short int block_align;
                short int bits_per_sample;
                char subchunk2_id[4];
                int frames;
                unsigned char DATA[];

        } file;
        typedef struct snd_file* snd_p;

        struct frq_rmst{ // Structure for the frequency and rms amplitude values.
           float freq;
           double rms;
           // Comparison operator.
           bool operator!=(const frq_rmst &frq_rmst1){
                bool co = false;
                if (freq != frq_rmst1.freq){
                    co = true;
                }else if (rms != frq_rmst1.rms){
                    co = true;
                }
                return co;           
           }                  
        };
        std::vector<std::complex <double> > FFT(std::vector<std::complex <double> >ar) // A recursive implementation of \
        the 1D Cooley-Turkey FFT algorithm.
        {
            int n = ar.size();
            if (n == 1){
                return ar;

            }else{
                std::vector<std::complex <double> > fev;
                std::vector<std::complex <double> > fod;
                std::vector<std::complex <double> > h;
                h.resize(n/2);
                std::vector<std::complex <double> > h1;
                h1.resize(n/2);

                // Split the data into two parts to reduce the exponent (size n) of e.
                int i;
                int pos = 0;
                for (i = 0; i < n/2; i++){
                     h[i] = ar[pos];
                     pos = pos+2;
                }
                pos = 1;
                for (i = 0; i < n/2; i++){
                     h1[i] = ar[pos];
                     pos = pos+2;
                }

                //Recursively go through this function while slowly reducing the exponent (size n) of e,\
                  which in turn reduces the computation time.
                std::vector<std::complex <double> > fff = FFT(h);
                fev = fff;
                std::vector<std::complex <double> > ff = FFT(h1);
                fod = ff;

                // Calculate the exponent product of e, and then combine it with the two sets of recursed data.
                std::vector<std::complex <double> > comb;
                comb.resize(n);
                for (i = 0; i < n; i++){
                     comb[i] = 0;
                }
                double pi = 3.141592653589793;
                for (int m = 0; m < n/2; m++){
                     std::complex <double> om = ((2.0 * pi * 1j * -m) / n);
                     std::complex <double> co = std::exp(om);
                     std::complex <double> com = fev[m]+co*fod[m];
                     comb[m] = com;
                     comb[m+(n/2)] = fev[m] - co * fod[m];
                }
                return comb;

            }

        }
        std::vector<Frequency::frq_rmst>  freq_rms(const char *filename) // Calculate frequency and rms amplitude values.
        {
            FILE * infile = fopen(filename,"rb");
            snd_p meta = (snd_p)malloc(sizeof(file));
            fread(meta, 1, (sizeof(file)), infile);

            //Read the raw frame data.
            int l = meta->frames;
            infile = fopen(filename,"rb");
            meta = (snd_p)malloc(sizeof(file)+l);
            fread(meta, 1, (sizeof(file)+l), infile);

            //Calculate the blackman window function.
            int srate = meta->samplerate;
            float chunk = 2048; //Increase chunk size for sensitivity.
            float pi = 3.141592653589793;
            std::vector<double> window;
            window.resize(chunk);
            float N = chunk;
            for (int n = 0; n < N; n++){
                 double w = 0.42-0.5*std::cos(2.0*pi*n/(N-1))+0.08*std::cos(4.0*pi*n/(N-1));
                 window[n] = w;
            }

            //Loop through the frames by chunk size.
            int poss = 0;
            int ch = l/(chunk*2);
            std::vector<Frequency::frq_rmst> frqs_rms;
            for (int ii = 0; ii < ch; ii++){
                 //Unpack the raw data and times it by the hamming window.
                 std::vector<std::complex <double> > ind;
                 ind.resize(chunk);
                 int i;
                 double on = 1;
                 double nn = 32768;
                 std::complex<double> short_normalize = on/nn;
                 std::complex<double> sum_squares = 0;
                 std::complex<double> chunkk = chunk;
                 for (i = 0; i < chunk; i++){
                      std::complex <double> sample = meta->DATA[poss]-meta->DATA[poss+1]+-1;
                      ind[i] = sample*window[i];
                      std::complex<double> n = sample*short_normalize;
                      sum_squares+=n*n;
                      poss = poss+2;
                 }
                 std::complex<double> rms = std::sqrt(sum_squares/chunkk);
                 frqs_rms.push_back(Frequency::frq_rmst());
                 frqs_rms[ii].rms = rms.real();
                 
                 //Put the data through the FFT algorithm and then square the absolute value.
                 std::vector<std::complex <double> > dftd = FFT(ind);
                 std::vector <double> dft;
                 dft.resize(chunk/2);
                 for (i = 0; i < chunk/2; i++){
                      std::complex <double> ab = std::pow(std::abs(dftd[i]),2);
                      dft[i] = ab.real();

                 }
                 //Find the maximum and its index.
                 dft[0] = 0;
                 double which = *std::max_element(dft.begin(), dft.end());
                 int pos = 0;
                 int sz = dft.size();
                 for (i = 0; i < sz; i++){
                      if (which == dft[i]){
                          which = pos;
                          break;

                      }else{
                          pos = pos+1;
                      }
                 }
                 
                 //Do quadratic interpolation around the maximum and then output the frequency.
                 float frq;
                 if (which != dft.size()){
                     double l0 = std::log(dft[pos-1]);
                     double l1 = std::log(dft[pos]);
                     double l2 = std::log(dft[pos+1]);
                     int x1 = ((l2-l1)*0.5)/(2*(l1)-(l2)-(l0));
                     frq = (which+x1)*srate/chunk;
                     frqs_rms[ii].freq = frq;

                 }else{
                     frq = which*srate/chunk;
                     frqs_rms[ii].freq = frq;
                 }

            }
            return frqs_rms;
        }



};

class Recognition{   
   public:
       int euclideanDistance(Frequency::frq_rmst x, Frequency::frq_rmst y) // A function to calculate the euclidean distance.
       {
           int sum = 0;
           int a = x.freq;
           int b = y.freq;
           sum+=((a-b)*(a-b)); // The sum of the squared difference between the two frequency values. 
           a = x.rms;
           b = y.rms;
           sum+=((a-b)*(a-b)); // The sum of the squared difference between the two rms amplitude values.
          
           int euclid = std::sqrt(sum);
           return euclid;
       }

       std::vector<std::vector<Frequency::frq_rmst> > partition(std::vector<Frequency::frq_rmst> points, int k, std::vector<Frequency::frq_rmst> means) // A \
       function to partition the data around the cluster centers.
       {
           
           std::vector<std::vector<Frequency::frq_rmst> > thePartition; // A two dimensional vector to store the clusters.
           thePartition.resize(k); // Resize the two dimensional vector by the number of cluster centers.
           
           
           for (int i = 0; i < points.size(); i++){
                std::vector<int> indices;
                indices.resize(k);
                Frequency::frq_rmst x = points[i];
                for (int ii = 0; ii < k; ii++){
                     Frequency::frq_rmst y = means[ii]; // Pick a cluster center at index ii.
                     int index = euclideanDistance(x, y); // Calculate the euclidean distance between the data points and the cluster centers.
                     indices[ii] = index;
                }
                // Find the lowest data point to cluster center and check which cluster center it belongs to.
                int closestindex = *std::min_element(indices.begin(), indices.end());
                int pos = 0;
                for (int in = 0; in < k; in++){
                     if (indices[in] == closestindex){
                         break;
                     
                     }else{
                         pos = pos+1;
                     }
                }              
                thePartition[pos].push_back(x); // Store the data point at the cluster center posistion.
           }     
           return thePartition;

       }

       Frequency::frq_rmst mean(std::vector<Frequency::frq_rmst> points) //Calculate the mean for the cluster centers.
       {
           int n = points.size();
           double sum = 0;
           for (int i = 0; i < n; i++){
                sum+=points[i].rms;
                
           }    
           int sum1 = 0;
           for (int i = 0; i < n; i++){
                sum1+=points[i].freq;
           }
           Frequency::frq_rmst m;
           m.freq = sum1/n; // The mean for the frequency values.
           m.rms = sum/n; // The mean for the rms amplitude values.
           return m;          

       }

       std::vector<std::vector<Frequency::frq_rmst> > kmeans(std::vector<Frequency::frq_rmst> points, int k, \
       std::vector<Frequency::frq_rmst> initialmeans) //K-means clustering algorithm.
       {
           std::vector<std::vector<Frequency::frq_rmst> > oldpartition; //initialise a two dimensional vector to store the old clusters. 
           std::vector<std::vector<Frequency::frq_rmst> > newpartition = partition(points, k, initialmeans); //Calculate the new clusters.
           
           //Loop until convergence.
           bool comp = true;
           while (comp){
                  oldpartition = newpartition;
                  std::vector<Frequency::frq_rmst> newmeans;
                  newmeans.resize(k);
                  //Calculate new cluster centers for the new clusters.
                  for (int i = 0; i < k; i++){
                       newmeans[i] = mean(oldpartition[i]);
                  }
                  newpartition = partition(points, k, newmeans);
                  
                  //Check for convergence.
                  bool tar = true;
                  for (int i = 0; i < k; i++){
                       for (int ii = 0; ii < oldpartition[i].size(); ii++){
                            Frequency::frq_rmst c = oldpartition[i][ii];
                            Frequency::frq_rmst c1 = newpartition[i][ii];
                            if (c != c1){
                                tar = false;
                                break;
                             }
                        }    
                  }
                  if (tar){
                      comp = false;
                  }
           }     
           return newpartition;

       }


       struct fb_ret{
              std::vector<std::map<const char *, float> > fwd;
              std::vector<std::map<const char *, float> > bkw;
              std::map<const char *, std::vector<float> > post;
       };       

       fb_ret fwd_bkw(std::vector<int> ob, std::vector<const char *> st, std::map<const char *, float> a_0, std::map<const char *, std::map<const char *, float> >a,\
       std::map<const char *, std::map<int, float> >e, const char * end_st) //The Forward-Backward algorithm.
       {
           int sz = ob.size();
           std::vector<std::map<const char *, float> > fwd;
           fwd.resize(sz);
           std::map<const char *, float> f_prv;
           int stz = st.size();
           //Itterate through the vectors for the forward part of the algorithm.
           for (int i = 0; i < sz; i++){
                std::map<const char *, float> f_curr;
                int x_i = ob[i];
                for (int ii = 0; ii < stz; ii++){
                     float prevs = 0;
                     const char * stt = st[ii];
                     if (i == 0){
                         //Base case of the forward part of the algorithm.
                         prevs = a_0[stt];

                     }else{  
                         for (int in = 0; in < stz; in++){
                              const char * k = st[in];
                              prevs+=f_prv[k]*a[k][stt];
                         }
                     }
                     f_curr[stt] = e[stt][x_i] * prevs;
                }
                fwd[i] = f_curr;
                f_prv = f_curr;
           }
           float p_fwd = 0;
           for (int i = 0; i < stz; i++){
                const char * k = st[i];
                p_fwd+=f_prv[k] * a[k][end_st];
           }

           std::vector<std::map<const char *, float> > bkw;
           bkw.resize(sz+1);
           std::map<const char *, float> b_prv;
           //Itterate through the vectors in reverse order for the backward part of the algorithm.
           for (int i = sz; i > 0; i--){
                int x_plus = ob[i];
                std::map<const char *, float> b_curr;
                for (int ii = 0; ii < stz; ii++){
                     const char * stt = st[ii];
                     if (i == sz){
                         //Base case of the backward part of the algorithm.
                         b_curr[stt] = a[stt][end_st];

                     }else{
                         float sum = 0;
                         for (int in = 0; in < stz; in++){
                              const char * l = st[in];
                              sum+=a[stt][l]*e[l][x_plus]*b_prv[l];
                         }
                         b_curr[stt] = sum;
                     }
                }
                bkw[i-1] = b_curr;
                b_prv = b_curr;
           }
           //Merge both the forwad part and the backward part.
           float p_bkw = 0;
           for (int i = 0; i < stz; i++){
                const char * l = st[i];
                p_bkw+=a_0[l] * e[l][ob[0]] * b_prv[l];
           }
           std::map<const char *, std::vector<float> > post;
           for (int i = 0; i < stz; i++){
                const char * stt = st[i];
                std::vector<float> postv;
                postv.resize(sz);
                for (int ii = 0; ii < sz; ii++){
                     postv[ii] = (fwd[ii][stt]*bkw[ii][stt])/p_fwd;
                }
                post[stt] = postv;
           }   
           
           fb_ret ret;
           ret.fwd = fwd;
           ret.bkw = bkw;
           ret.post = post;
           return ret;
       }



};
  
