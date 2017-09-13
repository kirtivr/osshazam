#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>

#ifndef NOCILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#endif

#include <iostream>

#define MAGIC_NUMBER 1024
#define UINT unsigned int

#define FREQ_LOWER 30.0
#define FREQ_UPPER 630.0
#define FREQ_INTERVAL 50.0
#define RANGE_MAX_INDEX 10
#define FALSE_SONG_LENGTH 66 // We are increasing signal length to closest power of two. Puts the perceived song length out of whack
#define TRUE_SONG_LENGTH 53
#define SAMPLE_RATE 16000
#define NUM_SAMPLES 848000
#include "fft.h"
#include "string.h"
#include "songcompare.h"

using namespace std;

double TIME_PER_CHUNK;
map <double,vector<double>> time2dominantfreq_st[2];
map <double,vector<double>> time2dominantfreq_ta[2];

// Good excuse to use an inline function. Copies from source to destination in chunks of N complex_t structs.
static inline void _complex_copy(complex_t* dest,complex_t* source,UINT N){
  memmove(dest,source,N * sizeof(source[0]));
}
  
// Copy the values obtained in y to a . Divide and conquer !
void _copy(complex_t* y, complex_t*a, UINT N) {
  if(N <= MAGIC_NUMBER) {
    _complex_copy(a,y,N);
    return;
  }

  cilk_spawn _copy(y, a, N/2);
  cilk_spawn _copy(y+N/2, a+N/2, N/2);
  cilk_sync;
}


// If the number of elements passed into parallel_combine is lte 1024, we just serially combine.
// N is the number of elements.
// y_alt is y+N/2. It has been passed as a function parameter because accessing the address sequentially gives better locality, letting the TLB do its thing.
// start_index keeps track of where we are in the recursion tree.
// num_elem is the number of elements in the subtree corresponding to _fft_conquer
void _fft_conquer(complex_t *y, complex_t *y_alt, UINT start_index, UINT num_elem, UINT N) {
  // N <= 1024(MAGIC_NUMBER) has been reached empirically. Performance is worse when the threads are too small than when you just compute serially.
  if(num_elem <= MAGIC_NUMBER) {
    complex_t y1,y2,omega;
    // omega is re-calculated every time we pass through the loop. This gives a more accurate value as compared to multiplying and obtaining omega.
    // twiddle is the twiddle factor
    // y1 is omega * y[k+N/2].
    // y2 is the value to be put at y[k+N/2]. I save it for the time being and copy it into y[k+N/2] once y[k] is set. 
    for(UINT k = 0; k < num_elem; k++) {
      COMPLEX_ROOT(omega,N,start_index+k);
      COMPLEX_MULT(y1,omega,y_alt[k]);
      COMPLEX_SUB(y2,y[k],y1);
      COMPLEX_ADD(y[k],y[k],y1);
      _complex_copy(&(y_alt[k]),&y2,1);
    }
    return;
  } else {
    cilk_spawn _fft_conquer(y,y_alt,start_index,num_elem/2,N);
    cilk_spawn _fft_conquer(y+num_elem/2,y_alt+num_elem/2,start_index+num_elem/2,num_elem/2,N);
    cilk_sync;
    return;
  }
}

// a : input of size n
// y : output of size n
// N : number of elements in the subtree corresponding to _fft
// step_a: takes values of 1,2,4,8.. . Combined with the fact that we divide _fft into two parallel threads with a and a+step_a as input, this ensures
//         that we are able to populate y with the (odd and even indexes of a) at every branch of _fft.
void _fft(complex_t *a, complex_t *y, UINT N, UINT step_a) {
  UINT step = 1;
  // copy value at a to corresponding index in y.
  // what essentially happens is : odd indices of a are stored from y[N/2]...y[N]. Within y[N/2:N] even indices of odd a' are stored from y[N/2:N/4].. and so on
  if(N==1) {
    _complex_copy(y,a,1);
    return;
  }

  cilk_spawn _fft(a, y, N/2, 2*step_a);
  cilk_spawn _fft(a+step_a, y+N/2, N/2, 2*step_a);
  cilk_sync;
  
  cilk_spawn _fft_conquer(y,y+N/2,0,N/2,N);
  cilk_sync;
}

/* Compute the FFT of an array a of length N. */
void FFT(complex_t *a, int N) {
  complex_t* y = (complex_t*) malloc(N * sizeof(complex_t));

  if(y != NULL) {
    cilk_spawn _fft(a,y,N,1);
    cilk_sync;
  
    cilk_spawn _copy(y,a,N);
    cilk_sync;

    free(y);
  } else {
    printf("malloc failed! This should not happen. Perhaps there is a heap corruption?\n");
  }
}


void matchSongs( map<double,vector<double>>time2dominantfreq_st[], map<double,vector<double>>time2dominantfreq_ta[],int num_window_types, int num_chunks) {
  int point_matched[num_window_types]; 
  int tolerance = 25.0;
 
  map<double,vector<double>>::iterator it_ta;
  map<double,vector<double>>::iterator it_st;
  for(int k=0; k<num_window_types; k++) {
     point_matched[k] = 0;
    it_st = (time2dominantfreq_st[k]).begin();
    it_ta = (time2dominantfreq_ta[k]).begin();

    int max_points_matched =0;

    while(it_st != (time2dominantfreq_st[k]).end() && it_ta != (time2dominantfreq_ta[k]).end()){
      double timest = it_st->first;
      double timeta = it_ta->first;

      if(timest > 53.0 || timeta > 53.0)
	break;

      vector<double> dominantfreq_st = it_st->second;
      vector<double> dominantfreq_ta = it_ta->second;
      int matched = 5;
      
      // we know both sizes are the same
      for(int i=0; i< dominantfreq_st.size(); i++) {
	if(fabs(dominantfreq_st[i]) > fabs( dominantfreq_ta[i]) &&  fabs(dominantfreq_st[i]) - fabs(dominantfreq_ta[i]) > tolerance) {
	  matched--;
	  if(matched == 0)
	    break;
	}

	if(fabs(dominantfreq_ta[i]) > fabs( dominantfreq_st[i]) &&  fabs(dominantfreq_ta[i]) - fabs(dominantfreq_st[i]) > tolerance) {
	  matched--;
	  if(matched == 0)
	    break;
	}
      }
      
      if(matched){
	max_points_matched++;
	if(max_points_matched >   point_matched[k]) {
	  point_matched[k] = max_points_matched;
	}
      } else {
	  max_points_matched =0 ;
      }

      it_st++;
      it_ta++;
    }
    cout<<endl<<"For window set "<<k+1<<": "<<point_matched[k]<<" consecutive chunks were very similar"<<endl;    
  }//for
  int greater = point_matched[0] >= point_matched[1] ?  point_matched[0]: point_matched[1];
  double ratio = (double)greater/num_chunks;
  cout<<endl<<"Similarity ratio is "<<ratio<<endl;
}
// Calculate the intensity of the frequency at any chunk location.
void getPowerSpectrum(complex_t* incoming_chunk, int windowSize, map<double,double>&freq2power) {
  double power;
  double frequency;
  
  for(int i=0; i<windowSize; i++ ) {   
    power = 1+ (incoming_chunk[i]).real * (incoming_chunk[i]).real +  (incoming_chunk[i]).imag * (incoming_chunk[i]).imag;
    power = log(power);
    //printf( "pwrspc[%d] = %f\n",
    //	    i,  powerSpectrum[i]);
    // Get frequency at this position
    // Size of FFT = windowSize
    // Sample rate is constant SAMPLE_RATE
    // what is the frequency ? i*SAMPLE_RATE/windowSize
    frequency = i*SAMPLE_RATE/windowSize;
    freq2power.insert(make_pair(frequency,power));
  }
}

int getRangeForFrequency(double frequency) {
  for(int i=0; i<=RANGE_MAX_INDEX; i++) {
    double lower = (double)(FREQ_LOWER + i*FREQ_INTERVAL);
    double upper = (double)(FREQ_LOWER + (i+1)*FREQ_INTERVAL);
    if(frequency >= lower && frequency <= upper)
      return i;
  }
  return -1;
}

void doTheRealStuff(vector<complex_t*>st_windows[],vector<complex_t*>ta_windows[],int num_window_types,int num_chunks,int windowSize) {
  int i,j,k;
  // Maps highest intensity frequencies
  map<double,double>freq2powerst[num_window_types][windowSize];
  map<double,double>freq2powerta[num_window_types][windowSize];

  for(k=0; k<num_window_types; k++) {
    for(j=0; j<num_chunks; j++) {     
      complex_t* stairway2hvn = (st_windows[k])[j];
      complex_t* taurusv = (ta_windows[k])[j];

      getPowerSpectrum(stairway2hvn,windowSize,freq2powerst[k][j]);
      getPowerSpectrum(taurusv,windowSize,freq2powerta[k][j]);
    }
      
  }
  // iterate through each chunk of each window, and find the most intense frequency in a given range
  // Take intervals of frequencies 30-80, 80-130, 130-180,180-230,230-280

  map<int,double>range2powst;
  map<int,double>range2powta;
  
  map<int,double>range2freqst;
  map<int,double>range2freqta;
  
  for(k=0; k<num_window_types; k++) {  
    for(j=0; j<num_chunks; j++) {
      
      for(i=0; i<=RANGE_MAX_INDEX; i++) {
	range2powst[i] = 0.0;
	range2powst[i] = 0.0;
      }
      
      // iterate through power2freq map
      
      for (map<double,double>::iterator it=(freq2powerst[k][j]).begin(); it!=(freq2powerst[k][j]).end(); ++it) {
	double frequency = it->first;
	double power = it->second;
	int range_no = getRangeForFrequency(frequency);
	//	if (j==3)
	//printf("stairway .. frequ %f pow %f\n",frequency,power);
	if(range_no != -1) {
	  double relevant_max_pow = range2powst[range_no];
	  // -1 means out of range
	  if(power > relevant_max_pow) {
	    range2powst[range_no] = power;
	    range2freqst[range_no] = frequency;
	  }//if
	}//if
      }

      
      for(i=0; i<=RANGE_MAX_INDEX; i++) {
	range2powta[i] = 0.0;
	range2powta[i] = 0.0;
      }
      
      // iterate through power2freq map
      for (map<double,double>::iterator it2=(freq2powerta[k][j]).begin(); it2!=(freq2powerta[k][j]).end(); ++it2) {
	double frequency = it2->first;
	double power = it2->second;
	int range_no = getRangeForFrequency(frequency);
	//	if (j==3)
	//printf("taurus .. frequ %f pow %f\n",frequency,power);
	if(range_no != -1) {
	  double relevant_max_pow = range2powta[range_no];
	  // -1 means out of range
	  if(power > relevant_max_pow) {
	    range2powta[range_no] = power;
	    range2freqta[range_no] = frequency;
	  }//if
	}//if
	
      } // we are done with a chunk, push most dominant frequency to the great hashmap
      /*
	for(i=0; i<=RANGE_MAX_INDEX; i++) {
	printf("stairway,, dominant powe in range i=%d is %f\n",i,range2powst[i]);
	printf("taurus,, dominant powee in range i=%d is %f\n",i,range2powta[i]);
	}
      */
      vector<double> dominantfreq_st;
      vector<double> dominantfreq_ta;
      
      for(i=0; i<=RANGE_MAX_INDEX; i++) {
	dominantfreq_st.push_back(range2freqst[i]);
	dominantfreq_ta.push_back(range2freqta[i]);
      }

      double windowOffset = ((double)k/num_window_types) * TIME_PER_CHUNK;
      double timestamp = windowOffset + TIME_PER_CHUNK * (double)j;

      //  printf("window offset = %f time stamp %f TIME_PER_CHUNK = %f\n",windowOffset,timestamp,TIME_PER_CHUNK);

      time2dominantfreq_st[k].insert(make_pair(timestamp,dominantfreq_st));
      time2dominantfreq_ta[k].insert(make_pair(timestamp,dominantfreq_ta));
    }// j chunk
  }// window

  cout<<"For stairway to heaven "<<endl<<endl;
    
  for(k=0; k<num_window_types; k++) {
    for (map<double,vector<double>>::iterator it=(time2dominantfreq_st[k]).begin(); it!=(time2dominantfreq_st[k]).end(); ++it) {
      double timest = it->first;
      if(timest <= TRUE_SONG_LENGTH){
	printf("Timestamp: %f\t",timest);
	vector<double>dom = it->second;

	for(i=0; i<dom.size(); i++) {
	  cout<<dom[i]<<", ";
	}
	cout<<endl;
      }
    }
  }

  cout<<endl<<endl<<"For taurus  "<<endl<<endl;
  for(k=0; k<num_window_types; k++) {
    for (map<double,vector<double>>::iterator it=(time2dominantfreq_ta[k]).begin(); it!=(time2dominantfreq_ta[k]).end(); ++it) {
      double timeta = it->first;
      if(timeta <= TRUE_SONG_LENGTH){
	printf("Timestamp: %f\t",timeta);
	vector<double>dom = it->second;

	for(i=0; i<dom.size(); i++) {
	  cout<<dom[i]<<", ";
	}
	cout<<endl;
      }
    }
  }

  matchSongs(time2dominantfreq_st,time2dominantfreq_ta,num_window_types,num_chunks);
}


// Create a hamming window of windowLength samples in buffer
void hamming(int windowLength, double *buffer) {
  
  for(int i = 0; i < windowLength; i++) {
    buffer[i] = 0.54 - (0.46 * cos( 2 * M_PI * (i / ((windowLength - 1) * 1.0))));
    //printf("value at i=%d is %f\n",i,buffer[i]);
  }
}

void STFFT(complex_t *stairway, complex_t *taurus, int signalLength, int windowSize, int hopSize) {
  int i;
  complex_t* st_chunk = (complex_t *) malloc( sizeof(complex_t) * windowSize);
  complex_t* ta_chunk = (complex_t *) malloc( sizeof(complex_t) * windowSize);
  complex_t* st_freq_chunk = (complex_t *) malloc( sizeof(complex_t) * windowSize);
  complex_t* ta_freq_chunk = (complex_t *) malloc( sizeof(complex_t) * windowSize);

  int num_chunks = signalLength/windowSize;
  
  int num_window_types = windowSize/hopSize;
  vector<complex_t*>st_windows[num_window_types];
  vector<complex_t*>ta_windows[num_window_types];

  //cout<<"going for ham"<<endl;
  // Create a hamming window of appropriate length
  double window[windowSize];
  hamming(windowSize, window);
 
  int chunkPosition = 0;
  int readIndex;
  // Should we stop reading in chunks? 
  int bStop = 0;
  int numChunks = 0;
  //cout<<endl<<"given input"<<endl;
  /*for(i=0; i<signalLength; i++){
    cout<<"stairway["<<i<<"]"<<"real "<< stairway[i].real<<"imag: "<<stairway[i].imag<<endl;
    cout<<"taurus["<<i<<"]"<<"real "<< taurus[i].real<<"imag: "<<taurus[i].imag<<endl;
    }*/
  
  //cout<<endl<<endl;
  // Process each chunk of the signal
  while(chunkPosition < signalLength && !bStop) {
    // Copy the chunk into our buffer
    for(i=0; i<windowSize; i++) {
      readIndex = chunkPosition + i;
      if(readIndex < signalLength) {
      // Note the windowing! 
      (st_chunk[i]).real = (stairway[readIndex]).real * window[i]; 
      (st_chunk[i]).imag = 0.0; 
      (ta_chunk[i]).real = (taurus[readIndex]).real * window[i];  
      (ta_chunk[i]).imag = 0.0; 
      } else {     
      // we have read beyond the signal, so zzero paf
      (st_chunk[i]).real = 0.0; 
      (st_chunk[i]).imag = 0.0; 
      (ta_chunk[i]).real = 0.0; 
      (ta_chunk[i]).imag = 0.0; 
      bStop = 1;
      }
    }// for windowsize
  
    FFT(st_chunk,windowSize);
    FFT(ta_chunk,windowSize);
    // Copy the first (windowSize/2 + 1) data points into spectrogram.
    // We do this because the FFT output is mirrored about the nyquist 
    // frequency, so the second half of the data is redundant.
    int window_num = numChunks%num_window_types;
    
    complex_t* temp_st_chunk = (complex_t *) malloc( sizeof(complex_t) * windowSize/2);
    complex_t* temp_ta_chunk = (complex_t *) malloc( sizeof(complex_t) * windowSize/2);

    _complex_copy(temp_st_chunk,st_chunk,windowSize/2);
    _complex_copy(temp_ta_chunk,ta_chunk,windowSize/2);
    
    (st_windows[window_num]).push_back(temp_st_chunk);
    (ta_windows[window_num]).push_back(temp_ta_chunk);
    
    chunkPosition += hopSize;
    numChunks++;
    
    // free dummy objs
    //free(temp_st_chunk);
    //free(temp_ta_chunk);
  }// while chunkpos
 
  doTheRealStuff(st_windows,ta_windows,num_window_types,num_chunks,windowSize/2);
}

void analyse(complex_t *stairway, complex_t *taurus, UINT signalLength) {
  int windowLength = 16384;
  int hopSize = windowLength/2;
  double chunks = (signalLength/(double)windowLength);
  TIME_PER_CHUNK = FALSE_SONG_LENGTH/chunks;
  STFFT(stairway,taurus,signalLength,windowLength,hopSize);
}

unsigned int closestPowOf2(unsigned int v) {
  const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
  const unsigned int S[] = {1, 2, 4, 8, 16};
  int i;

  register unsigned int r = 0; // result of log2(v) will go here
  for (i = 4; i >= 0; i--) // unroll for speed...
    {
      if (v & b[i])
	{
	  v >>= S[i];
	  r |= S[i];
	} 
    }
  return 1<<(r+1);
}

int main(int argc, char * argv[]) {
  complex_t *stairway_d;
  complex_t *taurus_d;
  unsigned int signalLength = closestPowOf2(NUM_SAMPLES);
  //unsigned int signalLength = 848000;
  //cout<<"signal length is "<<signalLength<<endl;
  //sleep(2);
  stairway_d = (complex_t*) malloc(sizeof(complex_t) * signalLength);
  taurus_d = (complex_t*) malloc(sizeof(complex_t) * signalLength);
  
  fstream st_stairway("./stairway.txt", ios_base::in);
  //fstream st_stairway("./u2_worwithoutyou.txt", ios_base::in);
  //fstream st_taurus("./africa.txt", ios_base::in);
  fstream st_taurus("./taurus.txt", ios_base::in);
  int i=0;

  string x,y;  
  double a,b;
  
  while (st_stairway >> x && st_taurus >> y)
    {
      // Audacity output has [-inf] values which need to be parsed out.
      if(x.compare("[-inf]") == 0) {
	if(i != 0) {
	  (stairway_d[i]).real = (stairway_d[i-1]).real;
	} else {
	  (stairway_d[i]).real = 0.0;
	}
      } else {
	a = atof(x.c_str());
	(stairway_d[i]).real = a;
      }

      
      if(y.compare("[-inf]") == 0) {
	if(i!=0){
	  (taurus_d[i]).real = (taurus_d[i-1]).real;
	} else {
	  (taurus_d[i]).real = 0.0;
	}
      } else {
	b = atof(y.c_str());
	(taurus_d[i]).real = b;
      }

      (stairway_d[i]).imag = 0.0;
      (taurus_d[i]).imag = 0.0;
      
      //printf("%f,%f read\n",a,b);
      i++;
    }
  
  // We want the signal to be a power of 2, so for good karma we pad with zeros
  /* while (i<signalLength){
    (stairway_d[i]).real = 0.0;
    (taurus_d[i]).real = 0.0;
    (stairway_d[i]).imag = 0.0;
    (taurus_d[i]).imag = 0.0;
    i++;
    }*/
  /*
    cout << "Stairway window.. " << endl;
    for(i=0; i<signalLength; i++ ) {
    printf( "input[%d] = { %2.2f, %2.2f }\n",
    i, (stairway_d[i]).real,(stairway_d[i].imag));
    }

    
    //input for the FFT calculation
    cout << "Taurus window.. " << endl;
    for(i=0; i<signalLength; i++ ) {
    printf( "input[%d] = { %2.2f, %2.2f }\n",
    i, (taurus_d[i]).real,(taurus_d[i].imag));
    }
  */
  //cout<<"calling analyse\n"<<endl;
  analyse(stairway_d,taurus_d,signalLength);
  
  return 0;  
} 
