// new code 2/25/2020
#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "string.h"
#include "math.h"
#include "signal.h"


#define FFTLEN 4096
#define FFTLEN2 2048
#define MAXBIN 1024 // about 12KHz top bin for harmonic analysis
#define MAXBLOCK 1024
//#define CONLY 1

void fft(int);
void hannCompute(void);
void hannMult(void);

void init_fft_twiddles(void);

#ifdef CONLY
void parseWavHeader();
float get_wav16(void);
int processAudioData(float *,int,float);
FILE *wav_infile; 
FILE *infile;
FILE *outfile;

char inbuf[300];
char fname[60];
char fname_r[40];
char fname_l[40];
unsigned char ch;
int bits_per_sample,num_samples;
int num_channels,block_align;
int block_GL;
#endif

int num_points;
int numcycles;


float fs_wav;
float re[FFTLEN] = {0.0};
float im[FFTLEN] = {0.0};
float twid_re[FFTLEN2] = {0.0};
float twid_im[FFTLEN2] = {0.0};

float fftmagdb[MAXBIN] = {0.0};
float fftpower_GL[MAXBIN] = {0.0};


float maxampl_db_GL = 0.0;





float win[FFTLEN] = {0.0};
float PI = 3.141592653589793;


float fbin = 1.0;

int init = 1;
int block_count_GL;

int peakList[MAXBIN] = {0};
float amplList[MAXBIN] = {0};
int peakList1[MAXBIN] = {0};
float amplList1[MAXBIN] = {0};


#ifdef CONLY

int main() 
{
int i,k,kk;
int num_return_samples = 1024;
float wavin[1000000] = {0.0};
int offseti;
float bufferIn[FFTLEN] = {0.0};

int numblocks,block,count,tempi;
float left_in,right_in,wav_in;

outfile = fopen("freqout.txt","w");
char buff[] = "./wavs/runyon_zz2.wav";
if( (wav_infile = fopen(buff,"rb")) == NULL) {
	printf("cant find it\n");
	exit(0);
	}
parseWavHeader(); /** sets global num_samples, num_channels, bits_per_sample **/
printf("parsed wav file, num_samples = %d\n",num_samples);

/*********************start main routine ***************/
count = 0;
block_count_GL = 0;
hannCompute();
for(i=0;i < num_samples;i++) {

	if(num_channels == 1) { // mono
		wavin[i] = get_wav16();
	} else { // stereo
		wavin[i] = get_wav16();
	    wavin[i] = get_wav16();
	}

	

  } /** end of while loop ***/

numblocks = (int)(num_samples/FFTLEN);
printf("numblocks = %d\n",numblocks);

// start main calling routine 
for(k=0;k < numblocks;k++) {
	block_GL = k;
	offseti = k*FFTLEN;
	for(kk = 0;kk < FFTLEN;kk++) {
		bufferIn[kk] = wavin[kk+offseti];
	}

	processAudioData(bufferIn,FFTLEN,fs_wav);	
}


fclose(outfile);
}
#endif



/****************************************************/

#ifndef CONLY
extern "C" {
#endif

int processAudioData(float * wav_in, int num_samples, float fs_wav) // returns a pointer to an array
{

int i,k,kk,count;
int numPeaks,peakListMaxIndex=0;
int ptr=0,ptr2=0,yy,numList1;
int bin3,bin4,bin5,bin6,testSigBin,testBetweenBin,firstPeaki;
float temp=0.0,spectMax=0.0;
float Beta,Alpha,Gamma,frac_bin,amplMax=0.0,freqMaxBin,low_limit;
float spectMaxiInterp,minpeak,maxpeak,snr,changeRatio;
float pwr,rms,gain;
float firstPeakInterp = 0.0, firstPeakInterp_div2 = 0.0;
float ratioMF=0.0,ratioMF_zm1=0.0,divby=1.0,fundInterp;
float snrGate = 1.0,noisePow = 0.0,sigPow = 0.0;
static float freq=0.0,freqZm1=0.0,freqGate=0.0;
static float freqRawZm1 =0.0,freqRaw=0.0,freqRatio=0.0;
static int waitCount;




if(init==1) {
	init_fft_twiddles();
	init = 0;
	waitCount = 0;
	printf("pitch3.cpp Hello\n");
}


//fs_wav = 44.1e3; // temporary, will need a way to pass this in from the JS side
//fs_wav = 48e3;


/*********************start main routine ***************/
count = 0;
block_count_GL = 0;
fbin = fs_wav/(float)FFTLEN;
//printf("num_samples = %d\n",num_samples);
hannCompute();



// load fft buffer fromm input buffer
for(k=0; k < FFTLEN;k++) {re[k]=wav_in[k];im[k] = 0.0;}
pwr = 0.0;
for(k=0; k < FFTLEN;k++) {pwr = pwr + re[k]*re[k];}
rms = sqrt(pwr/(float)FFTLEN);
gain = 1.0/rms;
for(k=0; k < FFTLEN;k++) {re[k]=re[k]*gain;}

hannMult(); // in place hann window on global re
fft(FFTLEN); // in-place fft on global buffer variables
for(k = 0;k < MAXBIN;k++) { // only look up to 5KHz
		fftpower_GL[k] = re[k]*re[k] + im[k]*im[k];
	}
	// find peaks, fill peaklist1,amplist1 arrays
	ptr = 0;
	for(k=2;k < MAXBIN-5;k++) {
		if( (fftpower_GL[k] > fftpower_GL[k-1]) & (fftpower_GL[k] > fftpower_GL[k+1]) \
			& (fftpower_GL[k] > 2.5*fftpower_GL[k-2]) & (fftpower_GL[k] > 2.5*fftpower_GL[k+2]) ) {
				peakList1[ptr] = k;
				amplList1[ptr] = fftpower_GL[k];
				ptr++;
				if(fftpower_GL[k] > spectMax) spectMax = fftpower_GL[k]; // only consider peaks when finding the max
			}
	}
	 // ptr holds the number of elements in peaklist1
	numList1 = ptr;
	// now narrow down the peaks to those that are > -40dB of spectMax
	low_limit = spectMax*0.0001;// 40 dB is 80 dB in power = 1/10000
	ptr2 = 0;
	amplMax = 0.0;
	for(k=0;k < numList1;k++) {
		if(amplList1[k] > low_limit) {
			amplList[ptr2] = amplList1[k];
			peakList[ptr2] = peakList1[k];
			if(amplList[ptr2] > amplMax) {
				amplMax = amplList[ptr2];
				peakListMaxIndex = ptr2;
			}
			ptr2++;
		}
		
	}
	
	numPeaks = ptr2;

	// now find the first peak that exceeds some fraction of amplMax
	for(k=0;k < numPeaks;k++) {
		if(amplList[k] > 0.01*amplMax) { // -20 dB thresh
			firstPeaki = k;
			break;
		}
	}
	



	// now interpolate the first peak

	//for(k=0;k < numPeaks;k++) {
	Beta = sqrt(fftpower_GL[peakList[firstPeaki]]);
    Alpha = sqrt(fftpower_GL[peakList[firstPeaki]-1]);
    Gamma = sqrt(fftpower_GL[peakList[firstPeaki]+1]);
    frac_bin = 0.5*(Alpha-Gamma)/(Alpha-2.0*Beta + Gamma);
    firstPeakInterp = peakList[firstPeaki] + frac_bin;
    firstPeakInterp_div2 = firstPeakInterp/2.0; 


 	bin3 = (int)floor(3.0*firstPeakInterp_div2 + 0.5); // bin corresponding to 3rd harm of 1/2 firstPeakInterp
    bin5 = (int)floor(5.0*firstPeakInterp_div2 + 0.5); // bin corresponding to 5th harm of 1/2 firstPeakInterp
    bin4 = (int)floor(4.0*firstPeakInterp_div2 + 0.5); // bin corresponding to 4th harm of 1/2 firstPeakInterp
    bin6 = (int)floor(6.0*firstPeakInterp_div2 + 0.5); // bin corresponding to 6th harm of 1/2 firstPeakInterp
    ratioMF_zm1 = ratioMF;
    ratioMF = (fftpower_GL[bin3] + fftpower_GL[bin5])/(fftpower_GL[bin4] + fftpower_GL[bin6]);

 // hysterisis and persistance in divby decision

    if(ratioMF > 3.0 & ratioMF_zm1 > 3.0) {
        divby = 2.0;
    }
    if(ratioMF < 0.25 & ratioMF_zm1 < 0.25) {
        divby = 1.0;
    }
    fundInterp = firstPeakInterp/divby;


// try to estimate snr by looking in-between the low harmonic bins
    // Matlab
    // testBetween = fundInterp*(1.5:1:4.5);
    // testBetweenBins = floor(testBetween + 0.5) + 1;
    // testSig = fundInterp*(1.0:1.0:4.0);
    // testSigBins = floor(testSig + 0.5) + 1;
    noisePow = 0.0;
    sigPow = 0.0;
    for(kk=0;kk < 4;kk++) {
    	testBetweenBin = (int)floor( fundInterp*(1.5+(float)kk) + 0.5);
    	testSigBin = (int)floor( fundInterp*(1.0+(float)kk) + 0.5);
    	noisePow = noisePow + fftpower_GL[testBetweenBin];
    	sigPow = sigPow + fftpower_GL[testSigBin];

    }
    snr = sqrt(sigPow/noisePow);
  // create snr gate with hysterisis
    if(snr < 6.0) {
        snrGate = 0.0;
    }
    if(snr >= 20) {
        snrGate = 1.0;
    }

    if(snrGate == 1.0) {
    freqRawZm1 = freqRaw;
    freqRaw = fbin*fundInterp;
    freqRatio = freqRawZm1/freqRaw;
    if(freqRatio > 1.7 || freqRatio < 0.6) {
    	waitCount = 3;
    } 

    if(waitCount >= 1) waitCount = waitCount-1;
 

    if(waitCount > 0) {
    	freq = freqZm1; // hold before making a large jump
    } else {
    	freqZm1 = freq;
    	freq = freqRaw;
	}
	}



    
  
    freqGate = freq*snrGate; // this is what gets returned






// ************ This is the return to JS side!! ******************

wav_in[0] = freqGate; // return pitch using unused portion of input array
wav_in[1] = 1.0; // just to test to comms

//printf("block, freq = %d %f\n",block_GL,freq);

return(1);

}

#ifndef CONLY
}
#endif




/**************************************************/
void hannCompute(void) {

  	int k;
  	for(k = 0;k < FFTLEN;k++) {
  		win[k] = ( 0.5 * (1.0 - cos (2.0*PI*(float)k/(float)(FFTLEN-1))) );
  	}
}


/**************************************************/
void hannMult(void) {

  	int k;
  	for(k = 0;k < FFTLEN;k++) {
  		re[k] = re[k]*win[k];
  	}
}




void init_fft_twiddles(void)
{
int k;
	
for(k=0;k < FFTLEN2;k++) {
	twid_re[k] = cos(2.0*PI * (float)(k) / (float)FFTLEN);
	twid_im[k] = sin(2.0*PI * (float)(k) / (float)FFTLEN);

}


}
/**************************************************/
void fft(int fn)
//int fn, ff;
/*** re and im arrays are global, so they can be huge!! **/

/* fft(re, im fn, ff) performs Fast Fourier Transform.        */
/*  re[] and im[] contain the real and imaginary part of data  */
/*   ( After transform they contain the transformed data)        */
/*  fn is the number of point. ( must be power of 2 )            */
/*  ff = -1 for forward transform.       			 */
/*     =  1 for inverse transform.       			 */

{
	float tempr, tempi, fwr, fwi, ftheta;
	//float testfwr,testfwi;
	int fj, fi, fm, mmax, tmax, istep;


	/*  do Bit-Reversal */
	fj = 1;
	for ( fi =1; fi <=fn; fi++)
	{
		if(fi < fj)
		{
			tempr = re[fj-1];
			tempi = im[fj-1];
			re[fj-1] = re[fi-1];
			im[fj-1] = im[fi-1];
			re[fi-1] = tempr;
			im[fi-1] = tempi;
		}
		fm = fn / 2.0 ;
		while ( fj > fm )
		{
			fj = fj - fm;
			fm = (fm + 1 ) / 2;
		}
		fj = fj + fm;
	}

	/* Radix two frequency decimation algorithm  */
	mmax = 1;
	tmax = FFTLEN2;
		while ( mmax < fn )
	{
		istep = 2 * mmax;
		for ( fm = 1; fm <= mmax; fm++ )
		{

			fwr = twid_re[(fm-1)*tmax];
			fwi = -twid_im[(fm-1)*tmax];
			//ftheta = PI * (float)(-(fm - 1)) / (float)mmax;
			//fwr = cos(ftheta);
			//fwi = sin(ftheta);
			//printf("%d %f %f\n",(fm-1)*tmax,fwi,testfwi);
			for(fi = fm; fi <= fn; fi = fi + istep )
			{
				fj = fi + mmax;
				tempr = fwr * re[fj-1] - fwi * im[fj-1];
				tempi = fwr * im[fj-1] + fwi * re[fj-1];
				re[fj-1] = re[fi-1] - tempr;
				im[fj-1] = im[fi-1] - tempi;
				re[fi-1] = re[fi-1] + tempr;
				im[fi-1] = im[fi-1] + tempi;
			}
		}
		mmax = istep;
		tmax = tmax/2;
	}
	
}


/**************************************************/
#ifdef CONLY

float get_wav16()
{
int  i;
unsigned char templo,temphi;
//int fill;
int sum;
sum=0;
/**** 16-bit ************/

		templo = fgetc(wav_infile); // low byte
		if(feof(wav_infile)) printf("End 2\n");
		temphi = fgetc(wav_infile); // hi byte
		if(feof(wav_infile)) printf("End 3\n");
		//fill = 0;
		//if((temphi & 0x80) == 0x80) fill = 0xffff0000;
		//sum = ((temphi << 8) | templo) & 0x0000ffff;
		//sum = ((temphi << 8) | templo);
		temphi = (temphi + 128) & 0xffff; // invert sign bit
		sum = 256*(long)temphi + (long)templo;
		//sum = sum - 32768;
		//sum = sum | fill;
return((float)sum - 32768.0);
}
#endif





/**************************************************/
#ifdef CONLY

void parseWavHeader()
{
	/** note, format is LS byte followed by MS byte **/
int i;
int b0,b1,b2,b3,fs0,fs1,fs2,fs3,num_bytes;
for(i=1;i <= 44;i++) { 
	ch = fgetc(wav_infile); // 8 bits
	//printf("wavfile header %d %d\n",i,ch);

	if(i == 25) {fs0 = (int )ch;}
	if(i == 26) {fs1 = (int )ch;}
	if(i == 27) {fs2 = (int )ch;}
	if(i == 28) {fs3 = (int )ch;}


	if(i == 41) {b0 = (int )ch;}
	if(i == 42) {b1 = (int )ch;}
	if(i == 43) {b2 = (int )ch;}
	if(i == 44) {b3 = (int )ch;}

	if(i==23) {
		if((int)ch == 1) {
			num_channels = 1;
		}
		else {
			num_channels = 2;
		}
	}

	if(i== 35) {
		bits_per_sample = (int )ch;
	}

	}
	num_bytes = ((b3 << 24) | (b2 << 16) | (b1 << 8) | b0)/num_channels;
	if(bits_per_sample == 16) {
		num_samples = (int)(num_bytes/(2*num_channels));
	} else { // 24 bit assume
		num_samples = (int)(num_bytes/(3*num_channels));
		} 

	fs_wav = (float)(((fs3 << 24) | (fs2 << 16) | (fs1 << 8) | fs0));

	printf("fs_wav = %f\n",fs_wav);
	printf("Number of bytes in data section = %d\n",num_bytes);
	printf("Number of samples in data section = %d\n",num_samples);
	printf("num_chan= %d\n",num_channels);
	printf("bits-per-sample = %d\n",bits_per_sample);
}
#endif






