// new code 2/25/2020
#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "string.h"
#include "math.h"
#include "signal.h"


#define FFTLEN 4096
#define FFTLEN2 2048
#define MAXBIN 1024 // about 5KHz top bin for harmonic analysis
#define MAXBLOCK 1024
#define CONLY 1

void fft(int);
void ifft(int);
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

float peakListInterp[MAXBIN] = {0.0};
float deltaPeakList[MAXBIN]= {0.0};
int numMatch[6] = {0};
int waitCount = 0;


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
int ptr=0,ptr2=0,yy,divBy=1,numList1;
float temp=0.0,spectMax=0.0;
float Beta,Alpha,Gamma,frac_bin,amplMax=0.0,freqMaxBin,low_limit;
float spectMaxiInterp,minpeak,maxpeak,snr,changeRatio;
static float freq=0.0,freqZm1=0.0;
static float freqRawZm1 =0.0,freqRaw=0.0;



if(init==1) {
	init_fft_twiddles();
	init = 0;
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

	// now interpolate everything in the peaklist

	for(k=0;k < numPeaks;k++) {
		Beta = sqrt(fftpower_GL[peakList[k]]);
        Alpha = sqrt(fftpower_GL[peakList[k]-1]);
        Gamma = sqrt(fftpower_GL[peakList[k]+1]);
        frac_bin = 0.5*(Alpha-Gamma)/(Alpha-2.0*Beta + Gamma);
        peakListInterp[k] = (float)peakList[k] + frac_bin;
        if(k > 0) {
        	deltaPeakList[k-1] = peakListInterp[k]-peakListInterp[k-1];
        }

	}
	// find the interpolated index of the max peak. This will be the starting point of guessing the frequency
	spectMaxiInterp = peakListInterp[peakListMaxIndex];
	freqMaxBin = fbin*spectMaxiInterp;

	/*** MATLAB
	numMatch = zeros(6,1);
    for kk = 1:6 % test for guess that the answer is spectMaxiInterp/kk
        if freqMaxBin > 70*kk
            gate = (abs(deltaPeakList - spectMaxiInterp/kk) < 0.5); % within 1/2 bin
            numMatch(kk) = sum(gate); % number of deltas that are close to the assumed fundamental
        end
    end
    ***/

    for(k=0;k < 6;k++) {
    	numMatch[k] = 0;
		if(freqMaxBin > (float)(70*(k+1))) {
			for(kk=0;kk < numPeaks-1;kk++) { // there are numPeaks-1 deltas
				if(fabs(deltaPeakList[kk] - spectMaxiInterp/(float)(k+1)) < 0.5) {
					numMatch[k] = numMatch[k]+1;
				}

			} 

		}

    }

    /*** Matlab
    % deal with ties in numMatch
    for kk = 6:-1:2
        xx = find(numMatch(1:kk-1)==numMatch(kk));
        if ~isnan(xx)
            numMatch(xx) = 0;
        end
    end
    ***/

    // deal with ties in numMatch
    for(k=5;k > 0;k=k-1) {
    	for(kk=0;kk < k;kk++) {
    		if(numMatch[kk] == numMatch[k]) {
    			numMatch[kk]=0;
    		}
    	}
 	}

 	/*** Matlab
 	[yy,divby] = max(numMatch);
    if yy==0 % single tone??
        divby=1;
    end
    ***/
    yy = 0;
    divBy = 1;
    for(k=0;k < 6;k++) {
    	if(numMatch[k] > yy) { // find max if numMatch, this should be what we divide by
    		yy = numMatch[k];
    		divBy = k+1;
    	}
    }
    if(yy==0) divBy = 1; // single tone, no deltas

    /*** Matlab

    freqRaw(k) = fbin*spectMaxiInterp/divby;
    changeRatio = freqRaw(k)/(freqRaw(k-1) + 1e-6);
    if ((changeRatio > 1.7) || (changeRatio < 0.6))  % allows re-trigger
    %if ((changeRatio > 1.3) || (changeRatio < 0.76)) && (waitCount==0) % no re-trigger
        waitCount = 1;
    end
    if waitCount > 0
        freq(k) = freq(k-1); % hold
    else
        freq(k) = freqRaw(k); % update
    end

    ***/
	// skip for now
  	freqRawZm1 = freqRaw;
    freqRaw = fbin*spectMaxiInterp/(float)divBy;
    changeRatio = freqRaw/(freqRawZm1 + 1e-6);
    if ((changeRatio > 1.7) | (changeRatio < 0.6)) waitCount=1; //% allows re-trigger
    if(waitCount > 0) {
        freq = freqZm1; 
        }//% hold
    else {
    	freqZm1 = freq;
        freq = freqRaw; //% update
    }

    if(waitCount > 0) {
    	waitCount=waitCount-1;
	}

	/**** Matlab
	snr(k) = max(amplList1)/min(amplList1);
    if(snr(k) < 300)
        freq(k) = 0.0; % pause display
    end

    ***/
    minpeak = 1e12;
    maxpeak = 0.0;
    for(k=0;k < numList1;k++) {
    	if(amplList1[k] < minpeak) minpeak = amplList1[k];
    	if(amplList1[k] > maxpeak) maxpeak = amplList1[k];

    } 
	snr = sqrt(maxpeak/minpeak);

	if(snr < 300.0) {
		freq = 0.0; // signal to freeze display
	}
		



// ************ This is the return to JS side!! ******************

wav_in[0] = freq; // return pitch using unused portion of input array
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
void ifft(int fn)
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
	int fj, fi, fm, mmax, istep;


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
		while ( mmax < fn )
	{
		istep = 2 * mmax;
		for ( fm = 1; fm <= mmax; fm++ )
		{
			ftheta = PI * (float)((fm - 1)) / (float)mmax;
			fwr = cos(ftheta);
			fwi = sin(ftheta);
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
	}
	// scaling for inverse transform only, could skip??
	for ( fi = 0; fi < fn; fi++)
	{
		re[fi] = re[fi] / (float)fn;
		im[fi] = im[fi] / (float)fn;
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






