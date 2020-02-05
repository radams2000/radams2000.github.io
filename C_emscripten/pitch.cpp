#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "string.h"
#include "math.h"
#include "signal.h"


#define FFTLEN 2048
#define FFTLEN2 1024
#define HOP 128
//#define CONLY 1
void fft();
void fft(int, int);
void hannCompute(void);
void hannMult(void);
int fftmax(void) ;
float compute_fftmag(void);
void compute_fftphs(void);
void compute_fftphs_zm1(void);
int compute_freqof_peaks(void);
//void Array_sort(float *, int);
void Array_sort(float array[] , int n);





#ifdef CONLY
void parseWavHeader();
float get_wav16(void);
int processAudioData(float *,int);
FILE *wav_infile; 
FILE *infile;

char inbuf[300];
char fname[60];
char fname_r[40];
char fname_l[40];
unsigned char ch;
int bits_per_sample,num_samples;
int num_channels,block_align;
#endif

int num_points;
int numcycles;

float fs_wav;
float re[FFTLEN] = {0.0};
float im[FFTLEN] = {0.0};
float fftmag[FFTLEN] = {0.0};
float fftphs[FFTLEN] = {0.0};
float fftphs_zm1[FFTLEN] = {0.0};
float peakfreqs_GL[FFTLEN2] = {0.0};
float peakfreqs_delta_GL[FFTLEN2] = {0.0};
float max_ampl_GL = 0.0;
int numpeaks_GL = 0;
float truefreq_GL=0;
float FS = 48e3;
float rms_sig = 0.0;



float win[FFTLEN] = {0.0};
int fftargmax[10000] = {0};
float PI = 3.141592653589793;
int testme[] = {1,2,3,4,5,6,7,8,9};
int testme2[] = {1,2,3,4,5,6,7,8,9};



#ifdef CONLY

int main() 
{
int i,k;
float wavin[1000000] = {0.0};

int numblocks,block,count,block_count,tempi;
float left_in,right_in,wav_in;

char buff[] = "./bob/rwa_runyon_a.wav";
if( (wav_infile = fopen(buff,"rb")) == NULL) {
	printf("cant find it\n");
	exit(0);
	}
parseWavHeader(); /** sets global num_samples, num_channels, bits_per_sample **/
printf("parsed wav file, num_samples = %d\n",num_samples);

/*********************start main routine ***************/
count = 0;
block_count = 0;
hannCompute();
for(i=0;i < num_samples;i++) {

	if(num_channels == 1) { // mono
		wavin[i] = get_wav16();
	} else { // stereo
		wavin[i] = get_wav16();
	    wavin[i] = get_wav16();
	}

	

  } /** end of while loop ***/

tempi = processAudioData(wavin,num_samples);	
}
#endif



/****************************************************/

#ifndef CONLY
extern "C" {
#endif

int processAudioData(float * wav_in, int num_samples,int num_return_samples) // returns a pointer to an array
{

int i,k,numblocks,count,block_count;

int num_samples_wav = num_samples-num_return_samples;

//int num_samples_wav = num_samples;


/*********************start main routine ***************/
count = 0;
block_count = 0;
printf("num_samples = %d\n",num_samples);
hannCompute();

for(i=num_samples_wav;i < num_samples;i++) { // zero out the return buffer
	wav_in[i] = 0.0;
}

for(i=0;i < num_samples_wav-FFTLEN;i++) {rms_sig = rms_sig + wav_in[i]*wav_in[i];}

rms_sig = sqrt(rms_sig/(float)num_samples);

for(i=0;i < num_samples_wav-FFTLEN;i++) {

	

	 if(count == (FFTLEN-1)) { // non-overlapping FFT's
	 	// load fft buffer fromm input buffer
	 	for(k=0; k < FFTLEN;k++) {re[k]=wav_in[i-FFTLEN+k];im[k] = 0.0;}
		hannMult(); // in place hann window on global re
		fft(FFTLEN,-1); // in-place fft on global buffer variables
		max_ampl_GL = compute_fftmag();
		numpeaks_GL = 0;
		if(max_ampl_GL > 50.0) {
			compute_fftphs_zm1();
	    	for(k=0; k < FFTLEN;k++) {re[k]=wav_in[i-FFTLEN+k+1];im[k] = 0.0;}
	    	hannMult();
			fft(FFTLEN,-1); 
			compute_fftphs();
			numpeaks_GL = compute_freqof_peaks();
		}
		fftargmax[block_count] = fftmax();
		wav_in[num_samples_wav+block_count] = (float)fftmax(); // return pitch using unused portion of input array
		//wav_in[num_samples_wav+block_count] = truefreq_GL;
		//wav_in[num_samples_wav+block_count] = (float)block_count; // debug
		printf("truefreq = %lf\n",truefreq_GL);
		//wav_in[num_samples_wav+block_count] = (float)numpeaks_GL;

	 	count = 0;
	 	block_count = block_count + 1;

	} else {
		count = count + 1;
	}
	//wav_in[i] = 1.0; // debug to see if I can access this on the js side

  } /** end of for loop ***/


//for(k=0:k < 16;k++) wav_in[k] = 1;
printf("done FFT loop, max block count =   %d\n",block_count-1);
//printf("wav_in[10000] =   %f\n",wav_in[10000]);

//printf("pointer to testme2 =   %p\n",&testme);
//for(k=num_samples_wav;k < num_samples;k++) wav_in[k] = 1.0;

//return(fftargmax[block_count-10]); // right now just returns the pitch of the 20th block
return block_count; 

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


/**************************************************/
float compute_fftmag(void) { // operates on global FFT registers
	int k;
	float temp;
	for(k = 0;k < FFTLEN2;k++) { // only look at 1st half of FFT buffer
		temp = re[k]*re[k] + im[k]*im[k];
		fftmag[k] = sqrt(temp);
		if(fftmag[k] > max_ampl_GL) {
			max_ampl_GL = fftmag[k];
		}
		
	}
	return(max_ampl_GL);
}
	


/**************************************************/
void compute_fftphs_zm1(void) { // operates on global FFT registers
	int k;
	for(k = 0;k < FFTLEN2;k++) { // only look at 1st half of FFT buffer
		fftphs_zm1[k] = atan2(im[k],re[k]);
	}
	return;
}

void compute_fftphs(void) { // operates on global FFT registers
	int k;
	for(k = 0;k < FFTLEN2;k++) { // only look at 1st half of FFT buffer
		fftphs[k] = atan2(im[k],re[k]);
	}
	return;
}

int compute_freqof_peaks(void) {
	int j,k;
	int numpeaks = 1; // manually set the 1st peak to bin 0
	peakfreqs_GL[0] = 0.0; // always pretend 0 is a peak so you count the fundamental as the 1st difference
	for(k = 1;k < FFTLEN2;k++) { 
		if((fftmag[k] > fftmag[k-1]) && (fftmag[k] > fftmag[k+1]) && (fftmag[k] > 0.03*max_ampl_GL)) {
			peakfreqs_GL[numpeaks] = (fftphs[k]-fftphs_zm1[k])*FS/(2.0*PI);
			if(peakfreqs_GL[numpeaks] < 0.0) {peakfreqs_GL[numpeaks]+=FS;}
			peakfreqs_delta_GL[numpeaks-1] = peakfreqs_GL[numpeaks] - peakfreqs_GL[numpeaks-1];
			numpeaks++;
		}

	}
	if(numpeaks > 4) {

		Array_sort(peakfreqs_delta_GL,numpeaks-1);
		truefreq_GL = peakfreqs_delta_GL[(numpeaks-1)/2]; // median of the delta's


	}
	return(numpeaks);
}



/**************************************************/
int fftmax(void) { // operates on global FFT registers
	int k;
	float maxval = 0;
	int argmax = 0;
	float temp;
	for(k = 0;k < FFTLEN2;k++) { // only look at 1st half of FFT buffer
		temp = re[k]*re[k] + im[k]*im[k];
		if(temp > maxval) {
			maxval = temp;
			argmax = k;
		}
	}
	return(argmax);
}


/**************************************************/
void fft(int fn, int ff)
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
			ftheta = PI * (float)(ff * (fm - 1)) / (float)mmax;
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
	if ( ff > 0 )
	{
		for ( fi = 0; fi < fn; fi++)
		{
			re[fi] = re[fi] / (float)fn;
			im[fi] = im[fi] / (float)fn;
		}
	}
}




// function to sort the array in ascending order
void Array_sort(float array[] , int n)
{ 
    // declare some local variables
    int i=0 , j=0 ;
    float temp=0.0;

    for(i=0 ; i<n ; i++)
    {
        for(j=0 ; j<n-1 ; j++)
        {
            if(array[j]>array[j+1])
            {
                temp        = array[j];
                array[j]    = array[j+1];
                array[j+1]  = temp;
            }
        }
    }

    // printf("\nThe array after sorting is..\n");
    // for(i=0 ; i<n ; i++)
    // {
    //     printf("\narray_1[%d] : %d",i,array[i]);
    // }
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






