#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "string.h"
#include "math.h"
#include "signal.h"
#define FFTLEN 2048
#define FFTLEN2 1024
#define HOP 128
#define CONLY 1
void fft();
void fft(int, int);
void hannCompute(void);
void hannMult(void);
int fftmax(void) ;

#ifdef CONLY
void parseWavHeader();
float get_wav16(void);
int processAudioData(float *,int);
//int processAudioData(float * wav_in, int num_samples);
FILE *wav_infile; 
FILE *infile;
//FILE *debugfile;
//FILE *debug2file;

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
float win[FFTLEN] = {0.0};
int fftargmax[10000] = {0};
float PI = 3.141592653589793;



// int main()
// {


// // int k;
// // for(k=0;k < 10;k++) {
// // 	printf("%d\n",k);
// // }


// int i,ii,j,jj,k;

// int numblocks,block,count,block_count;
// float left_in,right_in,wav_in;


// char buff[] = "./bob/rwa_runyon_a.wav";
// if( (wav_infile = fopen(buff,"rb")) == NULL) {
// 	printf("cant find it\n");
// 	exit(0);
// 	}
// //num_samples = 100000;
// parseWavHeader(); /** sets global num_samples, num_channels, bits_per_sample **/
// printf("parsed wav file, num_samples = %d\n",num_samples);

// /*********************start main routine ***************/
// count = 0;
// block_count = 0;
// hannCompute();
// for(i=0;i < num_samples;i++) {

// 	if(num_channels == 1) { // mono
// 		wav_in = get_wav16();
// 	} else { // stereo
// 		wav_in = get_wav16();
// 		//wav_in = 0.5*(get_wav16() + get_wav16());
// 	}

// 	re[count] = wav_in;
// 	im[count] = 0.0;

// 	 if(count == (FFTLEN-1)) {
// 		hannMult(); // in place hann window on global re
// 		fft(FFTLEN,1); // in-place fft on global buffer variables
// 		fftargmax[block_count] = fftmax();
// 	 	count = 0;
// 	 	if(block_count == 90) {
// 	 		printf("%d\n",fftargmax[block_count]);
// 	 	}
	 	
// 	 	block_count = block_count + 1;

// 	} else {
// 		count = count + 1;
// 	}

//   } /** end of while loop ***/
// printf("done FFT loop\n");
// // fclose(wav_infile);


// return(1);
// }

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

#ifndef CONLY
extern "C" {
#endif

int processAudioData(float * wav_in, int num_samples)
{


int i;

int numblocks,count,block_count;
// float left_in,right_in,wav_in;


// char buff[] = "./bob/rwa_runyon_a.wav";
// if( (wav_infile = fopen(buff,"rb")) == NULL) {
// 	printf("cant find it\n");
// 	exit(0);
// 	}
// //num_samples = 100000;
// parseWavHeader(); /** sets global num_samples, num_channels, bits_per_sample **/
// printf("parsed wav file, num_samples = %d\n",num_samples);

/*********************start main routine ***************/
count = 0;
block_count = 0;
printf("num_samples = %d\n",num_samples);
hannCompute();
for(i=0;i < num_samples;i++) {

	// if(num_channels == 1) { // mono
	// 	wav_in = get_wav16();
	// } else { // stereo
	// 	wav_in = get_wav16();
	// 	//wav_in = 0.5*(get_wav16() + get_wav16());
	// }

	re[count] = wav_in[i];
	im[count] = 0.0;


	 if(count == (FFTLEN-1)) {
		hannMult(); // in place hann window on global re
		fft(FFTLEN,1); // in-place fft on global buffer variables
		fftargmax[block_count] = fftmax();
	 	count = 0;
	 	// if(block_count == 90) {
	 	// 	printf("%d\n",fftargmax[block_count]);
	 	// }
	 	
	 	block_count = block_count + 1;

	} else {
		count = count + 1;
	}

  } /** end of while loop ***/
printf("done FFT loop, max block count =   %d\n",block_count-1);
// fclose(wav_infile);


return(fftargmax[block_count-10]); // right now just returns the pitch of the 20th block
}

#ifndef CONLY
}
#endif



/***************************************************/
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



#ifdef CONLY
/*******************************/
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





void hannCompute(void) {

  	int k;
  	for(k = 0;k < FFTLEN;k++) {
  		win[k] = ( 0.5 * (1.0 - cos (2.0*PI*(float)k/(float)(FFTLEN-1))) );
  	}
}

void hannMult(void) {

  	int k;
  	for(k = 0;k < FFTLEN;k++) {
  		re[k] = re[k]*win[k];
  	}
}

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





