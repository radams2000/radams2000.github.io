#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "string.h"
#include "math.h"
#include "signal.h"


#define FFTLEN 2048
#define FFTLEN2 1024
#define HOP 1024
#define MAXBIN 216 // about 5KHz top bin for harmonic analysis
#define MAXBLOCK 1024
//#define CONLY 1

void fft(int);
void ifft(int);
void hannCompute(void);
void hannMult(void);
int fftmax(void) ;
float compute_fftmagdb(void);
void Array_sort(float *, int);
float find_freq(void);
void init_fft_twiddles(void);

#ifdef CONLY
void parseWavHeader();
float get_wav16(void);
int processAudioData(float *,int, int);
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
#endif

int num_points;
int numcycles;
int numpeaks_GL = 0;

float fs_wav;
float re[FFTLEN] = {0.0};
float im[FFTLEN] = {0.0};
float twid_re[FFTLEN2] = {0.0};
float twid_im[FFTLEN2] = {0.0};

float fftmagdb[FFTLEN] = {0.0};

float maxampl_db_GL = 0.0;
float freqs[MAXBLOCK] = {0.0};

float win[FFTLEN] = {0.0};
float PI = 3.141592653589793;
int testme[] = {1,2,3,4,5,6,7,8,9};
int testme2[] = {1,2,3,4,5,6,7,8,9};

float peaks_interp[MAXBIN+1] = {0.0};
float delta_peaks[MAXBIN] = {0.0};
float scale_bins2freq_GL = 1.0;
int init = 1;



#ifdef CONLY

int main() 
{
int i,k;
int num_return_samples = 1024;
float wavin[1000000] = {0.0};

int numblocks,block,count,block_count,tempi;
float left_in,right_in,wav_in;

outfile = fopen("freqout.txt","w");
char buff[] = "./wavs/rwa_runyon_a.wav";
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

numblocks = processAudioData(wavin,num_samples,num_return_samples);	

for(k=1;k < numblocks;k++) {
	fprintf(outfile,"%lf\n",freqs[k]);

}
fclose(outfile);
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

if(init==1) {
	init_fft_twiddles();
	init = 0;
}


fs_wav = 44.1e3; // temporary, will need a way to pass this in from the JS side

//int num_samples_wav = num_samples;
if(num_return_samples > MAXBLOCK) {
	printf("error, num return samples is < MAXBLOCK\n");
	exit(0);

}

/*********************start main routine ***************/
count = 0;
block_count = 0;
scale_bins2freq_GL = fs_wav/(float)FFTLEN;
printf("num_samples = %d\n",num_samples);
hannCompute();

for(i=num_samples_wav;i < num_samples;i++) { // zero out the return buffer
	wav_in[i] = 0.0;
}



for(i=0;i < num_samples_wav-FFTLEN;i++) {

	

	 if(count == (FFTLEN-1)) { // non-overlapping FFT's
	 	// load fft buffer fromm input buffer
	 	for(k=0; k < FFTLEN;k++) {re[k]=wav_in[i-FFTLEN+k];im[k] = 0.0;}
		hannMult(); // in place hann window on global re
		fft(FFTLEN); // in-place fft on global buffer variables
		maxampl_db_GL = compute_fftmagdb(); // fills array fftmagdb and reurns the max peak ampl in dB
		
		freqs[block_count] = find_freq();


		// ************ This is the return to JS side!! ******************

		wav_in[num_samples_wav+block_count] = freqs[block_count]; // return pitch using unused portion of input array
	 	
	 	// ***********************************************************


	 	count = 0;
	 	block_count = block_count + 1;

	} else {
		count = count + 1;
	}
	//wav_in[i] = 1.0; // debug to see if I can access this on the js side

  } /** end of for loop ***/


//for(k=0:k < 16;k++) wav_in[k] = 1;
printf("done FFT loop, max block count =   %d\n",block_count-1);

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
float compute_fftmagdb(void) { // operates on global FFT registers
	int k;
	float temp;
	for(k = 0;k < MAXBIN;k++) { // only look up to 5KHz
		temp = re[k]*re[k] + im[k]*im[k];
		fftmagdb[k] = 10.0*log10(temp + 1e-12);
		if(fftmagdb[k] > maxampl_db_GL) {
			maxampl_db_GL = fftmagdb[k];
		}
		
	}
	return(maxampl_db_GL);
}
	


/**************************************************/
float find_freq(void) {
	int k,kk,cond1,cond2,cond3,cond4,cond5,pass, halfbuff,num_keepers;; 
	float Beta,Alpha,Gamma,frac_bin,delta_peaks_median,av_keepers;
	numpeaks_GL=1;// manually stuff peaks_interp[0]
	peaks_interp[0] = 0.0; // for pure sine wave, generates a delta that can be used
	for(k=3;k < MAXBIN;k++) {
		cond1 = (fftmagdb[k] > fftmagdb[k+1]);
		cond2 = (fftmagdb[k] > fftmagdb[k-1]);
		cond3 = (fftmagdb[k] > (fftmagdb[k+2]+4.0));
		cond4 = (fftmagdb[k] > (fftmagdb[k-2]+4.0));
		cond5 = (fftmagdb[k] > (maxampl_db_GL-32.0));
		pass = (cond1 & cond2 & cond3 & cond4 & cond5);
		if(pass==1) { // interpolate
			Beta = fftmagdb[k];
        	Alpha = fftmagdb[k-1];
        	Gamma = fftmagdb[k+1];
        	frac_bin = 0.5*(Alpha-Gamma)/(Alpha-2.0*Beta + Gamma);
        	peaks_interp[numpeaks_GL] = (float)k + frac_bin;
        	//printf("%f\n",peaks_interp[numpeaks_GL]);
			numpeaks_GL++;
		}
		

	}
	//numpeaks_GL=numpeaks_GL-1;
	// now generate delta-peaks
	for(k=1;k <= numpeaks_GL;k++) {
		delta_peaks[k-1] = peaks_interp[k]-peaks_interp[k-1];

	}
	// now sort the delta-peaks array
	Array_sort(delta_peaks,numpeaks_GL);
	halfbuff = numpeaks_GL/2;
	delta_peaks_median = delta_peaks[halfbuff];

	// now reject outliers
	num_keepers = 0;
	av_keepers = 0.0;
	for(kk=0;kk < numpeaks_GL-1;kk++) { // note there is 1 less delta_peak than there are numpeaks
		if((delta_peaks[kk] > 0.95*delta_peaks_median) & (delta_peaks[kk] < 1.05*delta_peaks_median ) ) {
			num_keepers++;
			av_keepers+= delta_peaks[kk];
		}

	}
	av_keepers = av_keepers/(float)num_keepers;



	return(scale_bins2freq_GL*delta_peaks_median); // absolute frequency of the median delta-freq
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









// function to sort the array in ascending order
void Array_sort(float *array , int n)
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






